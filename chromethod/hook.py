"""
Hook 验证层 — 物理化学约束规则
每条规则返回 (status, message)
  PASS   = 通过
  WARN   = 建议改进
  REJECT = 否决
"""
from dataclasses import dataclass, field
from typing import Optional

@dataclass
class HookResult:
    status: str      # PASS | WARN | REJECT
    rule: str        # 规则名称
    message: str     # 中文解释

@dataclass
class HookReport:
    results: list = field(default_factory=list)

    @property
    def rejected(self) -> bool:
        return any(r.status == "REJECT" for r in self.results)

    @property
    def warns(self) -> list:
        return [r for r in self.results if r.status == "WARN"]

    @property
    def errors(self) -> list:
        return [r for r in self.results if r.status == "REJECT"]

    def summary(self) -> str:
        lines = []
        for r in self.results:
            icon = {"PASS":"✅","WARN":"⚠️","REJECT":"❌"}[r.status]
            lines.append(f"  {icon} [{r.rule}] {r.message}")
        return "\n".join(lines)


def _col_polarity(col: dict) -> float:
    """从柱信息估算极性分数"""
    from .features import phase_to_polarity
    return phase_to_polarity(col.get("brand",""), col.get("phase",""))


def _max_oven_temp(col: dict) -> float:
    prog = col.get("oven_program", [])
    temps = []
    for s in prog:
        t = s.get("temp_c", 50) or 50
        temps.append(t)
        if s.get("target_c"):
            temps.append(s["target_c"])
    return max(temps) if temps else 200


def _col_limit_temp(col: dict) -> float:
    """从JSON字段读取柱温上限（优先程升限，其次通用限）"""
    return col.get("temp_limit_c") or col.get("temp_limit_isothermal_c") or 300


def _has_fg(substance: dict, *keywords) -> bool:
    fgs = [f.lower() for f in substance.get("functional_groups", [])]
    return any(k.lower() in " ".join(fgs) for k in keywords)


def _has_any_fg(substances: list, *keywords) -> bool:
    return any(_has_fg(s, *keywords) for s in substances)


def _max_bp(substances: list) -> float:
    return max((s.get("bp_c", 100) for s in substances), default=100)


def _detector_type(method: dict) -> str:
    return method["detectors"][0]["type"]


# ════════════════════════════════════════════════════════
# 规则 1: 极性-柱匹配
# ════════════════════════════════════════════════════════
def r_polarity_column(substances: list, method: dict) -> HookResult:
    col = method["columns"][0]
    col_pol = _col_polarity(col)
    max_pol = max((s.get("polarity", 3) for s in substances), default=3)

    if max_pol >= 5:  # 强极性
        if col_pol < 0.3:
            return HookResult("REJECT", "极性-柱匹配",
                f"强极性样品(max={max_pol})不能用非极性柱({col['brand']} 极性={col_pol:.2f})。需 PEG/WAX 或 FFAP 柱")
        if col_pol < 0.7:
            return HookResult("WARN", "极性-柱匹配",
                f"强极性样品建议用极性柱(PEG/WAX/FFAP)，当前{col['brand']}极性={col_pol:.2f}")
        return HookResult("PASS", "极性-柱匹配", f"极性匹配: 样品{max_pol}级 ↔ {col['brand']}")

    if max_pol >= 4:
        if col_pol < 0.3:
            return HookResult("WARN", "极性-柱匹配",
                f"较强极性样品不宜用非极性柱({col['brand']})，考虑 DB-1701/DB-624 或更极性柱")

    return HookResult("PASS", "极性-柱匹配", f"极性匹配: 样品{max_pol}级 ↔ {col['brand']}")


# ════════════════════════════════════════════════════════
# 规则 2: 柱温上限约束
# ════════════════════════════════════════════════════════
def r_temperature_ceiling(substances: list, method: dict) -> HookResult:
    col = method["columns"][0]
    max_oven = _max_oven_temp(col)
    col_limit = _col_limit_temp(col)
    max_bp = max((s.get("bp_c", 100) for s in substances), default=100)

    # 温度余量判断: 判定方法是程升还是等温
    prog = col.get("oven_program", [])
    is_isothermal = len(prog) == 1 and prog[0]["type"] == "hold"

    # 程升: 留10°C余量; 等温: 留20°C余量
    margin = 20 if is_isothermal else 10
    if max_oven > col_limit - margin:
        return HookResult("REJECT", "温度上限",
            f"柱箱温度{max_oven}°C超过{col['brand']}安全上限{col_limit}°C"
            f"({'等温' if is_isothermal else '程升'}留{margin}°C余量)")

    # 柱温必须 > 最高沸点 (保证出峰)
    if max_oven < max_bp:
        return HookResult("WARN", "温度上限",
            f"最高柱温{max_oven}°C < 样品最高沸点{max_bp}°C，高沸点组分可能不出峰")

    # 等温方法检查
    if is_isothermal:
        if (max_bp - prog[0].get("temp_c", 50)) > 30:
            return HookResult("WARN", "温度上限",
                f"等温方法(温度{prog[0]['temp_c']}°C)不适合沸点跨度大的样品(最高{max_bp}°C)")

    return HookResult("PASS", "温度上限", f"柱箱温度范围合适(最大{max_oven}°C)")


# ════════════════════════════════════════════════════════
# 规则 3: 检测器-元素路由
# ════════════════════════════════════════════════════════
def r_detector_element_routing(substances: list, method: dict) -> HookResult:
    dt = _detector_type(method)

    has_p = _has_any_fg(substances, "磷酸酯", "有机磷", "P=O")
    has_s = _has_any_fg(substances, "硫醇", "硫醚", "-S-", "-SH", "H₂S")
    has_cl = _has_any_fg(substances, "-Cl", "有机氯")
    has_br = _has_any_fg(substances, "-Br", "有机溴")
    has_n = _has_any_fg(substances, "-NH2", "氮杂环", "硝基", "嘌呤", "吡啶")

    # 含磷 → 必须 FPD 或 NPD，不能用 ECD
    if has_p:
        if dt == "ECD":
            return HookResult("REJECT", "检测器路由",
                "含磷化合物不能用ECD检测！需FPD(推荐)或NPD")
        if dt not in ("FPD", "PFPD", "NPD"):
            return HookResult("WARN", "检测器路由",
                f"含磷化合物建议用FPD/NPD，当前检测器{dt}可能灵敏度不足")

    # 含硫 → 推荐 SCD 或 FPD
    if has_s:
        if dt not in ("SCD", "FPD", "PFPD"):
            return HookResult("WARN", "检测器路由",
                f"含硫化合物建议用SCD或FPD，当前检测器{dt}选择性/灵敏度可能不足")

    # 含卤素 → ECD 最佳
    if (has_cl or has_br) and dt == "ECD":
        return HookResult("PASS", "检测器路由", "含卤素化合物→ECD检测器，最佳匹配")

    # 含氮 → 推荐 NPD 或 MS
    if has_n and dt not in ("NPD", "MS", "FPD"):
        return HookResult("WARN", "检测器路由",
            f"含氮化合物建议用NPD或MS，当前检测器{dt}")

    return HookResult("PASS", "检测器路由", f"检测器{dt}与样品元素兼容")


# ════════════════════════════════════════════════════════
# 规则 4: 检出限-检测器匹配
# ════════════════════════════════════════════════════════
def r_detection_limit(substances: list, method: dict) -> HookResult:
    dt = _detector_type(method)
    min_det = min((s.get("detection_limit", "ppm") for s in substances),
                  key=lambda x: {"ppm":0,"ppb":1,"ppt":2}.get(x,0))

    if min_det == "ppt":
        if dt == "FID":
            return HookResult("WARN", "检出限",
                f"ppt级检出限要求，FID灵敏度可能不够。建议ECD(卤代)/MS(SIM)/FPD")
        if dt == "TCD":
            return HookResult("REJECT", "检出限",
                "TCD检出限仅~ppm级，无法满足ppt级要求")

    if min_det == "ppb":
        if dt == "TCD":
            return HookResult("REJECT", "检出限",
                "TCD检出限仅~ppm级，无法满足ppb级要求")

    return HookResult("PASS", "检出限", f"检测器{dt}满足{min_det}级检出限")


# ════════════════════════════════════════════════════════
# 规则 5: 载气-检测器配对
# ════════════════════════════════════════════════════════
def r_carrier_detector(substances: list, method: dict) -> HookResult:
    gas = method["carrier_gas"]["type"]
    dt = _detector_type(method)

    if gas == "N2" and dt in ("FID", "MS"):
        return HookResult("WARN", "载气-检测器",
            f"N2载气配{dt}检测器，柱效和灵敏度低于He/H2。He推荐用于MS，H2推荐用于FID")

    if gas == "H2" and dt == "MS":
        return HookResult("WARN", "载气-检测器",
            "H2载气用于MS需注意真空系统负荷和氢脆风险")

    return HookResult("PASS", "载气-检测器", f"载气{gas}与检测器{dt}兼容")


# ════════════════════════════════════════════════════════
# 规则 6: 多组分程升补偿
# ════════════════════════════════════════════════════════
def r_multicomponent_ramp(substances: list, method: dict) -> HookResult:
    n = len(substances)
    prog = method["columns"][0].get("oven_program", [])
    ramp_steps = [s for s in prog if s["type"] == "ramp"]

    if n > 10:
        # 检查是否有足够的升温段
        if len(ramp_steps) < 2:
            return HookResult("WARN", "多组分程升",
                f"{n}个组分建议≥2段升温程序以优化分离")
        # 检查升温速率
        fast_ramps = [s for s in ramp_steps if s.get("rate", 0) > 10]
        if fast_ramps:
            return HookResult("WARN", "多组分程升",
                f"{n}个组分建议升温速率≤10°C/min，当前{max(s['rate'] for s in fast_ramps)}°C/min偏快")

    if n > 5:
        for s in prog:
            if s["type"] == "hold" and s.get("time_min", 0) is not None and s["time_min"] < 2:
                return HookResult("WARN", "多组分程升",
                    f"{n}个组分建议每段保持≥2min，当前{s['temp_c']}°C段仅{s['time_min']}min")

    return HookResult("PASS", "多组分程升", f"{n}个组分，升温程序合适")


# ════════════════════════════════════════════════════════
# 规则 7: 进样模式检查
# ════════════════════════════════════════════════════════
def r_injection_mode(substances: list, method: dict) -> HookResult:
    inlet = method.get("inlet", {})
    mode = inlet.get("mode", "")
    max_bp = _max_bp(substances)
    prep = method["preparation"]["technique"]

    # 高沸点 + 不分流
    if max_bp > 250:
        if "不分流" not in str(mode) and "splitless" not in str(mode).lower():
            return HookResult("WARN", "进样模式",
                f"高沸点样品(max BP={max_bp}°C)，建议不分流进样以提高灵敏度")

    # 低沸点VOC → 顶空或吹扫（气体样品除外）
    has_gas = any(s.get("phase", "") == "气体" for s in substances)
    has_voc = any((s.get("bp_c", 100) or 100) < 100 for s in substances)
    is_direct = "直接" in prep and "顶空" not in prep and "吹扫" not in prep
    if has_voc and is_direct and not has_gas:
        return HookResult("WARN", "进样模式",
            "低沸点VOC样品，考虑顶空或吹扫捕集进样以减少基质干扰")

    return HookResult("PASS", "进样模式", f"进样模式{mode}适合当前样品")


# ════════════════════════════════════════════════════════
# 规则 8: 基质-前处理匹配
# ════════════════════════════════════════════════════════
def r_matrix_prep(substances: list, method: dict) -> HookResult:
    prep = method["preparation"]["technique"]
    # 从方法适用基质推断
    app = method["applicability"]
    matrix = app.get("matrix", "")

    has_aqueous = any(k in matrix for k in ("水", "aqueous", "blood", "血液"))
    has_solid = any(k in matrix for k in ("固体", "solid", "粉", "材料", "土壤", "tissue"))

    if has_aqueous:
        if "顶空" not in prep and "吹扫" not in prep and "Pur&Trap" not in prep:
            return HookResult("WARN", "基质-前处理",
                "水溶液基质建议顶空或吹扫捕集进样，直接进样会污染衬管和柱")

    if has_solid:
        if prep in ("直接进样", "无"):
            return HookResult("WARN", "基质-前处理",
                "固体基质建议顶空/HS-Trap/热脱附等前处理，直接进样会堵塞衬管")

    return HookResult("PASS", "基质-前处理", "前处理方式与基质兼容")


# ════════════════════════════════════════════════════════
# 规则 9: 酚类 → FFAP 柱
# ════════════════════════════════════════════════════════
def r_phenol_ffap(substances: list, method: dict) -> HookResult:
    col = method["columns"][0]
    col_pol = _col_polarity(col)
    is_phenol = _has_any_fg(substances, "酚羟基", "Ar-OH")

    if is_phenol:
        if col_pol < 0.6:
            return HookResult("WARN", "酚类-柱选择",
                f"酚类化合物在非极性柱上可能拖尾，建议WAX/FFAP柱(当前{col['brand']}极性={col_pol:.2f})")

    return HookResult("PASS", "酚类-柱选择", "")


# ════════════════════════════════════════════════════════
# 规则 10: 羧酸 → FFAP 或衍生化
# ════════════════════════════════════════════════════════
def r_acid_column(substances: list, method: dict) -> HookResult:
    col = method["columns"][0]
    col_pol = _col_polarity(col)
    is_acid = _has_any_fg(substances, "-COOH", "羧酸")

    if is_acid:
        if col_pol < 0.7:
            return HookResult("WARN", "羧酸-柱选择",
                f"羧酸类在非极性/中极性柱上易拖尾，建议FFAP柱(Stabilwax-DA)或甲酯化衍生化(当前{col['brand']})")

    return HookResult("PASS", "羧酸-柱选择", "")


# ════════════════════════════════════════════════════════
# 规则 11: 高热不稳定 → 冷柱头
# ════════════════════════════════════════════════════════
def r_thermal_stability(substances: list, method: dict) -> HookResult:
    inlet = method.get("inlet", {})
    inlet_temp = inlet.get("temp_c", 250) or 250
    max_bp = _max_bp(substances)

    # 简单代理: 沸点 > 350°C 的高分子量物质可能热不稳定
    if max_bp > 350:
        if inlet_temp > 280:
            return HookResult("WARN", "热稳定性",
                f"高沸点物质(max BP={max_bp}°C)进样口{inlet_temp}°C可能导致热分解，考虑冷柱头进样(COC)")

    return HookResult("PASS", "热稳定性", "")


# ════════════════════════════════════════════════════════
# 规则 12: 水溶液-柱匹配
# ════════════════════════════════════════════════════════
def r_aqueous_column(substances: list, method: dict) -> HookResult:
    col = method["columns"][0]
    col_pol = _col_polarity(col)
    app = method["applicability"]
    matrix = app.get("matrix", "")

    if any(k in matrix for k in ("水", "aqueous", "水溶液")):
        if col_pol < 0.2:
            return HookResult("REJECT", "水溶液-柱",
                "水溶液样品不能用100%二甲基柱(DB-1类)，水峰会严重拖尾掩盖分析物。需中极性以上柱")

    return HookResult("PASS", "水溶液-柱", "")


# ── 规则集 ──
RULES = [
    ("极性-柱匹配", r_polarity_column),
    ("温度上限", r_temperature_ceiling),
    ("检测器路由", r_detector_element_routing),
    ("检出限", r_detection_limit),
    ("载气-检测器", r_carrier_detector),
    ("多组分程升", r_multicomponent_ramp),
    ("进样模式", r_injection_mode),
    ("基质-前处理", r_matrix_prep),
    ("酚类-柱选择", r_phenol_ffap),
    ("羧酸-柱选择", r_acid_column),
    ("热稳定性", r_thermal_stability),
    ("水溶液-柱", r_aqueous_column),
]


def validate(substances: list[dict], method: dict) -> HookReport:
    """对推荐方法运行全部验证规则"""
    report = HookReport()
    for name, rule_fn in RULES:
        result = rule_fn(substances, method)
        # 跳过空 PASS (无消息的)
        if result.status == "PASS" and not result.message.strip():
            continue
        report.results.append(result)
    return report

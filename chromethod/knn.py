"""
k-NN 检索：余弦距离 + Top-K 排序
"""
import json
import glob
from .features import encode_method, encode_substance, weighted_cosine, FUNC_GROUPS


def load_methods(data_dir: str = "data") -> list[dict]:
    """加载所有标准方法 JSON"""
    methods = []
    for path in sorted(glob.glob(f"{data_dir}/GC-STD-*.json")):
        with open(path, encoding="utf-8") as f:
            m = json.load(f)
            m["_file"] = path
            methods.append(m)
    return methods


def build_index(methods: list[dict]) -> list[list[float]]:
    """构建方法特征向量索引"""
    return [encode_method(m) for m in methods]


def search(substance: dict, methods: list[dict],
           vectors: list[list[float]], k: int = 3) -> list[tuple[float, dict]]:
    """检索 Top-K 最相似方法

    Returns: [(similarity, method_dict), ...]
    """
    query_vec = encode_substance(substance)
    scores = []
    # 官能团维度索引 (7 ~ 7+17)
    FG_START, FG_END = 7, 7 + len(FUNC_GROUPS)

    skip_kw = ("空白", "污染验证")
    for i, (method, method_vec) in enumerate(zip(methods, vectors)):
        cat = method["applicability"]["analyte_category"]
        if any(k in cat for k in skip_kw):
            continue

        sim = weighted_cosine(query_vec, method_vec)

        # 官能团缺失惩罚: 查询有的FG, 方法没有 → 降权
        query_fg_count = sum(1 for j in range(FG_START, FG_END) if query_vec[j] > 0.5)
        if query_fg_count > 0:
            matched = sum(1 for j in range(FG_START, FG_END)
                          if query_vec[j] > 0.5 and method_vec[j] > 0.5)
            fg_ratio = matched / query_fg_count
            # fg_ratio=0 → penalty=0.6,  fg_ratio=1 → penalty=1.0
            penalty = 0.6 + 0.4 * fg_ratio
            sim *= penalty

        scores.append((sim, i))

    scores.sort(key=lambda x: -x[0])
    return [(s, methods[i]) for s, i in scores[:k]]


def format_result(sim: float, method: dict) -> str:
    """格式化一条推荐结果"""
    col = method["columns"][0]
    det = method["detectors"][0]
    cg = method["carrier_gas"]
    inlet = method.get("inlet", {})
    prog = col.get("oven_program", [])

    # 载气行: 种类 + 流速/线速度/柱前压
    gas_parts = [f"载气: {cg['type']}"]
    if cg.get("flow_rate_ml_min"):
        gas_parts.append(f"{cg['flow_rate_ml_min']} mL/min")
    if cg.get("linear_velocity_cm_s"):
        gas_parts.append(f"({cg['linear_velocity_cm_s']} cm/s)")
    if cg.get("pressure_psi"):
        gas_parts.append(f"柱前压 {cg['pressure_psi']} psi")
    gas_str = " ".join(gas_parts)

    msgs = [
        f"  相似度: {sim:.3f}",
        f"  方法: {method['source']['standard']} — {method['source']['name']}",
        f"  色谱柱: {col['brand']} {col['length_m']}m×{col['id_mm']}mm×{col['film_um']}µm ({col['phase']})",
        gas_str,
        f"  检测器: {det['type']}",
    ]
    # 进样参数
    if inlet:
        inlet_parts = []
        if inlet.get("mode"):
            inlet_parts.append(f"进样: {inlet['mode']}")
        if inlet.get("split_ratio"):
            inlet_parts.append(f"分流比 {inlet['split_ratio']}")
        if inlet.get("injection_volume_ul"):
            inlet_parts.append(f"进样量 {inlet['injection_volume_ul']} µL")
        if inlet.get("temp_c"):
            inlet_parts.append(f"进样口 {inlet['temp_c']}°C")
        if inlet.get("liner"):
            inlet_parts.append(f"衬管: {inlet['liner']}")
        msgs.append("  " + " | ".join(inlet_parts))
    # 检测器附加参数
    det_parts = []
    if det.get("solvent_delay_min"):
        det_parts.append(f"溶剂延迟 {det['solvent_delay_min']} min")
    if det.get("transfer_line_temp_c"):
        det_parts.append(f"传输线 {det['transfer_line_temp_c']}°C")
    if det.get("source_temp_c"):
        det_parts.append(f"离子源 {det['source_temp_c']}°C")
    if det.get("temp_c"):
        det_parts.append(f"检测器 {det['temp_c']}°C")
    if det_parts:
        msgs.append("  " + " | ".join(det_parts))
    msgs.append(f"  前处理: {method['preparation']['technique']}")
    for step in prog[:4]:
        if step["type"] == "hold":
            msgs.append(f"    {step['temp_c']}°C 保持 {step['time_min']}min")
        else:
            msgs.append(f"    → {step['target_c']}°C @ {step['rate']}°C/min")
    return "\n".join(msgs)

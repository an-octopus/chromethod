"""
特征编码 v2：化合物 + 方法 → 数值向量
共享特征空间，加权余弦相似度
"""
import math
import logging

logger = logging.getLogger(__name__)

# ── 官能团定义 (17类) ──
FUNC_GROUPS = [
    "-OH 脂肪醇羟基",         # 0
    "酚羟基 (Ar-OH)",         # 1
    "-COOH 羧酸",             # 2
    "-COO- 酯基",             # 3
    "-C=O 酮/醛羰基",         # 4
    "-NH2 伯/仲胺",           # 5
    "含氮杂环 (吡啶/嘌呤等)",   # 6
    "-NO2 硝基",              # 7
    "磷酸酯/有机磷 (P=O/S)",   # 8
    "-Cl 有机氯",             # 9
    "-Br 有机溴",             # 10
    "-F 有机氟",              # 11
    "硫醇 -SH (H₂S/硫醇)",    # 12
    "硫醚 -S-",               # 13
    "-O- 醚",                 # 14
    "aromatic 芳香环",         # 15
    "alkene 烯烃",             # 16
]

POLARITY_MAP = {1: 0.0, 2: 0.25, 3: 0.5, 4: 0.75, 5: 1.0}
DET_MAP = {"ppm": 0.0, "ppb": 0.5, "ppt": 1.0}


def encode_substance(substance: dict) -> list[float]:
    """化合物 → 特征向量 (26维)"""
    vec = []

    # 0-2: 物态 one-hot
    phase = substance.get("phase", "液体")
    vec.extend([1.0 if phase == p else 0.0 for p in ("气体", "液体", "固体")])

    # 3: 极性 (柱选择系数高)
    vec.append(POLARITY_MAP.get(substance.get("polarity", 3), 0.5))

    # 4: 沸点归一化 (温度选择系数高)
    bp = substance.get("bp_c", 100)
    vec.append(min(max((bp - 30) / 370.0, 0), 1.0))

    # 5: LogP 归一化
    logp = substance.get("logp", 0)
    vec.append(min(max((logp + 5) / 15.0, 0), 1.0))

    # 6: 碳数归一化
    vec.append(min(substance.get("carbon_count", 1) / 40.0, 1.0))

    # 7-18: 官能团 one-hot
    fgroups = set(substance.get("functional_groups", []))
    for fg in FUNC_GROUPS:
        vec.append(1.0 if fg in fgroups else 0.0)

    # 19-20: 氢键
    vec.append(min(substance.get("hbd", 0) / 10.0, 1.0))
    vec.append(min(substance.get("hba", 0) / 15.0, 1.0))

    # 21: 检出限
    vec.append(DET_MAP.get(substance.get("detection_limit", "ppm"), 0.0))

    # 22: 物质初始沸点代理 (低=易挥发)
    vec.append(1.0 - vec[4])  # 沸点低 → 值高 = 挥发性VOC

    # 23: 是否需要不分流 (高沸点/痕量→不分流)
    vec.append(1.0 if bp > 150 else 0.0)

    # 24: 含杂原子 (Cl/Br/F/S/N/P → 特殊检测器)
    hetero_kw = ("-Cl", "-Br", "-F", "硫醇", "硫醚", "-NH2", "-NO2", "磷酸酯")
    vec.append(1.0 if any(k in fg for fg in fgroups for k in hetero_kw) else 0.0)

    # 25: 需顶空 (低沸点/水溶液 → 顶空)
    vec.append(1.0 if bp < 150 else 0.0)

    # 26: π共轭体系大小 (0=无, 0.25=单环, 0.5=双环, 0.75=三环, 1.0=四环+)
    aromatic = "aromatic 芳香环" in fgroups
    c_count = substance.get("carbon_count", 1)
    if aromatic:
        if c_count >= 20:
            vec.append(1.0)   # 四环+ (芘、苯并芘等)
        elif c_count >= 14:
            vec.append(0.75)  # 三环 (蒽、菲)
        elif c_count >= 10:
            vec.append(0.5)   # 双环 (萘)
        else:
            vec.append(0.25)  # 单环 (苯、甲苯)
    else:
        vec.append(0.0)

    # 27: 长链酯 (FAMEs特征)
    is_ester = any(fg in fgroups for fg in ("-COO- 酯基", "-C=O 酮/醛羰基"))
    is_long = c_count >= 12
    vec.append(1.0 if (is_ester and is_long) else 0.0)

    # 28: 检测器路由 (0=FID通用, 1=P→FPD/NPD, 2=S→SCD/FPD, 3=卤素→ECD, 4=N→NPD, 5=MS高分辨)
    has_p = "磷酸酯/有机磷 (P=O/S)" in fgroups
    has_s = any(fg in fgroups for fg in ("硫醇 -SH (H₂S/硫醇)", "硫醚 -S-"))
    has_halogen = any(fg in fgroups for fg in ("-Cl 有机氯", "-Br 有机溴", "-F 有机氟"))
    has_n = any(fg in fgroups for fg in ("-NH2 伯/仲胺", "含氮杂环 (吡啶/嘌呤等)", "-NO2 硝基"))
    if has_p: vec.append(1.0)
    elif has_s: vec.append(2.0/5)
    elif has_halogen: vec.append(3.0/5)
    elif has_n: vec.append(4.0/5)
    elif bp > 300 and substance.get("detection_limit","ppm") == "ppt": vec.append(5.0/5)
    else: vec.append(0.0)

    # 29: 有机金属 (Sn/Hg/Pb/As)
    vec.append(0.0)  # 查询端暂不设

    # 30: 含氮杂环药物
    vec.append(1.0 if "含氮杂环 (吡啶/嘌呤等)" in fgroups else 0.0)

    # 31: 邻苯二甲酸酯特异性 (二酯+芳香+高沸点+不饱和)
    has_aromatic = "aromatic 芳香环" in fgroups
    vec.append(1.0 if (is_ester and has_aromatic and bp > 250) else 0.0)

    # 32: 高苯基柱亲和 (大π体系需要高苯基%)
    vec.append(1.0 if (has_aromatic and c_count >= 14) else 0.0)

    return vec
def detector_to_sensitivity(det_type: str) -> float:
    d = {"TCD": 0.0, "FID": 0.4, "NPD": 0.6, "MS": 0.8, "ECD": 1.0, "SCD": 0.9}
    return d.get(det_type, 0.4)


def encode_method(method: dict) -> list[float]:
    """标准方法 → 特征向量 (26维, 与 encode_substance 对应)"""
    col = method["columns"][0]
    det = method["detectors"][0]
    prep = method["preparation"]["technique"]
    app = method["applicability"]
    cat = app["analyte_category"]
    inlet = method.get("inlet", {})

    vec = []

    # 0-2: 物态 — 从基质推导
    matrix = app.get("matrix", "")
    if any(k in matrix for k in ("气体", "空气", "氦气")):
        vec.extend([1.0, 0.0, 0.0])
    elif any(k in matrix for k in ("水", "血液")):
        vec.extend([0.0, 1.0, 0.0])
    elif any(k in matrix for k in ("固体", "材料", "粉")):
        vec.extend([0.0, 0.0, 1.0])
    else:
        vec.extend([0.0, 1.0, 0.0])

    # 3: 极性 — 柱固定相推导
    vec.append(phase_to_polarity(col["brand"], col["phase"]))

    # 4: 最高柱温归一化
    oven = col.get("oven_program", [])
    temps = []
    for s in oven:
        temps.append(s.get("temp_c", 50) or 50)
        if s["type"] == "ramp":
            temps.append(s.get("target_c", 50) or 50)
    max_t = max(temps) if temps else 200
    vec.append(min(max((max_t - 30) / 370.0, 0), 1.0))

    # 5: LogP 代理
    vec.append(0.5 - vec[3] * 0.4)

    # 6: 碳数代理
    vec.append(max_t / 400.0)

    # 7-23: 官能团 one-hot 17
    analytes_str = " ".join(app.get("target_analytes", []))
    al = analytes_str.lower()
    fg_vec = [0.0] * 17

    fg_rules = [
        # 0: -OH 脂肪醇
        (0, ("methanol", "ethanol", "propanol", "butanol", "isopropanol",
             "glycol", "diol", "alcohol", "hydroxy", "sterol", "cholesterol",
             "甲醇", "乙醇", "丙醇", "丁醇", "异丙醇", "乙二醇", "丙二醇", "二甘醇",
             "醇"), []),
        # 1: 酚羟基 Ar-OH
        (1, ("phenol", "cresol", "phenolic", "bisphenol", "xylenol",
             "苯酚", "甲酚", "二甲苯酚", "酚类"), []),
        # 2: -COOH 羧酸
        (2, ("acetic acid", "有机酸", "acid", "羧酸", "甲酸",
             "乙酸", "酸性"), []),
        # 3: -COO- 酯基
        (3, ("ester", "phthalate", "FAME", "acrylate",
             "methyl ", "ethyl ",
             "malathion", "parathion",
             "邻苯二甲酸酯", "甲酯", "乙酯", "邻苯二甲酸"), []),
        # 4: -C=O 酮/醛羰基
        (4, ("ketone", "acetone", "aldehyde", "furfural", "hexanal",
             "caffeine", "xanthine",
             "丙酮", "酮", "醛", "糠醛", "己醛", "咖啡因"), []),
        # 5: -NH2 伯/仲胺 (非杂环)
        (5, ("amine", "麻醉", "aniline", "amino",
             "benzocaine", "prilocaine", "lidocaine", "procaine",
             "tetracaine", "bupivacaine",
             "胺", "苯佐卡因", "利多卡因", "普鲁卡因", "丁卡因",
             "azine", "三嗪", "triazine"), []),
        # 6: 含氮杂环
        (6, ("caffeine", "pyridine", "imidazole", "purine", "alkaloid",
             "xanthine", "nicotine",
             "咖啡因", "嘌呤", "吡啶", "杂环",
             "azine", "三嗪", "triazine", "triazole", "azole"), []),
        # 7: -NO2 硝基
        (7, ("nitro", "nitrosamine", "nitroso", "nitroaromatic",
             "硝基", "亚硝胺"), []),
        # 8: 磷酸酯/有机磷
        (8, ("organophosphate", "phosphorus", "phosphor",
             "malathion", "parathion", "dimethoate", "diazinon",
             "马拉硫磷", "乐果", "二嗪农", "磷", "磷酸", "有机磷"), []),
        # 9: -Cl 有机氯
        (9, ("chloro", "OCP", "PCB", "DDT", "BHC", "CFC",
             "dichloro", "trichloro", "tetrachloro",
             "chlordane", "endrin", "dieldrin", "heptachlor",
             "DDE", "DDD", "chlorobenzene",
             "有机氯", "氯", "六六六", "滴滴涕", "七氯", "艾氏剂",
             "狄氏剂", "异狄氏剂", "氯丹", "多氯联苯",
             "atrazine", "simazine", "莠去津", "西玛津"), []),
        # 10: -Br 有机溴
        (10, ("bromo", "BFR", "PBDE", "brominated", "bromine",
              "溴", "多溴"), []),
        # 11: -F 有机氟
        (11, ("fluoro", "fluoranthene", "CFC", "fluorocarbon",
              "氟", "氟碳"), []),
        # 12: 硫醇 -SH
        (12, ("thiol", "mercaptan", "hydrogen sulfide", "H2S",
              "硫醇", "硫化氢"), []),
        # 13: 硫醚 -S-
        (13, ("sulfide", "disulfide", "disulphide", "thiobis", "thioether",
              "sulfone", "sulfate", "sulfur", "thiophene", "sulphur",
              "硫醚", "硫", "二硫", "硫丹", "扑草净"), []),
        # 14: -O- 醚
        (14, ("ether", "dioxane", "epoxy", "ethoxy", "methoxy", "furan",
              "醚", "二氧六环", "环氧", "四氢呋喃", "呋喃"), []),
        # 15: aromatic 芳香环
        (15, ("benzene", "phenyl", "PAH", "phenol", "anthracene", "pyrene",
              "naphthalene", "fluorene", "chrysene", "PCB", "phthalate",
              "xylene", "toluene", "styrene", "benzyl", "cresol", "biphenyl",
              "caffeine", "benzo", "perylene", "acenaph", "fluoranthene",
              "苯", "甲苯", "二甲苯", "乙苯", "苯乙烯", "萘", "蒽", "菲",
              "芘", "荧蒽", "屈", "苊", "芴", "甾", "胆甾",
              "苯并", "茚并", "二苯并", "邻苯", "多环芳烃", "芳烃",
              "苯酚", "甲酚", "DDT", "BHC",
              "azine", "triazine", "azole"), []),
        # 16: alkene 烯烃
        (16, ("alkene", "oleate", "linole", "linolen", "butadiene",
              "pentene", "hexene", "octene", "FAME", "fatty acid",
              "oil", "TPH", "烯", "油酸", "亚油酸", "亚麻酸",
              "丁二烯", "碳数分布", "sterol", "cholesterol", "甾"), []),
    ]
    for idx, en_keywords, cn_keywords in fg_rules:
        # 同时检查分析物名称和类别 (英文+中文)
        text = (al + " " + cat).lower()
        all_kw = tuple(en_keywords) + tuple(cn_keywords)
        if any(k.lower() in text for k in all_kw):
            fg_vec[idx] = 1.0
    vec.extend(fg_vec)

    # 19-20: 氢键标记 (基于新官能团)
    hbd_match = any(fg_vec[i] for i in (0, 1, 2, 5, 6, 12))  # OH,酚,COOH,NH2,杂环,SH
    hba_match = any(fg_vec[i] for i in (0, 1, 2, 3, 4, 14))  # OH,酚,COOH,酯,酮,醚
    vec.append(1.0 if hbd_match else 0.0)
    vec.append(1.0 if hba_match else 0.0)

    # 21: 检出限 — 检测器推导
    vec.append(detector_to_sensitivity(det["type"]))

    # 22: 初始柱温 (低=VOC方法, 高=SVOC方法)
    min_t = min(temps) if temps else 35
    vec.append(1.0 - min(max((min_t - 30) / 200.0, 0), 1.0))

    # 23: 是否不分流
    mode = inlet.get("mode", "")
    vec.append(1.0 if "不分流" in str(mode) or "splitless" in str(mode).lower() else 0.0)

    # 24: 是否含杂原子方法
    hetero = ("Cl", "Br", "F", "S", "N", "P", "氯", "溴", "氟", "硫", "氮", "磷",
              "nitro", "amino")
    vec.append(1.0 if any(k in cat for k in hetero) else 0.0)

    # 25: 是否需要顶空前处理
    vec.append(1.0 if any(k in prep for k in ("顶空", "吹扫", "HS", "Pur&Trap")) else 0.0)

    # 26: π共轭体系大小 (从类别+分析物推导)
    pah_kw = ("PAH", "多环芳烃", "pyrene", "naphthalene", "anthracene", "fluorene",
              "chrysene", "phenanthrene", "benzo", "perylene", "fluoranthene", "acenaph",
              "芘", "萘", "蒽", "菲", "荧蒽", "屈", "苯并", "苊", "芴")
    pcb_kw = ("PCB", "多氯联苯", "biphenyl", "联苯")
    if any(k in cat for k in pah_kw):
        # 计数实际PAH分析物中的最大环数
        astr = " ".join(app.get("target_analytes", []))
        if any(k in astr for k in ("benzo[", "indeno", "dibenz", "perylene",
                                    "苯并", "茚并", "二苯并")):
            vec.append(1.0)  # 4-5环
        elif any(k in astr for k in ("pyrene", "anthracene", "phenanthrene",
                                      "chrysene", "fluoranthene",
                                      "芘", "蒽", "菲", "荧蒽", "屈")):
            vec.append(0.75)  # 3-4环
        elif any(k in astr for k in ("naphthalene", "acenaph", "fluorene",
                                      "萘", "苊", "芴")):
            vec.append(0.5)   # 2-3环
        else:
            vec.append(0.6)
    elif any(k in cat for k in pcb_kw):
        vec.append(0.4)  # PCBs 有联苯结构
    elif any(k in cat for k in ("芳", "苯", "酚", "邻苯", "phthalate", "phenyl", "xylene", "toluene", "cresol", "甾", "麻醉")):
        vec.append(0.2)  # 单环芳烃
    else:
        vec.append(0.0)  # 非芳香

    # 27: FAMEs/长链酯方法
    fame_kw = ("FAME", "fatty", "甲酯", "脂肪酸")
    if any(k in cat for k in fame_kw):
        vec.append(1.0)
    elif any(k in cat for k in ("phthalate", "邻苯")):
        vec.append(0.6)
    else:
        vec.append(0.0)

    # 28: 检测器路由 (与方法detector对应)
    dt = det["type"]
    if dt in ("FPD", "PFPD"):
        vec.append(1.0)  # P/S → FPD
    elif dt == "SCD":
        vec.append(2.0/5)  # S → SCD
    elif dt == "ECD":
        vec.append(3.0/5)  # 卤素 → ECD
    elif dt == "NPD":
        vec.append(4.0/5)  # N/P → NPD
    elif dt == "MS":
        vec.append(5.0/5)  # 高分辨 → MS
    else:
        vec.append(0.0)  # FID → 通用

    # 29: 有机金属
    organometal_kw = ("有机锡", "organotin", "有机汞", "organomercury")
    vec.append(1.0 if any(k in cat for k in organometal_kw) else 0.0)

    # 30: 含氮杂环药物
    drug_kw = ("药物", "麻醉", "drug", "pharmaceutical", "caffeine", "咖啡因", "alkaloid")
    vec.append(1.0 if any(k in cat for k in drug_kw) or fg_vec[5] or fg_vec[6] else 0.0)

    # 31: 邻苯二甲酸酯特异性
    phthalate_kw = ("phthalate", "邻苯二甲酸酯", "邻苯")
    vec.append(1.0 if any(k in cat for k in phthalate_kw) else 0.0)

    # 32: 高苯基柱亲和 (50%苯基或35%苯基)
    high_phenyl = phase_to_polarity(col["brand"], col["phase"]) >= 0.50
    vec.append(1.0 if high_phenyl else 0.0)

    return vec


def phase_to_polarity(brand: str, phase: str) -> float:
    """柱固定相 → 极性分数"""
    s = (brand + " " + phase).lower()
    if "2560" in s or "双氰丙基" in s:
        return 0.95
    if "wax" in s and "da" in s:
        return 0.85
    if "wax" in s:
        return 0.78
    if "50%" in s or "17sil" in s:
        return 0.60
    if "17" in s and "1701" not in s:
        return 0.60
    if "35%" in s or "35sil" in s:
        return 0.50
    if "6%" in s or "624" in s or "1301" in s:
        return 0.35
    if "14%" in s or "1701" in s:
        return 0.25
    if "5%" in s or "5ms" in s or "rtx-5" in s:
        return 0.12 if "sil" in s else 0.10
    if "xlb" in s:
        return 0.55
    if "100%" in s or "二甲基" in s or "1ms" in s:
        return 0.05
    if "分子筛" in s or "msieve" in s or "alumina" in s:
        return 0.0
    if "bac" in s:
        return 0.35
    logger.warning("未识别的柱固定相: brand=%r phase=%r → fallback=0.3，可能影响极性匹配精度", brand, phase)
    return 0.3


# ════════════════════════════════════════════════════════
# 柱固定相数据库 — 替代推荐
# ════════════════════════════════════════════════════════
# phase_key → 特征
_PHASE_DB = {
    # name, polarity, backbone, phenyl%, cyano%, temp_limit, 适用场景
    "dm100":  ("100%二甲基聚硅氧烷", 0.05, "pdms", 0, 0, 350, "非极性通用，烷烃/石化"),
    "ph5":    ("5%苯基-95%二甲基", 0.10, "pdms", 5, 0, 350, "半挥发物/芳烃/通用"),
    "ph5ms":  ("5%苯基 MS级(低流失)", 0.12, "pdms", 5, 0, 350, "MS检测器首选"),
    "ph35":   ("35%苯基-65%二甲基", 0.50, "pdms", 35, 0, 340, "PAHs/PCBs/农药"),
    "ph50":   ("50%苯基-50%二甲基", 0.60, "pdms", 50, 0, 340, "PAHs/药物/固醇类"),
    "cn6":    ("6%氰丙基苯基-94%二甲基", 0.35, "pdms", 6, 6, 280, "VOCs/溶剂/挥发性有机物"),
    "cn14":   ("14%氰丙基苯基-86%二甲基", 0.25, "pdms", 14, 14, 280, "农药/除草剂/Aroclor"),
    "xlb":    ("专有低流失中极性", 0.55, "pdms", 20, 0, 360, "PCBs/OCPs/溴代阻燃剂"),
    "peg":    ("聚乙二醇 PEG/WAX", 0.78, "peg", 0, 0, 260, "极性样品/醇/酯/香料"),
    "ffap":   ("酸改性聚乙二醇 FFAP", 0.85, "peg", 0, 0, 260, "有机酸/酚类/极性酸"),
    "cn90":   ("双氰丙基 高极性", 0.95, "pdms", 0, 90, 250, "FAMEs/顺反异构/二噁英"),
    "alumina":("氧化铝 PLOT", 0.0, "plot", 0, 0, 200, "C1-C5烃/炼厂气/气固"),
    "msieve": ("分子筛 5A PLOT", 0.0, "plot", 0, 0, 350, "永久气体 O₂/N₂/CO/CH₄"),
    "bac":    ("专有血液酒精固定相", 0.35, "proprietary", 0, 0, 260, "血液酒精/挥发性醇"),
}

# 品牌/型号 → phase_key 映射
_BRAND_MAP = {
    # 100% dimethyl
    "db-1": "dm100", "hp-1": "dm100", "rtx-1": "dm100",
    "rxi-1ms": "dm100", "zb-1": "dm100", "spb-1": "dm100",
    # 5% phenyl
    "db-5": "ph5", "hp-5": "ph5", "rtx-5": "ph5",
    "zb-5": "ph5", "spb-5": "ph5", "cp-sil 8": "ph5",
    # 5% phenyl MS
    "db-5ms": "ph5ms", "hp-5ms": "ph5ms", "rtx-5ms": "ph5ms",
    "rxi-5ms": "ph5ms", "rxi-svocms": "ph5ms", "zb-5ms": "ph5ms",
    "vf-5ms": "ph5ms", "rmx-5sil ms": "ph5ms",
    # 35% phenyl
    "db-35": "ph35", "hp-35": "ph35", "rtx-35": "ph35",
    "rxi-35sil ms": "ph35", "zb-35": "ph35", "spb-35": "ph35",
    # 50% phenyl
    "db-17": "ph50", "hp-50+": "ph50", "rtx-50": "ph50",
    "rxi-17": "ph50", "rxi-17sil ms": "ph50", "zb-50": "ph50",
    # 6% cyanopropyl
    "db-624": "cn6", "hp-624": "cn6", "rtx-624": "cn6",
    "rxi-624sil ms": "cn6", "zb-624": "cn6", "rtx-1301": "cn6",
    "vf-624ms": "cn6",
    # 14% cyanopropyl
    "db-1701": "cn14", "hp-1701": "cn14", "rtx-1701": "cn14",
    "zb-1701": "cn14", "cp-sil 19": "cn14",
    # WAX/PEG
    "db-wax": "peg", "hp-wax": "peg", "rtx-wax": "peg",
    "stabilwax": "peg", "zb-wax": "peg", "carbowax": "peg",
    "cp-wax": "peg", "supelcowax": "peg",
    # FFAP
    "db-ffap": "ffap", "hp-ffap": "ffap", "stabilwax-da": "ffap",
    "zb-ffap": "ffap", "cp-wax 58": "ffap", "nukol": "ffap",
    # cyano
    "rt-2560": "cn90", "sp-2560": "cn90",
    # PLOT
    "rt-alumina": "alumina", "hp-plot al2o3": "alumina",
    "gs-alumina": "alumina",
    "rt-msieve 5a": "msieve", "hp-molesieve": "msieve",
    "cp-molsieve": "msieve",
    # proprietary
    "rtx-bac1": "bac", "rtx-bac2": "bac",
    "rxi-xlb": "xlb",
    # fill column
    "res-sil c": "msieve",  # 气固色谱填充柱 → 归入PLOT类
}


def _classify_column(brand: str, phase: str = "") -> str:
    """品牌+固定相 → phase_key"""
    s = (brand + " " + phase).lower()
    # 精确匹配品牌名
    for key, phase_key in sorted(_BRAND_MAP.items(), key=lambda x: -len(x[0])):
        if key in s:
            return phase_key
    # 从固定相描述推断
    if "wax" in s and "da" in s:
        return "ffap"
    if "wax" in s or "peg" in s:
        return "peg"
    if "2560" in s or "双氰丙基" in s:
        return "cn90"
    if "alumina" in s or "氧化铝" in s:
        return "alumina"
    if "分子筛" in s or "msieve" in s or "molsieve" in s:
        return "msieve"
    if "50%" in s or "17sil" in s:
        return "ph50"
    if "17" in s and "1701" not in s:
        return "ph50"
    if "35%" in s or "35sil" in s:
        return "ph35"
    if "14%" in s or "1701" in s:
        return "cn14"
    if "6%" in s or "624" in s or "1301" in s:
        return "cn6"
    if "xlb" in s:
        return "xlb"
    if "100%" in s or "二甲基" in s or "1ms" in s:
        return "dm100"
    if "5%" in s or "5ms" in s or "rtx-5" in s:
        return "ph5ms" if "sil" in s else "ph5"
    if "bac" in s:
        return "bac"
    logger.warning("未识别的柱固定相: brand=%r phase=%r → fallback=ph5 用于替代匹配", brand, phase)
    return "ph5"


def column_substitute(recommended_brand: str, recommended_phase: str,
                      available_columns: list[str]) -> list[tuple[str, str, str]]:
    """
    推荐柱 vs 用户手头柱子 → 最佳替代建议

    available_columns: 用户手头有的品牌名列表, 如 ["DB-5","DB-1701","DB-WAX"]

    Returns: [(品牌, 匹配等级, 说明), ...]
      等级: "等效" | "接近" | "勉强" | "不推荐"
    """
    rec_key = _classify_column(recommended_brand, recommended_phase)
    rec = _PHASE_DB.get(rec_key)
    if not rec:
        return []

    # 分类每个手头柱子
    available = []
    for col_name in available_columns:
        avail_key = _classify_column(col_name)
        avail_info = _PHASE_DB.get(avail_key)
        if avail_info:
            available.append((col_name, avail_key, avail_info))

    if not available:
        return []

    results = []
    for col_name, avail_key, avail_info in available:
        name_r, pol_r, bb_r, ph_r, cn_r, tmp_r, use_r = rec
        name_a, pol_a, bb_a, ph_a, cn_a, tmp_a, use_a = avail_info

        if rec_key == avail_key:
            results.append((col_name, "等效", f"与{recommended_brand}固定相相同，可直接用"))
            continue

        if bb_r != bb_a:
            if bb_r == "plot" or bb_a == "plot":
                results.append((col_name, "不推荐",
                    f"{name_a}是气固色谱柱，不能替代气液色谱柱{recommended_brand}"))
            elif bb_r == "peg" and bb_a == "pdms":
                results.append((col_name, "不推荐",
                    f"{name_a}(聚硅氧烷)与{recommended_brand}(PEG)化学性质差异大"))
            elif bb_r == "pdms" and bb_a == "peg":
                results.append((col_name, "勉强",
                    f"{name_a}(PEG)可做极性更高的替代，但选择性完全不同，保留时间和峰序会变"))
            else:
                results.append((col_name, "勉强",
                    f"骨架类型不同({bb_r}↔{bb_a})，出峰顺序可能不同"))
            continue

        # 同骨架，按极性差 + 苯基/氰基差判断
        pol_diff = abs(pol_r - pol_a)
        phenyl_diff = abs(ph_r - ph_a)
        cyano_diff = abs(cn_r - cn_a)

        if pol_diff < 0.06 and phenyl_diff <= 2 and cyano_diff == 0:
            results.append((col_name, "等效",
                f"{name_a}与{recommended_brand}极性极接近，可直接替代"))
        elif pol_diff < 0.20:
            if phenyl_diff >= 3:
                results.append((col_name, "接近",
                    f"{name_a}苯基含量({ph_a}%)与推荐({ph_r}%)不同，保留时间可能偏移，可调温度补偿"))
            elif cyano_diff >= 3:
                results.append((col_name, "接近",
                    f"{name_a}含氰丙基({cn_a}%)与推荐不同，选择性有差异，建议标样验证"))
            else:
                results.append((col_name, "接近",
                    f"{name_a}极性({pol_a:.2f})与推荐({pol_r:.2f})有差异，建议先做标样确认"))
        elif pol_diff < 0.40:
            results.append((col_name, "勉强",
                f"极性差{pol_diff:.2f}较大，建议调温度程序或考虑分两次跑"))
        else:
            results.append((col_name, "不推荐",
                f"极性差{pol_diff:.2f}过大，替代后分离度无法保证"))

    # 按匹配等级排序: 等效 > 接近 > 勉强 > 不推荐
    rank = {"等效": 0, "接近": 1, "勉强": 2, "不推荐": 3}
    results.sort(key=lambda x: rank.get(x[1], 4))
    return results


# ── 权重 — 极性/沸点主导 ──
FEATURE_WEIGHTS = (
    [0.01] * 3 +    # 物态 (合计0.03)
    [0.25] +        # 极性 ★ 柱选择
    [0.20] +        # 沸点 ★ 温度选择
    [0.03] +        # LogP
    [0.01] +        # 碳数
    [0.007] * 17 +  # 官能团 (合计0.119)
    [0.01] +        # HBD
    [0.01] +        # HBA
    [0.04] +        # 检出限
    [0.04] +        # 初温
    [0.03] +        # 不分流
    [0.04] +        # 杂原子
    [0.04] +        # 顶空
    [0.08] +        # π共轭体系 ★
    [0.09] +        # 长链酯 ★
    [0.06] +        # 检测器路由 ★
    [0.03] +        # 有机金属
    [0.05] +        # N-杂环药物 ★
    [0.03] +        # 邻苯特异性
    [0.04]          # 高苯基柱亲和 ★
)


def weighted_cosine(a: list[float], b: list[float]) -> float:
    """加权余弦相似度"""
    dot = sum(w * x * y for w, x, y in zip(FEATURE_WEIGHTS, a, b))
    na = math.sqrt(sum(w * x * x for w, x in zip(FEATURE_WEIGHTS, a)))
    nb = math.sqrt(sum(w * y * y for w, y in zip(FEATURE_WEIGHTS, b)))
    return dot / (na * nb) if na and nb else 0


def cosine_similarity(a: list[float], b: list[float]) -> float:
    dot = sum(x * y for x, y in zip(a, b))
    na = math.sqrt(sum(x * x for x in a))
    nb = math.sqrt(sum(y * y for y in b))
    return dot / (na * nb) if na and nb else 0

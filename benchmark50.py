#!/usr/bin/env python3
"""50种物质评测 — 统计 + 逐条输出"""
import os
from chromethod.knn import load_methods, build_index, search
from chromethod.hook import validate

_dir = os.path.dirname(os.path.abspath(__file__))
methods = load_methods(os.path.join(_dir, "chromethod", "data"))
vectors = build_index(methods)

# (name, fg_list, polarity, bp_c, carbon, phase, logp, hbd, hba, detlim, expected_category)
BENCH = [
    # === VOCs/残留溶剂 (10) ===
    ("苯",    ["aromatic 芳香环"],1,80,6,"液体",2.13,0,0,"ppb","VOCs/残留溶剂/芳烃"),
    ("甲苯",  ["aromatic 芳香环"],1,111,7,"液体",2.73,0,0,"ppb","VOCs/残留溶剂/芳烃"),
    ("乙苯",  ["aromatic 芳香环"],1,136,8,"液体",3.15,0,0,"ppb","VOCs/芳烃"),
    ("邻二甲苯",["aromatic 芳香环"],1,144,8,"液体",3.12,0,0,"ppb","VOCs/芳烃"),
    ("苯乙烯",["aromatic 芳香环","alkene 烯烃"],1,145,8,"液体",2.95,0,0,"ppb","VOCs/芳烃"),
    ("二氯甲烷",["-Cl 有机氯"],2,40,1,"液体",1.25,0,0,"ppm","残留溶剂/卤代"),
    ("三氯甲烷",["-Cl 有机氯"],2,61,1,"液体",1.97,0,0,"ppb","残留溶剂/卤代"),
    ("四氯化碳",["-Cl 有机氯"],1,77,1,"液体",2.83,0,0,"ppb","卤代/OCI"),
    ("丙酮",  ["-C=O 酮/醛羰基"],3,56,3,"液体",-0.16,0,1,"ppm","VOCs/酮"),
    ("正己烷",[],1,69,6,"液体",3.9,0,0,"ppm","VOCs/烷烃"),

    # === 醇类 (3) ===
    ("甲醇",  ["-OH 脂肪醇羟基"],5,65,1,"液体",-0.77,1,1,"ppm","醇/残留溶剂"),
    ("乙醇",  ["-OH 脂肪醇羟基"],4,78,2,"液体",-0.31,1,1,"ppm","醇/残留溶剂"),
    ("异丙醇",["-OH 脂肪醇羟基"],3,82,3,"液体",0.05,1,1,"ppm","醇/VOCs"),

    # === 有机酸/酚 (5) ===
    ("乙酸",  ["-COOH 羧酸"],5,118,2,"液体",-0.17,1,2,"ppm","有机酸/FFAP"),
    ("苯酚",  ["酚羟基 (Ar-OH)","aromatic 芳香环"],3,182,6,"固体",1.46,1,1,"ppm","酚类"),
    ("邻甲酚",["酚羟基 (Ar-OH)","aromatic 芳香环"],2,191,7,"固体",1.95,1,1,"ppm","酚类"),
    ("2,4-二甲基苯酚",["酚羟基 (Ar-OH)","aromatic 芳香环"],2,211,8,"固体",2.3,1,1,"ppm","酚类"),
    ("五氯苯酚",["酚羟基 (Ar-OH)","aromatic 芳香环","-Cl 有机氯"],3,310,6,"固体",5.0,1,1,"ppb","酚类/OCPs"),

    # === PAHs (6) ===
    ("萘",    ["aromatic 芳香环"],1,218,10,"固体",3.3,0,0,"ppb","PAHs"),
    ("联苯",  ["aromatic 芳香环"],1,255,12,"固体",4.0,0,0,"ppb","PAHs/PCB"),
    ("苊",    ["aromatic 芳香环"],1,279,12,"固体",3.92,0,0,"ppb","PAHs"),
    ("菲",    ["aromatic 芳香环"],1,340,14,"固体",4.46,0,0,"ppb","PAHs"),
    ("芘",    ["aromatic 芳香环"],1,393,16,"固体",4.88,0,0,"ppt","PAHs"),
    ("苯并[a]芘",["aromatic 芳香环"],1,495,20,"固体",6.13,0,0,"ppt","PAHs"),

    # === 有机氯农药 (5) ===
    ("γ-BHC(林丹)",["-Cl 有机氯"],3,323,6,"固体",3.72,0,0,"ppt","OCPs/ECD"),
    ("4,4'-DDT",["-Cl 有机氯","aromatic 芳香环"],2,260,14,"固体",6.91,0,0,"ppt","OCPs/ECD"),
    ("4,4'-DDE",["-Cl 有机氯","aromatic 芳香环"],2,316,14,"固体",6.0,0,0,"ppt","OCPs/ECD"),
    ("狄氏剂",["-Cl 有机氯","-O- 醚"],3,385,12,"固体",5.2,0,1,"ppt","OCPs/ECD"),
    ("硫丹I",["-Cl 有机氯","硫醚 -S-","-O- 醚"],3,350,9,"固体",3.6,0,1,"ppt","OCPs/ECD"),

    # === 有机磷农药 (4) ===
    ("马拉硫磷",["-COO- 酯基","硫醚 -S-","磷酸酯/有机磷 (P=O/S)"],3,156,10,"液体",2.36,0,6,"ppb","有机磷/NPD/FPD"),
    ("对硫磷",["-NO2 硝基","磷酸酯/有机磷 (P=O/S)","硫醚 -S-"],3,375,10,"液体",3.83,0,4,"ppb","有机磷/NPD/FPD"),
    ("乐果",  ["磷酸酯/有机磷 (P=O/S)","-C=O 酮/醛羰基","硫醚 -S-"],3,117,5,"固体",0.78,0,3,"ppb","有机磷/NPD"),
    ("毒死蜱",["磷酸酯/有机磷 (P=O/S)","-Cl 有机氯","硫醚 -S-","aromatic 芳香环"],3,350,9,"固体",4.96,0,3,"ppb","有机磷/含氯"),

    # === 含氮农药/杂环 (3) ===
    ("莠去津",["-Cl 有机氯","含氮杂环 (吡啶/嘌呤等)","-NH2 伯/仲胺"],3,200,8,"固体",2.61,2,5,"ppb","三嗪/NPD/MS"),
    ("西玛津",["-Cl 有机氯","含氮杂环 (吡啶/嘌呤等)","-NH2 伯/仲胺"],3,225,7,"固体",2.18,2,5,"ppb","三嗪/NPD/MS"),
    ("咖啡因",["-C=O 酮/醛羰基","含氮杂环 (吡啶/嘌呤等)","alkene 烯烃"],3,178,8,"固体",-0.07,0,3,"ppb","药物/含氮"),

    # === 硫化物 (3) ===
    ("硫化氢",["硫醇 -SH (H₂S/硫醇)"],1,-60,0,"气体",0.23,0,0,"ppm","硫化物/SCD/FPD"),
    ("甲硫醇",["硫醇 -SH (H₂S/硫醇)"],1,6,1,"气体",0.65,0,0,"ppb","硫化物/SCD"),
    ("二甲基二硫醚",["硫醚 -S-"],2,109,2,"液体",1.77,0,0,"ppb","硫化物/SCD/FPD"),

    # === FAMEs/脂肪酸(4) ===
    ("油酸甲酯",["-COO- 酯基","alkene 烯烃"],3,350,19,"液体",7.5,0,2,"ppm","FAMEs"),
    ("棕榈酸甲酯",["-COO- 酯基"],2,330,17,"液体",6.4,0,2,"ppm","FAMEs"),
    ("亚油酸甲酯",["-COO- 酯基","alkene 烯烃"],3,365,19,"液体",7.0,0,2,"ppm","FAMEs"),
    ("硬脂酸甲酯",["-COO- 酯基"],2,360,19,"液体",7.2,0,2,"ppm","FAMEs"),

    # === 邻苯二甲酸酯(3) ===
    ("DEHP",["-COO- 酯基","aromatic 芳香环"],2,385,24,"液体",7.6,0,4,"ppb","邻苯二甲酸酯"),
    ("DBP",["-COO- 酯基","aromatic 芳香环"],2,340,16,"液体",4.5,0,4,"ppb","邻苯二甲酸酯"),
    ("BBP",["-COO- 酯基","aromatic 芳香环"],2,370,19,"液体",4.73,0,4,"ppb","邻苯二甲酸酯"),

    # === 气体(4) ===
    ("甲烷",  [],1,-161,1,"气体",1.09,0,0,"ppm","永久气体/TCD"),
    ("二氧化碳",[],1,-78,1,"气体",-0.5,0,2,"ppm","永久气体/TCD"),
    ("氮气",  [],1,-196,0,"气体",0.0,0,0,"ppm","永久气体/TCD"),
    ("乙烯",  ["alkene 烯烃"],1,-104,2,"气体",1.13,0,0,"ppm","气体/烯烃"),

    # === 其他(3) ===
    ("胆固醇",["-OH 脂肪醇羟基","alkene 烯烃"],2,360,27,"固体",8.7,1,1,"ppm","甾醇/衍生化"),
    ("1,4-二氧六环",["-O- 醚"],3,101,4,"液体",-0.27,0,2,"ppb","残留溶剂/醚"),
    ("四氢呋喃",["-O- 醚"],3,66,4,"液体",0.46,0,1,"ppb","残留溶剂/醚"),
]

# ── 评测 ──
stats = {"knn_correct":0, "hook_warn":0, "hook_reject":0, "total":len(BENCH)}

for name, fgs, pol, bp, c, ph, logp, hbd, hba, detlim, exp in BENCH:
    q = {"phase":ph,"polarity":pol,"bp_c":bp,"logp":logp,"carbon_count":c,
         "functional_groups":fgs,"hbd":hbd,"hba":hba,"detection_limit":detlim}
    res = search(q, methods, vectors, k=1)
    if not res:
        print(f"✗ {name:12s} 无匹配")
        continue
    sim, method = res[0]
    col = method["columns"][0]; det = method["detectors"][0]
    cat = method["applicability"]["analyte_category"]
    src = method["source"]["standard"]
    mname = method["source"]["name"][:38]

    # k-NN 命中判定
    kmatch = any(k in cat+src+mname for k in exp.split("/"))
    if kmatch: stats["knn_correct"] += 1
    knn_icon = "✓" if kmatch else "?"

    # Hook 验证
    report = validate([q], method)
    h_icon = ""
    if report.rejected:
        h_icon = "R"  # REJECT
        stats["hook_reject"] += 1
    elif report.warns:
        h_icon = "W"  # WARN
        stats["hook_warn"] += 1

    print(f"{knn_icon} {h_icon} {name:12s} s={sim:.3f} {col['brand']:15s} {det['type']:4s} {mname}")

print(f"\n{'='*60}")
print(f"总计: {stats['total']} 种物质")
print(f"k-NN 命中: {stats['knn_correct']}/{stats['total']} ({100*stats['knn_correct']/stats['total']:.0f}%)")
print(f"Hook WARN: {stats['hook_warn']} | Hook REJECT: {stats['hook_reject']}")
print(f"k-NN 命中 + Hook 通过 = 推荐可用")
print(f"k-NN 未命中 + Hook 标记 = 触发兜底/LLM")

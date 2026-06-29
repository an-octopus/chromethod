#!/usr/bin/env python3
"""20种物质评测"""
import os
from chromethod.knn import load_methods, build_index, search

_dir = os.path.dirname(os.path.abspath(__file__))
methods = load_methods(os.path.join(_dir, "chromethod", "data"))
vectors = build_index(methods)
print(f"数据库: {len(methods)} 条方法\n")

queries = [
    ("苯", {"phase":"液体","polarity":1,"bp_c":80.1,"logp":2.13,"carbon_count":6,
            "functional_groups":["aromatic (芳香环)"],"hbd":0,"hba":0,"detection_limit":"ppb"}),
    ("甲苯", {"phase":"液体","polarity":1,"bp_c":110.6,"logp":2.73,"carbon_count":7,
             "functional_groups":["aromatic (芳香环)"],"hbd":0,"hba":0,"detection_limit":"ppb"}),
    ("二氯甲烷", {"phase":"液体","polarity":2,"bp_c":40,"logp":1.25,"carbon_count":1,
                 "functional_groups":["-Cl (有机氯)"],"hbd":0,"hba":0,"detection_limit":"ppm"}),
    ("丙酮", {"phase":"液体","polarity":3,"bp_c":56,"logp":-0.16,"carbon_count":3,
             "functional_groups":["-C=O (羰基/酮/醛)"],"hbd":0,"hba":1,"detection_limit":"ppm"}),
    ("正己烷", {"phase":"液体","polarity":1,"bp_c":69,"logp":3.9,"carbon_count":6,
              "functional_groups":[],"hbd":0,"hba":0,"detection_limit":"ppm"}),
    ("萘", {"phase":"固体","polarity":1,"bp_c":218,"logp":3.3,"carbon_count":10,
           "functional_groups":["aromatic (芳香环)"],"hbd":0,"hba":0,"detection_limit":"ppb"}),
    ("苯酚", {"phase":"固体","polarity":3,"bp_c":182,"logp":1.46,"carbon_count":6,
             "functional_groups":["-OH (醇/酚羟基)","aromatic (芳香环)"],"hbd":1,"hba":1,"detection_limit":"ppm"}),
    ("苯并[a]芘", {"phase":"固体","polarity":1,"bp_c":495,"logp":6.13,"carbon_count":20,
                  "functional_groups":["aromatic (芳香环)"],"hbd":0,"hba":0,"detection_limit":"ppt"}),
    ("邻苯二甲酸二(2-乙基己基)酯 DEHP", {"phase":"液体","polarity":2,"bp_c":385,"logp":7.6,
              "carbon_count":24,"functional_groups":["-C=O (羰基/酮/醛)","aromatic (芳香环)"],"hbd":0,"hba":4,"detection_limit":"ppb"}),
    ("4,4'-DDT", {"phase":"固体","polarity":2,"bp_c":260,"logp":6.91,"carbon_count":14,
                 "functional_groups":["-Cl (有机氯)","aromatic (芳香环)"],"hbd":0,"hba":0,"detection_limit":"ppt"}),
    ("林丹 γ-BHC", {"phase":"固体","polarity":3,"bp_c":323,"logp":3.72,"carbon_count":6,
                   "functional_groups":["-Cl (有机氯)"],"hbd":0,"hba":0,"detection_limit":"ppt"}),
    ("马拉硫磷", {"phase":"液体","polarity":3,"bp_c":156,"logp":2.36,"carbon_count":10,
                 "functional_groups":["-C=O (羰基/酮/醛)","-S- (硫醚/硫醇)"],"hbd":0,"hba":6,"detection_limit":"ppb"}),
    ("甲醇", {"phase":"液体","polarity":5,"bp_c":65,"logp":-0.77,"carbon_count":1,
             "functional_groups":["-OH (醇/酚羟基)"],"hbd":1,"hba":1,"detection_limit":"ppm"}),
    ("乙醇", {"phase":"液体","polarity":4,"bp_c":78,"logp":-0.31,"carbon_count":2,
             "functional_groups":["-OH (醇/酚羟基)"],"hbd":1,"hba":1,"detection_limit":"ppm"}),
    ("乙酸", {"phase":"液体","polarity":5,"bp_c":118,"logp":-0.17,"carbon_count":2,
             "functional_groups":["-COOH (羧酸)"],"hbd":1,"hba":2,"detection_limit":"ppm"}),
    ("甲烷", {"phase":"气体","polarity":1,"bp_c":-161,"logp":1.09,"carbon_count":1,
             "functional_groups":[],"hbd":0,"hba":0,"detection_limit":"ppm"}),
    ("硫化氢", {"phase":"气体","polarity":1,"bp_c":-60,"logp":0.23,"carbon_count":0,
               "functional_groups":["-S- (硫醚/硫醇)"],"hbd":0,"hba":0,"detection_limit":"ppm"}),
    ("油酸甲酯", {"phase":"液体","polarity":3,"bp_c":350,"logp":7.5,"carbon_count":19,
                 "functional_groups":["-C=O (羰基/酮/醛)","alkene (烯烃)"],"hbd":0,"hba":2,"detection_limit":"ppm"}),
    ("棕榈酸甲酯", {"phase":"液体","polarity":2,"bp_c":330,"logp":6.4,"carbon_count":17,
                   "functional_groups":["-C=O (羰基/酮/醛)"],"hbd":0,"hba":2,"detection_limit":"ppm"}),
    ("咖啡因", {"phase":"固体","polarity":3,"bp_c":178,"logp":-0.07,"carbon_count":8,
               "functional_groups":["-C=O (羰基/酮/醛)","-NH2 (伯胺)","alkene (烯烃)"],"hbd":0,"hba":3,"detection_limit":"ppb"}),
]

ok = 0
for name, q in queries:
    results = search(q, methods, vectors, k=1)
    if results:
        sim, m = results[0]
        col = m["columns"][0]
        det = m["detectors"][0]
        src = m["source"]["standard"]
        mname = m["source"]["name"][:35]
        # 人工标注期望
        expected = {
            "苯": "VOCs/残留溶剂",
            "甲苯": "VOCs/残留溶剂",
            "二氯甲烷": "残留溶剂/VOCs",
            "丙酮": "残留溶剂/VOCs",
            "正己烷": "残留溶剂/VOCs/TPH",
            "萘": "PAHs",
            "苯酚": "酚类/酸性",
            "苯并[a]芘": "PAHs",
            "DEHP": "邻苯二甲酸酯",
            "DDT": "OCPs/ECD",
            "林丹": "OCPs/ECD",
            "马拉硫磷": "有机磷/NPD",
            "甲醇": "残留溶剂/醇类",
            "乙醇": "残留溶剂/醇类",
            "乙酸": "有机酸/FFAP",
            "甲烷": "永久气体/TCD",
            "硫化氢": "硫化物/SCD",
            "油酸甲酯": "FAMEs",
            "棕榈酸甲酯": "FAMEs",
            "咖啡因": "药物/含氮",
        }
        exp = expected.get(name, "")
        cat = m["applicability"]["analyte_category"]
        hit = "✓" if any(k in cat+src for k in exp.split("/")) else "?"
        print(f"{hit} {name:12s} sim={sim:.3f} | {col['brand']:15s} {det['type']:3s} | {mname}")
        if "✓" in hit: ok += 1
    else:
        print(f"✗ {name:12s} 无匹配")

print(f"\n命中: {ok}/{len(queries)}")

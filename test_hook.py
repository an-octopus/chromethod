#!/usr/bin/env python3
"""Hook 验证测试"""
import os
from chromethod.knn import load_methods, build_index, search
from chromethod.hook import validate

_dir = os.path.dirname(os.path.abspath(__file__))
methods = load_methods(os.path.join(_dir, "chromethod", "data"))
vectors = build_index(methods)

# 问题案例: 物质 + k-NN Top-1
cases = [
    ("马拉硫磷", {"phase":"液体","polarity":3,"bp_c":156,"logp":2.36,
     "carbon_count":10,"functional_groups":["-COO- 酯基","硫醚 -S-","磷酸酯/有机磷 (P=O/S)"],
     "hbd":0,"hba":6,"detection_limit":"ppb"}),
    ("硫化氢", {"phase":"气体","polarity":1,"bp_c":-60,"logp":0.23,
     "carbon_count":0,"functional_groups":["硫醇 -SH (H₂S/硫醇)"],
     "hbd":0,"hba":0,"detection_limit":"ppm"}),
    ("DEHP", {"phase":"液体","polarity":2,"bp_c":385,"logp":7.6,
     "carbon_count":24,"functional_groups":["-COO- 酯基","aromatic 芳香环"],
     "hbd":0,"hba":4,"detection_limit":"ppb"}),
    ("咖啡因", {"phase":"固体","polarity":3,"bp_c":178,"logp":-0.07,
     "carbon_count":8,"functional_groups":["-C=O 酮/醛羰基","含氮杂环 (吡啶/嘌呤等)","alkene 烯烃"],
     "hbd":0,"hba":3,"detection_limit":"ppb"}),
    ("乙酸", {"phase":"液体","polarity":5,"bp_c":118,"logp":-0.17,
     "carbon_count":2,"functional_groups":["-COOH 羧酸"],
     "hbd":1,"hba":2,"detection_limit":"ppm"}),
    ("苯酚", {"phase":"固体","polarity":3,"bp_c":182,"logp":1.46,
     "carbon_count":6,"functional_groups":["酚羟基 (Ar-OH)","aromatic 芳香环"],
     "hbd":1,"hba":1,"detection_limit":"ppm"}),
    ("甲醇", {"phase":"液体","polarity":5,"bp_c":65,"logp":-0.77,
     "carbon_count":1,"functional_groups":["-OH 脂肪醇羟基"],
     "hbd":1,"hba":1,"detection_limit":"ppm"}),
]

for name, q in cases:
    print(f"\n{'='*60}")
    print(f"  {name}")
    results = search(q, methods, vectors, k=1)
    if not results:
        continue
    sim, method = results[0]
    method_name = f"{method['source']['standard']} — {method['source']['name'][:35]}"
    print(f"  k-NN Top-1: {method_name} (sim={sim:.3f})")

    report = validate([q], method)
    print(f"  Hook 结果: {'❌ REJECT' if report.rejected else '⚠ WARN' if report.warns else '✅ ALL PASS'}")
    for r in report.results:
        icon = {"PASS":"✅","WARN":"⚠️","REJECT":"❌"}[r.status]
        if r.status != "PASS":
            print(f"    {icon} [{r.rule}] {r.message}")

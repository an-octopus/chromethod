#!/usr/bin/env python3
"""Debug: 查看每条匹配的逐维度贡献"""
import os
from chromethod.knn import load_methods, build_index
from chromethod.features import (encode_substance, encode_method,
                                  weighted_cosine, FUNC_GROUPS, FEATURE_WEIGHTS)

_dir = os.path.dirname(os.path.abspath(__file__))
methods = load_methods(os.path.join(_dir, "chromethod", "data"))
vectors = [encode_method(m) for m in methods]

# 维度标签
DIM_LABELS = (
    ["物态:气", "物态:液", "物态:固"] +
    ["极性"] + ["沸点"] + ["LogP"] + ["碳数"] +
    [f"FG:{f}" for f in FUNC_GROUPS] +
    ["HBD", "HBA", "检出限", "初温(挥发)", "不分流", "杂原子", "顶空",
     "π体系", "长链酯"]
)

def debug_search(substance, name, methods, vectors, k=3):
    qvec = encode_substance(substance)
    scores = []
    for i, (m, mvec) in enumerate(zip(methods, vectors)):
        s = weighted_cosine(qvec, mvec)
        scores.append((s, i))
    scores.sort(key=lambda x: -x[0])

    # 计算逐维贡献
    ranked = []
    for sim, idx in scores[:k]:
        mvec = vectors[idx]
        dim_contrib = []
        for d in range(len(qvec)):
            w = FEATURE_WEIGHTS[d]
            # 该维度对余弦分子的贡献
            contrib = w * qvec[d] * mvec[d]
            dim_contrib.append((contrib, DIM_LABELS[d], qvec[d], mvec[d]))
        dim_contrib.sort(key=lambda x: -abs(x[0]))
        ranked.append((sim, idx, dim_contrib))
    return ranked

def run_debug(substance, name, methods, vectors):
    print(f"\n{'='*70}")
    print(f"DEBUG: {name}")
    qvec = encode_substance(substance)
    print(f"查询特征 ({len(qvec)}维):")
    nonzero = [(DIM_LABELS[i], qvec[i]) for i in range(len(qvec)) if qvec[i] > 0]
    for label, val in nonzero:
        print(f"  {label}: {val:.2f}")

    ranked = debug_search(substance, name, methods, vectors, k=3)
    for rank, (sim, idx, contribs) in enumerate(ranked):
        m = methods[idx]
        col = m["columns"][0]
        det = m["detectors"][0]
        print(f"\n  #{rank+1} {m['method_id']} sim={sim:.3f}")
        print(f"     {col['brand']} {det['type']} | {m['source']['name'][:50]}")

        # 查看方法关键特征
        mvec = vectors[idx]
        print(f"    方法非零特征:")
        for i in range(len(mvec)):
            if mvec[i] > 0:
                print(f"      {DIM_LABELS[i]}: {mvec[i]:.2f}")

        # Top-5 正贡献和负贡献
        print(f"    Top匹配维度:")
        for contrib, label, qv, mv in contribs[:8]:
            arrow = "↑↑" if qv > 0.5 and mv > 0.5 else ("↑" if qv > 0 and mv > 0 else " ")
            print(f"      {label:20s} q={qv:.2f} m={mv:.2f} → {contrib:+.4f} {arrow}")
        print(f"    最大分歧维度:")
        neg = [(c, l, qv, mv) for c, l, qv, mv in contribs if qv > 0.3 and mv < 0.1]
        for contrib, label, qv, mv in neg[:3]:
            print(f"      {label:20s} q={qv:.2f} m={mv:.2f} → 漏匹配")

# ── 逐个调试 ──
tests = [
    ("马拉硫磷", {"phase":"液体","polarity":3,"bp_c":156,"logp":2.36,
     "carbon_count":10,"functional_groups":["-COO- 酯基","硫醚 -S-","磷酸酯/有机磷 (P=O/S)"],
     "hbd":0,"hba":6,"detection_limit":"ppb"}),
    ("硫化氢", {"phase":"气体","polarity":1,"bp_c":-60,"logp":0.23,
     "carbon_count":0,"functional_groups":["硫醇 -SH (H₂S/硫醇)"],
     "hbd":0,"hba":0,"detection_limit":"ppm"}),
    ("棕榈酸甲酯", {"phase":"液体","polarity":2,"bp_c":330,"logp":6.4,
     "carbon_count":17,"functional_groups":["-COO- 酯基"],
     "hbd":0,"hba":2,"detection_limit":"ppm"}),
    ("咖啡因", {"phase":"固体","polarity":3,"bp_c":178,"logp":-0.07,
     "carbon_count":8,"functional_groups":["-C=O 酮/醛羰基","含氮杂环 (吡啶/嘌呤等)","alkene 烯烃"],
     "hbd":0,"hba":3,"detection_limit":"ppb"}),
]

for name, q in tests:
    run_debug(q, name, methods, vectors)

# 检查关键方法有哪些官能团
print(f"\n{'='*70}")
print("关键方法官能团检查")
for mid in ["GC-STD-017", "GC-STD-013", "GC-STD-023", "GC-STD-024", "GC-STD-027"]:
    for m in methods:
        if m["method_id"] == mid:
            mvec = encode_method(m)
            fgs = [(i, FUNC_GROUPS[i], mvec[7+i]) for i in range(17) if mvec[7+i] > 0]
            print(f"  {mid} ({m['source']['name'][:35]}):")
            for idx, fg, val in fgs:
                print(f"    [{idx}] {fg}: {val}")
            break

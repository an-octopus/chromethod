#!/usr/bin/env python3
"""AIGC Demo: k-NN 检索端到端"""
import json, os
from chromethod.knn import load_methods, build_index, search, format_result

_dir = os.path.dirname(os.path.abspath(__file__))
methods = load_methods(os.path.join(_dir, "chromethod", "data"))
vectors = build_index(methods)
print(f"已加载 {len(methods)} 条标准方法\n")

# ── 测试查询 1: 苯 (非极性 VOC) ──
print("=" * 60)
print("查询: 苯 — 非极性, 沸点80°C, 芳香, ppb级")
benzene = {
    "phase": "液体",
    "polarity": 1,
    "bp_c": 80.1,
    "logp": 2.13,
    "carbon_count": 6,
    "functional_groups": ["aromatic (芳香环)"],
    "hbd": 0,
    "hba": 0,
    "detection_limit": "ppb",
}
results = search(benzene, methods, vectors, k=3)
for sim, m in results:
    print(format_result(sim, m))
    print()

# ── 测试查询 2: 甲醇 (极性小分子) ──
print("=" * 60)
print("查询: 甲醇 — 强极性, 沸点65°C, -OH, ppm级")
methanol = {
    "phase": "液体",
    "polarity": 5,
    "bp_c": 65,
    "logp": -0.77,
    "carbon_count": 1,
    "functional_groups": ["-OH (醇/酚羟基)"],
    "hbd": 1,
    "hba": 1,
    "detection_limit": "ppm",
}
results = search(methanol, methods, vectors, k=3)
for sim, m in results:
    print(format_result(sim, m))
    print()

# ── 测试查询 3: DDT (有机氯, 半挥发性) ──
print("=" * 60)
print("查询: 4,4'-DDT — 弱极性, 沸点260°C, 含氯, ppt级")
ddt = {
    "phase": "固体",
    "polarity": 2,
    "bp_c": 260,
    "logp": 6.91,
    "carbon_count": 14,
    "functional_groups": ["-Cl (有机氯)", "aromatic (芳香环)"],
    "hbd": 0,
    "hba": 0,
    "detection_limit": "ppt",
}
results = search(ddt, methods, vectors, k=3)
for sim, m in results:
    print(format_result(sim, m))
    print()

# ── 测试查询 4: 萘 (PAH) ──
print("=" * 60)
print("查询: 萘 — 非极性, 沸点218°C, 多环芳烃, ppb级")
naphthalene = {
    "phase": "固体",
    "polarity": 1,
    "bp_c": 218,
    "logp": 3.3,
    "carbon_count": 10,
    "functional_groups": ["aromatic (芳香环)"],
    "hbd": 0,
    "hba": 0,
    "detection_limit": "ppb",
}
results = search(naphthalene, methods, vectors, k=3)
for sim, m in results:
    print(format_result(sim, m))
    print()

# ── 测试查询 5: 油酸甲酯 (FAME) ──
print("=" * 60)
print("查询: 油酸甲酯 — 中等极性, 沸点~350°C, 酯, ppm级")
oleate = {
    "phase": "液体",
    "polarity": 3,
    "bp_c": 350,
    "logp": 7.5,
    "carbon_count": 19,
    "functional_groups": ["-C=O (羰基/酮/醛)", "alkene (烯烃)"],
    "hbd": 0,
    "hba": 2,
    "detection_limit": "ppm",
}
results = search(oleate, methods, vectors, k=3)
for sim, m in results:
    print(format_result(sim, m))
    print()

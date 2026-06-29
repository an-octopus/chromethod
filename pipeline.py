#!/usr/bin/env python3
"""
AIGC 端到端流水线（LLM 集成版）
用法: python pipeline.py "茶叶中测咖啡因和茶碱，ppb级"
       python pipeline.py                     # 交互模式
"""
import os, sys, json

from chromethod.knn import load_methods, build_index, search
from chromethod.hook import validate
from chromethod.features import FUNC_GROUPS
from llm_agent import parse_user_input, pubchem_lookup, generate_report

_dir = os.path.dirname(os.path.abspath(__file__))
methods = load_methods(os.path.join(_dir, "chromethod", "data"))
vectors = build_index(methods)

# ── 官能团估算（简化版，不调LLM）──
def estimate_functional_groups(smiles: str, name_cn: str) -> list:
    """从SMILES和中文名估算官能团"""
    fgs = []
    s = smiles.lower()
    n = name_cn
    # 芳环
    if "c" in s and ("cc" in s or "c1" in s or "c(" in s):
        fgs.append("aromatic 芳香环")
    # 烯烃
    if "=" in s and "c=c" not in s.replace("c1", "").replace("c2", ""):
        if "c=c" not in s or s.count("c") > s.count("="):
            pass
    if "C=C" in s.upper() and "O=" not in s.upper():
        fgs.append("alkene 烯烃")
    # -OH 醇
    if ("O" in s or "o" in s) and "C-O" not in s and ("[OH]" not in s and "O=" not in s):
        if "醇" in n or "alcohol" in s:
            fgs.append("-OH 脂肪醇羟基")
    elif "酚" in n or "phenol" in s:
        fgs.append("酚羟基 (Ar-OH)")
    # -COOH
    if "C(=O)O" in s or "carboxylic" in s or "酸" in n:
        fgs.append("-COOH 羧酸")
    # 酯
    if "C(=O)OC" in s or "C(=O)O[C]" in s or "酯" in n or "FAME" in s.upper():
        fgs.append("-COO- 酯基")
    # 酮/醛
    elif "C(=O)" in s or "C=O" in s or "酮" in n or "醛" in n or "one" in s:
        fgs.append("-C=O 酮/醛羰基")
    # 胺
    if "N" in s and "amine" in s or "胺" in n:
        fgs.append("-NH2 伯/仲胺")
    # 氮杂环
    if any(k in n for k in ("吡啶","嘌呤","咖啡因","茶碱","三嗪","咪唑","唑")):
        fgs.append("含氮杂环 (吡啶/嘌呤等)")
    # 硝基
    if "[N+](=O)[O-]" in s or "硝基" in n or "nitro" in s:
        fgs.append("-NO2 硝基")
    # 磷酸酯
    if "P(=O)" in s or "P=S" in s or "磷" in n or "phosphate" in s:
        fgs.append("磷酸酯/有机磷 (P=O/S)")
    # 卤素
    if "Cl" in s or "氯" in n: fgs.append("-Cl 有机氯")
    if "Br" in s or "溴" in n: fgs.append("-Br 有机溴")
    if "F" in s and "Fluor" in s or "氟" in n: fgs.append("-F 有机氟")
    # 硫
    if "S" in s and "sulfide" in s or "硫醇" in n or "SH" in s or "H2S" in s.upper():
        fgs.append("硫醇 -SH (H₂S/硫醇)")
    elif "S" in s and ("sulfide" in s or "硫醚" in n or "thio" in s or "disulfide" in s):
        fgs.append("硫醚 -S-")
    elif "S" in s or "硫" in n:
        fgs.append("硫醚 -S-")
    # 醚
    if "C-O-C" in s or "醚" in n or "ether" in s or "epoxy" in s:
        fgs.append("-O- 醚")
    return fgs


def infer_polarity(fgs: list, logp: float | None) -> int:
    """从官能团和LogP推断极性等级 1-5"""
    score = 0
    if any(k in fgs for k in ("-COOH 羧酸", "磷酸酯/有机磷 (P=O/S)")):
        score += 2
    if any(k in fgs for k in ("-OH 脂肪醇羟基", "酚羟基 (Ar-OH)", "-NH2 伯/仲胺")):
        score += 1
    if any(k in fgs for k in ("-COO- 酯基", "-C=O 酮/醛羰基", "-O- 醚", "含氮杂环 (吡啶/嘌呤等)")):
        score += 0.5
    if logp is not None and logp < 0:
        score += 1
    if logp is not None and logp > 3:
        score -= 1
    return max(1, min(5, round(2 + score)))


def estimate_carbon_count(smiles: str) -> int:
    """从SMILES估算碳数"""
    import re
    carbons = re.findall(r'C(?![a-z])', smiles)
    return len(carbons) if carbons else 5


def run(user_input: str, verbose: bool = True):
    """端到端流水线"""
    if verbose:
        print("=" * 60)
        print(f"  AIGC 色谱方法推荐")
        print(f"  输入: {user_input}")
        print("=" * 60)

    # 1. LLM 解析
    if verbose: print("\n[1/5] LLM 解析输入...")
    parsed = parse_user_input(user_input)
    substances = parsed.get("substances", [])
    matrix = parsed.get("matrix", "")
    det_limit = parsed.get("detection_limit", "ppm")
    if verbose:
        print(f"  识别物质: {[s['name_cn'] for s in substances]}")
        print(f"  基质: {matrix}, 检出限: {det_limit}")

    # 2. PubChem 补属性
    if verbose: print("\n[2/5] PubChem 查询属性...")
    enriched = []
    for sub in substances:
        pub = pubchem_lookup(sub["name_en"])
        bp = pub.get("bp_c") or 100
        logp = pub.get("logp") or 0
        smiles = pub.get("smiles", "")
        fgs = estimate_functional_groups(smiles, sub["name_cn"])
        carbon = estimate_carbon_count(smiles)
        polarity = infer_polarity(fgs, logp)
        phase = "气体" if (bp is not None and bp < 25) else "液体" if (bp is not None and bp < 300) else "固体"
        sub_q = {
            "phase": phase,
            "polarity": polarity,
            "bp_c": bp,
            "logp": logp,
            "carbon_count": carbon,
            "functional_groups": fgs,
            "hbd": pub.get("hbd", 0) or 0,
            "hba": pub.get("hba", 0) or 0,
            "detection_limit": det_limit,
        }
        enriched.append(sub_q)
        if verbose:
            print(f"  {sub['name_cn']}: bp={bp}°C logP={logp} 极性={polarity} fgs={fgs}")

    # 3. 用沸点最高物质跑 k-NN
    if verbose: print("\n[3/5] k-NN 检索...")
    main = max(enriched, key=lambda s: s["bp_c"])
    top3 = search(main, methods, vectors, k=3)
    if verbose:
        for sim, m in top3:
            print(f"  {m['method_id']} sim={sim:.3f} {m['source']['name'][:40]}")

    # 4. Hook 验证所有物质
    if verbose: print("\n[4/5] Hook 验证...")
    validated = []
    for sim, method in top3:
        report = validate(enriched, method)
        validated.append((sim, method, report))
        if verbose:
            status = "❌ REJECT" if report.rejected else "⚠ WARN" if report.warns else "✅ PASS"
            print(f"  {method['method_id']}: {status}")

    # 5. LLM 生成报告
    if verbose: print("\n[5/5] 生成推荐报告...")
    report = generate_report(user_input, validated, [r for _, _, r in validated])
    if verbose:
        print("\n" + "=" * 60)
        print(report)
        print("=" * 60)
    return report


# ── 入口 ──
if __name__ == "__main__":
    if len(sys.argv) > 1:
        text = " ".join(sys.argv[1:])
    else:
        text = input("请输入分析需求: ")
    run(text)

"""
端到端流水线（无 LLM 依赖）
输入: 结构化物质列表 → 输出: Top-3 方法 + Hook 报告
LLM 解析和报告生成由 Claude Code 调用方处理
"""
import os, json
from .knn import load_methods, build_index, search
from .hook import validate
from .features import column_substitute

_pkg_dir = os.path.dirname(os.path.abspath(__file__))
_methods_cache = None
_vectors_cache = None


def _init():
    global _methods_cache, _vectors_cache
    if _methods_cache is None:
        _methods_cache = load_methods(os.path.join(_pkg_dir, "data"))
        _vectors_cache = build_index(_methods_cache)
    return _methods_cache, _vectors_cache


def recommend(substances: list[dict], top_k: int = 3,
              available_columns: list[str] | None = None) -> list[dict]:
    """
    输入: [{"name":"苯","bp_c":80,"polarity":1,"functional_groups":[...],...}, ...]
    available_columns: 用户手头的柱品牌列表, 如 ["DB-5","DB-1701","DB-WAX"]
    输出: Top-K 方法 + Hook 验证结果 + 柱替代建议(如有)
    策略: 沸点最高物质主导，被Hook否决则回退到下一高沸点物质
    """
    methods, vectors = _init()

    # 按沸点降序排列，用作回退链
    sorted_subs = sorted(substances,
                         key=lambda s: s.get("bp_c", 100) or 100,
                         reverse=True)
    seen_ids = set()
    output = []

    for anchor_sub in sorted_subs:
        candidates = search(anchor_sub, methods, vectors, k=top_k)

        for sim, method in candidates:
            mid = method["method_id"]
            if mid in seen_ids:
                continue
            seen_ids.add(mid)

            col = method["columns"][0]
            det = method["detectors"][0]
            hook_rpt = validate(substances, method)
            output.append({
                "method_id": mid,
                "similarity": round(sim, 3),
                "source": method["source"],
                "column": {
                    "brand": col["brand"],
                    "phase": col["phase"],
                    "specs": f"{col['length_m']}m×{col['id_mm']}mm×{col['film_um']}µm",
                    "temp_limit_c": col.get("temp_limit_c"),
                    "oven_program": col["oven_program"],
                },
                "carrier_gas": method["carrier_gas"],
                "detector": det,
                "preparation": method["preparation"],
                "inlet": method.get("inlet", {}),
                "applicability": method["applicability"],
                "hook": {
                    "rejected": hook_rpt.rejected,
                    "warnings": [{"rule": r.rule, "message": r.message} for r in hook_rpt.warns],
                    "errors": [{"rule": r.rule, "message": r.message} for r in hook_rpt.errors],
                    "all_pass": not hook_rpt.rejected and not hook_rpt.warns,
                },
                "column_alternatives": None,
            })

            # 柱替代: 用户手头没有推荐柱时找最接近的
            if available_columns:
                subs = column_substitute(
                    col["brand"], col["phase"], available_columns)
                if subs:
                    output[-1]["column_alternatives"] = [
                        {"brand": b, "grade": g, "note": n} for b, g, n in subs
                    ]

            if len(output) >= top_k:
                break

        if len(output) >= top_k:
            break

    # 混合样兼容性检查
    if len(substances) > 1 and output:
        _check_mixture_compat(substances, output)

    return output


def _check_mixture_compat(substances: list, output: list):
    """检查推荐方法是否覆盖所有物质的官能团，不覆盖则加警告"""
    from .features import encode_method, FUNC_GROUPS
    for r in output:
        for m in _methods_cache:
            if m["method_id"] == r["method_id"]:
                mvec = encode_method(m)
                uncovered = []
                for sub in substances:
                    fgs = sub.get("functional_groups", [])
                    missed = [fg for fg in fgs if fg in FUNC_GROUPS and mvec[7 + FUNC_GROUPS.index(fg)] < 0.5]
                    if missed:
                        uncovered.append((sub.get("name", "?"), missed))
                if len(uncovered) > len(substances) * 0.3:
                    names = ", ".join([f"{n}({','.join(f)})" for n, f in uncovered[:3]])
                    r["hook"]["warnings"].append({
                        "rule": "混合兼容性",
                        "message": f"该方法未覆盖{len(uncovered)}/{len(substances)}个组分的官能团: {names}。可能需分开分析"
                    })
                break


def list_methods() -> list:
    """列出所有已索引方法"""
    methods, _ = _init()
    return [{
        "id": m["method_id"],
        "source": m["source"]["standard"],
        "name": m["source"]["name"],
        "column": m["columns"][0]["brand"],
        "detector": m["detectors"][0]["type"],
        "category": m["applicability"]["analyte_category"],
    } for m in methods]

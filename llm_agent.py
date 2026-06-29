"""
LLM 代理：前端解析 + PubChem 查属性 + 后端解释
"""
import json, os, urllib.request, urllib.error
from anthropic import Anthropic

# 尝试多种方式获取 API key
_api_key = os.environ.get("ANTHROPIC_API_KEY") or os.environ.get("ANTHROPIC_AUTH_TOKEN") or ""
_base_url = os.environ.get("ANTHROPIC_BASE_URL") or "https://api.anthropic.com"

# fallback: 从 .claude 配置中读取
if not _api_key:
    creds_file = os.path.expanduser("~/.claude/credentials.json")
    if os.path.exists(creds_file):
        try:
            with open(creds_file) as f:
                creds = json.load(f)
            _api_key = creds.get("apiKey_hidden", "") or creds.get("api_key", "")
        except: pass

_client = Anthropic(api_key=_api_key, base_url=_base_url) if _api_key else None

MODEL = "claude-sonnet-4-6"


# ═══════════════════════════════════════
# PubChem 查询
# ═══════════════════════════════════════
def pubchem_lookup(name_en: str) -> dict:
    """查询 PubChem 获取化合物属性"""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name_en}/property/MolecularWeight,XLogP,HBondDonorCount,HBondAcceptorCount,BoilingPoint,CanonicalSMILES/JSON"
    try:
        with urllib.request.urlopen(url, timeout=10) as resp:
            data = json.loads(resp.read())
            props = data["PropertyTable"]["Properties"][0]
            return {
                "name_en": name_en,
                "mw": props.get("MolecularWeight"),
                "logp": props.get("XLogP"),
                "hbd": props.get("HBondDonorCount", 0) or 0,
                "hba": props.get("HBondAcceptorCount", 0) or 0,
                "bp_c": props.get("BoilingPoint"),
                "smiles": props.get("CanonicalSMILES", ""),
            }
    except Exception:
        return {"name_en": name_en, "mw": None, "logp": None, "hbd": 0, "hba": 0, "bp_c": None, "smiles": ""}


# ═══════════════════════════════════════
# 前端：自然语言→结构化
# ═══════════════════════════════════════
PARSE_PROMPT = """你是分析化学专家。从用户输入的色谱分析需求中提取结构化信息。

输出JSON，不要加任何解释:
{{
  "substances": [{{"name_cn": "中文名", "name_en": "English name"}}],
  "matrix": "水溶液/有机溶剂/气体/固体",
  "detection_limit": "ppm/ppb/ppt",
  "note": "补充说明"
}}

用户输入: {user_input}
"""


def parse_user_input(text: str) -> dict:
    """LLM 解析用户自然语言输入"""
    if _client is None:
        raise RuntimeError("No API key found. Set ANTHROPIC_API_KEY or place credentials in ~/.claude/credentials.json")
    msg = _client.messages.create(
        model=MODEL,
        max_tokens=500,
        temperature=0,
        messages=[{"role": "user", "content": PARSE_PROMPT.format(user_input=text)}],
    )
    raw = msg.content[0].text.strip()
    # 去掉可能的 markdown code block
    if raw.startswith("```"):
        raw = raw.split("\n", 1)[1]
        if raw.endswith("```"):
            raw = raw[:-3]
    return json.loads(raw)


# ═══════════════════════════════════════
# 后端：k-NN+Hook 结果→人话报告
# ═══════════════════════════════════════
REPORT_PROMPT = """你是分析化学专家。根据以下信息给用户写一份色谱方法推荐报告。

## 用户需求
{{user_input}}

## 推荐方法 (k-NN 检索 Top-3)
{{top3_summary}}

## 验证结果
{{hook_summary}}

## 要求
1. 用中文写，面向实验室分析人员
2. 第一段：一句话总结推荐哪个方法
3. 第二段：解释为什么推荐这个方法（和样品特征的匹配点）
4. 第三段：如果有 WARN/REJECT，解释风险和建议
5. 如果所有方法都被 Hook 否决，诚实告诉用户，给出分两次跑或其他替代建议
6. 不要编造方法参数。参数从上面的推荐方法中引用
7. 300字以内"""


def generate_report(user_input: str, top3: list, hook_reports: list) -> str:
    """生成面向用户的推荐报告"""
    top3_text = ""
    for i, (sim, method, hook_rpt) in enumerate(top3):
        col = method["columns"][0]
        det = method["detectors"][0]
        top3_text += f"\n#{i+1} [{method['method_id']}] {method['source']['standard']} — {method['source']['name']}\n"
        top3_text += f"  相似度: {sim:.3f}\n"
        top3_text += f"  柱: {col['brand']} {col['length_m']}m×{col['id_mm']}mm×{col['film_um']}µm\n"
        top3_text += f"  检测器: {det['type']}\n"
        top3_text += f"  前处理: {method['preparation']['technique']}\n"
        if hook_rpt:
            top3_text += f"  Hook: {hook_rpt.summary()}\n"

    hook_text = "\n".join([r.summary() for r in hook_reports]) if hook_reports else "全部通过"

    if _client is None:
        return "[API key not configured]"
    prompt = REPORT_PROMPT.replace("{{user_input}}", user_input)
    prompt = prompt.replace("{{top3_summary}}", top3_text)
    prompt = prompt.replace("{{hook_summary}}", hook_text)
    msg = _client.messages.create(
        model=MODEL,
        max_tokens=600,
        temperature=0.2,
        messages=[{"role": "user", "content": prompt}],
    )
    return msg.content[0].text.strip()

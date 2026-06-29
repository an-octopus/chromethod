# AIGC — 气相色谱方法起点推荐

输入样品特征，推荐最接近的标准 GC 方法。

## 快速开始

Python 3.11+，零外部依赖。

```bash
pip install git+https://github.com/an-octopus/chromethod.git
```

告诉它你的组分及其浓度，如果是纯物质会提供溶剂选择，然后返回推荐的方法：

```python
from chromethod import recommend

# 测水中苯，ppb 级
results = recommend([
    {"bp_c": 80, "polarity": 1, "functional_groups": ["aromatic 芳香环"],
     "phase": "液体", "hbd": 0, "hba": 0, "detection_limit": "ppb"},
])
# → 返回标准方法，每条带完整色谱条件 + Hook 验证结果
```

手头柱子不全？传 `available_columns`，自动匹配最佳替代。

## 架构

```
样品特征 → k-NN 检索 (余弦距离) → Hook 验证 → 推荐报告
              ↑ 标准方法库
```

- **k-NN**：找固定相、温度程序、检测器都最接近的已知方法，每条可溯源到 EPA/USP/Restek 原文
- **Hook**：在检索外面作为规则
- **柱替代**：推荐柱你没有的时候，可以自动从手头柱子里找最接近的

## 数据

所有方法来自 Restek 官网色谱图库、EPA、USP 及知网文献。覆盖 FID、MS、ECD、TCD、SCD、NPD、FPD。

## 路线图

- **RDKit 可选依赖**：有 RDKit 时用 Morgan 指纹 + 子结构匹配替代正则解析 SMILES，提升官能团识别精度和同分异构体区分能力；没有则 fallback 当前零依赖方案
- **负向推荐 / 故障归因**：从"该用什么方法"扩展到"用了之后哪里最容易翻车 + 症状 + 怎么修"，纯经验规则 Hook，零物理计算

## 局限

- 库很小，很多场景没覆盖
- 当前 SMILES 解析基于正则匹配，官能团识别精度有限，同分异构体暂无法区分（待 RDKit 可选依赖解决）
- 只做推荐，不负责安全。实验前请自行验证方法可行性
- 结果仅供参考，不能替代实验验证

## 许可

MIT

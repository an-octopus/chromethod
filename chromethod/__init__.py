"""
AIGC — AI 气相色谱方法推荐引擎
k-NN 检索 + Hook 物理约束
"""
from .pipeline import recommend, list_methods
from .hook import validate, HookReport, HookResult

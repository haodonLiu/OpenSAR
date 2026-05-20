"""CSAR - 化合物结构-活性关系分析工具包.

本包提供了用于分析化合物分子结构与其生物活性之间关系的完整工具链，
包括分子读取、聚类分析、最大公共子结构(MCS)查找、骨架聚类和SAR可视化等功能。

主要模块:
    - io: 分子文件读写模块
    - clustering: 分子聚类 (指纹相似度 + Bemis-Murcko 骨架)
    - mcs: 最大公共子结构查找模块
    - sar: 活性预处理 + 结构-活性关系分析
    - visualization: SAR结果可视化模块

示例:
    >>> from src import MoleculeReader, MolecularClusterer, ScaffoldClusterer
    >>> from src import ActivityPreprocessor, SARAnalyzer
    >>> reader = MoleculeReader()
    >>> molecules = reader.read("data.sdf")
"""

__version__ = "0.2.0"  # 版本号 (Phase 2 升级)

# 从各子模块导入主要类
from .io.reader import MoleculeReader  # 分子读取器
from .io.writer import MoleculeWriter  # 分子写入器
from .clustering.fingerprinter import MolecularFingerprinter  # 分子指纹计算器
from .clustering.cluster import MolecularClusterer  # 分子聚类器 (基于指纹)
from .clustering.scaffold import ScaffoldClusterer  # 骨架聚类器 (Bemis-Murcko)
from .mcs.finder import MCSFinder  # 最大公共子结构查找器
from .sar.preprocessor import ActivityPreprocessor  # 活性数据预处理器
from .sar.analyzer import SARAnalyzer  # SAR分析器
from .visualization.renderer import SARRenderer  # SAR可视化渲染器

# 公开API列表
__all__ = [
    "MoleculeReader",  # 分子读取器
    "MoleculeWriter",  # 分子写入器
    "MolecularFingerprinter",  # 分子指纹计算器
    "MolecularClusterer",  # 分子聚类器 (Tanimoto/Butina)
    "ScaffoldClusterer",  # 骨架聚类器 (Bemis-Murcko)
    "MCSFinder",  # 最大公共子结构查找器
    "ActivityPreprocessor",  # 活性数据预处理器
    "SARAnalyzer",  # SAR分析器
    "SARRenderer",  # SAR可视化渲染器
]

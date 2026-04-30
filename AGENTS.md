<!-- AGENTS.md for CSAR (Compound Structure-Activity Relationship) -->

# CSAR - 化合物结构-活性关系分析工作流

## 项目概述

CSAR (Compound Structure-Activity Relationship) 是一个用于分析分子活性数据的化学信息学 Python 包。核心工作流：

1. **读取分子**：从输入文件（SDF、Excel、CSV）读取分子结构和活性数据（如 IC50）
2. **分子聚类**：基于分子指纹相似性（Tanimoto 或 Butina 算法）对分子聚类
3. **查找 MCS**：为每个聚类查找最大公共子结构（Maximum Common Substructure）
4. **SAR 分析**：计算每个聚类的活性统计，识别取代基对活性的贡献
5. **生成可视化**：SAR 汇总图、相似度矩阵热图、SAR 表格图像

项目包含两套可视化系统：
- **主工作流**（`src/main.py`）：集成聚类、MCS、SAR 分析和基础可视化
- **SAR Image Generator**（`src/scripts/sar_image_generator.py`）：独立的增强版 SAR 图像生成器，支持层次聚类和更专业的 SAR 报告图像

---

## 技术栈

- **Python 3.10+**
- **RDKit** (>=2024.3.0) - 分子指纹、MCS、结构渲染、SMILES 解析
- **Pandas** (>=2.0.0) - Excel/CSV 数据处理
- **NumPy** (>=1.24.0) - 数值运算、相似度矩阵
- **NetworkX** (>=3.0) - 图操作
- **Matplotlib** (>=3.7.0) / **Seaborn** (>=0.12.0) - 统计图表
- **Pillow (PIL)** - SAR 表格图像处理

---

## 项目结构

```
openSAR/
├── src/                          # 源代码 (v0.2.0)
│   ├── __init__.py               # 包导出 (v0.2.0 新增 ScaffoldClusterer, ActivityPreprocessor)
│   ├── main.py                   # CLI 入口和完整工作流编排
│   ├── io/                       # 输入/输出
│   │   ├── __init__.py
│   │   ├── reader.py             # MoleculeReader - 读取 SDF/Excel/CSV
│   │   └── writer.py             # MoleculeWriter - 写入 SDF/CSV/Excel（含描述符）
│   ├── clustering/               # 聚类算法
│   │   ├── __init__.py
│   │   ├── fingerprinter.py      # MolecularFingerprinter (Morgan/MACCS/RDKit)
│   │   ├── cluster.py            # MolecularClusterer (Tanimoto/Butina)
│   │   └── scaffold.py           # [v0.2.0] ScaffoldClusterer - Bemis-Murcko 骨架聚类
│   ├── mcs/                      # 最大公共子结构
│   │   ├── __init__.py
│   │   └── finder.py             # MCSFinder, MCSResult, ScaffoldInfo, SubstituentInfo
│   ├── sar/                      # SAR 分析 + 预处理
│   │   ├── __init__.py
│   │   ├── preprocessor.py       # [v0.2.0] ActivityPreprocessor - 单位转换/pActivity/异常值/删失数据
│   │   └── analyzer.py           # SARAnalyzer, SARResult
│   ├── scripts/                  # [已整理] 独立命令行工具
│   │   ├── __init__.py
│   │   ├── sar_image_generator.py  # 独立 SAR 图像生成器（原根目录）
│   │   └── proc.py                 # Excel 去重脚本（原根目录）
│   └── visualization/            # 可视化
│       ├── __init__.py
│       ├── utils.py              # [v0.2.0] 统一工具函数
│       ├── renderer.py           # SARRenderer, PlotSettings
│       └── sar_table.py          # SAR 表格（委托到 utils.py）
├── tests/                        # 单元测试
│   ├── __init__.py
│   ├── test_reader.py            # MoleculeReader 测试
│   ├── test_fingerprinter.py     # 指纹生成和相似度测试
│   ├── test_clustering.py        # 聚类算法测试 (Tanimoto/Butina)
│   ├── test_scaffold.py          # [v0.2.0] ScaffoldClusterer / Murcko 骨架测试
│   ├── test_mcs.py               # MCS 查找测试
│   ├── test_preprocessor.py      # [v0.2.0] ActivityPreprocessor 测试
│   ├── test_sar.py               # SAR 分析测试
│   └── test_render.py            # 快速渲染测试脚本（原根目录）
├── README.md                     # 项目介绍（中英双语）
├── sar_image_generator.py        # [已迁移至 src/scripts/]
├── test_render.py                # [已迁移至 tests/]
├── proc.py                       # [已迁移至 src/scripts/]
├── USAGE_SAR_IMAGE.md            # SAR 图像生成器使用说明（中文）
├── pyproject.toml                # 项目配置和依赖
├── setup.py                      # setuptools 配置
├── data/                         # 输入数据目录
└── output/                       # 输出目录
```

---

## 构建、测试和运行命令

### 安装

```bash
# 开发安装
pip install -e ".[dev]"

# 完整安装
pip install -e ".[all]"
```

### 运行测试

```bash
pytest                          # 运行所有测试
pytest tests/test_reader.py     # 单个测试文件
pytest tests/test_reader.py::test_read_sdf_file -v  # 单个测试函数
pytest --cov=src --cov-report=html  # 带覆盖率
pytest -k "clustering"          # 按模式匹配
```

### 代码检查与格式化

```bash
ruff check src/                 # 代码检查
ruff format src/                # 代码格式化
mypy src/                       # 类型检查

# 全部检查（模拟 pre-commit）
ruff check src/ && ruff format --check src/ && mypy src/
```

### 运行项目

**主 CLI 入口**（`setup.py` 定义 `csar=src.main:main`）：

```bash
# 安装后
pip install -e .

# 命令行使用
csar input.xlsx -o output/
python -m src.main input.xlsx -o output/

# 完整选项
python -m src.main input.xlsx \
    -o output/ \
    --smiles-column SMILES \
    --name-column Name \
    --activity-column Activity \
    --cluster-threshold 0.7 \
    --cluster-method tanimoto \
    --fp-type Morgan \
    --mcs-timeout 30 \
    --ic50-nm-column IC50_nM \
    --ic50-um-column IC50_uM
```

**SAR Image Generator**（独立脚本）：

```bash
# 生成所有聚类的 SAR 图像
python -m src.scripts.sar_image_generator input.csv output/sar/

# 只生成指定聚类
python -m src.scripts.sar_image_generator input.csv output/sar/ --cluster 1

# 自定义参数
python -m src.scripts.sar_image_generator input.csv output/sar/ --threshold 0.6 --max 8 --min-subs 6
```

**编程方式调用**：

```python
from src.main import run_workflow
from pathlib import Path

run_workflow(
    input_file=Path("data.xlsx"),
    output_dir=Path("output"),
    smiles_column="SMILES",
    activity_column="IC50_nM",
    cluster_threshold=0.7,
    skip_visualization=False,
)
```

### 构建分发包

```bash
python -m build                 # 构建分发包
pip install -e .                # 本地安装
```

---

## 代码风格指南

### 类型提示

所有函数签名必须有类型提示。使用 `from __future__ import annotations` 处理前向引用。

```python
from __future__ import annotations
from typing import Optional, List, Dict, Tuple, Union

def cluster_molecules(
    molecules: List[Dict[str, Any]],
    threshold: float = 0.7,
    method: str = "tanimoto"
) -> List[ClusterResult]:
    ...
```

使用 `Optional[X]` 而不是 `X | None`。使用 `List`、`Dict`、`Tuple` 而不是内置泛型。

### 导入组织

导入分三个部分，用空行分隔：

1. 标准库
2. 第三方包（RDKit、pandas 等）
3. 本地应用导入

```python
from __future__ import annotations

from pathlib import Path
from typing import List, Dict, Optional, Tuple

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

from src.clustering.fingerprinter import MolecularFingerprinter
```

避免通配符导入（`from rdkit import *`）。导入具体的类/函数。

### 命名规范

| 元素 | 规范 | 示例 |
|------|------|------|
| 模块 | 小写，snake_case | `cluster_engine.py` |
| 类 | PascalCase | `MolecularCluster` |
| 函数 | snake_case | `calculate_tanimoto()` |
| 变量 | snake_case | `mol_list`, `threshold` |
| 常量 | UPPER_SNAKE | `DEFAULT_THRESHOLD` |

### 模块结构（`__init__.py`）

```python
"""Short description of module."""

from .submodule import (
    MyClass,
    my_function,
)

__all__ = ["MyClass", "my_function"]
```

### 异常处理

每个模块有自己的异常类：

```python
class MoleculeReadError(Exception):
    """Failed to read molecule from file."""
    pass

class ClusteringError(Exception):
    """Clustering operation failed."""
    pass

class MCSError(Exception):
    """Maximum common substructure finding failed."""
    pass
```

优先捕获具体异常：

```python
try:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
except ValueError as e:
    raise MoleculeReadError(f"Cannot parse SMILES '{smiles}': {e}") from e
```

### 文档字符串

使用 Google 风格文档字符串。代码注释和文档字符串使用**中文**，类名/函数名/变量名使用英文。

```python
def find_mcs(molecules: List[Dict[str, Any]], timeout: int = 5) -> Tuple[Mol, float]:
    """查找分子间的最大公共子结构.

    使用 RDKit 的 MCS 算法识别所有输入分子共有的最大子结构。

    Args:
        molecules: 分子字典列表（至少 2 个）.
        timeout: 每个聚类的 MCS 搜索最大时间（秒）.

    Returns:
        (MCS 分子对象, MCS 得分) 元组.

    Raises:
        MCSError: 分子少于 2 个或 MCS 搜索超时.
    """
```

### RDKit 约定

- 总是检查 `mol is None` 后再使用
- SMARTS 分子转换为 SMILES 再转回普通分子用于渲染
- 使用 `rdMolDraw2D.MolDraw2DCairo` 生成高质量 PNG
- 使用 `Chem.MolFragmentToSmiles()` 提取取代基（不是 `PathToSubmol`）

### 分子字典约定

分子以字典形式传递，包含以下键：

```python
{
    "mol": Chem.Mol,           # RDKit 分子对象
    "smiles": str,             # SMILES 字符串
    "name": str,               # 分子名称
    "activity": float,         # 活性值（可选）
    "activity_raw": str,       # 原始活性值字符串（可选）
}
```

### 日志

每个模块使用 `logging.getLogger(__name__)`。日志消息使用中文。

```python
import logging

logger = logging.getLogger(__name__)
```

---

## 测试策略

### 测试文件

| 文件 | 覆盖内容 |
|------|----------|
| `test_reader.py` | MoleculeReader 初始化、方法存在性、不支持的格式错误 |
| `test_fingerprinter.py` | Morgan/MACCS/RDKit 指纹生成、相似度计算、成对矩阵、无效 SMILES/None 处理 |
| `test_clustering.py` | Tanimoto/Butina 聚类、空列表、结果结构、不同阈值 |
| `test_scaffold.py` | **[v0.2.0]** Murcko 骨架提取、ScaffoldClusterer 聚类/多样性统计/过滤 |
| `test_mcs.py` | MCS 查找（2 分子、相同分子、分子不足）、结果属性 |
| `test_preprocessor.py` | **[v0.2.0]** 单位检测转换、删失解析、pActivity 转换、IQR 异常值检测 |
| `test_sar.py` | SARAnalyzer 初始化、带/不带活性的聚类分析、总体统计 |

### 测试模式

- 使用 `pytest` 和 fixture
- 分子 fixture 返回 `dict`（含 `mol`, `smiles`, `name`, `activity` 键）或 `Chem.Mol`
- 验证正常路径和错误路径（无效 SMILES、None 分子、分子不足）
- 浮点数比较使用 `pytest.approx`

### 运行配置（`pyproject.toml`）

```toml
[tool.pytest.ini_options]
testpaths = ["tests"]
python_files = ["test_*.py"]
python_classes = ["Test*"]
python_functions = ["test_*"]
```

---

## 工具配置

### Ruff（`pyproject.toml`）

```toml
[tool.ruff]
line-length = 88
target-version = "py310"

[tool.ruff.lint]
select = ["E", "F", "W", "I", "N", "UP", "B", "C4", "SIM"]
ignore = ["E501"]
```

### MyPy（`pyproject.toml`）

```toml
[tool.mypy]
python_version = "3.10"
strict = true
warn_return_any = true
warn_unused_configs = true
disallow_untyped_defs = true
```

---

## 关键数据类

| 数据类 | 说明 | 关键字段 |
|--------|------|----------|
| `FingerprintResult` | 指纹结果 | `fingerprint`, `fp_type`, `n_bits` |
| `ClusterResult` | 聚类结果 (Tanimoto/Butina) | `cluster_id`, `molecules`, `representative_idx`, `size`, `avg_similarity` |
| `ScaffoldClusterResult` | **[v0.2.0]** 骨架聚类结果 | `scaffold_smiles`, `scaffold_mol`, `molecules`, `num_molecules`, `activities`, `mean_activity` |
| `MCSResult` | MCS 结果 | `mcs_mol`, `smiles`, `num_atoms`, `num_bonds`, `score` |
| `ScaffoldInfo` | MCS 骨架信息 | `scaffold_mol`, `scaffold_smiles`, `r_positions`, `num_r_groups` |
| `SubstituentInfo` | 取代基信息 | `position_idx`, `substituent_smiles`, `substituent_mol`, `activity` |
| `PreprocessingResult` | **[v0.2.0]** 预处理结果 | `values_original`, `values_processed`, `unit_detected`, `censoring_map`, `outlier_indices` |
| `SARResult` | SAR 分析结果 | `cluster_id`, `mcs`, `activities`, `mean_activity`, `contributions` |
| `PlotSettings` | 绘图设置 | `figure_size`, `dpi`, `color_scheme` |

---

## 开发工作流

### 提交前检查

1. `ruff check src/`
2. `ruff format src/`
3. `mypy src/`
4. `pytest`

### 添加新功能

1. 先写测试（TDD）
2. 在 `src/` 中实现
3. 更新 `__init__.py` 导出
4. 更新文档

---

## 注意事项

- **双语代码**：文档字符串和注释使用中文；变量/类名使用英文
- **SMARTS 渲染**：MCS 返回的是 SMARTS 模式分子，渲染前需转换为 SMILES 再转回普通分子
- **聚类大小限制**：主工作流中 MCS 查找跳过超过 30 个分子的聚类（`max_mcs_cluster_size = 30`）
- **SAR Image Generator 的层次聚类**：两层聚类（粗聚类 threshold-0.15，细聚类 threshold），每簇最多 10 个分子
- **R-group 动态数量** (`--r-groups N`)：SAR Image Generator 支持配置 MCS 最大 R 基团数（1-6，默认 2）
- **字体依赖**：可视化模块尝试加载系统字体（`DejaVuSans.ttf`、`arial.ttf`），失败时回退到默认字体
- **统一表格图像**：`create_sar_table_image()` 已从 `sar_table.py` 和 `renderer.py` 提取到 `utils.py`，消除代码重复

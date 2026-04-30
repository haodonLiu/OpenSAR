# OpenSAR

**化合物结构-活性关系分析工具包** | *Compound Structure-Activity Relationship Toolkit*

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue)](https://python.org)
[![RDKit](https://img.shields.io/badge/RDKit-2024.03%2B-green)](https://rdkit.org)
[![License](https://img.shields.io/badge/License-MIT-yellow)](LICENSE)

OpenSAR 是一个面向药物化学研究的 Python 工具包，提供从分子读取、聚类、骨架分析到 SAR 可视化的完整工作流。

---

## 📦 功能特性

| 模块 | 功能 |
|------|------|
| **分子读取** | 支持 SDF / Excel / CSV 格式，自动解析 SMILES 和活性数据 |
| **指纹聚类** | Morgan/MACCS/RDKit 指纹 + Tanimoto/Butina 聚类 |
| **骨架聚类** | Bemis-Murcko 通用骨架提取与聚类 (v0.2.0 新增) |
| **MCS 分析** | 最大公共子结构查找、R 基团提取 |
| **活性预处理** | 单位自动检测与转换、删失数据解析、pActivity 转换、异常值检测 (v0.2.0 新增) |
| **SAR 分析** | 聚类活性统计、取代基贡献分析 |
| **可视化** | SAR 汇总图、相似度矩阵热图、带取代基图像的 SAR 表格 |

---

## 🚀 快速开始

### 安装

```bash
# 克隆仓库
git clone https://github.com/haodonLiu/OpenSAR.git
cd OpenSAR

# 开发安装（推荐）
pip install -e ".[dev]"
```

### 基本使用

```python
from src import MoleculeReader, MolecularClusterer, MCSFinder, SARAnalyzer, SARRenderer

# 1. 读取分子
reader = MoleculeReader()
molecules = reader.read("data.xlsx", smiles_col="SMILES", activity_col="IC50_nM")

# 2. 指纹聚类
clusterer = MolecularClusterer(threshold=0.7)
clusters = clusterer.cluster(molecules)

# 3. 查找 MCS
finder = MCSFinder()
mcs_results = finder.find_all(clusters)

# 4. SAR 分析
analyzer = SARAnalyzer()
sar_results = analyzer.analyze_clusters(clusters, mcs_results)

# 5. 可视化
renderer = SARRenderer()
renderer.render_summary(clusters, mcs_results, sar_results, output_dir="output/")
```

### 骨架聚类

```python
from src import ScaffoldClusterer

scaffold_clusterer = ScaffoldClusterer()
results = scaffold_clusterer.cluster(molecules)
print(f"Unique scaffolds: {len(results)}")

# 骨架多样性分析
diversity = scaffold_clusterer.get_scaffold_diversity(molecules)
print(f"Diversity index: {diversity['scaffold_diversity_index']}")
```

### 活性数据预处理

```python
from src.sar import ActivityPreprocessor

prep = ActivityPreprocessor()
preprocessed = prep.preprocess(
    ["150", "3.2", ">10000", "12.5", "<1.0"],
    raw=True,
)
print(preprocessed.values_processed)  # [150.0, 3.2, 10000.0, 12.5, 1.0]
print(preprocessed.unit_detected)     # "nM"
print(preprocessed.p_activity)        # [-log10] 转换值
```

### 命令行 SAR 图像生成

```bash
# 从 CSV 生成所有聚类的 SAR 图片
python -m src.scripts.sar_image_generator data.csv output/sar/

# 自定义参数
python -m src.scripts.sar_image_generator data.csv output/sar/ \
    --threshold 0.6 --max 8 --r-groups 4

# 仅处理指定聚类
python -m src.scripts.sar_image_generator data.csv output/sar/ --cluster 1
```

### 数据去重

```bash
python -m src.scripts.proc input.xlsx output.xlsx
```

---

## 📁 项目结构

```
OpenSAR/
├── src/                          # 核心源代码
│   ├── __init__.py               # 包导出入口
│   ├── main.py                   # 主工作流 CLI
│   ├── io/                       # 分子文件读写
│   │   ├── reader.py             # MoleculeReader
│   │   └── writer.py             # MoleculeWriter
│   ├── clustering/               # 聚类算法
│   │   ├── fingerprinter.py      # 分子指纹
│   │   ├── cluster.py            # Tanimoto/Butina 聚类
│   │   └── scaffold.py           # Bemis-Murcko 骨架聚类
│   ├── mcs/                      # 最大公共子结构
│   │   └── finder.py             # MCSFinder
│   ├── sar/                      # SAR 分析
│   │   ├── preprocessor.py       # 活性数据预处理
│   │   └── analyzer.py           # SARAnalyzer
│   ├── visualization/            # 可视化
│   │   ├── renderer.py           # SARRenderer
│   │   ├── sar_table.py          # SAR 表格
│   │   └── utils.py              # 统一工具函数
│   └── scripts/                  # 命令行工具
│       ├── sar_image_generator.py  # SAR 图像生成器
│       └── proc.py                 # Excel 去重工具
├── tests/                        # 单元测试
├── docs/                         # 文档
├── data/                         # 输入数据 (git 忽略)
├── output/                       # 输出结果 (git 忽略)
├── AGENTS.md                     # AI 代理工作指南
├── USAGE_SAR_IMAGE.md            # SAR 图像生成器详细说明
├── setup.py                      # 安装配置
└── pyproject.toml                # 项目配置
```

---

## 🧪 运行测试

```bash
# 运行所有测试
pytest

# 运行特定测试
pytest tests/test_scaffold.py -v

# 带覆盖率
pytest --cov=src --cov-report=html
```

---

## 🛠 技术栈

- **Python 3.10+**
- **RDKit** — 分子指纹、MCS、结构渲染、SMILES 解析
- **Pandas / NumPy** — 数据处理
- **NetworkX** — 图操作
- **Matplotlib / Seaborn** — 统计图表
- **Pillow** — 表格图像处理

---

## 📄 许可证

[MIT License](LICENSE)

---

## 👤 作者

[haodonLiu](https://github.com/haodonLiu)

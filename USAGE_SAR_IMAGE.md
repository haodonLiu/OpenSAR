# SAR图片生成器使用说明

## 功能特性

✅ **层次聚类**
- 两层嵌套聚类：大类 → 小类
- 第一层粗聚类（较低相似度阈值）
- 第二层细聚类（较高相似度阈值）

✅ **MCS限制**
- 每个MCS最多2个取代基（R基团）
- 超过限制时自动过滤，保留取代基最多的2个

✅ **SAR图片**
- 顶部：MCS骨架结构（带R1, R2标签）
- 底部：SAR表格（取代基结构 + 波浪线连接符 + 活性数据 + 理化性质）
- 颜色编码：绿色=高活性，红色=低活性

## 使用方法

### 基本用法

```bash
# 生成所有簇的SAR图片
python sar_image_generator.py input.csv output/sar/

# 生成特定簇的SAR图片
python sar_image_generator.py input.csv output/sar/ --cluster 1

# 自定义聚类参数
python sar_image_generator.py input.csv output/sar/ --threshold 0.6 --max 8
```

### 参数说明

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `input.csv` | 输入CSV文件路径 | 必需 |
| `output/sar/` | 输出目录 | `output/sar/` |
| `--cluster ID` | 只生成指定簇（1-based） | 所有簇 |
| `--threshold N` | 聚类相似度阈值（0-1） | 0.7 |
| `--max N` | 每簇最大分子数 | 10 |

### 输出示例

```
SAR Image Generation
  Input: data.csv
  Output: output/sar/
  Clustering threshold: 0.7
  Max molecules per cluster: 8

2026-04-16 - INFO - Reading CSV: data.csv
2026-04-16 - INFO - Parsed 244 valid molecules
2026-04-16 - INFO - Clustering 244 molecules...
2026-04-16 - INFO - Level 1: Coarse clustering (threshold=0.55)...
2026-04-16 - INFO - Created 64 coarse clusters
2026-04-16 - INFO - Level 2: Fine clustering...
2026-04-16 - INFO - Final: 100 clusters
2026-04-16 - INFO - === Processing Cluster 1/100 (2 molecules) ===
2026-04-16 - INFO -   MCS has 2 R-group(s)
2026-04-16 - INFO - Saved: output/sar/sar_cluster_1.png

Generated 85 SAR image(s):
  - output/sar/sar_cluster_1.png
  - output/sar/sar_cluster_2.png
  ...
```

## CSV文件格式

```csv
Compound_ID,SMILES,IC50_nM
NJTQYSY202203005,Cc1ccccc1,150
NJTQYSY202203009,CC(C1=CN(C)N=C1)...,15630
```

**必需列：**
- `SMILES`: 分子SMILES字符串

**可选列：**
- `Compound_ID`: 分子名称
- `IC50_nM`: 活性值（nM）

## 生成的图片结构

```
┌─────────────────────────────────────────────┐
│   SAR Analysis - Cluster 1 (5 compounds)    │
├─────────────────────────────────────────────┤
│                                             │
│         [MCS Scaffold with R labels]        │
│              R1        R2                   │
│               \        /                    │
│                [骨架结构]                     │
│                                             │
├─────────────────────────────────────────────┤
│  R1     │  R2     │ IC50(nM) │ MW  │ LogP  │
├─────────┼─────────┼──────────┼─────┼───────┤
│ ~[结构] │   H     │  21.0    │ 450 │ 3.2   │  ← 高活性（绿色）
│ ~[结构] │ ~[结构] │ 150.5    │ 480 │ 3.5   │  ← 中活性
│ ~[结构] │   H     │ 15630    │ 520 │ 4.1   │  ← 低活性（红色）
└─────────────────────────────────────────────┘
```

## 聚类算法

### 层次聚类流程

```
244个分子
    ↓
Level 1: 粗聚类 (threshold=0.55)
    ↓
64个大类
    ↓
Level 2: 细聚类 (threshold=0.7, max=8)
    ↓
100个小类（每类≤8个分子）
    ↓
每个类生成一张SAR图片
```

### MCS取代基限制

```
找到MCS → 提取所有R基团
    ↓
R基团数量 > 2？
    ├─ 是 → 保留取代基最多的2个
    └─ 否 → 使用所有R基团
    ↓
生成SAR图片
```

## 注意事项

1. **性能优化**
   - 大数据集会进行层次聚类
   - 每个簇最多10个分子（可调整）
   - MCS查找使用8个分子随机采样

2. **R基团限制**
   - 最多2个R基团/每个MCS
   - 确保SAR表格清晰易读
   - 适合科研报告展示

3. **输出文件**
   - 每个簇一个PNG文件
   - 文件名：`sar_cluster_1.png`, `sar_cluster_2.png`, ...
   - 高分辨率（300 DPI）

## 常见问题

**Q: 为什么有些簇没有生成图片？**
A: 簇中分子数<2，或无法找到有效的MCS结构。

**Q: 如何调整聚类数量？**
A: 修改`--threshold`参数：
   - 较低阈值（0.5-0.6）→ 更多簇
   - 较高阈值（0.8-0.9）→ 更少簇

**Q: 如何增加每簇的分子数？**
A: 使用`--max N`参数，例如`--max 15`。

**Q: 为什么MCS只有1个R基团？**
A: 该簇分子结构差异较大，只在1个位置有共同取代基。

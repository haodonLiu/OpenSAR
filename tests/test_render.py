"""Quick test to verify molecular structure image rendering with CSV data."""
from __future__ import annotations

import sys
from pathlib import Path
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import Draw
from PIL import Image
import io
import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

output_dir = Path("output")
output_dir.mkdir(exist_ok=True)

# ====== 步骤1: 读取CSV数据 ======
logger.info("步骤1: 读取CSV数据")
csv_file = Path("MRGPRX2_activity_data_deduplicated.csv")
df = pd.read_csv(csv_file)
logger.info(f"读取到 {len(df)} 条分子数据")

# 过滤掉无效活性值，只保留数值型
df["activity"] = pd.to_numeric(df["IC50_nM"], errors="coerce")
df_valid = df.dropna(subset=["activity"])
df_valid = df_valid[df_valid["activity"] > 0]
logger.info(f"有效活性分子: {len(df_valid)} 条")

# 只取前20个分子用于快速测试
max_test = 20
df_test = df_valid.head(max_test)
logger.info(f"测试使用: {len(df_test)} 条分子")

# ====== 步骤2: 解析SMILES ======
logger.info("步骤2: 解析SMILES")
molecules = []
for _, row in df_test.iterrows():
    smi = row["SMILES"]
    mol = Chem.MolFromSmiles(smi)
    if mol is not None:
        molecules.append({
            "mol": mol,
            "name": row["Compound_ID"],
            "smiles": smi,
            "activity": row["activity"]
        })
    else:
        logger.warning(f"无效SMILES: {smi[:50]}...")

logger.info(f"成功解析: {len(molecules)} 个分子")

if len(molecules) < 2:
    logger.error("至少需要2个有效分子")
    sys.exit(1)

# ====== 步骤3: 渲染单个分子 ======
logger.info("步骤3: 渲染单个分子")
for i, mol_data in enumerate(molecules[:5]):
    mol = mol_data["mol"]
    AllChem.Compute2DCoords(mol)
    d = rdMolDraw2D.MolDraw2DCairo(400, 300)
    d.drawOptions().addStereoAnnotation = True
    d.DrawMolecule(mol)
    d.FinishDrawing()
    png_bytes = d.GetDrawingText()
    img = Image.open(io.BytesIO(png_bytes))
    img.save(output_dir / f"mol_{i}_{mol_data['name']}.png")
    logger.info(f"  已保存 mol_{i}_{mol_data['name']}.png ({len(png_bytes)} bytes, IC50={mol_data['activity']:.1f} nM)")

# ====== 步骤4: 渲染MCS ======
logger.info("步骤4: 查找并渲染MCS")
from rdkit.Chem import rdFMCS

# 取前6个分子找MCS
test_mols = [m["mol"] for m in molecules[:6]]
if len(test_mols) >= 2:
    mcs_result = rdFMCS.FindMCS(
        test_mols,
        ringMatchesRingOnly=True,
        completeRingsOnly=True,
        timeout=30
    )
    mcs_smarts = mcs_result.smartsString
    logger.info(f"  MCS SMARTS: {mcs_smarts}")
    logger.info(f"  MCS 原子数: {mcs_result.numAtoms}")

    # 将SMARTS转换为普通分子用于渲染
    mcs_mol = Chem.MolFromSmarts(mcs_smarts)
    if mcs_mol is not None:
        try:
            # 转换为SMILES再转回普通分子（修复渲染问题）
            mcs_smiles = Chem.MolToSmiles(mcs_mol)
            proper_mol = Chem.MolFromSmiles(mcs_smiles)
            if proper_mol is not None:
                AllChem.Compute2DCoords(proper_mol)
                d = rdMolDraw2D.MolDraw2DCairo(500, 400)
                d.drawOptions().addStereoAnnotation = True
                d.DrawMolecule(proper_mol)
                d.FinishDrawing()
                png_bytes = d.GetDrawingText()
                img = Image.open(io.BytesIO(png_bytes))
                img.save(output_dir / "mcs_structure.png")
                logger.info(f"  已保存 mcs_structure.png (SMILES: {mcs_smiles})")
        except Exception as e:
            logger.error(f"  MCS渲染失败: {e}")

# ====== 步骤5: 测试取代基提取（修复的核心Bug） ======
logger.info("步骤5: 测试取代基提取")
if len(molecules) >= 2:
    # 用第一个分子测试取代基提取
    mol = molecules[0]["mol"]
    # 简单测试：提取第一个非氢原子作为取代基
    test_atom_idx = 0
    sub_atoms = {test_atom_idx}
    
    # 新方法：MolFragmentToSmiles
    try:
        sub_smiles = Chem.MolFragmentToSmiles(mol, atomsToUse=list(sub_atoms))
        sub_mol = Chem.MolFromSmiles(sub_smiles)
        if sub_mol is not None:
            AllChem.Compute2DCoords(sub_mol)
            d = rdMolDraw2D.MolDraw2DCairo(200, 150)
            d.DrawMolecule(sub_mol)
            d.FinishDrawing()
            img = Image.open(io.BytesIO(d.GetDrawingText()))
            img.save(output_dir / "test_substituent.png")
            logger.info(f"  已保存 test_substituent.png (SMILES: {sub_smiles})")
    except Exception as e:
        logger.error(f"  取代基提取失败: {e}")

logger.info("\n=== 所有测试完成！===")
logger.info(f"输出目录: {output_dir.absolute()}")
logger.info("请检查 output/ 目录中的图片文件")

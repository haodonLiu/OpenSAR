"""ScaffoldClusterer 和 Murcko 骨架提取的单元测试."""

from __future__ import annotations

import pytest
from rdkit import Chem

from src.clustering.scaffold import (
    ScaffoldClusterer,
    ScaffoldClusteringError,
    ScaffoldClusterResult,
    get_murcko_scaffold,
    get_murcko_scaffold_smiles,
)


# ===== Fixtures =====


@pytest.fixture
def aspirin() -> Chem.Mol:
    """阿司匹林: CC(=O)Oc1ccccc1C(=O)O -> 骨架: c1ccccc1 (苯基)"""
    return Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")


@pytest.fixture
def ibuprofen() -> Chem.Mol:
    """布洛芬: CC(C)Cc1ccc(cc1)C(C)C(=O)O -> 骨架: c1ccccc1"""
    return Chem.MolFromSmiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O")


@pytest.fixture
def paracetamol() -> Chem.Mol:
    """对乙酰氨基酚: CC(=O)Nc1ccc(O)cc1 -> 骨架: c1ccccc1"""
    return Chem.MolFromSmiles("CC(=O)Nc1ccc(O)cc1")


@pytest.fixture
def caffeine() -> Chem.Mol:
    """咖啡因: 含双环 (嘌呤) 骨架"""
    return Chem.MolFromSmiles("Cn1cnc2c1c(=O)n(c(=O)n2C)C")


@pytest.fixture
def benzene_derivatives_molecules() -> list[dict]:
    """苯衍生物分子列表 (相同骨架)."""
    smiles_list = [
        "CC(=O)Oc1ccccc1C(=O)O",   # aspirin
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # ibuprofen
        "CC(=O)Nc1ccc(O)cc1",       # paracetamol
        "COc1ccc(CC(=O)O)cc1",      # 4-methoxyphenylacetic acid
    ]
    molecules = []
    for i, smi in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            molecules.append({
                "mol": mol,
                "smiles": smi,
                "name": f"benzene_der_{i+1}",
                "activity": float((i + 1) * 100),
            })
    return molecules


@pytest.fixture
def mixed_scaffold_molecules() -> list[dict]:
    """混合骨架分子列表."""
    data = [
        ("CC(=O)Oc1ccccc1C(=O)O", "asp", 10.0),     # 苯基
        ("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "ibu", 50.0),  # 苯基
        ("Cn1cnc2c1c(=O)n(c(=O)n2C)C", "caff", 500.0),  # 嘌呤
        ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "theo", 200.0),  # 嘌呤
        ("c1ccccc1", "benzene", 1000.0),  # 苯基本身
        ("CCO", "ethanol", 9999.0),  # 无环 (无骨架)
    ]
    molecules = []
    for smi, name, act in data:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            molecules.append({
                "mol": mol,
                "smiles": smi,
                "name": name,
                "activity": act,
            })
    return molecules


# ===== get_murcko_scaffold 测试 =====


class TestGetMurckoScaffold:
    """get_murcko_scaffold 函数测试."""

    def test_benzene_derivative(self, aspirin):
        """苯衍生物应返回苯基骨架."""
        scaffold = get_murcko_scaffold(aspirin)
        assert scaffold is not None
        smiles = Chem.MolToSmiles(scaffold)
        assert smiles == "c1ccccc1"

    def test_none_input(self):
        """None 输入应返回 None."""
        assert get_murcko_scaffold(None) is None

    def test_returns_mol_object(self, aspirin):
        """应返回有效的 Mol 对象."""
        scaffold = get_murcko_scaffold(aspirin)
        assert scaffold is not None
        assert isinstance(scaffold, Chem.Mol)
        assert scaffold.GetNumAtoms() > 0

    def test_purine_scaffold(self, caffeine):
        """嘌呤类化合物应返回双环骨架."""
        scaffold = get_murcko_scaffold(caffeine)
        assert scaffold is not None
        # 嘌呤骨架应有至少 5 个原子 (两个稠合环)
        assert scaffold.GetNumAtoms() >= 5

    def test_no_rings_returns_none(self):
        """无环分子应无法提取骨架 (或返回极小/空)."""
        mol = Chem.MolFromSmiles("CCO")  # ethanol - 无环
        scaffold = get_murcko_scaffold(mol)
        # 无环分子的 Murcko 骨架为空或非常小
        if scaffold is not None:
            assert scaffold.GetNumAtoms() == 0 or Chem.MolToSmiles(scaffold) == ""


# ===== get_murcko_scaffold_smiles 测试 =====


class TestGetMurckoScaffoldSmiles:
    """get_murcko_scaffold_smiles 函数测试."""

    def test_returns_string(self, aspirin):
        """应返回非空字符串."""
        result = get_murcko_scaffold_smiles(aspirin)
        assert isinstance(result, str)
        assert len(result) > 0

    def test_none_input_returns_empty(self):
        """None 输入应返回空字符串."""
        assert get_murcko_scaffold_smiles(None) == ""

    def test_correct_smiles(self, paracetamol):
        """对乙酰氨基酚应返回苯基 SMILES."""
        result = get_murcko_scaffold_smiles(paracetamol)
        assert result == "c1ccccc1"


# ===== ScaffoldClusterer 初始化测试 =====


class TestScaffoldClustererInit:
    """ScaffoldClusterer 初始化测试."""

    def test_default_init(self):
        """默认参数初始化."""
        clusterer = ScaffoldClusterer()
        assert clusterer.include_chirality is False
        assert clusterer.merge_aromatic_kekulize is True
        assert clusterer.min_cluster_size == 1

    def test_custom_params(self):
        """自定义参数初始化."""
        clusterer = ScaffoldClusterer(
            include_chirality=True,
            merge_aromatic_kekulize=False,
            min_cluster_size=3,
        )
        assert clusterer.include_chirality is True
        assert clusterer.merge_aromatic_kekulize is False
        assert clusterer.min_cluster_size == 3


# ===== ScaffoldClusterer.cluster 测试 =====


class TestScaffoldClustererCluster:
    """ScaffoldClusterer.cluster 方法测试."""

    def test_same_scaffold_groups_together(self, benzene_derivatives_molecules):
        """相同骨架的分子应归为一簇."""
        clusterer = ScaffoldClusterer()
        results = clusterer.cluster(benzene_derivatives_molecules)

        # 所有苯衍生物应归为一个骨架簇
        assert len(results) >= 1
        # 总分子数应为输入数量
        total_in_results = sum(r.num_molecules for r in results)
        assert total_in_results == len(benzene_derivatives_molecules)

    def test_different_scaffolds_separated(self, mixed_scaffold_molecules):
        """不同骨架应分到不同簇."""
        clusterer = ScaffoldClusterer()
        results = clusterer.cluster(mixed_scaffold_molecules)

        # 至少有苯基和嘌呤两种骨架 (乙醇无环可能被过滤)
        scaffold_smiles_set = {r.scaffold_smiles for r in results}
        assert len(scaffold_smiles_set) >= 2

    def test_empty_raises_error(self):
        """空列表应抛出异常."""
        clusterer = ScaffoldClusterer()
        with pytest.raises(ScaffoldClusteringError):
            clusterer.cluster([])

    def test_result_ordering_by_size(self, mixed_scaffold_molecules):
        """结果应按分子数降序排列."""
        clusterer = ScaffoldClusterer()
        results = clusterer.cluster(mixed_scaffold_molecules)

        if len(results) >= 2:
            sizes = [r.num_molecules for r in results]
            assert sizes == sorted(sizes, reverse=True)

    def test_min_cluster_size_filter(self, benzene_derivatives_molecules):
        """min_cluster_size 应正确过滤小簇."""
        clusterer = ScaffoldClusterer(min_cluster_size=10)
        results = clusterer.cluster(benzene_derivatives_molecules)

        # 苯衍生物只有 4 个，min_cluster_size=10 应该返回空列表
        assert len(results) == 0

    def test_result_attributes(self, benzene_derivatives_molecules):
        """结果对象应包含正确的属性."""
        clusterer = ScaffoldClusterer()
        results = clusterer.cluster(benzene_derivatives_molecules)

        assert len(results) >= 1
        result = results[0]
        assert isinstance(result, ScaffoldClusterResult)
        assert isinstance(result.scaffold_smiles, str)
        assert result.num_molecules > 0
        assert isinstance(result.molecules, list)
        assert len(result.activities) > 0
        assert result.mean_activity > 0

    def test_activity_statistics(self, mixed_scaffold_molecules):
        """活性统计计算应正确."""
        clusterer = ScaffoldClusterer()
        results = clusterer.cluster(mixed_scaffold_molecules)

        for result in results:
            if result.activities:
                assert result.min_activity == min(result.activities)
                assert result.max_activity == max(result.activities)
                expected_mean = sum(result.activities) / len(result.activities)
                assert abs(result.mean_activity - expected_mean) < 0.01

    def test_molecules_without_activity(self):
        """没有活性的分子也能正常聚类."""
        mols_data = [
            ("CC(=O)Oc1ccccc1C(=O)O", "a"),
            ("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "b"),
        ]
        molecules = []
        for smi, name in mols_data:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                molecules.append({"mol": mol, "smiles": smi, "name": name})

        clusterer = ScaffoldClusterer()
        results = clusterer.cluster(molecules)
        assert len(results) >= 1
        # 没有活性时统计值应为 0
        assert results[0].mean_activity == 0.0

    def test_skip_invalid_molecules(self):
        """跳过无效分子 (None mol)."""
        valid_mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        molecules = [
            {"mol": None, "name": "invalid"},
            {"mol": valid_mol, "smiles": "CC(=O)Oc1ccccc1C(=O)O", "name": "valid"},
        ]

        clusterer = ScaffoldClusterer()
        results = clusterer.cluster(molecules)
        # 只有有效分子被聚类
        assert len(results) >= 1
        assert sum(r.num_molecules for r in results) == 1


# ===== get_scaffold_diversity 测试 =====


class TestScaffoldDiversity:
    """骨架多样性分析测试."""

    def test_diversity_stats(self, mixed_scaffold_molecules):
        """多样性统计应包含所有必需字段."""
        clusterer = ScaffoldClusterer()
        stats = clusterer.get_scaffold_diversity(mixed_scaffold_molecules)

        expected_keys = {
            "total_molecules",
            "unique_scaffolds",
            "singleton_scaffolds",
            "largest_cluster_size",
            "scaffold_diversity_index",
        }
        assert set(stats.keys()) == expected_keys
        assert stats["total_molecules"] == len(mixed_scaffold_molecules)
        assert 0 <= stats["scaffold_diversity_index"] <= 1

    def test_uniform_scaffold_low_diversity(self, benzene_derivatives_molecules):
        """所有分子同骨架 -> 低多样性指数."""
        clusterer = ScaffoldClusterer()
        stats = clusterer.get_scaffold_diversity(benzene_derivatives_molecules)

        # 全部是苯基骨架，unique scaffolds 应为 1
        assert stats["unique_scaffolds"] == 1
        assert stats["scaffold_diversity_index"] == 1.0 / len(benzene_derivatives_molecules)


# ===== ScaffoldClusterResult repr 测试 =====


class TestScaffoldClusterResultRepr:
    """ScaffoldClusterResult __repr__ 测试."""

    def test_repr_format(self, benzene_derivatives_molecules):
        """repr 应包含关键信息."""
        clusterer = ScaffoldClusterer()
        results = clusterer.cluster(benzene_derivatives_molecules)
        if results:
            r = repr(results[0])
            assert "ScaffoldClusterResult" in r
            assert "n=" in r

"""Tests for MolecularClusterer."""

from __future__ import annotations

import pytest
import numpy as np
from rdkit import Chem

from src.clustering.cluster import MolecularClusterer, ClusteringError, ClusterResult
from src.clustering.fingerprinter import MolecularFingerprinter


class TestMolecularClusterer:
    """Tests for MolecularClusterer class."""

    @pytest.fixture
    def benzene(self) -> Chem.Mol:
        """Benzene molecule fixture."""
        return Chem.MolFromSmiles("c1ccccc1")

    @pytest.fixture
    def toluene(self) -> Chem.Mol:
        """Toluene molecule fixture."""
        return Chem.MolFromSmiles("Cc1ccccc1")

    @pytest.fixture
    def naphthalene(self) -> Chem.Mol:
        """Naphthalene molecule fixture."""
        return Chem.MolFromSmiles("c1ccc2ccccc2c1")

    @pytest.fixture
    def anthracene(self) -> Chem.Mol:
        """Anthracene molecule fixture."""
        return Chem.MolFromSmiles("c1ccc2ccccc2c2ccc3ccccc3c2c1")

    @pytest.fixture
    def molecules_benzene_toluene(
        self, benzene: Chem.Mol, toluene: Chem.Mol
    ) -> list[dict]:
        """Molecule list with benzene and toluene."""
        return [
            {"mol": benzene, "name": "benzene", "smiles": "c1ccccc1"},
            {"mol": toluene, "name": "toluene", "smiles": "Cc1ccccc1"},
        ]

    @pytest.fixture
    def molecules_mixed(
        self,
        benzene: Chem.Mol,
        toluene: Chem.Mol,
        naphthalene: Chem.Mol,
        anthracene: Chem.Mol,
    ) -> list[dict]:
        """Mixed molecule list."""
        return [
            {"mol": benzene, "name": "benzene", "smiles": "c1ccccc1"},
            {"mol": toluene, "name": "toluene", "smiles": "Cc1ccccc1"},
            {"mol": naphthalene, "name": "naphthalene", "smiles": "c1ccc2ccccc2c1"},
            {
                "mol": anthracene,
                "name": "anthracene",
                "smiles": "c1ccc2ccccc2c2ccc3ccccc3c2c1",
            },
        ]

    @pytest.fixture
    def clusterer(self) -> MolecularClusterer:
        """Clusterer fixture."""
        return MolecularClusterer(threshold=0.7, method="tanimoto")

    def test_cluster_two_similar_molecules(
        self, clusterer: MolecularClusterer, molecules_benzene_toluene: list[dict]
    ) -> None:
        """Test clustering two similar molecules (using lower threshold for these dissimilar mols)."""
        low_threshold_clusterer = MolecularClusterer(threshold=0.2, method="tanimoto")
        results, sim_matrix = low_threshold_clusterer.cluster(molecules_benzene_toluene)

        assert len(results) == 1
        assert results[0].size == 2
        assert results[0].avg_similarity > 0.2
        assert sim_matrix.shape == (2, 2)

    def test_cluster_dissimilar_molecules(
        self, clusterer: MolecularClusterer, molecules_mixed: list[dict]
    ) -> None:
        """Test clustering dissimilar molecules."""
        results, sim_matrix = clusterer.cluster(molecules_mixed)

        assert len(results) >= 1
        for result in results:
            assert result.size >= 1

    def test_cluster_empty_list(self, clusterer: MolecularClusterer) -> None:
        """Test clustering empty list."""
        results, sim_matrix = clusterer.cluster([])
        assert len(results) == 0
        assert sim_matrix.size == 0

    def test_cluster_results_have_valid_structure(
        self, clusterer: MolecularClusterer, molecules_mixed: list[dict]
    ) -> None:
        """Test cluster results have valid structure."""
        results, sim_matrix = clusterer.cluster(molecules_mixed)

        assert sim_matrix.shape[0] == len(molecules_mixed)
        for result in results:
            assert isinstance(result, ClusterResult)
            assert result.cluster_id >= 0
            assert result.size == len(result.molecules)
            assert result.representative_idx >= 0
            assert 0.0 <= result.avg_similarity <= 1.0

    def test_cluster_with_different_thresholds(
        self, molecules_mixed: list[dict]
    ) -> None:
        """Test clustering with different thresholds."""
        high_thresholder = MolecularClusterer(threshold=0.9)
        low_thresholder = MolecularClusterer(threshold=0.4)

        results_high, _ = high_thresholder.cluster(molecules_mixed)
        results_low, _ = low_thresholder.cluster(molecules_mixed)

        assert len(results_high) >= len(results_low)

    def test_butina_clustering(self, molecules_mixed: list[dict]) -> None:
        """Test Butina clustering method."""
        clusterer = MolecularClusterer(method="butina")
        results, sim_matrix = clusterer.cluster(molecules_mixed)

        assert len(results) >= 1
        total_mols = sum(r.size for r in results)
        assert total_mols == len(molecules_mixed)

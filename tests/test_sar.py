"""Tests for SARAnalyzer."""

from __future__ import annotations

import pytest
from rdkit import Chem

from src.sar.analyzer import SARAnalyzer, SARError, SARResult


class TestSARAnalyzer:
    """Tests for SARAnalyzer class."""

    @pytest.fixture
    def benzene(self) -> dict:
        """Benzene molecule fixture."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        return {"mol": mol, "smiles": "c1ccccc1", "name": "benzene", "activity": 5.0}

    @pytest.fixture
    def toluene(self) -> dict:
        """Toluene molecule fixture."""
        mol = Chem.MolFromSmiles("Cc1ccccc1")
        return {"mol": mol, "smiles": "Cc1ccccc1", "name": "toluene", "activity": 7.5}

    @pytest.fixture
    def xylene(self) -> dict:
        """Xylene molecule fixture."""
        mol = Chem.MolFromSmiles("Cc1ccc(C)cc1")
        return {"mol": mol, "smiles": "Cc1ccc(C)cc1", "name": "xylene", "activity": 8.0}

    @pytest.fixture
    def analyzer(self) -> SARAnalyzer:
        """SARAnalyzer fixture."""
        return SARAnalyzer()

    def test_analyzer_initialization(self, analyzer: SARAnalyzer) -> None:
        """Test analyzer initializes correctly."""
        assert analyzer.activity_threshold is None
        assert analyzer.use_mcs is True

    def test_analyze_cluster_with_activities(
        self, analyzer: SARAnalyzer, benzene: dict, toluene: dict, xylene: dict
    ) -> None:
        """Test cluster analysis with activity data."""
        molecules = [benzene, toluene, xylene]
        result = analyzer.analyze_cluster(0, molecules)

        assert isinstance(result, SARResult)
        assert result.cluster_id == 0
        assert result.num_compounds == 3
        assert len(result.activities) == 3
        assert result.mean_activity == pytest.approx(6.833, rel=0.01)
        assert result.max_activity == 8.0
        assert result.min_activity == 5.0

    def test_analyze_cluster_without_activities(self, analyzer: SARAnalyzer) -> None:
        """Test cluster analysis without activity data."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        molecules = [{"mol": mol, "smiles": "c1ccccc1", "name": "benzene"}]
        result = analyzer.analyze_cluster(0, molecules)

        assert isinstance(result, SARResult)
        assert result.num_compounds == 1
        assert len(result.activities) == 0
        assert result.mean_activity == 0.0

    def test_get_activity_stats(
        self, analyzer: SARAnalyzer, benzene: dict, toluene: dict, xylene: dict
    ) -> None:
        """Test overall activity statistics."""
        molecules = [benzene, toluene, xylene]
        result = analyzer.analyze_cluster(0, molecules)

        stats = analyzer.get_activity_stats({0: result})

        assert stats["total_compounds"] == 3
        assert stats["mean"] == pytest.approx(6.833, rel=0.01)
        assert stats["max"] == 8.0
        assert stats["min"] == 5.0

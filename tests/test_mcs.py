"""Tests for MCSFinder."""

from __future__ import annotations

import pytest
from rdkit import Chem

from src.mcs.finder import MCSFinder, MCSError, MCSResult


class TestMCSFinder:
    """Tests for MCSFinder class."""

    @pytest.fixture
    def benzene(self) -> dict:
        """Benzene molecule fixture."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        return {"mol": mol, "smiles": "c1ccccc1", "name": "benzene"}

    @pytest.fixture
    def toluene(self) -> dict:
        """Toluene molecule fixture."""
        mol = Chem.MolFromSmiles("Cc1ccccc1")
        return {"mol": mol, "smiles": "Cc1ccccc1", "name": "toluene"}

    @pytest.fixture
    def naphthalene(self) -> dict:
        """Naphthalene molecule fixture."""
        mol = Chem.MolFromSmiles("c1ccc2ccccc2c1")
        return {"mol": mol, "smiles": "c1ccc2ccccc2c1", "name": "naphthalene"}

    @pytest.fixture
    def anthracene(self) -> dict:
        """Anthracene molecule fixture."""
        mol = Chem.MolFromSmiles("c1ccc2ccccc2c2ccc3ccccc3c2c1")
        return {
            "mol": mol,
            "smiles": "c1ccc2ccccc2c2ccc3ccccc3c2c1",
            "name": "anthracene",
        }

    @pytest.fixture
    def mcs_finder(self) -> MCSFinder:
        """MCSFinder fixture."""
        return MCSFinder(timeout=10)

    def test_find_mcs_two_molecules(
        self, mcs_finder: MCSFinder, benzene: dict, toluene: dict
    ) -> None:
        """Test MCS finding between two molecules."""
        result = mcs_finder.find_mcs([benzene, toluene])

        assert result is not None
        assert isinstance(result, MCSResult)
        assert result.num_atoms > 0
        assert result.smiles_a == benzene["smiles"]
        assert result.smiles_b == toluene["smiles"]

    def test_find_mcs_identical_molecules(
        self, mcs_finder: MCSFinder, benzene: dict
    ) -> None:
        """Test MCS finding with identical molecules."""
        result = mcs_finder.find_mcs([benzene, benzene])

        assert result is not None
        assert result.num_atoms == benzene["mol"].GetNumAtoms()

    def test_find_mcs_insufficient_molecules(
        self, mcs_finder: MCSFinder, benzene: dict
    ) -> None:
        """Test MCS raises error with fewer than 2 molecules."""
        with pytest.raises(MCSError, match="At least 2 molecules"):
            mcs_finder.find_mcs([benzene])

    def test_find_mcs_similar_molecules(
        self, mcs_finder: MCSFinder, benzene: dict, toluene: dict
    ) -> None:
        """Test MCS finds meaningful MCS for similar molecules."""
        result = mcs_finder.find_mcs([benzene, toluene])

        assert result is not None
        assert result.num_atoms >= 6

    def test_find_mcs_dissimilar_molecules(
        self, mcs_finder: MCSFinder, benzene: dict, naphthalene: dict
    ) -> None:
        """Test MCS for dissimilar molecules."""
        result = mcs_finder.find_mcs([benzene, naphthalene])

        assert result is not None
        assert result.num_atoms >= 6

    def test_mcs_result_attributes(
        self, mcs_finder: MCSFinder, benzene: dict, toluene: dict
    ) -> None:
        """Test MCS result has all required attributes."""
        result = mcs_finder.find_mcs([benzene, toluene])

        assert result is not None
        assert hasattr(result, "mcs_mol")
        assert hasattr(result, "smiles")
        assert hasattr(result, "num_atoms")
        assert hasattr(result, "num_bonds")
        assert hasattr(result, "smiles_a")
        assert hasattr(result, "smiles_b")
        assert hasattr(result, "bond_matches")
        assert hasattr(result, "atom_matches")
        assert hasattr(result, "score")

"""Tests for MolecularFingerprinter."""

from __future__ import annotations

import pytest
import numpy as np
from rdkit import Chem

from src.clustering.fingerprinter import MolecularFingerprinter, FingerprintResult


class TestMolecularFingerprinter:
    """Tests for MolecularFingerprinter class."""

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
    def fingerprinter(self) -> MolecularFingerprinter:
        """Fingerprinter fixture."""
        return MolecularFingerprinter(fp_type="Morgan", radius=2, n_bits=2048)

    def test_fingerprint_generation_morgan(
        self, fingerprinter: MolecularFingerprinter, benzene: Chem.Mol
    ) -> None:
        """Test Morgan fingerprint generation."""
        result = fingerprinter.fingerprint(benzene)

        assert isinstance(result, FingerprintResult)
        assert result.fingerprint.shape == (2048,)
        assert result.fp_type == "Morgan"
        assert result.n_bits == 2048
        assert result.fingerprint.dtype == np.uint8

    def test_fingerprint_from_smiles(
        self, fingerprinter: MolecularFingerprinter
    ) -> None:
        """Test fingerprint generation from SMILES string."""
        result = fingerprinter.fingerprint("c1ccccc1")

        assert isinstance(result, FingerprintResult)
        assert result.fingerprint.shape == (2048,)

    def test_fingerprint_invalid_smiles(
        self, fingerprinter: MolecularFingerprinter
    ) -> None:
        """Test fingerprint raises error for invalid SMILES."""
        with pytest.raises(ValueError, match="Invalid SMILES"):
            fingerprinter.fingerprint("invalid")

    def test_fingerprint_none_mol(self, fingerprinter: MolecularFingerprinter) -> None:
        """Test fingerprint raises error for None molecule."""
        with pytest.raises(ValueError, match="Molecule cannot be None"):
            fingerprinter.fingerprint(None)

    def test_similarity_benzene_itself(
        self, fingerprinter: MolecularFingerprinter, benzene: Chem.Mol
    ) -> None:
        """Test similarity of molecule to itself is 1.0."""
        sim = fingerprinter.similarity(benzene, benzene)
        assert sim == 1.0

    def test_similarity_benzene_toluene(
        self,
        fingerprinter: MolecularFingerprinter,
        benzene: Chem.Mol,
        toluene: Chem.Mol,
    ) -> None:
        """Test similarity between benzene and toluene."""
        sim = fingerprinter.similarity(benzene, toluene)
        assert 0.0 < sim < 1.0

    def test_similarity_benzene_naphthalene(
        self,
        fingerprinter: MolecularFingerprinter,
        benzene: Chem.Mol,
        naphthalene: Chem.Mol,
    ) -> None:
        """Test similarity between benzene and naphthalene."""
        sim = fingerprinter.similarity(benzene, naphthalene)
        assert 0.0 < sim < 1.0

    def test_similarity_from_smiles(
        self, fingerprinter: MolecularFingerprinter
    ) -> None:
        """Test similarity calculation from SMILES strings."""
        sim = fingerprinter.similarity("c1ccccc1", "Cc1ccccc1")
        assert 0.0 < sim < 1.0

    def test_pairwise_similarity_matrix(
        self,
        fingerprinter: MolecularFingerprinter,
        benzene: Chem.Mol,
        toluene: Chem.Mol,
        naphthalene: Chem.Mol,
    ) -> None:
        """Test pairwise similarity matrix computation."""
        mols = [benzene, toluene, naphthalene]
        sim_matrix = fingerprinter.pairwise_similarity(mols)

        assert sim_matrix.shape == (3, 3)
        assert np.allclose(np.diag(sim_matrix), 1.0)
        assert sim_matrix[0, 1] == sim_matrix[1, 0]
        assert sim_matrix[0, 2] == sim_matrix[2, 0]
        assert np.all(sim_matrix >= 0.0)
        assert np.all(sim_matrix <= 1.0)

    def test_different_fingerprint_types(self) -> None:
        """Test different fingerprint types."""
        mol = Chem.MolFromSmiles("c1ccccc1")

        for fp_type in ["Morgan", "MACCS", "RDKit"]:
            fp = MolecularFingerprinter(fp_type=fp_type)
            result = fp.fingerprint(mol)
            assert result.fp_type == fp_type

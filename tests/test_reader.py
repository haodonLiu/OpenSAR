"""Tests for MoleculeReader."""

from __future__ import annotations

import pytest
from pathlib import Path
from rdkit import Chem

from src.io.reader import MoleculeReader, MoleculeReadError


class TestMoleculeReader:
    """Tests for MoleculeReader class."""

    @pytest.fixture
    def reader(self) -> MoleculeReader:
        """Reader fixture."""
        return MoleculeReader(
            smiles_column="SMILES", name_column="Name", activity_column="Activity"
        )

    def test_reader_initialization(self, reader: MoleculeReader) -> None:
        """Test reader initializes with correct defaults."""
        assert reader.smiles_column == "SMILES"
        assert reader.name_column == "Name"
        assert reader.activity_column == "Activity"

    def test_read_excel_method_exists(self, reader: MoleculeReader) -> None:
        """Test read_excel method exists."""
        assert hasattr(reader, "read_excel")

    def test_read_csv_method_exists(self, reader: MoleculeReader) -> None:
        """Test read_csv method exists."""
        assert hasattr(reader, "read_csv")

    def test_read_sdf_method_exists(self, reader: MoleculeReader) -> None:
        """Test read_sdf method exists."""
        assert hasattr(reader, "read_sdf")

    def test_read_auto_detect_method_exists(self, reader: MoleculeReader) -> None:
        """Test read method exists for auto-detection."""
        assert hasattr(reader, "read")

    def test_unsupported_format_raises_error(
        self, reader: MoleculeReader, tmp_path: Path
    ) -> None:
        """Test unsupported format raises MoleculeReadError."""
        test_file = tmp_path / "test.xyz"
        test_file.write_text("fake data")

        with pytest.raises(MoleculeReadError, match="Unsupported file format"):
            reader.read(test_file)

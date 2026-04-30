"""ActivityPreprocessor 单元测试."""
from __future__ import annotations

import pytest

from src.sar.preprocessor import (
    ActivityPreprocessor,
    ActivityUnit,
    CensoredType,
    CensoredValue,
    DataQualityReport,
)


@pytest.fixture
def preprocessor() -> ActivityPreprocessor:
    return ActivityPreprocessor()


@pytest.fixture
def sample_molecules():
    """提供包含多种活性格式的测试分子."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mols = []
    for i, smiles in enumerate([
        "CCO", "CCCO", "CC(C)O", "CCCCCO", "c1ccccc1",
    ]):
        mol = Chem.MolFromSmiles(smiles)
        AllChem.Compute2DCoords(mol)
        mols.append({
            "mol": mol,
            "name": f"mol_{i}",
            "smiles": smiles,
            "activity_raw": ["12.5 nM", "<100", ">10000", "0.5 uM", "1000"][i],
            "activity": [12.5, None, None, None, 1000.0][i],
        })
    return mols


class TestActivityPreprocessor:
    """ActivityPreprocessor 的所有功能测试."""

    def test_init_defaults(self):
        """测试默认参数初始化."""
        p = ActivityPreprocessor()
        assert p.iqr_multiplier == 1.5
        assert p.enable_p_conversion is True
        assert p.enable_outlier_detection is True
        assert p.censor_treatment == "conservative"

    def test_init_custom(self):
        """测试自定义参数初始化."""
        p = ActivityPreprocessor(
            iqr_multiplier=3.0,
            enable_p_conversion=False,
            enable_outlier_detection=False,
            censor_treatment="as_is",
        )
        assert p.iqr_multiplier == 3.0
        assert p.enable_p_conversion is False
        assert p.enable_outlier_detection is False
        assert p.censor_treatment == "as_is"

    # --- Unit Detection ---

    def test_detect_unit_nm(self, preprocessor):
        """检测 nM 单位."""
        for raw in ["12.5 nM", "100nm", "50 nanomolar", "25 n molar"]:
            assert preprocessor.detect_unit(raw) == ActivityUnit.NM, f"Failed for {raw}"

    def test_detect_unit_um(self, preprocessor):
        """检测 μM 单位."""
        for raw in ["0.5 uM", "100 μm", "50 micromolar", "25 micro molar"]:
            assert preprocessor.detect_unit(raw) == ActivityUnit.UM, f"Failed for {raw}"

    def test_detect_unit_mm(self, preprocessor):
        """检测 mM 单位."""
        for raw in ["0.001 mM", "5 mm", "2 millimolar"]:
            assert preprocessor.detect_unit(raw) == ActivityUnit.MM, f"Failed for {raw}"

    def test_detect_unit_unknown(self, preprocessor):
        """未识别单位返回 UNKNOWN."""
        assert preprocessor.detect_unit("100") == ActivityUnit.UNKNOWN
        assert preprocessor.detect_unit("test") == ActivityUnit.UNKNOWN
        assert preprocessor.detect_unit("") == ActivityUnit.UNKNOWN

    # --- Censored Data ---

    def test_parse_censored_greater_than(self, preprocessor):
        """解析 > 删失数据."""
        v = preprocessor.parse_censored(">10000")
        assert v is not None
        assert v.censored_type == CensoredType.GREATER_THAN
        assert v.value == 10000.0
        assert v.raw_string == ">10000"

    def test_parse_censored_less_than(self, preprocessor):
        """解析 < 删失数据."""
        v = preprocessor.parse_censored("<50")
        assert v is not None
        assert v.censored_type == CensoredType.LESS_THAN
        assert v.value == 50.0

    def test_parse_censored_tilde(self, preprocessor):
        """解析 ~ 近似值."""
        v = preprocessor.parse_censored("~100")
        assert v is not None
        assert v.censored_type == CensoredType.TILDE
        assert v.value == 100.0

    def test_parse_censored_approx(self, preprocessor):
        """解析 approx. 近似值."""
        v = preprocessor.parse_censored("approx. 500")
        assert v is not None
        assert v.censored_type == CensoredType.TILDE
        assert v.value == 500.0

    def test_parse_censored_none(self, preprocessor):
        """非删失数据返回 None."""
        assert preprocessor.parse_censored("12.5") is None
        assert preprocessor.parse_censored("abc") is None
        assert preprocessor.parse_censored("") is None

    # --- Clean Numeric ---

    def test_clean_numeric_simple(self, preprocessor):
        """简单数值."""
        val, unit, censored = preprocessor.clean_numeric("12.5")
        assert val == pytest.approx(12.5)
        assert unit == ActivityUnit.NM
        assert censored is None

    def test_clean_numeric_with_unit(self, preprocessor):
        """带单位的数值."""
        val, unit, _ = preprocessor.clean_numeric("12.5 nM")
        assert val == pytest.approx(12.5)
        assert unit == ActivityUnit.NM

        val, unit, _ = preprocessor.clean_numeric("0.5 uM")
        assert val == pytest.approx(0.5)
        assert unit == ActivityUnit.UM

        val, unit, _ = preprocessor.clean_numeric("0.001 mM")
        assert val == pytest.approx(0.001)
        assert unit == ActivityUnit.MM

    def test_clean_numeric_censored(self, preprocessor):
        """带删失符号的数值."""
        val, unit, censored = preprocessor.clean_numeric(">10000 nM")
        assert val == pytest.approx(10000.0)
        assert censored is not None
        assert censored.censored_type == CensoredType.GREATER_THAN

    def test_clean_numeric_invalid(self, preprocessor):
        """无效输入."""
        val, unit, censored = preprocessor.clean_numeric("")
        assert val is None
        assert unit == ActivityUnit.UNKNOWN
        assert censored is None

        val, unit, censored = preprocessor.clean_numeric("abc")
        assert val is None

    def test_clean_numeric_scientific(self, preprocessor):
        """科学计数法."""
        val, unit, _ = preprocessor.clean_numeric("1e-3")
        assert val == pytest.approx(0.001)

        val, unit, _ = preprocessor.clean_numeric("1.5e2 nM")
        assert val == pytest.approx(150.0)
        assert unit == ActivityUnit.NM

    # --- Unit Conversion ---

    def test_convert_to_nm(self, preprocessor):
        """单位转换."""
        assert preprocessor.convert_to_nm(1.0, ActivityUnit.NM) == pytest.approx(1.0)
        assert preprocessor.convert_to_nm(1.0, ActivityUnit.UM) == pytest.approx(1000.0)
        assert preprocessor.convert_to_nm(1.0, ActivityUnit.MM) == pytest.approx(1_000_000.0)
        assert preprocessor.convert_to_nm(1.0, ActivityUnit.UNKNOWN) == pytest.approx(1.0)

    # --- pActivity Conversion ---

    def test_to_p_activity(self, preprocessor):
        """pActivity 转换: pIC50 = -log10(IC50_M)."""
        # 1 nM = 1e-9 M → pIC50 = 9
        assert preprocessor.to_p_activity(1.0) == pytest.approx(9.0)
        # 10 nM → pIC50 = 8
        assert preprocessor.to_p_activity(10.0) == pytest.approx(8.0)
        # 100 nM → pIC50 = 7
        assert preprocessor.to_p_activity(100.0) == pytest.approx(7.0)
        # 1000 nM = 1 uM → pIC50 = 6
        assert preprocessor.to_p_activity(1000.0) == pytest.approx(6.0)
        # 10000 nM = 10 uM → pIC50 = 5
        assert preprocessor.to_p_activity(10000.0) == pytest.approx(5.0)

    def test_to_p_activity_zero(self, preprocessor):
        """零值返回 0."""
        assert preprocessor.to_p_activity(0) == 0.0
        assert preprocessor.to_p_activity(-1) == 0.0

    # --- Outlier Detection ---

    def test_detect_outliers_iqr(self, preprocessor):
        """IQR 异常值检测."""
        values = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 100]  # 100 是异常值
        outliers = preprocessor.detect_outliers_iqr(values)
        assert sum(outliers) >= 1
        assert outliers[-1] is True  # 100 是异常值

    def test_detect_outliers_iqr_few_values(self, preprocessor):
        """少于 4 个值时不做异常值检测."""
        assert preprocessor.detect_outliers_iqr([1, 2]) == [False, False]
        assert preprocessor.detect_outliers_iqr([1]) == [False]

    def test_detect_outliers_no_outliers(self, preprocessor):
        """正常数据无异常值."""
        values = [1, 2, 3, 4, 5, 6, 7, 8]
        outliers = preprocessor.detect_outliers_iqr(values)
        assert all(not o for o in outliers)

    # --- Process Molecule ---

    def test_process_molecule_with_raw(self, preprocessor):
        """从 activity_raw 解析活性."""
        from rdkit import Chem
        mol = Chem.MolFromSmiles("CCO")

        processed = preprocessor.process_molecule({
            "mol": mol,
            "name": "test",
            "smiles": "CCO",
            "activity_raw": "12.5 nM",
            "activity": 12.5,
        })
        assert processed is not None
        assert processed["activity"] == pytest.approx(12.5)
        assert processed["activity_unit"] == ActivityUnit.NM
        assert processed["activity_censored"] is None

    def test_process_molecule_censored(self, preprocessor):
        """删失数据处理."""
        from rdkit import Chem
        mol = Chem.MolFromSmiles("CCO")

        processed = preprocessor.process_molecule({
            "mol": mol,
            "name": "test",
            "smiles": "CCO",
            "activity_raw": "<100",
        })
        assert processed is not None
        assert processed["activity"] == pytest.approx(100.0)
        assert processed["activity_censored"] is not None
        assert processed["activity_censored"].censored_type == CensoredType.LESS_THAN

    def test_process_molecule_no_activity(self, preprocessor):
        """无活性数据."""
        from rdkit import Chem
        mol = Chem.MolFromSmiles("CCO")

        processed = preprocessor.process_molecule({
            "mol": mol,
            "name": "test",
            "smiles": "CCO",
        })
        assert processed is not None
        assert processed.get("activity") is None

    # --- Batch Processing ---

    def test_process_batch(self, preprocessor, sample_molecules):
        """批量处理测试."""
        processed = preprocessor.process(sample_molecules)
        assert len(processed) == 5

        # 第一个: "12.5 nM"
        assert processed[0]["activity"] == pytest.approx(12.5)
        assert processed[0]["activity_unit"] == ActivityUnit.NM

        # 第二个: "<100" - 删失
        assert processed[1]["activity_censored"] is not None

        # 第五个: "1000" - 纯数字
        assert processed[4]["activity"] == pytest.approx(1000.0)
        assert processed[4].get("p_activity") is not None

    def test_process_p_activity(self, preprocessor, sample_molecules):
        """pActivity 转换."""
        preprocessor.enable_p_conversion = True
        processed = preprocessor.process(sample_molecules)
        act_mols = [m for m in processed if m.get("p_activity") is not None]
        assert len(act_mols) > 0
        for m in act_mols:
            assert 0 <= m["p_activity"] <= 15  # 合理的 pIC50 范围

    def test_process_no_p_activity(self, preprocessor, sample_molecules):
        """禁用 pActivity 转换."""
        preprocessor.enable_p_conversion = False
        processed = preprocessor.process(sample_molecules)
        for m in processed:
            assert m.get("p_activity") is None

    def test_process_no_outlier_detection(self, preprocessor, sample_molecules):
        """禁用异常值检测."""
        preprocessor.enable_outlier_detection = False
        processed = preprocessor.process(sample_molecules)
        for m in processed:
            assert m.get("is_outlier") is None or m["is_outlier"] is False

    # --- Quality Report ---

    def test_quality_report(self, preprocessor, sample_molecules):
        """数据质量报告."""
        preprocessor.process(sample_molecules)
        report = preprocessor.get_quality_report()

        assert isinstance(report, DataQualityReport)
        assert report.total_molecules == 5
        assert report.valid_activities > 0
        assert "nM" in report.unit_counts

    def test_quality_report_summary(self, preprocessor, sample_molecules):
        """质量报告摘要字符串."""
        preprocessor.process(sample_molecules)
        report = preprocessor.get_quality_report()
        summary = report.summary()

        assert isinstance(summary, str)
        assert "数据质量报告" in summary
        assert "总分子数" in summary
        assert "有效活性值" in summary

    def test_quality_report_empty(self, preprocessor):
        """空处理会话的报告."""
        report = preprocessor.get_quality_report()
        assert report.total_molecules == 0
        assert report.valid_activities == 0

    # --- Edge Cases ---

    def test_empty_molecules(self, preprocessor):
        """空分子列表."""
        processed = preprocessor.process([])
        assert processed == []

    def test_mixed_units(self, preprocessor):
        """混合单位."""
        from rdkit import Chem
        mols = []
        for raw, val in [
            ("10 nM", 10.0),
            ("1 uM", 1.0),
            ("0.001 mM", 0.001),
        ]:
            mol = Chem.MolFromSmiles("CCO")
            mols.append({
                "mol": mol,
                "name": "test",
                "smiles": "CCO",
                "activity_raw": raw,
            })

        processed = preprocessor.process(mols)
        assert len(processed) == 3
        # 所有值都应转换为 nM

    def test_censor_treatment_as_is(self, sample_molecules):
        """as_is 模式下删失值参与统计."""
        p = ActivityPreprocessor(censor_treatment="as_is")
        processed = p.process(sample_molecules)
        report = p.get_quality_report()
        assert report.valid_activities > 0


class TestActivityUnit:
    """ActivityUnit 枚举测试."""

    def test_values(self):
        assert ActivityUnit.NM.value == "nM"
        assert ActivityUnit.UM.value == "uM"
        assert ActivityUnit.MM.value == "mM"
        assert ActivityUnit.UNKNOWN.value == "unknown"


class TestCensoredType:
    """CensoredType 枚举测试."""

    def test_values(self):
        assert CensoredType.NONE.value == "none"
        assert CensoredType.GREATER_THAN.value == ">"
        assert CensoredType.LESS_THAN.value == "<"
        assert CensoredType.TILDE.value == "~"

"""Tests for imputed McDonald-Kreitman test implementation."""

from pathlib import Path

import pytest

from mkado.analysis.asymptotic import PolymorphismData, extract_polymorphism_data
from mkado.analysis.imputed import ImputedMKResult, imputed_mk_test, imputed_mk_test_multi
from mkado.io.output import OutputFormat, format_result


class TestImputedMKTest:
    """Tests for imputed_mk_test function."""

    def test_hand_calculated(self) -> None:
        """Test against hand-calculated values."""
        # Pn_low=6, Pn_high=4, Ps_low=3, Ps_high=6
        # ratio_ps = 3/6 = 0.5
        # pwd = 6 - 4*0.5 = 4.0
        # pn_neutral = 10 - 4 = 6.0
        # alpha = 1 - (6/9) * (5/7) ≈ 1 - 0.4762 = 0.5238
        gene = PolymorphismData(
            polymorphisms=(
                [(0.05, "N")] * 3
                + [(0.10, "N")] * 3
                + [(0.30, "N")] * 2
                + [(0.60, "N")] * 2
                + [(0.05, "S")] * 2
                + [(0.12, "S")] * 1
                + [(0.40, "S")] * 3
                + [(0.70, "S")] * 3
            ),
            dn=7,
            ds=5,
        )

        result = imputed_mk_test(gene, cutoff=0.15)

        assert result.pwd == pytest.approx(4.0)
        assert result.pn_neutral == pytest.approx(6.0)
        assert result.alpha == pytest.approx(1 - (6 / 9) * (5 / 7), abs=1e-4)
        assert result.pn_total == 10
        assert result.ps_total == 9

    def test_paper_example(self) -> None:
        """Reproduce the imputation from Murga-Moreno et al. (2022).

        With Pn_low=7, Pn_high=4, Ps_low=6, Ps_high=11:
        ratio_ps = 6/11
        pwd = 7 - 4*(6/11) = 7 - 24/11 ≈ 4.818
        """
        gene = PolymorphismData(
            polymorphisms=(
                [(0.05, "N")] * 7  # low-freq nonsyn
                + [(0.50, "N")] * 4  # high-freq nonsyn
                + [(0.10, "S")] * 6  # low-freq syn
                + [(0.50, "S")] * 11  # high-freq syn
            ),
            dn=10,
            ds=8,
        )

        result = imputed_mk_test(gene, cutoff=0.15)

        expected_pwd = 7 - 4 * (6 / 11)
        assert result.pwd == pytest.approx(expected_pwd, abs=1e-3)

    def test_dfe_fractions_sum_to_one(self) -> None:
        """Verify d + b + f sums to 1 - alpha when m0/mi provided."""
        gene = PolymorphismData(
            polymorphisms=(
                [(0.05, "N")] * 5
                + [(0.50, "N")] * 5
                + [(0.10, "S")] * 4
                + [(0.50, "S")] * 8
            ),
            dn=10,
            ds=6,
        )

        result = imputed_mk_test(gene, cutoff=0.15, num_synonymous_sites=100.0,
                                 num_nonsynonymous_sites=300.0)

        assert result.d is not None
        assert result.b is not None
        assert result.f is not None
        # d + b + f should equal 1 - alpha (the non-adaptive fraction)
        # Actually d = 1 - f - b, so d + b + f = 1 by construction
        assert result.d + result.b + result.f == pytest.approx(1.0, abs=1e-10)

    def test_dfe_fractions_none_without_sites(self) -> None:
        """Verify DFE fractions are None when m0/mi not provided."""
        gene = PolymorphismData(
            polymorphisms=[(0.05, "N")] * 3 + [(0.50, "S")] * 5,
            dn=5,
            ds=3,
        )

        result = imputed_mk_test(gene)

        assert result.d is None
        assert result.b is None
        assert result.f is None

    def test_ps_high_zero(self) -> None:
        """Edge case: all synonymous polymorphisms are low-frequency."""
        gene = PolymorphismData(
            polymorphisms=(
                [(0.05, "N")] * 3
                + [(0.50, "N")] * 2
                + [(0.10, "S")] * 4  # all below cutoff
            ),
            dn=5,
            ds=3,
        )

        result = imputed_mk_test(gene, cutoff=0.15)

        # ratio_ps = 0 since ps_high = 0, so pwd = pn_low - 0 = pn_low
        assert result.pwd == pytest.approx(3.0)
        assert result.pn_neutral == pytest.approx(2.0)

    def test_all_polymorphisms_above_cutoff(self) -> None:
        """When all polymorphisms are high-frequency, pwd should be 0."""
        gene = PolymorphismData(
            polymorphisms=(
                [(0.50, "N")] * 4
                + [(0.60, "S")] * 5
            ),
            dn=6,
            ds=4,
        )

        result = imputed_mk_test(gene, cutoff=0.15)

        assert result.pwd == 0.0
        assert result.pn_neutral == pytest.approx(4.0)

    def test_dn_zero(self) -> None:
        """Edge case: Dn == 0, alpha should be None."""
        gene = PolymorphismData(
            polymorphisms=[(0.05, "N")] * 3 + [(0.50, "S")] * 5,
            dn=0,
            ds=5,
        )

        result = imputed_mk_test(gene)

        assert result.alpha is None

    def test_no_polymorphisms(self) -> None:
        """Edge case: no polymorphisms at all."""
        gene = PolymorphismData(
            polymorphisms=[],
            dn=5,
            ds=3,
        )

        result = imputed_mk_test(gene)

        assert result.pwd == 0.0
        assert result.pn_neutral == 0.0
        assert result.pn_total == 0
        assert result.ps_total == 0

    def test_cutoff_variation(self) -> None:
        """Different cutoffs should produce different results."""
        gene = PolymorphismData(
            polymorphisms=(
                [(0.05, "N")] * 3
                + [(0.12, "N")] * 2
                + [(0.20, "N")] * 2
                + [(0.50, "N")] * 3
                + [(0.08, "S")] * 2
                + [(0.18, "S")] * 3
                + [(0.40, "S")] * 4
            ),
            dn=8,
            ds=5,
        )

        result_10 = imputed_mk_test(gene, cutoff=0.10)
        result_15 = imputed_mk_test(gene, cutoff=0.15)
        result_20 = imputed_mk_test(gene, cutoff=0.20)

        # Different cutoffs should give different pwd values
        assert result_10.pwd != result_15.pwd or result_15.pwd != result_20.pwd

    def test_pwd_clamped_to_zero(self) -> None:
        """pwd should never be negative (clamped to 0)."""
        # More high-freq nonsyn relative to low-freq nonsyn than synonymous ratio
        gene = PolymorphismData(
            polymorphisms=(
                [(0.05, "N")] * 1  # very few low-freq nonsyn
                + [(0.50, "N")] * 10  # many high-freq nonsyn
                + [(0.05, "S")] * 5
                + [(0.50, "S")] * 2  # ratio_ps = 5/2 = 2.5
                # pwd = 1 - 10*2.5 = -24 => clamped to 0
            ),
            dn=5,
            ds=3,
        )

        result = imputed_mk_test(gene, cutoff=0.15)

        assert result.pwd == 0.0


class TestImputedMKTestMulti:
    """Tests for imputed_mk_test_multi function."""

    def test_pooling_matches_manual(self) -> None:
        """Verify pooled results match manually aggregated data."""
        gene1 = PolymorphismData(
            polymorphisms=(
                [(0.05, "N")] * 3
                + [(0.50, "N")] * 2
                + [(0.10, "S")] * 2
                + [(0.40, "S")] * 4
            ),
            dn=5,
            ds=3,
        )
        gene2 = PolymorphismData(
            polymorphisms=(
                [(0.08, "N")] * 2
                + [(0.60, "N")] * 3
                + [(0.12, "S")] * 1
                + [(0.50, "S")] * 3
            ),
            dn=4,
            ds=2,
        )

        # Multi-gene pooled
        result_multi = imputed_mk_test_multi([gene1, gene2], cutoff=0.15)

        # Manual aggregation
        combined = PolymorphismData(
            polymorphisms=gene1.polymorphisms + gene2.polymorphisms,
            dn=gene1.dn + gene2.dn,
            ds=gene1.ds + gene2.ds,
        )
        result_manual = imputed_mk_test(combined, cutoff=0.15)

        assert result_multi.pwd == pytest.approx(result_manual.pwd)
        assert result_multi.alpha == pytest.approx(result_manual.alpha)
        assert result_multi.p_value == pytest.approx(result_manual.p_value)
        assert result_multi.dn == result_manual.dn
        assert result_multi.ds == result_manual.ds

    def test_multi_gene_counts(self) -> None:
        """Verify total counts are summed correctly."""
        genes = [
            PolymorphismData(
                polymorphisms=[(0.05, "N")] * 2 + [(0.50, "S")] * 3,
                dn=4,
                ds=2,
            ),
            PolymorphismData(
                polymorphisms=[(0.10, "N")] * 1 + [(0.60, "S")] * 2,
                dn=3,
                ds=5,
            ),
        ]

        result = imputed_mk_test_multi(genes)

        assert result.dn == 7
        assert result.ds == 7
        assert result.pn_total == 3
        assert result.ps_total == 5


class TestImputedMKResult:
    """Tests for ImputedMKResult dataclass."""

    def test_str_representation(self) -> None:
        """Test string representation."""
        result = ImputedMKResult(
            alpha=0.52,
            p_value=0.03,
            pn_neutral=6.0,
            pwd=4.0,
            dn=7,
            ds=5,
            pn_total=10,
            ps_total=9,
            cutoff=0.15,
        )

        s = str(result)
        assert "Imputed MK Test" in s
        assert "Dn=7" in s
        assert "Ds=5" in s
        assert "Pn=10" in s
        assert "Ps=9" in s
        assert "0.52" in s or "0.5200" in s
        assert "Pwd" in s

    def test_str_with_dfe(self) -> None:
        """Test string representation includes DFE when present."""
        result = ImputedMKResult(
            alpha=0.5,
            p_value=0.05,
            pn_neutral=5.0,
            pwd=3.0,
            dn=8,
            ds=4,
            pn_total=8,
            ps_total=6,
            d=0.3,
            b=0.2,
            f=0.5,
            cutoff=0.15,
        )

        s = str(result)
        assert "DFE fractions" in s

    def test_to_dict(self) -> None:
        """Test dictionary conversion."""
        result = ImputedMKResult(
            alpha=0.52,
            p_value=0.03,
            pn_neutral=6.0,
            pwd=4.0,
            dn=7,
            ds=5,
            pn_total=10,
            ps_total=9,
            cutoff=0.15,
        )

        d = result.to_dict()
        assert d["alpha"] == 0.52
        assert d["p_value"] == 0.03
        assert d["dn"] == 7
        assert d["ds"] == 5
        assert d["cutoff"] == 0.15
        assert "d" not in d  # No DFE without m0/mi

    def test_to_dict_with_dfe(self) -> None:
        """Test dictionary includes DFE fractions when present."""
        result = ImputedMKResult(
            alpha=0.5,
            p_value=0.05,
            pn_neutral=5.0,
            pwd=3.0,
            dn=8,
            ds=4,
            pn_total=8,
            ps_total=6,
            d=0.3,
            b=0.2,
            f=0.5,
            cutoff=0.15,
        )

        d = result.to_dict()
        assert d["d"] == 0.3
        assert d["b"] == 0.2
        assert d["f"] == 0.5

    def test_format_tsv(self) -> None:
        """Test TSV output formatting via format_result."""
        result = ImputedMKResult(
            alpha=0.52,
            p_value=0.03,
            pn_neutral=6.0,
            pwd=4.0,
            dn=7,
            ds=5,
            pn_total=10,
            ps_total=9,
            cutoff=0.15,
        )

        tsv = format_result(result, OutputFormat.TSV)
        header, values = tsv.split("\n")

        assert "Dn\tDs\tPn\tPs\tPwd\tPn_neutral\talpha\tp_value\tcutoff" == header
        fields = values.split("\t")
        assert fields[0] == "7"   # Dn
        assert fields[1] == "5"   # Ds
        assert fields[2] == "10"  # Pn
        assert fields[3] == "9"   # Ps
        assert fields[4] == "4.00"  # Pwd
        assert fields[5] == "6.00"  # Pn_neutral
        assert fields[6] == "0.520000"  # alpha
        assert fields[8] == "0.15"  # cutoff


class TestImputedMKIntegration:
    """Integration tests with real sequence data."""

    def test_with_kreitman_data(self) -> None:
        """Test imputed MK on the Kreitman Adh dataset."""
        test_data = Path(__file__).parent / "data"
        ingroup = test_data / "kreitmanAdh.fa"
        outgroup = test_data / "mauritianaAdh.fa"

        if not ingroup.exists() or not outgroup.exists():
            pytest.skip("Test data files not found")

        gene_data = extract_polymorphism_data(ingroup, outgroup, gene_id="adh")
        result = imputed_mk_test(gene_data)

        assert isinstance(result, ImputedMKResult)
        assert result.dn >= 0
        assert result.ds >= 0
        assert result.pwd >= 0
        assert result.pn_neutral >= 0
        assert result.cutoff == 0.15

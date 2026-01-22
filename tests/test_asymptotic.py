"""Tests for asymptotic MK test."""

from pathlib import Path

import numpy as np
import pytest

from mkado.analysis.asymptotic import (
    AggregatedSFS,
    AsymptoticMKResult,
    PolymorphismData,
    _compute_ci_monte_carlo,
    _exponential_model,
    _linear_model,
    aggregate_polymorphism_data,
    asymptotic_mk_test,
    asymptotic_mk_test_aggregated,
    extract_polymorphism_data,
)


class TestAsymptoticMKResult:
    """Tests for AsymptoticMKResult class."""

    def test_result_to_dict(self) -> None:
        """Test converting result to dictionary."""
        result = AsymptoticMKResult(
            frequency_bins=[0.1, 0.3, 0.5],
            alpha_by_freq=[0.1, 0.2, 0.3],
            alpha_asymptotic=0.4,
            ci_low=0.2,
            ci_high=0.6,
            fit_a=0.4,
            fit_b=0.3,
            fit_c=2.0,
            dn=10,
            ds=20,
        )

        d = result.to_dict()

        assert d["alpha_asymptotic"] == 0.4
        assert d["ci_low"] == 0.2
        assert d["ci_high"] == 0.6
        assert d["dn"] == 10
        assert d["ds"] == 20
        assert "fit_parameters" in d

    def test_result_str(self) -> None:
        """Test string representation."""
        result = AsymptoticMKResult(
            alpha_asymptotic=0.4,
            ci_low=0.2,
            ci_high=0.6,
            dn=10,
            ds=20,
        )

        s = str(result)
        assert "0.4" in s or "0.40" in s
        assert "Dn=10" in s


class TestAsymptoticMKTest:
    """Integration tests for asymptotic MK test."""

    def test_asymptotic_mk_test_basic(self, tmp_path: Path) -> None:
        """Test basic asymptotic MK test."""
        # Create test data with some variation
        ingroup_fa = tmp_path / "ingroup.fa"
        ingroup_fa.write_text(""">seq1
ATGATGATGATGATGATG
>seq2
ATGCTGATGATGATGATG
>seq3
ATGATGATGCTGATGATG
>seq4
ATGATGATGATGCTGATG
>seq5
ATGATGATGATGATGATG
""")

        outgroup_fa = tmp_path / "outgroup.fa"
        outgroup_fa.write_text(""">out1
ATGGTGATGGTGATGGTG
""")

        result = asymptotic_mk_test(ingroup_fa, outgroup_fa, num_bins=5)

        # Should complete without errors
        assert isinstance(result, AsymptoticMKResult)
        assert result.dn >= 0
        assert result.ds >= 0

    def test_asymptotic_mk_test_insufficient_data(self, tmp_path: Path) -> None:
        """Test asymptotic MK test with insufficient data."""
        # Create minimal test data
        ingroup_fa = tmp_path / "ingroup.fa"
        ingroup_fa.write_text(""">seq1
ATGATG
>seq2
ATGATG
""")

        outgroup_fa = tmp_path / "outgroup.fa"
        outgroup_fa.write_text(""">out1
ATGATG
""")

        result = asymptotic_mk_test(ingroup_fa, outgroup_fa)

        # Should handle gracefully
        assert isinstance(result, AsymptoticMKResult)

    def test_asymptotic_mk_test_with_kreitman_data(self) -> None:
        """Test asymptotic MK test with Kreitman data."""
        test_data = Path(__file__).parent / "data"
        ingroup = test_data / "kreitmanAdh.fa"
        outgroup = test_data / "mauritianaAdh.fa"

        if not ingroup.exists() or not outgroup.exists():
            pytest.skip("Test data files not found")

        result = asymptotic_mk_test(ingroup, outgroup, num_bins=5)

        assert isinstance(result, AsymptoticMKResult)
        # With only 4 ingroup sequences, might not have enough frequency data
        # but should still complete


class TestPolymorphismData:
    """Tests for PolymorphismData dataclass."""

    def test_creation(self) -> None:
        """Test creating PolymorphismData."""
        data = PolymorphismData(
            polymorphisms=[(0.2, "N"), (0.5, "S")],
            dn=5,
            ds=10,
            gene_id="gene1",
        )
        assert len(data.polymorphisms) == 2
        assert data.dn == 5
        assert data.ds == 10
        assert data.gene_id == "gene1"

    def test_empty_creation(self) -> None:
        """Test creating empty PolymorphismData."""
        data = PolymorphismData()
        assert len(data.polymorphisms) == 0
        assert data.dn == 0
        assert data.ds == 0


class TestAggregatedSFS:
    """Tests for AggregatedSFS dataclass."""

    def test_creation(self) -> None:
        """Test creating AggregatedSFS."""
        sfs = AggregatedSFS(
            bin_edges=np.linspace(0, 1, 11),
            pn_counts=np.array([1, 2, 3, 4, 5, 5, 4, 3, 2, 1]),
            ps_counts=np.array([2, 3, 4, 5, 6, 6, 5, 4, 3, 2]),
            dn_total=100,
            ds_total=200,
            num_genes=50,
        )
        assert sfs.dn_total == 100
        assert sfs.ds_total == 200
        assert sfs.num_genes == 50
        assert len(sfs.bin_edges) == 11


class TestExtractPolymorphismData:
    """Tests for extract_polymorphism_data function."""

    def test_extract_basic(self, tmp_path: Path) -> None:
        """Test basic extraction."""
        ingroup_fa = tmp_path / "ingroup.fa"
        ingroup_fa.write_text(""">seq1
ATGATGATGATGATGATG
>seq2
ATGCTGATGATGATGATG
>seq3
ATGATGATGCTGATGATG
>seq4
ATGATGATGATGCTGATG
>seq5
ATGATGATGATGATGATG
""")

        outgroup_fa = tmp_path / "outgroup.fa"
        outgroup_fa.write_text(""">out1
ATGGTGATGGTGATGGTG
""")

        data = extract_polymorphism_data(ingroup_fa, outgroup_fa, gene_id="test_gene")

        assert isinstance(data, PolymorphismData)
        assert data.gene_id == "test_gene"
        assert data.dn >= 0
        assert data.ds >= 0

    def test_extract_matches_standard_mk_counts(self, tmp_path: Path) -> None:
        """Test that extracted counts match standard MK test."""
        from mkado.analysis.mk_test import mk_test

        ingroup_fa = tmp_path / "ingroup.fa"
        ingroup_fa.write_text(""">seq1
ATGATGATGATGATGATG
>seq2
ATGCTGATGATGATGATG
>seq3
ATGATGATGCTGATGATG
>seq4
ATGATGATGATGCTGATG
>seq5
ATGATGATGATGATGATG
""")

        outgroup_fa = tmp_path / "outgroup.fa"
        outgroup_fa.write_text(""">out1
ATGGTGATGGTGATGGTG
""")

        data = extract_polymorphism_data(ingroup_fa, outgroup_fa)
        mk_result = mk_test(ingroup_fa, outgroup_fa)

        # Divergence counts should match
        assert data.dn == mk_result.dn
        assert data.ds == mk_result.ds


class TestAggregatePolymorphismData:
    """Tests for aggregate_polymorphism_data function."""

    def test_aggregate_single_gene(self) -> None:
        """Test aggregation with a single gene."""
        gene_data = [
            PolymorphismData(
                polymorphisms=[(0.2, "N"), (0.5, "S"), (0.8, "N")],
                dn=5,
                ds=10,
                gene_id="gene1",
            )
        ]

        agg = aggregate_polymorphism_data(gene_data, num_bins=10)

        assert agg.dn_total == 5
        assert agg.ds_total == 10
        assert agg.num_genes == 1
        assert np.sum(agg.pn_counts) == 2  # Two N polymorphisms
        assert np.sum(agg.ps_counts) == 1  # One S polymorphism

    def test_aggregate_multiple_genes(self) -> None:
        """Test aggregation with multiple genes."""
        gene_data = [
            PolymorphismData(
                polymorphisms=[(0.2, "N"), (0.5, "S")],
                dn=5,
                ds=10,
                gene_id="gene1",
            ),
            PolymorphismData(
                polymorphisms=[(0.3, "N"), (0.7, "S"), (0.8, "N")],
                dn=3,
                ds=7,
                gene_id="gene2",
            ),
        ]

        agg = aggregate_polymorphism_data(gene_data, num_bins=10)

        assert agg.dn_total == 8  # 5 + 3
        assert agg.ds_total == 17  # 10 + 7
        assert agg.num_genes == 2
        assert np.sum(agg.pn_counts) == 3  # Three N polymorphisms total
        assert np.sum(agg.ps_counts) == 2  # Two S polymorphisms total

    def test_aggregate_preserves_totals(self) -> None:
        """Test that aggregation preserves total counts."""
        gene_data = [
            PolymorphismData(
                polymorphisms=[(0.1, "N")] * 10 + [(0.5, "S")] * 20,
                dn=15,
                ds=25,
            ),
            PolymorphismData(
                polymorphisms=[(0.2, "N")] * 5 + [(0.8, "S")] * 15,
                dn=10,
                ds=20,
            ),
        ]

        agg = aggregate_polymorphism_data(gene_data, num_bins=20)

        assert agg.dn_total == 25
        assert agg.ds_total == 45
        assert int(np.sum(agg.pn_counts)) == 15
        assert int(np.sum(agg.ps_counts)) == 35


class TestLinearModel:
    """Tests for _linear_model function."""

    def test_linear_model(self) -> None:
        """Test linear model evaluation."""
        x = np.array([0.0, 0.5, 1.0])
        result = _linear_model(x, a=0.1, b=0.8)

        np.testing.assert_array_almost_equal(result, [0.1, 0.5, 0.9])


class TestMonteCarloCi:
    """Tests for Monte Carlo CI calculation."""

    def test_monte_carlo_ci_linear(self) -> None:
        """Test Monte Carlo CI with linear model."""
        # Well-determined parameters
        popt = np.array([0.1, 0.8])
        pcov = np.array([[0.001, 0.0], [0.0, 0.001]])

        ci_low, ci_high = _compute_ci_monte_carlo(
            popt, pcov, _linear_model, n_sim=5000
        )

        # CI should bracket the point estimate at x=1
        alpha_at_1 = 0.1 + 0.8  # = 0.9
        assert ci_low < alpha_at_1
        assert ci_high > alpha_at_1
        assert ci_high - ci_low < 0.5  # CI should be reasonably tight

    def test_monte_carlo_ci_exponential(self) -> None:
        """Test Monte Carlo CI with exponential model."""
        # Model: α(x) = a + b * exp(-c * x)
        # For typical positive selection signal, b is negative (alpha increases with x)
        popt = np.array([0.5, -0.3, 2.0])
        pcov = np.array([
            [0.001, 0.0, 0.0],
            [0.0, 0.001, 0.0],
            [0.0, 0.0, 0.01],
        ])

        ci_low, ci_high = _compute_ci_monte_carlo(
            popt, pcov, _exponential_model, n_sim=5000
        )

        # CI should bracket the point estimate
        # α(1) = a + b * exp(-c) = 0.5 + (-0.3) * exp(-2.0)
        alpha_at_1 = 0.5 + (-0.3) * np.exp(-2.0)
        assert ci_low < alpha_at_1
        assert ci_high > alpha_at_1


class TestAsymptoticMKTestAggregated:
    """Tests for asymptotic_mk_test_aggregated function."""

    def test_aggregated_basic(self) -> None:
        """Test basic aggregated analysis."""
        # Create synthetic gene data with clear signal
        gene_data = []
        for i in range(10):
            # Create polymorphisms across frequency spectrum
            poly = []
            for freq in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]:
                poly.extend([(freq, "N")] * 2)
                poly.extend([(freq, "S")] * 3)
            gene_data.append(
                PolymorphismData(
                    polymorphisms=poly,
                    dn=10,
                    ds=15,
                    gene_id=f"gene{i}",
                )
            )

        result = asymptotic_mk_test_aggregated(gene_data, num_bins=10)

        assert isinstance(result, AsymptoticMKResult)
        assert result.num_genes == 10
        assert result.dn == 100  # 10 genes * 10 Dn each
        assert result.ds == 150  # 10 genes * 15 Ds each
        assert result.pn_total > 0
        assert result.ps_total > 0

    def test_aggregated_with_kreitman_data(self) -> None:
        """Test aggregated analysis with real data (single gene as sanity check)."""
        test_data = Path(__file__).parent / "data"
        ingroup = test_data / "kreitmanAdh.fa"
        outgroup = test_data / "mauritianaAdh.fa"

        if not ingroup.exists() or not outgroup.exists():
            pytest.skip("Test data files not found")

        gene_data = [
            extract_polymorphism_data(ingroup, outgroup, gene_id="adh")
        ]

        result = asymptotic_mk_test_aggregated(gene_data, num_bins=5)

        assert isinstance(result, AsymptoticMKResult)
        assert result.num_genes == 1

    def test_aggregated_result_str(self) -> None:
        """Test string representation of aggregated result."""
        result = AsymptoticMKResult(
            alpha_asymptotic=0.35,
            ci_low=0.25,
            ci_high=0.45,
            dn=100,
            ds=200,
            num_genes=50,
            model_type="exponential",
            pn_total=500,
            ps_total=800,
        )

        s = str(result)
        assert "Genes aggregated: 50" in s
        assert "Pn=500" in s
        assert "Ps=800" in s

    def test_aggregated_result_to_dict(self) -> None:
        """Test dictionary conversion of aggregated result."""
        result = AsymptoticMKResult(
            alpha_asymptotic=0.35,
            ci_low=0.25,
            ci_high=0.45,
            dn=100,
            ds=200,
            num_genes=50,
            model_type="linear",
            pn_total=500,
            ps_total=800,
            fit_a=0.1,
            fit_b=0.25,
        )

        d = result.to_dict()

        assert d["num_genes"] == 50
        assert d["pn_total"] == 500
        assert d["ps_total"] == 800
        assert d["model_type"] == "linear"
        assert "c" not in d["fit_parameters"]  # Linear model has no c

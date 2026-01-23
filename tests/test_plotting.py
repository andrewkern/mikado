"""Tests for plotting functionality."""

from pathlib import Path
from tempfile import TemporaryDirectory

import pytest

from mkado.analysis.asymptotic import AsymptoticMKResult
from mkado.analysis.mk_test import MKResult
from mkado.analysis.polarized import PolarizedMKResult
from mkado.io.plotting import create_asymptotic_plot, create_volcano_plot


class TestVolcanoPlot:
    """Tests for volcano plot generation."""

    def test_basic_volcano_plot(self) -> None:
        """Test basic volcano plot generation with MKResult."""
        results = [
            ("gene1", MKResult(dn=10, ds=20, pn=5, ps=10, p_value=0.5, ni=1.0, alpha=0.0)),
            ("gene2", MKResult(dn=15, ds=10, pn=3, ps=8, p_value=0.01, ni=0.5, alpha=0.5)),
            ("gene3", MKResult(dn=8, ds=25, pn=6, ps=12, p_value=0.001, ni=1.5, alpha=-0.5)),
        ]

        with TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "volcano.png"
            create_volcano_plot(results, output_path)
            assert output_path.exists()
            assert output_path.stat().st_size > 0

    def test_volcano_plot_with_polarized_results(self) -> None:
        """Test volcano plot with PolarizedMKResult."""
        results = [
            (
                "gene1",
                PolarizedMKResult(
                    dn_ingroup=10,
                    ds_ingroup=20,
                    pn_ingroup=5,
                    ps_ingroup=10,
                    dn_outgroup=3,
                    ds_outgroup=5,
                    dn_unpolarized=2,
                    ds_unpolarized=3,
                    pn_unpolarized=1,
                    ps_unpolarized=2,
                    p_value_ingroup=0.5,
                    ni_ingroup=1.0,
                    alpha_ingroup=0.0,
                ),
            ),
            (
                "gene2",
                PolarizedMKResult(
                    dn_ingroup=15,
                    ds_ingroup=10,
                    pn_ingroup=3,
                    ps_ingroup=8,
                    dn_outgroup=5,
                    ds_outgroup=7,
                    dn_unpolarized=1,
                    ds_unpolarized=2,
                    pn_unpolarized=0,
                    ps_unpolarized=1,
                    p_value_ingroup=0.01,
                    ni_ingroup=0.5,
                    alpha_ingroup=0.5,
                ),
            ),
        ]

        with TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "volcano_polarized.png"
            create_volcano_plot(results, output_path)
            assert output_path.exists()
            assert output_path.stat().st_size > 0

    def test_volcano_plot_skips_invalid_ni(self) -> None:
        """Test that genes with NI=0 or NI=None are skipped."""
        results = [
            ("gene1", MKResult(dn=10, ds=20, pn=5, ps=10, p_value=0.5, ni=1.0, alpha=0.0)),
            ("gene2", MKResult(dn=15, ds=10, pn=0, ps=8, p_value=1.0, ni=0.0, alpha=1.0)),
            ("gene3", MKResult(dn=8, ds=25, pn=6, ps=0, p_value=1.0, ni=None, alpha=None)),
        ]

        with TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "volcano.png"
            create_volcano_plot(results, output_path)
            assert output_path.exists()

    def test_volcano_plot_no_valid_data(self) -> None:
        """Test that ValueError is raised when no valid data."""
        results = [
            ("gene1", MKResult(dn=10, ds=20, pn=0, ps=10, p_value=1.0, ni=0.0, alpha=1.0)),
            ("gene2", MKResult(dn=15, ds=10, pn=0, ps=8, p_value=1.0, ni=None, alpha=None)),
        ]

        with TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "volcano.png"
            with pytest.raises(ValueError, match="No valid NI/p-value pairs"):
                create_volcano_plot(results, output_path)

    def test_volcano_plot_different_formats(self) -> None:
        """Test saving volcano plot in different formats."""
        results = [
            ("gene1", MKResult(dn=10, ds=20, pn=5, ps=10, p_value=0.5, ni=1.0, alpha=0.0)),
            ("gene2", MKResult(dn=15, ds=10, pn=3, ps=8, p_value=0.01, ni=0.5, alpha=0.5)),
        ]

        with TemporaryDirectory() as tmpdir:
            for ext in ["png", "pdf", "svg"]:
                output_path = Path(tmpdir) / f"volcano.{ext}"
                create_volcano_plot(results, output_path)
                assert output_path.exists()
                assert output_path.stat().st_size > 0


class TestAsymptoticPlot:
    """Tests for asymptotic MK test plot generation."""

    def test_basic_asymptotic_plot_exponential(self) -> None:
        """Test asymptotic plot with exponential fit."""
        result = AsymptoticMKResult(
            frequency_bins=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
            alpha_by_freq=[-0.2, 0.0, 0.1, 0.15, 0.2, 0.25, 0.28, 0.3, 0.32],
            alpha_x_values=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
            alpha_asymptotic=0.35,
            ci_low=0.30,
            ci_high=0.40,
            fit_a=0.35,
            fit_b=-0.55,
            fit_c=5.0,
            dn=100,
            ds=200,
            num_genes=50,
            model_type="exponential",
            pn_total=150,
            ps_total=300,
        )

        with TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "asymptotic.png"
            create_asymptotic_plot(result, output_path)
            assert output_path.exists()
            assert output_path.stat().st_size > 0

    def test_asymptotic_plot_linear(self) -> None:
        """Test asymptotic plot with linear fit."""
        result = AsymptoticMKResult(
            frequency_bins=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
            alpha_by_freq=[0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5],
            alpha_x_values=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
            alpha_asymptotic=0.55,
            ci_low=0.50,
            ci_high=0.60,
            fit_a=0.05,
            fit_b=0.5,
            fit_c=0.0,
            dn=80,
            ds=150,
            num_genes=0,  # Single gene, not aggregated
            model_type="linear",
            pn_total=0,
            ps_total=0,
        )

        with TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "asymptotic_linear.png"
            create_asymptotic_plot(result, output_path)
            assert output_path.exists()
            assert output_path.stat().st_size > 0

    def test_asymptotic_plot_no_data(self) -> None:
        """Test that ValueError is raised when no frequency data."""
        result = AsymptoticMKResult(
            frequency_bins=[],
            alpha_by_freq=[],
            alpha_x_values=[],
            alpha_asymptotic=0.35,
            ci_low=0.30,
            ci_high=0.40,
            fit_a=0.35,
            fit_b=-0.55,
            fit_c=5.0,
            dn=100,
            ds=200,
            model_type="exponential",
        )

        with TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "asymptotic.png"
            with pytest.raises(ValueError, match="No frequency bin data"):
                create_asymptotic_plot(result, output_path)

    def test_asymptotic_plot_different_formats(self) -> None:
        """Test saving asymptotic plot in different formats."""
        result = AsymptoticMKResult(
            frequency_bins=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
            alpha_by_freq=[-0.2, 0.0, 0.1, 0.15, 0.2, 0.25, 0.28, 0.3, 0.32],
            alpha_x_values=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
            alpha_asymptotic=0.35,
            ci_low=0.30,
            ci_high=0.40,
            fit_a=0.35,
            fit_b=-0.55,
            fit_c=5.0,
            dn=100,
            ds=200,
            model_type="exponential",
        )

        with TemporaryDirectory() as tmpdir:
            for ext in ["png", "pdf", "svg"]:
                output_path = Path(tmpdir) / f"asymptotic.{ext}"
                create_asymptotic_plot(result, output_path)
                assert output_path.exists()
                assert output_path.stat().st_size > 0

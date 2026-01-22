"""Plotting functions for MK test results."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from mkado.analysis.asymptotic import AsymptoticMKResult
    from mkado.analysis.mk_test import MKResult
    from mkado.analysis.polarized import PolarizedMKResult


def create_volcano_plot(
    results: list[tuple[str, MKResult | PolarizedMKResult]],
    output_path: Path,
    alpha: float = 0.05,
) -> None:
    """Create a volcano plot from batch MK test results.

    The volcano plot shows -log10(NI) on the X-axis and -log10(p-value) on the Y-axis.
    A horizontal line indicates the Bonferroni-corrected significance threshold.

    Args:
        results: List of (gene_name, result) tuples from batch MK tests
        output_path: Path to save the plot (PNG, PDF, or SVG)
        alpha: Significance level for Bonferroni correction (default 0.05)
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    from mkado.analysis.mk_test import MKResult
    from mkado.analysis.polarized import PolarizedMKResult

    # Set seaborn dark grid style for modern look
    sns.set_theme(style="darkgrid")

    # Extract NI and p-values
    ni_values = []
    p_values = []
    gene_names = []

    for name, result in results:
        if isinstance(result, MKResult):
            ni = result.ni
            pval = result.p_value
        elif isinstance(result, PolarizedMKResult):
            ni = result.ni_ingroup
            pval = result.p_value_ingroup
        else:
            continue

        # Skip if NI is None or invalid
        if ni is None or ni <= 0 or pval <= 0:
            continue

        ni_values.append(ni)
        p_values.append(pval)
        gene_names.append(name)

    if not ni_values:
        raise ValueError("No valid NI/p-value pairs found in results")

    # Convert to numpy arrays
    ni_arr = np.array(ni_values)
    p_arr = np.array(p_values)

    # Calculate -log10 values
    neg_log10_ni = -np.log10(ni_arr)
    neg_log10_p = -np.log10(p_arr)

    # Bonferroni correction threshold
    n_tests = len(results)
    bonferroni_threshold = alpha / n_tests
    neg_log10_threshold = -np.log10(bonferroni_threshold)

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))

    # Determine point colors based on significance
    significant = p_arr < bonferroni_threshold
    colors = np.where(significant, "#e74c3c", "#3498db")  # Red for significant, blue otherwise

    # Create scatter plot
    scatter = ax.scatter(
        neg_log10_ni,
        neg_log10_p,
        c=colors,
        alpha=0.7,
        edgecolors="white",
        linewidth=0.5,
        s=60,
    )

    # Add Bonferroni threshold line
    ax.axhline(
        y=neg_log10_threshold,
        color="#e74c3c",
        linestyle="--",
        linewidth=1.5,
        label=f"Bonferroni threshold (p = {bonferroni_threshold:.2e})",
    )

    # Add vertical line at NI = 1 (neutral expectation)
    ax.axvline(
        x=0,  # -log10(1) = 0
        color="#95a5a6",
        linestyle=":",
        linewidth=1,
        alpha=0.7,
        label="NI = 1 (neutral)",
    )

    # Labels and title
    ax.set_xlabel("-log$_{10}$(NI)", fontsize=12)
    ax.set_ylabel("-log$_{10}$(p-value)", fontsize=12)
    ax.set_title("Volcano Plot: McDonald-Kreitman Test Results", fontsize=14, fontweight="bold")

    # Add legend
    ax.legend(loc="upper right", framealpha=0.9)

    # Add annotation for interpretation
    xlim = ax.get_xlim()
    ax.text(
        xlim[0] + 0.02 * (xlim[1] - xlim[0]),
        neg_log10_threshold + 0.3,
        f"n = {len(ni_values)} genes, {sum(significant)} significant",
        fontsize=9,
        color="#7f8c8d",
    )

    # Tight layout and save
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def create_asymptotic_plot(
    result: AsymptoticMKResult,
    output_path: Path,
) -> None:
    """Create an asymptotic MK test plot showing alpha(x) vs derived allele frequency.

    This plot follows the style of Messer & Petrov (2013), showing:
    - Scatter points of alpha at each frequency bin
    - The fitted curve (exponential or linear)
    - A horizontal line at alpha_asymptotic with confidence interval band

    Args:
        result: AsymptoticMKResult from asymptotic MK test
        output_path: Path to save the plot (PNG, PDF, or SVG)
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    # Set seaborn dark grid style for modern look
    sns.set_theme(style="darkgrid")

    # Extract data - use alpha_x_values if available, otherwise fall back to frequency_bins
    if result.alpha_x_values and len(result.alpha_x_values) == len(result.alpha_by_freq):
        x_data = np.array(result.alpha_x_values)
    else:
        x_data = np.array(result.frequency_bins)
    y_data = np.array(result.alpha_by_freq)

    if len(x_data) == 0 or len(y_data) == 0:
        raise ValueError("No frequency bin data available for plotting")

    if len(x_data) != len(y_data):
        raise ValueError(
            f"Mismatched data: {len(x_data)} x values vs {len(y_data)} alpha values"
        )

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot scatter points
    ax.scatter(
        x_data,
        y_data,
        c="#2c3e50",
        alpha=0.8,
        edgecolors="white",
        linewidth=0.5,
        s=80,
        zorder=3,
        label="Observed α(x)",
    )

    # Plot confidence interval band for alpha_asymptotic
    ax.fill_between(
        [0, 1],
        [result.ci_low, result.ci_low],
        [result.ci_high, result.ci_high],
        color="#95a5a6",
        alpha=0.3,
        zorder=1,
        label=f"95% CI ({result.ci_low:.2f} - {result.ci_high:.2f})",
    )

    # Plot alpha_asymptotic horizontal line
    ax.axhline(
        y=result.alpha_asymptotic,
        color="#e74c3c",
        linestyle="--",
        linewidth=2,
        zorder=2,
        label=f"α$_{{asym}}$ = {result.alpha_asymptotic:.2f}",
    )

    # Plot fitted curve
    x_curve = np.linspace(0.05, 1.0, 100)
    if result.model_type == "exponential":
        y_curve = result.fit_a + result.fit_b * np.exp(-result.fit_c * x_curve)
        fit_label = f"Fit: a + b·e$^{{-cx}}$"
    else:
        y_curve = result.fit_a + result.fit_b * x_curve
        fit_label = f"Fit: a + b·x"

    ax.plot(
        x_curve,
        y_curve,
        color="#e74c3c",
        linewidth=2,
        zorder=2,
        label=fit_label,
    )

    # Labels and title
    ax.set_xlabel("Derived allele frequency x", fontsize=12)
    ax.set_ylabel("MK α(x)", fontsize=12)
    ax.set_title("Asymptotic MK Test: α(x) vs Frequency", fontsize=14, fontweight="bold")

    # Set axis limits
    ax.set_xlim(0, 1.05)

    # Add legend
    ax.legend(loc="lower right", framealpha=0.9)

    # Add annotation with fit parameters
    if result.model_type == "exponential":
        param_text = f"a = {result.fit_a:.3f}, b = {result.fit_b:.3f}, c = {result.fit_c:.3f}"
    else:
        param_text = f"a = {result.fit_a:.3f}, b = {result.fit_b:.3f}"

    # Add gene count if aggregated
    if result.num_genes > 0:
        param_text += f"\nn = {result.num_genes} genes"

    ax.text(
        0.02,
        0.98,
        param_text,
        transform=ax.transAxes,
        fontsize=9,
        verticalalignment="top",
        color="#7f8c8d",
        family="monospace",
    )

    # Tight layout and save
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

"""Asymptotic McDonald-Kreitman test implementation.

Based on Messer & Petrov (2013) PNAS.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Callable

import numpy as np
from scipy import optimize

from mkado.core.alignment import AlignedPair
from mkado.core.codons import DEFAULT_CODE, GeneticCode
from mkado.core.sequences import SequenceSet

if TYPE_CHECKING:
    pass


@dataclass
class PolymorphismData:
    """Intermediate data from a single gene for aggregation."""

    polymorphisms: list[tuple[float, str]] = field(default_factory=list)
    """List of (derived_frequency, type) where type is 'N' or 'S'."""
    dn: int = 0
    ds: int = 0
    gene_id: str = ""


@dataclass
class AggregatedSFS:
    """Aggregated site frequency spectrum across genes."""

    bin_edges: np.ndarray = field(default_factory=lambda: np.array([]))
    pn_counts: np.ndarray = field(default_factory=lambda: np.array([]))
    """Pn per bin."""
    ps_counts: np.ndarray = field(default_factory=lambda: np.array([]))
    """Ps per bin."""
    dn_total: int = 0
    ds_total: int = 0
    num_genes: int = 0


@dataclass
class AsymptoticMKResult:
    """Results from an asymptotic McDonald-Kreitman test."""

    # Frequency bins and alpha estimates
    frequency_bins: list[float] = field(default_factory=list)
    alpha_by_freq: list[float] = field(default_factory=list)

    # Asymptotic alpha estimate (extrapolated to x=1)
    alpha_asymptotic: float = 0.0

    # 95% confidence interval
    ci_low: float = 0.0
    ci_high: float = 0.0

    # Fit parameters (alpha = a + b*exp(-c*x) for exponential, a + b*x for linear)
    fit_a: float = 0.0
    fit_b: float = 0.0
    fit_c: float = 0.0

    # Counts used in calculation
    dn: int = 0
    ds: int = 0

    # Aggregation info (0 = single gene, >0 = aggregated)
    num_genes: int = 0

    # Model type ("exponential" or "linear")
    model_type: str = "exponential"

    # Total polymorphism counts (useful for aggregated results)
    pn_total: int = 0
    ps_total: int = 0

    def __str__(self) -> str:
        lines = [
            "Asymptotic MK Test Results:",
            f"  Asymptotic α: {self.alpha_asymptotic:.4f} "
            f"(95% CI: {self.ci_low:.4f} - {self.ci_high:.4f})",
            f"  Divergence: Dn={self.dn}, Ds={self.ds}",
        ]
        if self.pn_total > 0 or self.ps_total > 0:
            lines.append(f"  Polymorphism: Pn={self.pn_total}, Ps={self.ps_total}")
        if self.num_genes > 0:
            lines.append(f"  Genes aggregated: {self.num_genes}")
        if self.model_type == "exponential":
            lines.append(
                f"  Fit ({self.model_type}): α(x) = {self.fit_a:.4f} + "
                f"({self.fit_b:.4f}) * exp(-{self.fit_c:.4f} * x)"
            )
        else:
            lines.append(
                f"  Fit ({self.model_type}): α(x) = {self.fit_a:.4f} + "
                f"{self.fit_b:.4f} * x"
            )
        return "\n".join(lines)

    def to_dict(self) -> dict:
        """Convert results to a dictionary."""
        result = {
            "alpha_asymptotic": self.alpha_asymptotic,
            "ci_low": self.ci_low,
            "ci_high": self.ci_high,
            "dn": self.dn,
            "ds": self.ds,
            "model_type": self.model_type,
            "fit_parameters": {
                "a": self.fit_a,
                "b": self.fit_b,
            },
            "frequency_bins": self.frequency_bins,
            "alpha_by_freq": self.alpha_by_freq,
        }
        if self.model_type == "exponential":
            result["fit_parameters"]["c"] = self.fit_c
        if self.num_genes > 0:
            result["num_genes"] = self.num_genes
        if self.pn_total > 0 or self.ps_total > 0:
            result["pn_total"] = self.pn_total
            result["ps_total"] = self.ps_total
        return result


def _exponential_model(x: np.ndarray, a: float, b: float, c: float) -> np.ndarray:
    """Exponential model: α(x) = a + b * exp(-c * x)

    Following the asymptoticMK R package convention.
    For positive selection with deleterious mutations, b is typically negative
    (alpha increases with frequency).
    """
    return a + b * np.exp(-c * x)


def _linear_model(x: np.ndarray, a: float, b: float) -> np.ndarray:
    """Linear model: α(x) = a + b * x"""
    return a + b * x


def _compute_ci_monte_carlo(
    popt: np.ndarray,
    pcov: np.ndarray,
    model_func: Callable[..., np.ndarray],
    n_sim: int = 10000,
    seed: int = 42,
) -> tuple[float, float]:
    """Compute CI via Monte Carlo simulation from covariance matrix.

    Args:
        popt: Optimal parameter values from curve_fit
        pcov: Covariance matrix from curve_fit
        model_func: Model function (_exponential_model or _linear_model)
        n_sim: Number of Monte Carlo simulations
        seed: Random seed for reproducibility

    Returns:
        Tuple of (ci_low, ci_high) for alpha at x=1
    """
    rng = np.random.default_rng(seed)

    # Handle case where covariance contains inf/nan
    if not np.all(np.isfinite(pcov)):
        alpha_at_1 = float(model_func(np.array([1.0]), *popt)[0])
        return (alpha_at_1, alpha_at_1)

    # Sample parameters from multivariate normal
    try:
        param_samples = rng.multivariate_normal(popt, pcov, size=n_sim)
    except np.linalg.LinAlgError:
        alpha_at_1 = float(model_func(np.array([1.0]), *popt)[0])
        return (alpha_at_1, alpha_at_1)

    # Evaluate model at x=1 for each sample
    alpha_samples = []
    for params in param_samples:
        try:
            alpha = float(model_func(np.array([1.0]), *params)[0])
            if np.isfinite(alpha):
                alpha_samples.append(alpha)
        except (ValueError, RuntimeWarning):
            continue

    if len(alpha_samples) < 100:
        alpha_at_1 = float(model_func(np.array([1.0]), *popt)[0])
        return (alpha_at_1, alpha_at_1)

    ci_low = float(np.percentile(alpha_samples, 2.5))
    ci_high = float(np.percentile(alpha_samples, 97.5))
    return (ci_low, ci_high)


def _compute_aic(n: int, rss: float, k: int) -> float:
    """Compute Akaike Information Criterion.

    Args:
        n: Number of data points
        rss: Residual sum of squares
        k: Number of parameters

    Returns:
        AIC value
    """
    if rss <= 0 or n <= k:
        return float("inf")
    return n * np.log(rss / n) + 2 * k


def extract_polymorphism_data(
    ingroup: SequenceSet | str | Path,
    outgroup: SequenceSet | str | Path,
    reading_frame: int = 1,
    genetic_code: GeneticCode | None = None,
    pool_polymorphisms: bool = False,
    gene_id: str = "",
) -> PolymorphismData:
    """Extract polymorphism and divergence data without curve fitting.

    This function extracts the raw data needed for asymptotic MK analysis,
    allowing multiple genes to be processed and aggregated before fitting.

    Args:
        ingroup: SequenceSet or path to FASTA file for ingroup sequences
        outgroup: SequenceSet or path to FASTA file for outgroup sequences
        reading_frame: Reading frame (1, 2, or 3)
        genetic_code: Genetic code for translation (uses standard if None)
        pool_polymorphisms: If True, consider sites polymorphic in either
            population (libsequence convention). Frequencies are still
            calculated from ingroup only.
        gene_id: Identifier for this gene

    Returns:
        PolymorphismData containing polymorphisms with frequencies and divergence counts
    """
    code = genetic_code or DEFAULT_CODE

    # Load sequences if paths provided
    if isinstance(ingroup, (str, Path)):
        ingroup = SequenceSet.from_fasta(ingroup, reading_frame, code)
    if isinstance(outgroup, (str, Path)):
        outgroup = SequenceSet.from_fasta(outgroup, reading_frame, code)

    # Create aligned pair
    pair = AlignedPair(ingroup=ingroup, outgroup=outgroup, genetic_code=code)

    # Count divergence
    dn = 0
    ds = 0

    for codon_idx in range(pair.num_codons):
        if pair.is_fixed_between(codon_idx):
            result = pair.classify_fixed_difference(codon_idx)
            if result is not None:
                nonsyn, syn = result
                dn += nonsyn
                ds += syn

    # Collect polymorphisms with their frequencies
    poly_data: list[tuple[float, str]] = []

    # Get polymorphic sites based on mode
    if pool_polymorphisms:
        poly_sites = pair.polymorphic_sites_pooled()
    else:
        poly_sites = pair.polymorphic_sites_ingroup()

    for codon_idx in poly_sites:
        codons = list(ingroup.codon_set_clean(codon_idx))
        if len(codons) < 2:
            continue

        # Get frequency spectrum
        freqs = ingroup.site_frequency_spectrum(codon_idx)
        if not freqs:
            continue

        # Get outgroup codon to determine ancestral state
        out_codons = outgroup.codon_set_clean(codon_idx)
        if not out_codons:
            continue

        # Find ancestral codon: must be shared between ingroup and outgroup
        ingroup_codons = set(freqs.keys())
        shared_codons = ingroup_codons & out_codons
        if not shared_codons:
            continue
        ancestral = max(shared_codons, key=lambda c: freqs.get(c, 0))

        # Calculate derived allele frequency
        derived_freq = 1.0 - freqs[ancestral]

        if derived_freq <= 0 or derived_freq >= 1:
            continue

        # Classify the polymorphism
        result = pair.classify_polymorphism(codon_idx)
        if result is not None:
            nonsyn, syn = result
            for _ in range(nonsyn):
                poly_data.append((derived_freq, "N"))
            for _ in range(syn):
                poly_data.append((derived_freq, "S"))

    return PolymorphismData(
        polymorphisms=poly_data,
        dn=dn,
        ds=ds,
        gene_id=gene_id,
    )


def aggregate_polymorphism_data(
    gene_data: list[PolymorphismData],
    num_bins: int = 20,
) -> AggregatedSFS:
    """Aggregate polymorphism data from multiple genes into SFS bins.

    Args:
        gene_data: List of PolymorphismData from individual genes
        num_bins: Number of frequency bins

    Returns:
        AggregatedSFS with per-bin Pn/Ps counts and total divergence
    """
    # Sum divergence across all genes
    dn_total = sum(g.dn for g in gene_data)
    ds_total = sum(g.ds for g in gene_data)

    # Combine all polymorphisms
    all_polymorphisms: list[tuple[float, str]] = []
    for g in gene_data:
        all_polymorphisms.extend(g.polymorphisms)

    # Create frequency bins
    bin_edges = np.linspace(0, 1, num_bins + 1)

    # Bin polymorphisms by derived allele frequency
    pn_counts = np.zeros(num_bins)
    ps_counts = np.zeros(num_bins)

    for freq, poly_type in all_polymorphisms:
        # Find which bin this frequency falls into
        bin_idx = np.searchsorted(bin_edges[1:], freq, side="right")
        bin_idx = min(bin_idx, num_bins - 1)

        if poly_type == "N":
            pn_counts[bin_idx] += 1
        else:
            ps_counts[bin_idx] += 1

    return AggregatedSFS(
        bin_edges=bin_edges,
        pn_counts=pn_counts,
        ps_counts=ps_counts,
        dn_total=dn_total,
        ds_total=ds_total,
        num_genes=len(gene_data),
    )


def asymptotic_mk_test_aggregated(
    gene_data: list[PolymorphismData],
    num_bins: int = 20,
    ci_replicates: int = 10000,
    frequency_cutoffs: tuple[float, float] = (0.1, 0.9),
) -> AsymptoticMKResult:
    """Genome-wide asymptotic MK test on aggregated data.

    This follows the approach of Messer & Petrov (2013) and the asymptoticMK
    R package: aggregate polymorphism SFS data and divergence counts from
    many genes, then fit a single exponential curve to the aggregated data.

    Args:
        gene_data: List of PolymorphismData from individual genes
        num_bins: Number of frequency bins (default 20)
        ci_replicates: Number of Monte Carlo replicates for CI (default 10000)
        frequency_cutoffs: (low, high) frequency range for fitting (default 0.1-0.9)

    Returns:
        AsymptoticMKResult with aggregated analysis
    """
    # Aggregate the data
    agg = aggregate_polymorphism_data(gene_data, num_bins)

    dn = agg.dn_total
    ds = agg.ds_total
    bin_edges = agg.bin_edges
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Calculate total polymorphisms
    pn_total = int(np.sum(agg.pn_counts))
    ps_total = int(np.sum(agg.ps_counts))

    # Calculate per-bin alpha for each frequency bin (matching asymptoticMK R package)
    # α(x) = 1 - (d0/d) * (p(x)/p0(x)) where p(x) and p0(x) are per-bin counts
    alpha_values = []
    valid_centers = []

    for i in range(num_bins):
        # Per-bin counts (not cumulative)
        pn_bin = int(agg.pn_counts[i])
        ps_bin = int(agg.ps_counts[i])

        if ds > 0 and dn > 0 and ps_bin > 0:
            alpha_x = 1.0 - (ds / dn) * (pn_bin / ps_bin)
            alpha_values.append(alpha_x)
            valid_centers.append(bin_centers[i])

    # Not enough data for curve fitting
    if len(valid_centers) < 3:
        simple_alpha = None
        if ds > 0 and dn > 0 and ps_total > 0:
            simple_alpha = 1.0 - (ds * pn_total) / (dn * ps_total)

        return AsymptoticMKResult(
            frequency_bins=list(bin_centers),
            alpha_by_freq=alpha_values,
            alpha_asymptotic=simple_alpha or 0.0,
            ci_low=simple_alpha or 0.0,
            ci_high=simple_alpha or 0.0,
            dn=dn,
            ds=ds,
            num_genes=agg.num_genes,
            pn_total=pn_total,
            ps_total=ps_total,
        )

    x_data = np.array(valid_centers)
    y_data = np.array(alpha_values)

    # Filter to frequency cutoffs for fitting
    low_cut, high_cut = frequency_cutoffs
    mask = (x_data >= low_cut) & (x_data <= high_cut)
    x_fit = x_data[mask]
    y_fit = y_data[mask]

    if len(x_fit) < 3:
        # Fall back to full range if cutoffs are too restrictive
        x_fit = x_data
        y_fit = y_data

    # Try exponential fit first
    exp_success = False
    exp_popt = None
    exp_pcov = None
    exp_alpha = 0.0
    exp_ci_width = float("inf")

    try:
        # Initial guesses: a is asymptotic value, b is negative for increasing alpha
        p0_exp = [y_fit[-1], y_fit[0] - y_fit[-1], 5.0]
        # b can be negative (typically is for increasing alpha with frequency)
        bounds_exp = ([-2.0, -2.0, 0.001], [2.0, 2.0, 100.0])

        exp_popt, exp_pcov = optimize.curve_fit(
            _exponential_model,
            x_fit,
            y_fit,
            p0=p0_exp,
            bounds=bounds_exp,
            maxfev=10000,
        )

        a, b, c = exp_popt
        # α(x=1) = a + b * exp(-c)
        exp_alpha = float(a + b * np.exp(-c))

        # Compute CI via Monte Carlo
        ci_low_exp, ci_high_exp = _compute_ci_monte_carlo(
            exp_popt, exp_pcov, _exponential_model, ci_replicates
        )
        exp_ci_width = ci_high_exp - ci_low_exp
        exp_success = True

    except (RuntimeError, ValueError):
        pass

    # Try linear fit
    lin_success = False
    lin_popt = None
    lin_pcov = None
    lin_alpha = 0.0

    try:
        p0_lin = [y_fit[0], y_fit[-1] - y_fit[0]]
        bounds_lin = ([-2.0, -2.0], [2.0, 2.0])

        lin_popt, lin_pcov = optimize.curve_fit(
            _linear_model,
            x_fit,
            y_fit,
            p0=p0_lin,
            bounds=bounds_lin,
            maxfev=10000,
        )

        a_lin, b_lin = lin_popt
        lin_alpha = float(a_lin + b_lin)

        ci_low_lin, ci_high_lin = _compute_ci_monte_carlo(
            lin_popt, lin_pcov, _linear_model, ci_replicates
        )
        lin_success = True

    except (RuntimeError, ValueError):
        pass

    # Select best model
    # Following asymptoticMK R package: prefer exponential unless CI > 100, then use linear
    # Also compare AIC when both succeed
    use_exponential = False

    if exp_success and lin_success:
        # Calculate AIC for both models
        exp_residuals = y_fit - _exponential_model(x_fit, *exp_popt)
        lin_residuals = y_fit - _linear_model(x_fit, *lin_popt)
        exp_rss = float(np.sum(exp_residuals**2))
        lin_rss = float(np.sum(lin_residuals**2))

        exp_aic = _compute_aic(len(x_fit), exp_rss, 3)
        lin_aic = _compute_aic(len(x_fit), lin_rss, 2)

        # Following asymptoticMK: use linear if exponential CI exceeds 100
        # or if linear AIC is better
        if exp_ci_width > 100:
            use_exponential = False
        elif lin_aic < exp_aic:
            use_exponential = False
        else:
            use_exponential = True
    elif exp_success:
        use_exponential = exp_ci_width <= 100
    elif lin_success:
        use_exponential = False
    else:
        # Both failed, use last value
        return AsymptoticMKResult(
            frequency_bins=list(bin_centers),
            alpha_by_freq=alpha_values,
            alpha_asymptotic=y_data[-1],
            ci_low=y_data[-1],
            ci_high=y_data[-1],
            dn=dn,
            ds=ds,
            num_genes=agg.num_genes,
            pn_total=pn_total,
            ps_total=ps_total,
        )

    if use_exponential:
        a, b, c = exp_popt
        return AsymptoticMKResult(
            frequency_bins=list(bin_centers),
            alpha_by_freq=alpha_values,
            alpha_asymptotic=exp_alpha,
            ci_low=ci_low_exp,
            ci_high=ci_high_exp,
            fit_a=float(a),
            fit_b=float(b),
            fit_c=float(c),
            dn=dn,
            ds=ds,
            num_genes=agg.num_genes,
            model_type="exponential",
            pn_total=pn_total,
            ps_total=ps_total,
        )
    else:
        a_lin, b_lin = lin_popt
        return AsymptoticMKResult(
            frequency_bins=list(bin_centers),
            alpha_by_freq=alpha_values,
            alpha_asymptotic=lin_alpha,
            ci_low=ci_low_lin,
            ci_high=ci_high_lin,
            fit_a=float(a_lin),
            fit_b=float(b_lin),
            fit_c=0.0,
            dn=dn,
            ds=ds,
            num_genes=agg.num_genes,
            model_type="linear",
            pn_total=pn_total,
            ps_total=ps_total,
        )


def asymptotic_mk_test(
    ingroup: SequenceSet | str | Path,
    outgroup: SequenceSet | str | Path,
    reading_frame: int = 1,
    genetic_code: GeneticCode | None = None,
    num_bins: int = 10,
    bootstrap_replicates: int = 100,
    pool_polymorphisms: bool = False,
) -> AsymptoticMKResult:
    """Perform the asymptotic McDonald-Kreitman test.

    This method accounts for slightly deleterious mutations by examining
    how alpha changes across the frequency spectrum and extrapolating
    to the asymptotic value at x=1.

    Based on Messer & Petrov (2013) PNAS.

    Args:
        ingroup: SequenceSet or path to FASTA file for ingroup sequences
        outgroup: SequenceSet or path to FASTA file for outgroup sequences
        reading_frame: Reading frame (1, 2, or 3)
        genetic_code: Genetic code for translation (uses standard if None)
        num_bins: Number of frequency bins (default 10)
        bootstrap_replicates: Number of bootstrap replicates for CI (default 100)
        pool_polymorphisms: If True, consider sites polymorphic in either
            population (libsequence convention). Frequencies are still
            calculated from ingroup only. If False (default), only consider
            ingroup polymorphisms (DnaSP/original MK convention).

    Returns:
        AsymptoticMKResult containing test statistics
    """
    code = genetic_code or DEFAULT_CODE

    # Load sequences if paths provided
    if isinstance(ingroup, (str, Path)):
        ingroup = SequenceSet.from_fasta(ingroup, reading_frame, code)
    if isinstance(outgroup, (str, Path)):
        outgroup = SequenceSet.from_fasta(outgroup, reading_frame, code)

    # Create aligned pair
    pair = AlignedPair(ingroup=ingroup, outgroup=outgroup, genetic_code=code)

    # Count divergence
    dn = 0
    ds = 0

    for codon_idx in range(pair.num_codons):
        if pair.is_fixed_between(codon_idx):
            result = pair.classify_fixed_difference(codon_idx)
            if result is not None:
                nonsyn, syn = result
                dn += nonsyn
                ds += syn

    # Collect polymorphisms with their frequencies
    poly_data: list[tuple[float, str]] = []  # (frequency, type: 'N' or 'S')

    # Get polymorphic sites based on mode
    if pool_polymorphisms:
        poly_sites = pair.polymorphic_sites_pooled()
    else:
        poly_sites = pair.polymorphic_sites_ingroup()

    for codon_idx in poly_sites:
        codons = list(ingroup.codon_set_clean(codon_idx))
        if len(codons) < 2:
            continue

        # Get frequency spectrum
        freqs = ingroup.site_frequency_spectrum(codon_idx)
        if not freqs:
            continue

        # Get outgroup codon to determine ancestral state
        out_codons = outgroup.codon_set_clean(codon_idx)
        if not out_codons:
            continue

        # Find ancestral codon: must be shared between ingroup and outgroup
        # to properly polarize the polymorphism
        ingroup_codons = set(freqs.keys())
        shared_codons = ingroup_codons & out_codons
        if not shared_codons:
            # No shared allele - can't determine ancestral state
            continue
        # Use the most frequent shared codon as ancestral
        ancestral = max(shared_codons, key=lambda c: freqs.get(c, 0))

        # Calculate derived allele frequency
        derived_freq = 1.0 - freqs[ancestral]

        if derived_freq <= 0 or derived_freq >= 1:
            continue

        # Classify the polymorphism
        result = pair.classify_polymorphism(codon_idx)
        if result is not None:
            nonsyn, syn = result
            for _ in range(nonsyn):
                poly_data.append((derived_freq, "N"))
            for _ in range(syn):
                poly_data.append((derived_freq, "S"))

    # Create frequency bins
    bin_edges = np.linspace(0, 1, num_bins + 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Bin polymorphisms by derived allele frequency
    pn_per_bin = [0] * num_bins
    ps_per_bin = [0] * num_bins

    for freq, poly_type in poly_data:
        # Find which bin this frequency falls into
        bin_idx = int(np.searchsorted(bin_edges[1:], freq, side="right"))
        bin_idx = min(bin_idx, num_bins - 1)
        if poly_type == "N":
            pn_per_bin[bin_idx] += 1
        else:
            ps_per_bin[bin_idx] += 1

    # Calculate per-bin alpha (matching asymptoticMK R package)
    # α(x) = 1 - (d0/d) * (p(x)/p0(x)) where p(x) and p0(x) are per-bin counts
    alpha_values = []
    valid_bins = []
    valid_centers = []

    for i in range(num_bins):
        pn_bin = pn_per_bin[i]
        ps_bin = ps_per_bin[i]

        if ds > 0 and dn > 0 and ps_bin > 0:
            alpha_x = 1.0 - (ds / dn) * (pn_bin / ps_bin)
            alpha_values.append(alpha_x)
            valid_bins.append(i)
            valid_centers.append(bin_centers[i])

    if len(valid_centers) < 3:
        # Not enough data points for curve fitting
        simple_alpha = None
        if ds > 0 and dn > 0:
            total_pn = sum(1 for _, t in poly_data if t == "N")
            total_ps = sum(1 for _, t in poly_data if t == "S")
            if total_ps > 0:
                simple_alpha = 1.0 - (ds * total_pn) / (dn * total_ps)

        return AsymptoticMKResult(
            frequency_bins=list(bin_centers),
            alpha_by_freq=alpha_values,
            alpha_asymptotic=simple_alpha or 0.0,
            ci_low=simple_alpha or 0.0,
            ci_high=simple_alpha or 0.0,
            dn=dn,
            ds=ds,
        )

    # Fit exponential model: α(x) = a + b * exp(-c * x)
    x_data = np.array(valid_centers)
    y_data = np.array(alpha_values)

    try:
        # Initial guesses: a is asymptotic value, b is negative for increasing alpha
        p0 = [y_data[-1], y_data[0] - y_data[-1], 5.0]

        # Bounds: b can be negative (typical for increasing alpha with frequency)
        bounds = ([-2.0, -2.0, 0.001], [2.0, 2.0, 100.0])

        popt, _ = optimize.curve_fit(
            _exponential_model,
            x_data,
            y_data,
            p0=p0,
            bounds=bounds,
            maxfev=10000,
        )

        a, b, c = popt
        # Evaluate the model at x=1: α(1) = a + b * exp(-c)
        alpha_asymp = float(a + b * np.exp(-c))

    except (RuntimeError, ValueError):
        # Curve fitting failed, use last value
        a, b, c = y_data[-1], 0.0, 1.0
        alpha_asymp = y_data[-1]

    # Bootstrap for confidence intervals
    bootstrap_alphas = []
    rng = np.random.default_rng(42)

    for _ in range(bootstrap_replicates):
        # Resample polymorphism data
        if len(poly_data) == 0:
            continue

        indices = rng.integers(0, len(poly_data), size=len(poly_data))
        boot_data = [poly_data[i] for i in indices]

        # Bin resampled polymorphisms
        boot_pn_per_bin = [0] * num_bins
        boot_ps_per_bin = [0] * num_bins

        for freq, poly_type in boot_data:
            bin_idx = int(np.searchsorted(bin_edges[1:], freq, side="right"))
            bin_idx = min(bin_idx, num_bins - 1)
            if poly_type == "N":
                boot_pn_per_bin[bin_idx] += 1
            else:
                boot_ps_per_bin[bin_idx] += 1

        # Recalculate per-bin alpha values
        boot_alpha_values = []
        boot_centers = []

        for i in range(num_bins):
            pn_bin = boot_pn_per_bin[i]
            ps_bin = boot_ps_per_bin[i]

            if ds > 0 and dn > 0 and ps_bin > 0:
                alpha_x = 1.0 - (ds / dn) * (pn_bin / ps_bin)
                boot_alpha_values.append(alpha_x)
                boot_centers.append(bin_centers[i])

        if len(boot_centers) >= 3:
            try:
                boot_x = np.array(boot_centers)
                boot_y = np.array(boot_alpha_values)
                popt_boot, _ = optimize.curve_fit(
                    _exponential_model,
                    boot_x,
                    boot_y,
                    p0=[a, b, c],
                    bounds=bounds,
                    maxfev=5000,
                )
                # Evaluate at x=1: α(1) = a + b * exp(-c)
                a_boot, b_boot, c_boot = popt_boot
                bootstrap_alphas.append(a_boot + b_boot * np.exp(-c_boot))
            except (RuntimeError, ValueError):
                bootstrap_alphas.append(boot_alpha_values[-1] if boot_alpha_values else alpha_asymp)

    # Calculate 95% CI
    if bootstrap_alphas:
        ci_low = float(np.percentile(bootstrap_alphas, 2.5))
        ci_high = float(np.percentile(bootstrap_alphas, 97.5))
    else:
        ci_low = alpha_asymp
        ci_high = alpha_asymp

    return AsymptoticMKResult(
        frequency_bins=list(bin_centers),
        alpha_by_freq=alpha_values,
        alpha_asymptotic=alpha_asymp,
        ci_low=ci_low,
        ci_high=ci_high,
        fit_a=float(a),
        fit_b=float(b),
        fit_c=float(c),
        dn=dn,
        ds=ds,
    )

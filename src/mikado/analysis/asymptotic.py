"""Asymptotic McDonald-Kreitman test implementation.

Based on Messer & Petrov (2013) PNAS.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from scipy import optimize

from mikado.core.alignment import AlignedPair
from mikado.core.codons import DEFAULT_CODE, GeneticCode
from mikado.core.sequences import SequenceSet

if TYPE_CHECKING:
    pass


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

    # Fit parameters (alpha = a - b*exp(-c*x))
    fit_a: float = 0.0
    fit_b: float = 0.0
    fit_c: float = 0.0

    # Counts used in calculation
    dn: int = 0
    ds: int = 0

    def __str__(self) -> str:
        return (
            f"Asymptotic MK Test Results:\n"
            f"  Asymptotic α: {self.alpha_asymptotic:.4f} "
            f"(95% CI: {self.ci_low:.4f} - {self.ci_high:.4f})\n"
            f"  Divergence: Dn={self.dn}, Ds={self.ds}\n"
            f"  Fit: α(x) = {self.fit_a:.4f} - {self.fit_b:.4f} * exp(-{self.fit_c:.4f} * x)"
        )

    def to_dict(self) -> dict:
        """Convert results to a dictionary."""
        return {
            "alpha_asymptotic": self.alpha_asymptotic,
            "ci_low": self.ci_low,
            "ci_high": self.ci_high,
            "dn": self.dn,
            "ds": self.ds,
            "fit_parameters": {
                "a": self.fit_a,
                "b": self.fit_b,
                "c": self.fit_c,
            },
            "frequency_bins": self.frequency_bins,
            "alpha_by_freq": self.alpha_by_freq,
        }


def _exponential_model(x: np.ndarray, a: float, b: float, c: float) -> np.ndarray:
    """Exponential model: α(x) = a - b * exp(-c * x)"""
    return a - b * np.exp(-c * x)


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

    # Calculate alpha for each frequency bin
    alpha_values = []
    valid_bins = []
    valid_centers = []

    for i in range(num_bins):
        high = bin_edges[i + 1]

        # Count polymorphisms up to this frequency
        pn_cumulative = sum(1 for f, t in poly_data if f <= high and t == "N")
        ps_cumulative = sum(1 for f, t in poly_data if f <= high and t == "S")

        if ds > 0 and dn > 0 and ps_cumulative > 0:
            # α(x) = 1 - (d0/d) * (p(x)/p0(x))
            # where d=Dn, d0=Ds, p(x)=Pn at freq≤x, p0(x)=Ps at freq≤x
            alpha_x = 1.0 - (ds / dn) * (pn_cumulative / ps_cumulative)
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

    # Fit exponential model
    x_data = np.array(valid_centers)
    y_data = np.array(alpha_values)

    try:
        # Initial parameter guesses
        p0 = [y_data[-1], y_data[-1] - y_data[0], 5.0]

        # Bounds to ensure reasonable fits
        bounds = ([-1.0, 0.0, 0.1], [1.0, 2.0, 50.0])

        popt, _ = optimize.curve_fit(
            _exponential_model,
            x_data,
            y_data,
            p0=p0,
            bounds=bounds,
            maxfev=10000,
        )

        a, b, c = popt
        # Evaluate the model at x=1 (not just 'a', since exp(-c) is not zero)
        alpha_asymp = float(a - b * np.exp(-c))

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

        # Recalculate alpha values
        boot_alpha_values = []
        boot_centers = []

        for i in range(num_bins):
            high = bin_edges[i + 1]
            pn_cum = sum(1 for f, t in boot_data if f <= high and t == "N")
            ps_cum = sum(1 for f, t in boot_data if f <= high and t == "S")

            if ds > 0 and dn > 0 and ps_cum > 0:
                alpha_x = 1.0 - (ds / dn) * (pn_cum / ps_cum)
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
                # Evaluate at x=1
                a_boot, b_boot, c_boot = popt_boot
                bootstrap_alphas.append(a_boot - b_boot * np.exp(-c_boot))
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

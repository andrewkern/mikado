"""DFE-based alpha estimation following GRAPES methodology.

Based on GRAPES implementation (Galtier 2016; Al-Saffar & Hahn 2022).
Reference: https://github.com/BioPP/grapes

This module provides methods for estimating alpha (proportion of adaptive
amino acid substitutions) using Distribution of Fitness Effects (DFE) models.
According to Al-Saffar & Hahn (2022), these methods significantly outperform
the asymptotic MK approach.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from mkado.analysis.asymptotic import (
    AggregatedSFS,
    PolymorphismData,
    extract_polymorphism_data,
    aggregate_polymorphism_data,
)
from mkado.analysis.dfe_models import (
    DFEInput,
    DFEResult,
    PrecomputedData,
    get_model,
    MODELS,
)
from mkado.core.codons import DEFAULT_CODE, GeneticCode
from mkado.core.sequences import SequenceSet

if TYPE_CHECKING:
    pass


__all__ = [
    "DFEInput",
    "DFEResult",
    "dfe_alpha",
    "dfe_alpha_aggregated",
    "polymorphism_data_to_dfe_input",
    "AVAILABLE_MODELS",
]


AVAILABLE_MODELS = list(MODELS.keys())


def _extract_sfs_from_polymorphisms(
    polymorphisms: list[tuple[float, str]],
    n_samples: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Extract folded SFS from polymorphism data.

    Args:
        polymorphisms: List of (frequency, type) tuples where type is 'N' or 'S'
        n_samples: Number of chromosomes sampled

    Returns:
        Tuple of (sfs_neutral, sfs_selected) as folded SFS arrays
    """
    # Create unfolded SFS bins (count at each frequency 1/n, 2/n, ..., (n-1)/n)
    sfs_neutral_unfolded = np.zeros(n_samples - 1)
    sfs_selected_unfolded = np.zeros(n_samples - 1)

    for freq, poly_type in polymorphisms:
        # Convert frequency to count
        count = int(round(freq * n_samples))
        if count < 1 or count >= n_samples:
            continue

        # Map to bin (0-indexed, so count 1 goes to bin 0)
        bin_idx = count - 1

        if poly_type == "S":
            sfs_neutral_unfolded[bin_idx] += 1
        else:  # "N"
            sfs_selected_unfolded[bin_idx] += 1

    # Fold the SFS
    folded_len = n_samples // 2
    sfs_neutral = np.zeros(folded_len)
    sfs_selected = np.zeros(folded_len)

    for j in range(1, folded_len + 1):
        complement = n_samples - j
        if j == complement:
            # At midpoint (n even)
            sfs_neutral[j - 1] = sfs_neutral_unfolded[j - 1]
            sfs_selected[j - 1] = sfs_selected_unfolded[j - 1]
        else:
            # Combine j and n-j
            sfs_neutral[j - 1] = (
                sfs_neutral_unfolded[j - 1] + sfs_neutral_unfolded[complement - 1]
            )
            sfs_selected[j - 1] = (
                sfs_selected_unfolded[j - 1] + sfs_selected_unfolded[complement - 1]
            )

    return sfs_neutral, sfs_selected


def polymorphism_data_to_dfe_input(
    poly_data: PolymorphismData,
    n_samples: int,
) -> DFEInput:
    """Convert PolymorphismData to DFE input format.

    Args:
        poly_data: Polymorphism data from extract_polymorphism_data()
        n_samples: Number of chromosomes sampled

    Returns:
        DFEInput ready for DFE model fitting
    """
    sfs_neutral, sfs_selected = _extract_sfs_from_polymorphisms(
        poly_data.polymorphisms, n_samples
    )

    return DFEInput(
        sfs_neutral=sfs_neutral,
        sfs_selected=sfs_selected,
        divergence_neutral=poly_data.ds,
        divergence_selected=poly_data.dn,
        n_samples=n_samples,
    )


def aggregated_sfs_to_dfe_input(
    agg_sfs: AggregatedSFS,
    n_samples: int,
) -> DFEInput:
    """Convert AggregatedSFS to DFE input format.

    Note: AggregatedSFS uses continuous frequency bins while DFE methods
    expect discrete count-based SFS. This conversion approximates by
    mapping bins to the nearest count class.

    Args:
        agg_sfs: Aggregated SFS from aggregate_polymorphism_data()
        n_samples: Number of chromosomes sampled

    Returns:
        DFEInput ready for DFE model fitting
    """
    # AggregatedSFS has continuous frequency bins
    # We need to convert to count-based SFS
    bin_edges = agg_sfs.bin_edges
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Create count-based SFS
    sfs_neutral_unfolded = np.zeros(n_samples - 1)
    sfs_selected_unfolded = np.zeros(n_samples - 1)

    for i, center in enumerate(bin_centers):
        # Map continuous frequency to discrete count
        count = int(round(center * n_samples))
        if count < 1 or count >= n_samples:
            continue

        bin_idx = count - 1
        sfs_neutral_unfolded[bin_idx] += agg_sfs.ps_counts[i]
        sfs_selected_unfolded[bin_idx] += agg_sfs.pn_counts[i]

    # Fold the SFS
    folded_len = n_samples // 2
    sfs_neutral = np.zeros(folded_len)
    sfs_selected = np.zeros(folded_len)

    for j in range(1, folded_len + 1):
        complement = n_samples - j
        if j == complement:
            sfs_neutral[j - 1] = sfs_neutral_unfolded[j - 1]
            sfs_selected[j - 1] = sfs_selected_unfolded[j - 1]
        else:
            sfs_neutral[j - 1] = (
                sfs_neutral_unfolded[j - 1] + sfs_neutral_unfolded[complement - 1]
            )
            sfs_selected[j - 1] = (
                sfs_selected_unfolded[j - 1] + sfs_selected_unfolded[complement - 1]
            )

    return DFEInput(
        sfs_neutral=sfs_neutral,
        sfs_selected=sfs_selected,
        divergence_neutral=agg_sfs.ds_total,
        divergence_selected=agg_sfs.dn_total,
        n_samples=n_samples,
    )


def dfe_alpha(
    ingroup: SequenceSet | str | Path,
    outgroup: SequenceSet | str | Path,
    reading_frame: int = 1,
    genetic_code: GeneticCode | None = None,
    model: str = "GammaExpo",
    pool_polymorphisms: bool = False,
) -> DFEResult:
    """Estimate alpha using DFE-based methods.

    This method fits a Distribution of Fitness Effects (DFE) model to the
    site frequency spectrum and uses it to estimate the proportion of
    adaptive amino acid substitutions (alpha).

    According to Al-Saffar & Hahn (2022), DFE-based methods significantly
    outperform the asymptotic MK approach, particularly the GammaExpo model.

    Args:
        ingroup: Ingroup sequences (SequenceSet or path to FASTA)
        outgroup: Outgroup sequences (SequenceSet or path to FASTA)
        reading_frame: Reading frame (1, 2, or 3)
        genetic_code: Genetic code for translation (uses standard if None)
        model: DFE model to use (GammaZero, GammaExpo, GammaGamma, DisplacedGamma)
        pool_polymorphisms: If True, consider sites polymorphic in either population

    Returns:
        DFEResult with alpha estimate and model parameters
    """
    code = genetic_code or DEFAULT_CODE

    # Load sequences if paths provided
    if isinstance(ingroup, (str, Path)):
        ingroup = SequenceSet.from_fasta(ingroup, reading_frame, code)
    if isinstance(outgroup, (str, Path)):
        outgroup = SequenceSet.from_fasta(outgroup, reading_frame, code)

    # Get sample size
    n_samples = len(ingroup)

    # Extract polymorphism data using existing function
    poly_data = extract_polymorphism_data(
        ingroup=ingroup,
        outgroup=outgroup,
        reading_frame=reading_frame,
        genetic_code=code,
        pool_polymorphisms=pool_polymorphisms,
    )

    # Convert to DFE input
    dfe_input = polymorphism_data_to_dfe_input(poly_data, n_samples)

    # Pre-compute integration grids
    precalc = PrecomputedData(n_samples)

    # Fit model
    dfe_model = get_model(model)
    result = dfe_model.fit(dfe_input, precalc)

    return result


def dfe_alpha_aggregated(
    gene_data: list[PolymorphismData],
    n_samples: int,
    model: str = "GammaExpo",
    num_bins: int = 20,
) -> DFEResult:
    """Estimate alpha from pre-aggregated multi-gene data.

    This function aggregates SFS data from multiple genes and then
    fits a DFE model to estimate alpha.

    Args:
        gene_data: List of PolymorphismData from individual genes
        n_samples: Number of chromosomes sampled (must be consistent across genes)
        model: DFE model to use
        num_bins: Number of frequency bins for aggregation

    Returns:
        DFEResult with alpha estimate and model parameters
    """
    # Aggregate polymorphism data
    agg = aggregate_polymorphism_data(gene_data, num_bins)

    # Convert to DFE input
    dfe_input = aggregated_sfs_to_dfe_input(agg, n_samples)

    # Pre-compute integration grids
    precalc = PrecomputedData(n_samples)

    # Fit model
    dfe_model = get_model(model)
    result = dfe_model.fit(dfe_input, precalc)

    return result


def dfe_alpha_from_sfs(
    sfs_neutral: np.ndarray,
    sfs_selected: np.ndarray,
    dn: int,
    ds: int,
    n_samples: int,
    model: str = "GammaExpo",
    n_sites_neutral: float | None = None,
    n_sites_selected: float | None = None,
    n_sites_div_neutral: float | None = None,
    n_sites_div_selected: float | None = None,
) -> DFEResult:
    """Estimate alpha from pre-computed SFS data.

    This function allows direct input of SFS data without going through
    sequence alignment, useful for data from external sources.

    Args:
        sfs_neutral: Folded SFS for synonymous sites
        sfs_selected: Folded SFS for nonsynonymous sites
        dn: Nonsynonymous divergence count
        ds: Synonymous divergence count
        n_samples: Number of chromosomes sampled
        model: DFE model to use
        n_sites_neutral: Number of synonymous sites for polymorphism (estimated if None)
        n_sites_selected: Number of nonsynonymous sites for polymorphism (estimated if None)
        n_sites_div_neutral: Number of synonymous sites for divergence (uses n_sites_neutral if None)
        n_sites_div_selected: Number of nonsynonymous sites for divergence (uses n_sites_selected if None)

    Returns:
        DFEResult with alpha estimate and model parameters
    """
    dfe_input = DFEInput(
        sfs_neutral=sfs_neutral,
        sfs_selected=sfs_selected,
        divergence_neutral=ds,
        divergence_selected=dn,
        n_samples=n_samples,
        n_sites_neutral=n_sites_neutral,
        n_sites_selected=n_sites_selected,
        n_sites_div_neutral=n_sites_div_neutral,
        n_sites_div_selected=n_sites_div_selected,
    )

    # Pre-compute integration grids
    precalc = PrecomputedData(n_samples)

    # Fit model
    dfe_model = get_model(model)
    result = dfe_model.fit(dfe_input, precalc)

    return result


def compare_models(
    ingroup: SequenceSet | str | Path,
    outgroup: SequenceSet | str | Path,
    reading_frame: int = 1,
    genetic_code: GeneticCode | None = None,
    models: list[str] | None = None,
    pool_polymorphisms: bool = False,
) -> list[DFEResult]:
    """Fit multiple DFE models and compare them by AIC.

    Args:
        ingroup: Ingroup sequences
        outgroup: Outgroup sequences
        reading_frame: Reading frame (1, 2, or 3)
        genetic_code: Genetic code for translation
        models: List of model names to compare (all if None)
        pool_polymorphisms: If True, consider sites polymorphic in either population

    Returns:
        List of DFEResult sorted by AIC (best first)
    """
    if models is None:
        models = list(MODELS.keys())

    code = genetic_code or DEFAULT_CODE

    # Load sequences if paths provided
    if isinstance(ingroup, (str, Path)):
        ingroup = SequenceSet.from_fasta(ingroup, reading_frame, code)
    if isinstance(outgroup, (str, Path)):
        outgroup = SequenceSet.from_fasta(outgroup, reading_frame, code)

    # Get sample size
    n_samples = len(ingroup)

    # Extract polymorphism data
    poly_data = extract_polymorphism_data(
        ingroup=ingroup,
        outgroup=outgroup,
        reading_frame=reading_frame,
        genetic_code=code,
        pool_polymorphisms=pool_polymorphisms,
    )

    # Convert to DFE input
    dfe_input = polymorphism_data_to_dfe_input(poly_data, n_samples)

    # Pre-compute integration grids (shared across models)
    precalc = PrecomputedData(n_samples)

    # Fit all models
    results = []
    for model_name in models:
        dfe_model = get_model(model_name)
        result = dfe_model.fit(dfe_input, precalc)
        results.append(result)

    # Sort by AIC
    results.sort(key=lambda r: r.aic)

    return results

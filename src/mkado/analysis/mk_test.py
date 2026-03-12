"""Standard McDonald-Kreitman test implementation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from mkado.analysis.statistics import alpha, dos, fishers_exact, neutrality_index
from mkado.core.alignment import AlignedPair
from mkado.core.codons import DEFAULT_CODE, GeneticCode
from mkado.core.sequences import SequenceSet

if TYPE_CHECKING:
    pass


@dataclass
class MKResult:
    """Results from a McDonald-Kreitman test."""

    dn: int  # Non-synonymous divergence (fixed differences)
    ds: int  # Synonymous divergence (fixed differences)
    pn: int  # Non-synonymous polymorphisms
    ps: int  # Synonymous polymorphisms
    p_value: float  # Fisher's exact test p-value
    ni: float | None  # Neutrality Index
    alpha: float | None  # Proportion of adaptive substitutions
    dos: float | None  # Direction of Selection

    def __str__(self) -> str:
        ni_str = f"{self.ni:.4f}" if self.ni is not None else "N/A"
        alpha_str = f"{self.alpha:.4f}" if self.alpha is not None else "N/A"
        dos_str = f"{self.dos:.4f}" if self.dos is not None else "N/A"
        return (
            f"MK Test Results:\n"
            f"  Divergence:    Dn={self.dn}, Ds={self.ds}\n"
            f"  Polymorphism:  Pn={self.pn}, Ps={self.ps}\n"
            f"  Fisher's exact p-value: {self.p_value:.4g}\n"
            f"  Neutrality Index (NI):  {ni_str}\n"
            f"  Alpha (α):              {alpha_str}\n"
            f"  DoS:                    {dos_str}"
        )

    def to_dict(self) -> dict:
        """Convert results to a dictionary."""
        return {
            "dn": self.dn,
            "ds": self.ds,
            "pn": self.pn,
            "ps": self.ps,
            "p_value": self.p_value,
            "ni": self.ni,
            "alpha": self.alpha,
            "dos": self.dos,
        }


def mk_test(
    ingroup: SequenceSet | str | Path,
    outgroup: SequenceSet | str | Path,
    reading_frame: int = 1,
    genetic_code: GeneticCode | None = None,
    pool_polymorphisms: bool = False,
    min_frequency: float = 0.0,
) -> MKResult:
    """Perform the standard McDonald-Kreitman test.

    The MK test compares the ratio of non-synonymous to synonymous changes
    within species (polymorphism) to the ratio between species (divergence).

    Under neutrality, these ratios should be equal. Deviations suggest selection:
    - Excess divergence relative to polymorphism suggests positive selection
    - Excess polymorphism relative to divergence suggests segregating weakly deleterious polymorphisms

    Args:
        ingroup: SequenceSet or path to FASTA file for ingroup sequences
        outgroup: SequenceSet or path to FASTA file for outgroup sequences
        reading_frame: Reading frame (1, 2, or 3)
        genetic_code: Genetic code for translation (uses standard if None)
        pool_polymorphisms: If True, count polymorphisms from both ingroup
            and outgroup (libsequence convention). If False (default), count
            only ingroup polymorphisms (DnaSP/original MK convention).
        min_frequency: Minimum derived allele frequency for a polymorphism
            to be counted (default 0.0, i.e., include all polymorphisms).
            Useful for filtering out rare variants that may be slightly
            deleterious.

    Returns:
        MKResult containing test statistics
    """
    code = genetic_code or DEFAULT_CODE

    # Load sequences if paths provided
    if isinstance(ingroup, (str, Path)):
        ingroup = SequenceSet.from_fasta(ingroup, reading_frame, code)
    if isinstance(outgroup, (str, Path)):
        outgroup = SequenceSet.from_fasta(outgroup, reading_frame, code)

    # Create aligned pair
    pair = AlignedPair(ingroup=ingroup, outgroup=outgroup, genetic_code=code)

    # Count divergence (fixed differences)
    dn = 0
    ds = 0
    processed_codons: set[int] = set()

    for codon_idx in range(pair.num_codons):
        if pair.is_fixed_between(codon_idx):
            if codon_idx in processed_codons:
                continue

            result = pair.classify_fixed_difference(codon_idx)
            if result is not None:
                nonsyn, syn = result
                dn += nonsyn
                ds += syn
                processed_codons.add(codon_idx)

    # Count polymorphisms
    pn = 0
    ps = 0

    if pool_polymorphisms:
        # Pooled mode: count polymorphisms from both populations (libsequence convention)
        poly_sites = pair.polymorphic_sites_pooled()
        classify_func = pair.classify_polymorphism_pooled
    else:
        # Standard mode: count only ingroup polymorphisms (DnaSP convention)
        poly_sites = pair.polymorphic_sites_ingroup()
        classify_func = pair.classify_polymorphism

    for codon_idx in poly_sites:
        # Skip if also a fixed difference (shouldn't happen, but be safe)
        if codon_idx in processed_codons:
            continue

        # Apply frequency filter if specified
        if min_frequency > 0:
            # Get frequency spectrum from ingroup
            freqs = ingroup.site_frequency_spectrum(codon_idx)
            if not freqs:
                continue

            # Get outgroup codons to determine ancestral state
            out_codons = outgroup.codon_set_clean(codon_idx)
            if not out_codons:
                continue

            # Find ancestral codon (shared between ingroup and outgroup)
            ingroup_codons = set(freqs.keys())
            shared_codons = ingroup_codons & out_codons
            if not shared_codons:
                continue

            # Use the most frequent shared codon as ancestral
            ancestral = max(shared_codons, key=lambda c: freqs.get(c, 0))

            # Calculate derived allele frequency
            derived_freq = 1.0 - freqs[ancestral]

            # Skip if below minimum frequency threshold
            if derived_freq < min_frequency + 1e-10:
                continue

        result = classify_func(codon_idx)
        if result is not None:
            nonsyn, syn = result
            pn += nonsyn
            ps += syn

    # Calculate statistics
    p_val = fishers_exact(dn, ds, pn, ps)
    ni = neutrality_index(dn, ds, pn, ps)
    a = alpha(dn, ds, pn, ps)
    d = dos(dn, ds, pn, ps)

    return MKResult(
        dn=dn,
        ds=ds,
        pn=pn,
        ps=ps,
        p_value=p_val,
        ni=ni,
        alpha=a,
        dos=d,
    )


def mk_test_from_counts(
    dn: int, ds: int, pn: int, ps: int
) -> MKResult:
    """Create MK test results from pre-computed counts.

    Useful for testing or when counts are already available.

    Args:
        dn: Non-synonymous divergence
        ds: Synonymous divergence
        pn: Non-synonymous polymorphisms
        ps: Synonymous polymorphisms

    Returns:
        MKResult with calculated statistics
    """
    p_val = fishers_exact(dn, ds, pn, ps)
    ni = neutrality_index(dn, ds, pn, ps)
    a = alpha(dn, ds, pn, ps)
    d = dos(dn, ds, pn, ps)

    return MKResult(
        dn=dn,
        ds=ds,
        pn=pn,
        ps=ps,
        p_value=p_val,
        ni=ni,
        alpha=a,
        dos=d,
    )

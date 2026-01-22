"""Polarized McDonald-Kreitman test implementation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from mkado.analysis.statistics import alpha, fishers_exact, neutrality_index
from mkado.core.alignment import PolarizedAlignedPair
from mkado.core.codons import DEFAULT_CODE, GeneticCode
from mkado.core.sequences import SequenceSet

if TYPE_CHECKING:
    pass


@dataclass
class PolarizedMKResult:
    """Results from a polarized McDonald-Kreitman test.

    The polarized test uses a second outgroup to assign mutations to
    specific lineages (ingroup or first outgroup).
    """

    # Ingroup lineage counts
    dn_ingroup: int
    ds_ingroup: int
    pn_ingroup: int
    ps_ingroup: int

    # Outgroup lineage counts
    dn_outgroup: int
    ds_outgroup: int

    # Unpolarized counts (could not determine lineage)
    dn_unpolarized: int
    ds_unpolarized: int

    # Statistics for ingroup lineage
    p_value_ingroup: float
    ni_ingroup: float | None
    alpha_ingroup: float | None

    def __str__(self) -> str:
        ni_str = f"{self.ni_ingroup:.4f}" if self.ni_ingroup is not None else "N/A"
        alpha_str = f"{self.alpha_ingroup:.4f}" if self.alpha_ingroup is not None else "N/A"
        return (
            f"Polarized MK Test Results:\n"
            f"  Ingroup Lineage:\n"
            f"    Divergence:    Dn={self.dn_ingroup}, Ds={self.ds_ingroup}\n"
            f"    Polymorphism:  Pn={self.pn_ingroup}, Ps={self.ps_ingroup}\n"
            f"    Fisher's p-value: {self.p_value_ingroup:.4g}\n"
            f"    NI: {ni_str}, Alpha: {alpha_str}\n"
            f"  Outgroup Lineage:\n"
            f"    Divergence:    Dn={self.dn_outgroup}, Ds={self.ds_outgroup}\n"
            f"  Unpolarized:\n"
            f"    Divergence:    Dn={self.dn_unpolarized}, Ds={self.ds_unpolarized}"
        )

    def to_dict(self) -> dict:
        """Convert results to a dictionary."""
        return {
            "ingroup": {
                "dn": self.dn_ingroup,
                "ds": self.ds_ingroup,
                "pn": self.pn_ingroup,
                "ps": self.ps_ingroup,
                "p_value": self.p_value_ingroup,
                "ni": self.ni_ingroup,
                "alpha": self.alpha_ingroup,
            },
            "outgroup": {
                "dn": self.dn_outgroup,
                "ds": self.ds_outgroup,
            },
            "unpolarized": {
                "dn": self.dn_unpolarized,
                "ds": self.ds_unpolarized,
            },
        }


def polarized_mk_test(
    ingroup: SequenceSet | str | Path,
    outgroup1: SequenceSet | str | Path,
    outgroup2: SequenceSet | str | Path,
    reading_frame: int = 1,
    genetic_code: GeneticCode | None = None,
    pool_polymorphisms: bool = False,
    min_frequency: float = 0.0,
) -> PolarizedMKResult:
    """Perform a polarized McDonald-Kreitman test.

    Uses a second outgroup to determine the ancestral state and assign
    mutations to specific lineages.

    Args:
        ingroup: SequenceSet or path to FASTA file for ingroup sequences
        outgroup1: SequenceSet or path to FASTA file for first outgroup
        outgroup2: SequenceSet or path to FASTA file for second outgroup
            (used for polarization)
        reading_frame: Reading frame (1, 2, or 3)
        genetic_code: Genetic code for translation (uses standard if None)
        pool_polymorphisms: If True, count polymorphisms from both ingroup
            and outgroup1 (libsequence convention). If False (default), count
            only ingroup polymorphisms (DnaSP/original MK convention).
        min_frequency: Minimum derived allele frequency for a polymorphism
            to be counted (default 0.0, i.e., include all polymorphisms).
            Useful for filtering out rare variants that may be slightly
            deleterious.

    Returns:
        PolarizedMKResult containing test statistics
    """
    code = genetic_code or DEFAULT_CODE

    # Load sequences if paths provided
    if isinstance(ingroup, (str, Path)):
        ingroup = SequenceSet.from_fasta(ingroup, reading_frame, code)
    if isinstance(outgroup1, (str, Path)):
        outgroup1 = SequenceSet.from_fasta(outgroup1, reading_frame, code)
    if isinstance(outgroup2, (str, Path)):
        outgroup2 = SequenceSet.from_fasta(outgroup2, reading_frame, code)

    # Create polarized aligned pair
    pair = PolarizedAlignedPair(
        ingroup=ingroup,
        outgroup=outgroup1,
        outgroup2=outgroup2,
        genetic_code=code,
    )

    # Count divergence with polarization
    dn_in = 0
    ds_in = 0
    dn_out = 0
    ds_out = 0
    dn_unpol = 0
    ds_unpol = 0

    for codon_idx in range(pair.num_codons):
        if pair.is_fixed_between(codon_idx):
            polarization = pair.polarize_fixed_difference(codon_idx)

            if polarization is None:
                # Cannot polarize - count as unpolarized
                result = pair.classify_fixed_difference(codon_idx)
                if result is not None:
                    nonsyn, syn = result
                    dn_unpol += nonsyn
                    ds_unpol += syn
            else:
                lineage, (nonsyn, syn) = polarization
                if lineage == "ingroup":
                    dn_in += nonsyn
                    ds_in += syn
                else:
                    dn_out += nonsyn
                    ds_out += syn

    # Count polymorphisms
    pn_in = 0
    ps_in = 0

    if pool_polymorphisms:
        # Pooled mode: count polymorphisms from both populations (libsequence convention)
        poly_sites = pair.polymorphic_sites_pooled()
        classify_func = pair.classify_polymorphism_pooled
    else:
        # Standard mode: count only ingroup polymorphisms (DnaSP convention)
        poly_sites = pair.polymorphic_sites_ingroup()
        classify_func = pair.classify_polymorphism

    for codon_idx in poly_sites:
        # Apply frequency filter if specified
        if min_frequency > 0:
            # Get frequency spectrum from ingroup
            freqs = ingroup.site_frequency_spectrum(codon_idx)
            if not freqs:
                continue

            # Get outgroup2 codons to determine ancestral state
            out_codons = outgroup2.codon_set_clean(codon_idx)
            if not out_codons:
                continue

            # Find ancestral codon (shared between ingroup and outgroup2)
            ingroup_codons = set(freqs.keys())
            shared_codons = ingroup_codons & out_codons
            if not shared_codons:
                continue

            # Use the most frequent shared codon as ancestral
            ancestral = max(shared_codons, key=lambda c: freqs.get(c, 0))

            # Calculate derived allele frequency
            derived_freq = 1.0 - freqs[ancestral]

            # Skip if below minimum frequency threshold
            if derived_freq < min_frequency:
                continue

        result = classify_func(codon_idx)
        if result is not None:
            nonsyn, syn = result
            pn_in += nonsyn
            ps_in += syn

    # Calculate statistics for ingroup lineage
    p_val = fishers_exact(dn_in, ds_in, pn_in, ps_in)
    ni = neutrality_index(dn_in, ds_in, pn_in, ps_in)
    a = alpha(dn_in, ds_in, pn_in, ps_in)

    return PolarizedMKResult(
        dn_ingroup=dn_in,
        ds_ingroup=ds_in,
        pn_ingroup=pn_in,
        ps_ingroup=ps_in,
        dn_outgroup=dn_out,
        ds_outgroup=ds_out,
        dn_unpolarized=dn_unpol,
        ds_unpolarized=ds_unpol,
        p_value_ingroup=p_val,
        ni_ingroup=ni,
        alpha_ingroup=a,
    )

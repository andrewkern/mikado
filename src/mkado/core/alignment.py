"""Alignment comparison utilities for MK test."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from mkado.core.codons import DEFAULT_CODE, GeneticCode
from mkado.core.sequences import SequenceSet

if TYPE_CHECKING:
    pass


@dataclass
class AlignedPair:
    """Represents a pair of aligned sequence sets (ingroup and outgroup)."""

    ingroup: SequenceSet
    outgroup: SequenceSet
    genetic_code: GeneticCode = field(default_factory=lambda: DEFAULT_CODE)

    def __post_init__(self) -> None:
        if self.ingroup.num_codons != self.outgroup.num_codons:
            raise ValueError(
                f"Alignment length mismatch: ingroup has {self.ingroup.num_codons} codons, "
                f"outgroup has {self.outgroup.num_codons} codons"
            )

    @property
    def num_codons(self) -> int:
        """Number of codons in the alignment."""
        return self.ingroup.num_codons

    def combined_codon_set(self, codon_index: int) -> set[str]:
        """Get all unique codons at a position from both groups.

        Args:
            codon_index: Zero-based codon index

        Returns:
            Set of unique codon strings from both ingroup and outgroup
        """
        return self.ingroup.codon_set(codon_index) | self.outgroup.codon_set(codon_index)

    def combined_codon_set_clean(self, codon_index: int) -> set[str]:
        """Get all clean unique codons at a position from both groups.

        Args:
            codon_index: Zero-based codon index

        Returns:
            Set of unique valid codon strings
        """
        return self.ingroup.codon_set_clean(codon_index) | self.outgroup.codon_set_clean(
            codon_index
        )

    def is_fixed_between(self, codon_index: int) -> bool:
        """Check if a codon position is fixed differently between groups.

        A position is a fixed difference if:
        - The ingroup has a single codon (or all share the same)
        - The outgroup has a single codon (or all share the same)
        - The ingroup and outgroup codons are different

        Args:
            codon_index: Zero-based codon index

        Returns:
            True if this is a fixed difference between groups
        """
        in_codons = self.ingroup.codon_set_clean(codon_index)
        out_codons = self.outgroup.codon_set_clean(codon_index)

        if not in_codons or not out_codons:
            return False

        # Fixed if each group has one codon and they're different
        if len(in_codons) == 1 and len(out_codons) == 1:
            return in_codons != out_codons

        return False

    def is_polymorphic_within_ingroup(self, codon_index: int) -> bool:
        """Check if a codon position is polymorphic within the ingroup only.

        Args:
            codon_index: Zero-based codon index

        Returns:
            True if polymorphic within ingroup but not a fixed difference
        """
        return self.ingroup.is_polymorphic(codon_index)

    def fixed_differences(self) -> list[int]:
        """Get all codon indices with fixed differences between groups.

        Returns:
            List of codon indices
        """
        return [i for i in range(self.num_codons) if self.is_fixed_between(i)]

    def polymorphic_sites_ingroup(self) -> list[int]:
        """Get all codon indices polymorphic within the ingroup.

        Returns:
            List of codon indices
        """
        return self.ingroup.polymorphic_codons()

    def polymorphic_sites_outgroup(self) -> list[int]:
        """Get all codon indices polymorphic within the outgroup.

        Returns:
            List of codon indices
        """
        return self.outgroup.polymorphic_codons()

    def polymorphic_sites_pooled(self) -> list[int]:
        """Get all codon indices polymorphic in either population (union).

        This follows the libsequence convention of pooling polymorphisms
        from both populations.

        Returns:
            List of unique codon indices (sorted)
        """
        ingroup_poly = set(self.ingroup.polymorphic_codons())
        outgroup_poly = set(self.outgroup.polymorphic_codons())
        return sorted(ingroup_poly | outgroup_poly)

    def classify_fixed_difference(
        self, codon_index: int
    ) -> tuple[int, int] | None:
        """Classify a fixed difference as synonymous/non-synonymous.

        Uses the shortest mutational path to count the minimum number of
        synonymous and non-synonymous changes.

        Args:
            codon_index: Zero-based codon index

        Returns:
            Tuple of (non_synonymous_count, synonymous_count), or None if
            not a valid fixed difference
        """
        in_codons = self.ingroup.codon_set_clean(codon_index)
        out_codons = self.outgroup.codon_set_clean(codon_index)

        if not in_codons or not out_codons:
            return None

        # Get representative codons
        in_codon = next(iter(in_codons))
        out_codon = next(iter(out_codons))

        if in_codon == out_codon:
            return None

        path = self.genetic_code.get_path(in_codon, out_codon)
        if not path:
            return None

        nonsyn = sum(1 for change_type, _ in path if change_type == "R")
        syn = sum(1 for change_type, _ in path if change_type == "S")

        return (nonsyn, syn)

    def classify_polymorphism(
        self, codon_index: int
    ) -> tuple[int, int] | None:
        """Classify a polymorphism as synonymous/non-synonymous.

        Args:
            codon_index: Zero-based codon index

        Returns:
            Tuple of (non_synonymous_count, synonymous_count), or None if
            not a valid polymorphism
        """
        codons = list(self.ingroup.codon_set_clean(codon_index))

        if len(codons) < 2:
            return None

        if len(codons) == 2:
            path = self.genetic_code.get_path(codons[0], codons[1])
            if not path:
                return None
            nonsyn = sum(1 for change_type, _ in path if change_type == "R")
            syn = sum(1 for change_type, _ in path if change_type == "S")
            return (nonsyn, syn)

        # For >2 codons, find shortest paths between all pairs
        total_nonsyn = 0
        total_syn = 0
        counted_positions: set[int] = set()

        # Use a simple approach: compare each codon to the most common one
        freqs = self.ingroup.site_frequency_spectrum(codon_index)
        if not freqs:
            return None

        major_codon = max(freqs.keys(), key=lambda c: freqs[c])
        for codon in codons:
            if codon == major_codon:
                continue
            path = self.genetic_code.get_path(major_codon, codon)
            if path:
                for change_type, pos in path:
                    if pos not in counted_positions:
                        if change_type == "R":
                            total_nonsyn += 1
                        else:
                            total_syn += 1
                        counted_positions.add(pos)

        return (total_nonsyn, total_syn)

    def classify_polymorphism_pooled(
        self, codon_index: int
    ) -> tuple[int, int] | None:
        """Classify a polymorphism using codons from both populations.

        This follows the libsequence convention of using all unique codons
        from both ingroup and outgroup when classifying polymorphisms.

        Args:
            codon_index: Zero-based codon index

        Returns:
            Tuple of (non_synonymous_count, synonymous_count), or None if
            not a valid polymorphism
        """
        # Get unique codons from both populations
        codons = list(self.combined_codon_set_clean(codon_index))

        if len(codons) < 2:
            return None

        if len(codons) == 2:
            path = self.genetic_code.get_path(codons[0], codons[1])
            if not path:
                return None
            nonsyn = sum(1 for change_type, _ in path if change_type == "R")
            syn = sum(1 for change_type, _ in path if change_type == "S")
            return (nonsyn, syn)

        # For >2 codons, find shortest paths from most common to others
        # Use combined frequency spectrum
        total_nonsyn = 0
        total_syn = 0
        counted_positions: set[int] = set()

        # Get frequencies from both populations
        in_freqs = self.ingroup.site_frequency_spectrum(codon_index)
        out_freqs = self.outgroup.site_frequency_spectrum(codon_index)

        # Merge frequencies
        combined_freqs: dict[str, int] = {}
        for codon, count in (in_freqs or {}).items():
            combined_freqs[codon] = combined_freqs.get(codon, 0) + count
        for codon, count in (out_freqs or {}).items():
            combined_freqs[codon] = combined_freqs.get(codon, 0) + count

        if not combined_freqs:
            return None

        major_codon = max(combined_freqs.keys(), key=lambda c: combined_freqs[c])
        for codon in codons:
            if codon == major_codon:
                continue
            path = self.genetic_code.get_path(major_codon, codon)
            if path:
                for change_type, pos in path:
                    if pos not in counted_positions:
                        if change_type == "R":
                            total_nonsyn += 1
                        else:
                            total_syn += 1
                        counted_positions.add(pos)

        return (total_nonsyn, total_syn)


@dataclass
class PolarizedAlignedPair(AlignedPair):
    """Aligned pair with a second outgroup for polarization."""

    outgroup2: SequenceSet | None = None

    def polarize_fixed_difference(
        self, codon_index: int
    ) -> tuple[str, tuple[int, int]] | None:
        """Polarize a fixed difference to determine which lineage changed.

        Uses the second outgroup to determine the ancestral state.

        Args:
            codon_index: Zero-based codon index

        Returns:
            Tuple of (lineage, (nonsyn, syn)) where lineage is 'ingroup' or
            'outgroup', or None if cannot be polarized
        """
        if self.outgroup2 is None:
            return None

        in_codons = self.ingroup.codon_set_clean(codon_index)
        out_codons = self.outgroup.codon_set_clean(codon_index)
        out2_codons = self.outgroup2.codon_set_clean(codon_index)

        if not in_codons or not out_codons or not out2_codons:
            return None

        in_codon = next(iter(in_codons))
        out_codon = next(iter(out_codons))
        out2_codon = next(iter(out2_codons))

        # Determine ancestral state
        if out_codon == out2_codon:
            # Outgroup agrees - ingroup changed
            ancestral = out_codon
            derived = in_codon
            lineage = "ingroup"
        elif in_codon == out2_codon:
            # Ingroup matches outgroup2 - outgroup1 changed
            ancestral = in_codon
            derived = out_codon
            lineage = "outgroup"
        else:
            # Cannot polarize
            return None

        if ancestral == derived:
            return None

        path = self.genetic_code.get_path(ancestral, derived)
        if not path:
            return None

        nonsyn = sum(1 for change_type, _ in path if change_type == "R")
        syn = sum(1 for change_type, _ in path if change_type == "S")

        return (lineage, (nonsyn, syn))

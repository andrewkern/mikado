"""Sequence and SequenceSet classes for handling aligned sequences."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from mkado.core.codons import DEFAULT_CODE, GeneticCode
from mkado.io.fasta import read_fasta

if TYPE_CHECKING:
    pass


@dataclass
class Sequence:
    """Represents a single biological sequence."""

    name: str
    sequence: str

    def __len__(self) -> int:
        return len(self.sequence)

    def __getitem__(self, index: int | slice) -> str:
        return self.sequence[index]

    def get_codon(self, codon_index: int, reading_frame: int = 1) -> str:
        """Get the codon at a given index.

        Args:
            codon_index: Zero-based codon index
            reading_frame: Reading frame (1, 2, or 3)

        Returns:
            Three-letter codon string
        """
        start = (reading_frame - 1) + (codon_index * 3)
        return self.sequence[start : start + 3]

    def num_codons(self, reading_frame: int = 1) -> int:
        """Get the number of complete codons in the sequence."""
        return (len(self.sequence) - (reading_frame - 1)) // 3


@dataclass
class SequenceSet:
    """Represents a set of aligned sequences (e.g., from one species)."""

    sequences: list[Sequence] = field(default_factory=list)
    reading_frame: int = 1
    genetic_code: GeneticCode = field(default_factory=lambda: DEFAULT_CODE)

    @classmethod
    def from_fasta(
        cls,
        path: str | Path,
        reading_frame: int = 1,
        genetic_code: GeneticCode | None = None,
    ) -> SequenceSet:
        """Load sequences from a FASTA file.

        Args:
            path: Path to FASTA file
            reading_frame: Reading frame (1, 2, or 3)
            genetic_code: Genetic code for translation

        Returns:
            SequenceSet containing all sequences from the file
        """
        sequences = [Sequence(name, seq) for name, seq in read_fasta(path)]
        return cls(
            sequences=sequences,
            reading_frame=reading_frame,
            genetic_code=genetic_code or DEFAULT_CODE,
        )

    def __len__(self) -> int:
        return len(self.sequences)

    def __getitem__(self, index: int) -> Sequence:
        return self.sequences[index]

    @property
    def num_codons(self) -> int:
        """Number of complete codons (based on first sequence)."""
        if not self.sequences:
            return 0
        return self.sequences[0].num_codons(self.reading_frame)

    @property
    def alignment_length(self) -> int:
        """Length of the alignment (based on first sequence)."""
        if not self.sequences:
            return 0
        return len(self.sequences[0])

    def codon_set(self, codon_index: int) -> set[str]:
        """Get the set of unique codons at a position across all sequences.

        Args:
            codon_index: Zero-based codon index

        Returns:
            Set of unique codon strings
        """
        codons = set()
        for seq in self.sequences:
            codon = seq.get_codon(codon_index, self.reading_frame)
            codons.add(codon)
        return codons

    def codon_set_clean(self, codon_index: int) -> set[str]:
        """Get unique codons at a position, excluding those with N or gaps.

        Args:
            codon_index: Zero-based codon index

        Returns:
            Set of unique valid codon strings
        """
        codons = self.codon_set(codon_index)
        return {c for c in codons if "N" not in c and "-" not in c}

    def amino_set(self, codon_index: int) -> set[str]:
        """Get the set of unique amino acids at a codon position.

        Args:
            codon_index: Zero-based codon index

        Returns:
            Set of unique amino acid characters
        """
        codons = self.codon_set(codon_index)
        return {self.genetic_code.translate(c) for c in codons}

    def amino_set_clean(self, codon_index: int) -> set[str]:
        """Get unique amino acids at a position, excluding ambiguous codons.

        Args:
            codon_index: Zero-based codon index

        Returns:
            Set of unique amino acid characters
        """
        codons = self.codon_set_clean(codon_index)
        return {self.genetic_code.translate(c) for c in codons}

    def is_polymorphic(self, codon_index: int) -> bool:
        """Check if a codon position is polymorphic (has >1 unique codon).

        Args:
            codon_index: Zero-based codon index

        Returns:
            True if polymorphic
        """
        codons = self.codon_set_clean(codon_index)
        return len(codons) > 1

    def site_set(self, site_index: int) -> set[str]:
        """Get the set of unique nucleotides at a site across all sequences.

        Args:
            site_index: Zero-based nucleotide position

        Returns:
            Set of unique nucleotide characters
        """
        return {seq.sequence[site_index] for seq in self.sequences if site_index < len(seq)}

    def site_set_clean(self, site_index: int) -> set[str]:
        """Get unique nucleotides at a site, excluding N and gaps.

        Args:
            site_index: Zero-based nucleotide position

        Returns:
            Set of unique valid nucleotide characters
        """
        sites = self.site_set(site_index)
        return {s for s in sites if s not in ("N", "-")}

    def codon_array(self) -> np.ndarray:
        """Get a 2D array of codons (n_sequences x n_codons).

        Returns:
            NumPy array of codon strings
        """
        n_seqs = len(self.sequences)
        n_codons = self.num_codons
        arr = np.empty((n_seqs, n_codons), dtype="U3")
        for i, seq in enumerate(self.sequences):
            for j in range(n_codons):
                arr[i, j] = seq.get_codon(j, self.reading_frame)
        return arr

    def polymorphic_codons(self) -> list[int]:
        """Get indices of polymorphic codon positions.

        Returns:
            List of codon indices that are polymorphic
        """
        return [i for i in range(self.num_codons) if self.is_polymorphic(i)]

    def site_frequency_spectrum(self, codon_index: int) -> dict[str, float]:
        """Calculate allele frequencies at a codon position.

        Args:
            codon_index: Zero-based codon index

        Returns:
            Dict mapping codon strings to their frequencies (0-1)
        """
        counts: dict[str, int] = {}
        total = 0
        for seq in self.sequences:
            codon = seq.get_codon(codon_index, self.reading_frame)
            if "N" not in codon and "-" not in codon:
                counts[codon] = counts.get(codon, 0) + 1
                total += 1

        if total == 0:
            return {}

        return {codon: count / total for codon, count in counts.items()}

    def derived_allele_frequency(
        self, codon_index: int, ancestral_codon: str
    ) -> float | None:
        """Calculate the derived allele frequency at a codon position.

        Args:
            codon_index: Zero-based codon index
            ancestral_codon: The ancestral (outgroup) codon

        Returns:
            Frequency of derived alleles (0-1), or None if not computable
        """
        freqs = self.site_frequency_spectrum(codon_index)
        if not freqs:
            return None

        ancestral_codon = ancestral_codon.upper()
        if ancestral_codon not in freqs:
            # Ancestral allele not present - all are derived
            return 1.0

        return 1.0 - freqs[ancestral_codon]

    def filter_by_name(self, pattern: str) -> SequenceSet:
        """Filter sequences by name pattern (substring match).

        Args:
            pattern: Substring to match in sequence names

        Returns:
            New SequenceSet containing only matching sequences
        """
        filtered = [seq for seq in self.sequences if pattern in seq.name]
        return SequenceSet(
            sequences=filtered,
            reading_frame=self.reading_frame,
            genetic_code=self.genetic_code,
        )

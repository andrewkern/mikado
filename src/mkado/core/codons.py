"""Genetic code and codon utilities."""

from __future__ import annotations

from functools import lru_cache
from typing import TYPE_CHECKING

from mkado.data.genetic_codes import (
    STANDARD_CODE,
    _build_code_table,
    _compute_codon_paths,
    get_codon_paths,
)

if TYPE_CHECKING:
    pass


class GeneticCode:
    """Represents a genetic code for translation and codon path computation."""

    def __init__(
        self,
        code: dict[str, str] | None = None,
        *,
        table_id: int | None = None,
    ):
        """Initialize with a codon-to-amino-acid mapping.

        Args:
            code: Dict mapping codons to single-letter amino acids.
                  Must be a complete mapping of all 64 codons.
                  Uses standard code if neither code nor table_id is provided.
            table_id: NCBI genetic code table ID (1-33). If provided, overrides code.
        """
        if table_id is not None:
            self.code = _build_code_table(table_id)
            self._table_id = table_id
        elif code is not None:
            self.code = code
            self._table_id = None
        else:
            self.code = STANDARD_CODE
            self._table_id = 1

        # Use cached paths when constructed by table_id; compute for custom dicts
        if self._table_id is not None:
            self._paths = get_codon_paths(self._table_id)
        else:
            self._paths = _compute_codon_paths(self.code)

        self._syn_sites_cache: dict[str, float] = {}

    @lru_cache(maxsize=4096)
    def translate(self, codon: str) -> str:
        """Translate a codon to its amino acid.

        Args:
            codon: Three-letter codon string

        Returns:
            Single-letter amino acid code, or 'X' for unknown
        """
        codon = codon.upper()
        return self.code.get(codon, "X")

    def translate_sequence(self, sequence: str, reading_frame: int = 1) -> str:
        """Translate a nucleotide sequence to amino acids.

        Args:
            sequence: Nucleotide sequence string
            reading_frame: Reading frame (1, 2, or 3)

        Returns:
            Amino acid sequence string
        """
        start = reading_frame - 1
        amino_acids = []
        for i in range(start, len(sequence) - 2, 3):
            codon = sequence[i : i + 3]
            amino_acids.append(self.translate(codon))
        return "".join(amino_acids)

    def get_path(self, codon1: str, codon2: str) -> list[tuple[str, int]]:
        """Get the shortest mutational path between two codons.

        Args:
            codon1: Starting codon
            codon2: Ending codon

        Returns:
            List of (change_type, position) tuples where change_type is
            'R' for replacement (non-synonymous) or 'S' for synonymous,
            and position is the codon position (0, 1, or 2).
        """
        codon1 = codon1.upper()
        codon2 = codon2.upper()
        return self._paths.get((codon1, codon2), [])

    def count_synonymous_sites(self, codon: str) -> float:
        """Count the number of synonymous sites in a codon.

        Uses Nei-Gojobori method: for each site, calculate fraction of
        possible changes that are synonymous.

        Args:
            codon: Three-letter codon string

        Returns:
            Number of synonymous sites (0-3)
        """
        codon = codon.upper()
        cached = self._syn_sites_cache.get(codon)
        if cached is not None:
            return cached

        if "N" in codon or "-" in codon:
            self._syn_sites_cache[codon] = 0.0
            return 0.0

        aa = self.translate(codon)
        if aa == "*" or aa == "X":
            self._syn_sites_cache[codon] = 0.0
            return 0.0

        syn_sites = 0.0
        nucleotides = ["A", "C", "G", "T"]

        for pos in range(3):
            syn_count = 0
            total_count = 0
            for nt in nucleotides:
                if nt == codon[pos]:
                    continue
                new_codon = codon[:pos] + nt + codon[pos + 1 :]
                new_aa = self.translate(new_codon)
                if new_aa != "*":
                    total_count += 1
                    if new_aa == aa:
                        syn_count += 1
            if total_count > 0:
                syn_sites += syn_count / total_count

        self._syn_sites_cache[codon] = syn_sites
        return syn_sites

    def is_synonymous_change(self, codon1: str, codon2: str) -> bool | None:
        """Check if a single-nucleotide change is synonymous.

        Args:
            codon1: Original codon
            codon2: Changed codon

        Returns:
            True if synonymous, False if non-synonymous, None if not
            a single-nucleotide change or invalid codons
        """
        codon1 = codon1.upper()
        codon2 = codon2.upper()

        # Check single nucleotide difference
        diffs = sum(1 for i in range(3) if codon1[i] != codon2[i])
        if diffs != 1:
            return None

        aa1 = self.translate(codon1)
        aa2 = self.translate(codon2)

        if aa1 == "X" or aa2 == "X":
            return None

        return aa1 == aa2


# Default genetic code instance
DEFAULT_CODE = GeneticCode()

"""Tests for genetic code and codon utilities."""

import pytest

from mkado.core.codons import DEFAULT_CODE, GeneticCode


class TestGeneticCode:
    """Tests for GeneticCode class."""

    def test_translate_basic(self) -> None:
        """Test basic codon translation."""
        code = GeneticCode()

        assert code.translate("ATG") == "M"  # Start codon
        assert code.translate("TAA") == "*"  # Stop codon
        assert code.translate("TGG") == "W"  # Tryptophan
        assert code.translate("GCT") == "A"  # Alanine

    def test_translate_lowercase(self) -> None:
        """Test that lowercase codons are handled."""
        code = GeneticCode()

        assert code.translate("atg") == "M"

    def test_translate_unknown(self) -> None:
        """Test translation of unknown codons."""
        code = GeneticCode()

        assert code.translate("NNN") == "X"
        assert code.translate("---") == "X"

    def test_translate_sequence(self) -> None:
        """Test translating a full sequence."""
        code = GeneticCode()

        # ATG = M, GCT = A, TAA = *
        seq = "ATGGCTTAA"
        assert code.translate_sequence(seq) == "MA*"

    def test_translate_sequence_reading_frame(self) -> None:
        """Test translation in different reading frames."""
        code = GeneticCode()

        # Frame 1: ATG GCT -> MA
        # Frame 2: TGG CT -> W (incomplete)
        # Frame 3: GGC T -> G (incomplete)
        seq = "ATGGCT"
        assert code.translate_sequence(seq, reading_frame=1) == "MA"
        assert code.translate_sequence(seq, reading_frame=2) == "W"
        assert code.translate_sequence(seq, reading_frame=3) == "G"

    def test_get_path_single_change(self) -> None:
        """Test path for single nucleotide change."""
        code = GeneticCode()

        # AAA (Lys) -> AAG (Lys) - synonymous at position 2
        path = code.get_path("AAA", "AAG")
        assert len(path) == 1
        assert path[0] == ("S", 2)

        # AAA (Lys) -> GAA (Glu) - replacement at position 0
        path = code.get_path("AAA", "GAA")
        assert len(path) == 1
        assert path[0] == ("R", 0)

    def test_get_path_two_changes(self) -> None:
        """Test path for two nucleotide changes."""
        code = GeneticCode()

        # AAA (Lys) -> AGA (Arg) - two changes at positions 1 and 2
        # A->G at position 1, A->A at position 2 (wait, that's one change)
        # Let's use AAA -> ACA (Thr) - one change at position 1
        # Or AAT (Asn) -> ACT (Thr) - one change at position 1
        # Actually, let's use a true 2-change case: AAA -> GAC
        # AAA (Lys) -> GAC (Asp) - changes at positions 0 and 2
        path = code.get_path("AAA", "GAC")
        assert len(path) == 2
        # Path should have some R and/or S changes

    def test_get_path_same_codon(self) -> None:
        """Test path for identical codons."""
        code = GeneticCode()

        path = code.get_path("ATG", "ATG")
        assert path == []

    def test_count_synonymous_sites(self) -> None:
        """Test counting synonymous sites."""
        code = GeneticCode()

        # ATG (Met) has no synonymous sites (only one codon for Met)
        assert code.count_synonymous_sites("ATG") == 0.0

        # Leucine codons have variable degeneracy
        # TTT (Phe) - position 2 can be C and still be Phe
        syn_sites = code.count_synonymous_sites("TTT")
        assert syn_sites > 0

    def test_count_synonymous_sites_ambiguous(self) -> None:
        """Test that ambiguous codons return 0 synonymous sites."""
        code = GeneticCode()

        assert code.count_synonymous_sites("NNN") == 0.0
        assert code.count_synonymous_sites("AT-") == 0.0

    def test_is_synonymous_change(self) -> None:
        """Test identifying synonymous changes."""
        code = GeneticCode()

        # TTT -> TTC (both Phe) - synonymous
        assert code.is_synonymous_change("TTT", "TTC") is True

        # TTT -> CTT (Phe -> Leu) - non-synonymous
        assert code.is_synonymous_change("TTT", "CTT") is False

    def test_is_synonymous_change_multiple_diffs(self) -> None:
        """Test that multiple changes return None."""
        code = GeneticCode()

        # Two differences
        assert code.is_synonymous_change("AAA", "GGG") is None

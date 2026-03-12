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


class TestGeneticCodeTables:
    """Tests for alternate genetic code table support."""

    def test_table_id_standard(self) -> None:
        """Standard code via table_id gives same results as default."""
        default = GeneticCode()
        standard = GeneticCode(table_id=1)
        assert default.translate("ATG") == standard.translate("ATG")
        assert default.translate("TGA") == standard.translate("TGA")

    def test_table_id_vertebrate_mito(self) -> None:
        """Vertebrate mitochondrial code has expected differences."""
        code = GeneticCode(table_id=2)
        assert code.translate("AGA") == "*"  # Stop, not Arg
        assert code.translate("AGG") == "*"  # Stop, not Arg
        assert code.translate("ATA") == "M"  # Met, not Ile
        assert code.translate("TGA") == "W"  # Trp, not Stop

    def test_table_id_invertebrate_mito(self) -> None:
        """Invertebrate mitochondrial code has expected differences."""
        code = GeneticCode(table_id=5)
        assert code.translate("AGA") == "S"  # Ser, not Arg
        assert code.translate("TGA") == "W"  # Trp, not Stop

    def test_paths_differ_between_codes(self) -> None:
        """Paths should differ when amino acid assignments differ."""
        std = GeneticCode()
        mito = GeneticCode(table_id=2)
        # ATG->AGA: in standard, AGA=R so this is R; in mito, AGA=* (stop)
        std_path = std.get_path("ATG", "AGA")
        mito_path = mito.get_path("ATG", "AGA")
        assert std_path != mito_path

    def test_codon_paths_cached(self) -> None:
        """Repeated construction with same table_id shares cached paths."""
        code1 = GeneticCode(table_id=5)
        code2 = GeneticCode(table_id=5)
        assert code1._paths is code2._paths

    def test_invalid_table_id(self) -> None:
        """Invalid table ID raises ValueError."""
        with pytest.raises(ValueError, match="Unknown genetic code"):
            GeneticCode(table_id=99)


class TestResolveCodeTable:
    """Tests for resolve_code_table."""

    def test_numeric_id(self) -> None:
        from mkado.data.genetic_codes import resolve_code_table

        assert resolve_code_table("1") == 1
        assert resolve_code_table("2") == 2
        assert resolve_code_table("5") == 5

    def test_name_alias(self) -> None:
        from mkado.data.genetic_codes import resolve_code_table

        assert resolve_code_table("standard") == 1
        assert resolve_code_table("vertebrate-mito") == 2
        assert resolve_code_table("invertebrate-mito") == 5

    def test_case_insensitive(self) -> None:
        from mkado.data.genetic_codes import resolve_code_table

        assert resolve_code_table("Vertebrate-Mito") == 2
        assert resolve_code_table("STANDARD") == 1

    def test_unknown_name_raises(self) -> None:
        from mkado.data.genetic_codes import resolve_code_table

        with pytest.raises(ValueError, match="Unknown genetic code"):
            resolve_code_table("not-a-code")

    def test_unknown_id_raises(self) -> None:
        from mkado.data.genetic_codes import resolve_code_table

        with pytest.raises(ValueError, match="Unknown genetic code"):
            resolve_code_table("99")

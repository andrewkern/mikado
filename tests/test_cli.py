"""Tests for CLI option validation."""

from pathlib import Path

import pytest
from typer.testing import CliRunner

from mkado.cli import app

runner = CliRunner()


class TestOptionValidation:
    """Tests for CLI option compatibility validation."""

    def test_asymptotic_with_min_freq_error(self, tmp_path: Path) -> None:
        """Test that --asymptotic and --min-freq cannot be used together."""
        # Create a minimal test file
        fasta = tmp_path / "test.fa"
        fasta.write_text(""">speciesA_1
ATGATGATG
>speciesB_1
ATGGTGATG
""")

        result = runner.invoke(
            app,
            [
                "test",
                str(fasta),
                "-i", "speciesA",
                "-o", "speciesB",
                "--asymptotic",
                "--min-freq", "0.1",
            ],
        )

        assert result.exit_code == 1
        assert "--min-freq cannot be used with --asymptotic" in result.output
        assert "--freq-cutoffs" in result.output

    def test_alpha_tg_with_min_freq_error(self, tmp_path: Path) -> None:
        """Test that --alpha-tg and --min-freq cannot be used together."""
        # Create a minimal test directory with alignment
        alignment_dir = tmp_path / "alignments"
        alignment_dir.mkdir()
        fasta = alignment_dir / "test.fa"
        fasta.write_text(""">speciesA_1
ATGATGATG
>speciesB_1
ATGGTGATG
""")

        result = runner.invoke(
            app,
            [
                "batch",
                str(alignment_dir),
                "-i", "speciesA",
                "-o", "speciesB",
                "--alpha-tg",
                "--min-freq", "0.1",
            ],
        )

        assert result.exit_code == 1
        assert "--min-freq cannot be used with --alpha-tg" in result.output

    def test_alpha_tg_with_asymptotic_error(self, tmp_path: Path) -> None:
        """Test that --alpha-tg and --asymptotic cannot be used together."""
        alignment_dir = tmp_path / "alignments"
        alignment_dir.mkdir()
        fasta = alignment_dir / "test.fa"
        fasta.write_text(""">speciesA_1
ATGATGATG
>speciesB_1
ATGGTGATG
""")

        result = runner.invoke(
            app,
            [
                "batch",
                str(alignment_dir),
                "-i", "speciesA",
                "-o", "speciesB",
                "--alpha-tg",
                "--asymptotic",
            ],
        )

        assert result.exit_code == 1
        assert "--alpha-tg and --asymptotic are mutually exclusive" in result.output

    def test_asymptotic_with_polarize_match_error(self, tmp_path: Path) -> None:
        """Test that --asymptotic and --polarize-match cannot be used together."""
        fasta = tmp_path / "test.fa"
        fasta.write_text(""">speciesA_1
ATGATGATG
>speciesB_1
ATGGTGATG
>speciesC_1
ATGATGATG
""")

        result = runner.invoke(
            app,
            [
                "test",
                str(fasta),
                "-i", "speciesA",
                "-o", "speciesB",
                "--asymptotic",
                "--polarize-match", "speciesC",
            ],
        )

        assert result.exit_code == 1
        assert "Polarized asymptotic test not supported" in result.output

    def test_asymptotic_without_min_freq_succeeds(self, tmp_path: Path) -> None:
        """Test that --asymptotic works without --min-freq."""
        fasta = tmp_path / "test.fa"
        fasta.write_text(""">speciesA_1
ATGATGATGATGATGATG
>speciesA_2
ATGCTGATGATGATGATG
>speciesB_1
ATGGTGATGATGATGATG
""")

        result = runner.invoke(
            app,
            [
                "test",
                str(fasta),
                "-i", "speciesA",
                "-o", "speciesB",
                "--asymptotic",
            ],
        )

        # Should succeed (exit code 0)
        assert result.exit_code == 0

    def test_min_freq_without_asymptotic_succeeds(self, tmp_path: Path) -> None:
        """Test that --min-freq works without --asymptotic (standard MK test)."""
        fasta = tmp_path / "test.fa"
        fasta.write_text(""">speciesA_1
ATGATGATG
>speciesA_2
ATGCTGATG
>speciesB_1
ATGGTGATG
""")

        result = runner.invoke(
            app,
            [
                "test",
                str(fasta),
                "-i", "speciesA",
                "-o", "speciesB",
                "--min-freq", "0.1",
            ],
        )

        # Should succeed
        assert result.exit_code == 0


class TestBatchOptionValidation:
    """Tests for batch command option validation."""

    def test_batch_asymptotic_with_polarize_match_error(self, tmp_path: Path) -> None:
        """Test that batch --asymptotic and --polarize-match cannot be used together."""
        alignment_dir = tmp_path / "alignments"
        alignment_dir.mkdir()
        fasta = alignment_dir / "test.fa"
        fasta.write_text(""">speciesA_1
ATGATGATG
>speciesB_1
ATGGTGATG
>speciesC_1
ATGATGATG
""")

        result = runner.invoke(
            app,
            [
                "batch",
                str(alignment_dir),
                "-i", "speciesA",
                "-o", "speciesB",
                "--asymptotic",
                "--polarize-match", "speciesC",
            ],
        )

        assert result.exit_code == 1
        assert "--asymptotic and --polarize-match are mutually exclusive" in result.output


class TestNoSingletonsOption:
    """Tests for --no-singletons option."""

    def test_no_singletons_with_min_freq_error(self, tmp_path: Path) -> None:
        """Test that --no-singletons and --min-freq cannot be used together."""
        fasta = tmp_path / "test.fa"
        fasta.write_text(""">speciesA_1
ATGATGATG
>speciesB_1
ATGGTGATG
""")

        result = runner.invoke(
            app,
            [
                "test",
                str(fasta),
                "-i", "speciesA",
                "-o", "speciesB",
                "--no-singletons",
                "--min-freq", "0.1",
            ],
        )

        assert result.exit_code == 1
        assert "--no-singletons and --min-freq cannot be used together" in result.output

    def test_no_singletons_with_asymptotic_error(self, tmp_path: Path) -> None:
        """Test that --no-singletons and --asymptotic cannot be used together."""
        fasta = tmp_path / "test.fa"
        fasta.write_text(""">speciesA_1
ATGATGATG
>speciesB_1
ATGGTGATG
""")

        result = runner.invoke(
            app,
            [
                "test",
                str(fasta),
                "-i", "speciesA",
                "-o", "speciesB",
                "--no-singletons",
                "--asymptotic",
            ],
        )

        assert result.exit_code == 1
        assert "--no-singletons cannot be used with --asymptotic" in result.output

    def test_no_singletons_succeeds(self, tmp_path: Path) -> None:
        """Test that --no-singletons works correctly."""
        fasta = tmp_path / "test.fa"
        fasta.write_text(""">speciesA_1
ATGATGATG
>speciesA_2
ATGCTGATG
>speciesB_1
ATGGTGATG
""")

        result = runner.invoke(
            app,
            [
                "test",
                str(fasta),
                "-i", "speciesA",
                "-o", "speciesB",
                "--no-singletons",
            ],
        )

        assert result.exit_code == 0
        # Should show the singleton exclusion message
        assert "Excluding singletons" in result.output
        assert "min frequency" in result.output

    def test_batch_no_singletons_with_alpha_tg_error(self, tmp_path: Path) -> None:
        """Test that --no-singletons and --alpha-tg cannot be used together."""
        alignment_dir = tmp_path / "alignments"
        alignment_dir.mkdir()
        fasta = alignment_dir / "test.fa"
        fasta.write_text(""">speciesA_1
ATGATGATG
>speciesB_1
ATGGTGATG
""")

        result = runner.invoke(
            app,
            [
                "batch",
                str(alignment_dir),
                "-i", "speciesA",
                "-o", "speciesB",
                "--no-singletons",
                "--alpha-tg",
            ],
        )

        assert result.exit_code == 1
        assert "--no-singletons cannot be used with --alpha-tg" in result.output

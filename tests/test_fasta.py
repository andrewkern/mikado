"""Tests for FASTA parser."""

from pathlib import Path

import pytest

from mkado.io.fasta import read_fasta, write_fasta


def test_read_fasta(tmp_path: Path) -> None:
    """Test reading a FASTA file."""
    fasta_content = """>seq1
ATGCATGC
ATGCATGC
>seq2
GGGGCCCC
"""
    fasta_file = tmp_path / "test.fa"
    fasta_file.write_text(fasta_content)

    sequences = list(read_fasta(fasta_file))

    assert len(sequences) == 2
    assert sequences[0] == ("seq1", "ATGCATGCATGCATGC")
    assert sequences[1] == ("seq2", "GGGGCCCC")


def test_read_fasta_with_header_info(tmp_path: Path) -> None:
    """Test that only the first word of the header is used as name."""
    fasta_content = """>seq1 description here
ATGCATGC
"""
    fasta_file = tmp_path / "test.fa"
    fasta_file.write_text(fasta_content)

    sequences = list(read_fasta(fasta_file))

    assert sequences[0][0] == "seq1"


def test_read_fasta_uppercase(tmp_path: Path) -> None:
    """Test that sequences are converted to uppercase."""
    fasta_content = """>seq1
atgcATGC
"""
    fasta_file = tmp_path / "test.fa"
    fasta_file.write_text(fasta_content)

    sequences = list(read_fasta(fasta_file))

    assert sequences[0][1] == "ATGCATGC"


def test_write_fasta(tmp_path: Path) -> None:
    """Test writing a FASTA file."""
    sequences = [("seq1", "ATGCATGC"), ("seq2", "GGGGCCCC")]
    fasta_file = tmp_path / "output.fa"

    write_fasta(sequences, fasta_file, line_width=4)

    content = fasta_file.read_text()
    lines = content.strip().split("\n")

    assert lines[0] == ">seq1"
    assert lines[1] == "ATGC"
    assert lines[2] == "ATGC"
    assert lines[3] == ">seq2"


def test_read_empty_fasta(tmp_path: Path) -> None:
    """Test reading an empty FASTA file."""
    fasta_file = tmp_path / "empty.fa"
    fasta_file.write_text("")

    sequences = list(read_fasta(fasta_file))

    assert len(sequences) == 0

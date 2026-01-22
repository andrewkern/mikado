"""FASTA file parser and writer."""

from __future__ import annotations

from pathlib import Path
from typing import Iterator


def read_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """Parse a FASTA file and yield (name, sequence) tuples.

    Args:
        path: Path to FASTA file

    Yields:
        Tuples of (sequence_name, sequence_string)
    """
    path = Path(path)
    name: str | None = None
    seq_parts: list[str] = []

    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq_parts)
                name = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.upper())

        if name is not None:
            yield name, "".join(seq_parts)


def write_fasta(
    sequences: list[tuple[str, str]],
    path: str | Path,
    line_width: int = 70,
) -> None:
    """Write sequences to a FASTA file.

    Args:
        sequences: List of (name, sequence) tuples
        path: Output file path
        line_width: Characters per line for sequence wrapping
    """
    path = Path(path)
    with path.open("w") as f:
        for name, seq in sequences:
            f.write(f">{name}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i : i + line_width] + "\n")

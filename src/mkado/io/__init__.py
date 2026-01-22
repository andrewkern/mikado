"""Input/output utilities."""

from mkado.io.fasta import read_fasta, write_fasta
from mkado.io.output import format_result, OutputFormat

__all__ = ["read_fasta", "write_fasta", "format_result", "OutputFormat"]

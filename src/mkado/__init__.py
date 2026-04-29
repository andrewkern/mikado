"""MKado: Modern Python implementation of the McDonald-Kreitman test toolkit."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("mkado")
except PackageNotFoundError:
    __version__ = "0.0.0+unknown"

from mkado.core.sequences import Sequence, SequenceSet
from mkado.core.codons import GeneticCode
from mkado.analysis.mk_test import MKResult, mk_test
from mkado.analysis.polarized import PolarizedMKResult, polarized_mk_test
from mkado.analysis.asymptotic import AsymptoticMKResult, asymptotic_mk_test

__all__ = [
    "Sequence",
    "SequenceSet",
    "GeneticCode",
    "MKResult",
    "mk_test",
    "PolarizedMKResult",
    "polarized_mk_test",
    "AsymptoticMKResult",
    "asymptotic_mk_test",
]

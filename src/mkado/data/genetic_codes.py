"""Genetic code tables and codon data."""

from __future__ import annotations

# Standard genetic code (NCBI table 1)
STANDARD_CODE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# All 64 codons in alphabetical order
CODONS = sorted(STANDARD_CODE.keys())

# Build codon lookup table
CODON_TABLE: dict[str, int] = {codon: i for i, codon in enumerate(CODONS)}

# Nucleotides
NUCLEOTIDES = ["A", "C", "G", "T"]


def _compute_codon_paths() -> dict[tuple[str, str], list[tuple[str, int]]]:
    """Pre-compute shortest paths between all codon pairs.

    Returns a dict mapping (codon1, codon2) to a list of (type, position) tuples,
    where type is 'R' for replacement or 'S' for synonymous, and position is
    the codon position (0, 1, or 2) where the change occurs.
    """
    paths: dict[tuple[str, str], list[tuple[str, int]]] = {}

    for c1 in CODONS:
        for c2 in CODONS:
            if c1 == c2:
                paths[(c1, c2)] = []
                continue

            # Find positions that differ
            diffs = [i for i in range(3) if c1[i] != c2[i]]

            if len(diffs) == 1:
                # Single nucleotide change
                pos = diffs[0]
                aa1 = STANDARD_CODE[c1]
                aa2 = STANDARD_CODE[c2]
                change_type = "S" if aa1 == aa2 else "R"
                paths[(c1, c2)] = [(change_type, pos)]

            elif len(diffs) == 2:
                # Two changes - find shortest path (minimize replacements)
                best_path: list[tuple[str, int]] = []
                best_replacements = 999

                for first_pos in diffs:
                    second_pos = [p for p in diffs if p != first_pos][0]
                    # Make intermediate codon
                    inter = list(c1)
                    inter[first_pos] = c2[first_pos]
                    inter_codon = "".join(inter)

                    if STANDARD_CODE.get(inter_codon) == "*":
                        continue  # Skip paths through stop codons

                    # Calculate path through this intermediate
                    aa1 = STANDARD_CODE[c1]
                    aa_inter = STANDARD_CODE[inter_codon]
                    aa2 = STANDARD_CODE[c2]

                    step1_type = "S" if aa1 == aa_inter else "R"
                    step2_type = "S" if aa_inter == aa2 else "R"

                    path = [(step1_type, first_pos), (step2_type, second_pos)]
                    n_replacements = sum(1 for t, _ in path if t == "R")

                    if n_replacements < best_replacements:
                        best_replacements = n_replacements
                        best_path = path

                paths[(c1, c2)] = best_path

            else:
                # Three changes - find shortest path
                best_path = []
                best_replacements = 999

                # Try all 6 orderings
                from itertools import permutations

                for order in permutations(diffs):
                    current = list(c1)
                    path = []
                    valid = True

                    for i, pos in enumerate(order):
                        old_codon = "".join(current)
                        current[pos] = c2[pos]
                        new_codon = "".join(current)

                        aa_old = STANDARD_CODE.get(old_codon)
                        aa_new = STANDARD_CODE.get(new_codon)

                        if aa_new == "*" and i < 2:  # Don't go through stop codons
                            valid = False
                            break

                        change_type = "S" if aa_old == aa_new else "R"
                        path.append((change_type, pos))

                    if valid:
                        n_replacements = sum(1 for t, _ in path if t == "R")
                        if n_replacements < best_replacements:
                            best_replacements = n_replacements
                            best_path = path

                paths[(c1, c2)] = best_path

    return paths


# Pre-computed codon paths
CODON_PATHS = _compute_codon_paths()

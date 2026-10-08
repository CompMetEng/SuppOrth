"""Expand known protein identifiers in each table row into unique cross-species combinations."""

from __future__ import annotations

import re
from collections.abc import Sequence
from itertools import combinations
from pathlib import Path
from typing import AbstractSet, Iterator

# Do not split on "|": NCBI identifiers such as gi|3844770|gb|AAC71393.1|
# would otherwise never match the FASTA ids.
_TOKEN_SPLIT = re.compile(r"[\t,;\s]+")


def tokenize(line: str) -> list[str]:
    return [token.strip() for token in _TOKEN_SPLIT.split(line.strip()) if token.strip()]


def pairs_from_table(
    path: Path,
    *id_sets: AbstractSet[str],
    max_row_pairs: int | None = None,
) -> Iterator[tuple[str, str]]:
    """Yield cross-species combinations in input order. An explicit row limit raises an error rather than dropping predictions."""
    if len(id_sets) < 2:
        raise ValueError("pairs_from_table needs at least two identifier sets")

    owner: dict[str, int] = {}
    for index, ids in enumerate(id_sets):
        for gene in ids:
            owner[gene] = index

    with Path(path).open(errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            tokens = dict.fromkeys(tokenize(line))
            if not tokens:
                continue

            by_species: dict[int, list[str]] = {}
            for token in tokens:
                index = owner.get(token)
                if index is None:
                    continue
                by_species.setdefault(index, []).append(token)
            if len(by_species) < 2:
                continue

            n_pairs = 0
            species = sorted(by_species)
            for left_i, right_i in combinations(species, 2):
                n_pairs += len(by_species[left_i]) * len(by_species[right_i])
            if max_row_pairs is not None and n_pairs > max_row_pairs:
                raise ValueError(
                    f"{path}:{line_number} contains {n_pairs} pairs, "
                    f"exceeding the requested limit of {max_row_pairs}"
                )

            for left_i, right_i in combinations(species, 2):
                for left in by_species[left_i]:
                    for right in by_species[right_i]:
                        yield left, right


def id_sets(proteomes: Sequence) -> tuple[frozenset[str], ...]:
    return tuple(proteome.ids for proteome in proteomes)

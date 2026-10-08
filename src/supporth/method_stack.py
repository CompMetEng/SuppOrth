"""Default method descriptions for the supported predictors. OMA is imported from an external run."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Stack:
    """One predictor's place in the complementary split."""

    search: str
    clustering: str
    alignment: str
    tree: str
    orthology: str


ORTHOFINDER = Stack(
    search="diamond",
    clustering="MCL",
    alignment="mafft",
    tree="iqtree",
    orthology="gene-tree reconciliation",
)

BROCCOLI = Stack(
    search="diamond",
    clustering="Broccoli graph/network",
    alignment="Broccoli internal",
    tree="FastTree",
    orthology="hybrid graph and phylogeny",
)

SONICPARANOID = Stack(
    search="mmseqs",
    clustering="MCL",
    alignment="none",
    tree="none",
    orthology="InParanoid-like graph with Pfam domain analysis",
)

OMA = Stack(
    search="OMA native Smith-Waterman",
    clustering="OMA groups / HOGs",
    alignment="OMA internal",
    tree="OMA species-tree framework",
    orthology="evolutionary-distance consistency",
)

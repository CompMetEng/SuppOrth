"""Import existing OMA pairwise tables. Protein identifiers are matched against the input FASTAs."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

from ..canonical import PredictorResult, build_result
from ..method_stack import OMA as OMA_STACK
from ..proteomes import Proteome
from .base import Adapter, AdapterError, MethodProfile, StageContext
from .tabular import id_sets, pairs_from_table


class OMA(Adapter):
    name = "oma"
    label = "O"
    repo = ""
    default_ref = ""
    entrypoint_candidates = ()
    binaries = ()
    install_notes = (
        "not installed by suppOrth; pass --oma to an existing OMA standalone run"
    )

    def install(self, **kwargs):  # noqa: ANN003
        raise AdapterError(
            "OMA is not installed by suppOrth. Run OMA standalone yourself and "
            "pass the result directory with: suppOrth run ... --oma PATH"
        )

    def check(self, options: list[str] | None = None) -> tuple[bool, str]:
        return True, "user-supplied output; nothing to install"

    def method_profile(self, context: StageContext) -> MethodProfile:
        return MethodProfile(
            search=OMA_STACK.search,
            clustering=OMA_STACK.clustering,
            alignment=OMA_STACK.alignment,
            tree=OMA_STACK.tree,
            orthology=OMA_STACK.orthology,
        )

    def run(self, context: StageContext) -> Path:
        raise AdapterError(
            "OMA is not run by suppOrth; pass --oma PATH to an existing Output/ directory"
        )

    def parse(self, native_root: Path, proteomes: Sequence[Proteome] | Proteome, *rest: Proteome) -> PredictorResult:
        loaded = (proteomes, *rest) if isinstance(proteomes, Proteome) else tuple(proteomes)
        tables = self.result_tables(native_root, proteomes=loaded)
        if not tables:
            raise AdapterError(
                f"no OMA pairwise tables under {native_root}. "
                f"Expected Output/PairwiseOrthologs/*.txt from an OMA standalone run."
            )
        raw: list[tuple[str, str]] = []
        for table in tables:
            raw.extend(pairs_from_table(table, *id_sets(loaded)))
        return build_result(self.name, self.label, raw, loaded, source_files=tables)

    def result_tables(
        self,
        native_root: Path,
        proteomes: Sequence[Proteome] | None = None,
        sp1: Proteome | None = None,
        sp2: Proteome | None = None,
    ) -> list[Path]:
        if proteomes is None and sp1 is not None and sp2 is not None:
            proteomes = (sp1, sp2)
        root = Path(native_root).expanduser()
        if root.is_file():
            return [root]
        pairwise = self._pairwise_tables(root)
        if pairwise:
            return pairwise
        return sorted(p for p in root.rglob("OrthologousGroups.txt") if p.is_file())

    def _pairwise_tables(self, root: Path) -> list[Path]:
        directories = [p for p in root.rglob("PairwiseOrthologs") if p.is_dir()]
        if root.name == "PairwiseOrthologs" and root.is_dir():
            directories.append(root)
        found: list[Path] = []
        seen: set[Path] = set()
        for directory in directories:
            for path in directory.iterdir():
                if path in seen:
                    continue
                if not path.is_file() or path.suffix != ".txt":
                    continue
                if "cutted" in path.name.lower():
                    continue
                seen.add(path)
                found.append(path)
        return sorted(found)

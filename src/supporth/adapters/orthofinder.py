"""Run OrthoFinder and read its pairwise orthologue tables. Orthogroup-only output is marked as a fallback."""

from __future__ import annotations

import csv
import shutil
from dataclasses import dataclass
from collections.abc import Sequence
from itertools import combinations
from pathlib import Path

from ..canonical import PredictorResult, build_result
from ..proteomes import Proteome
from .. import shell
from ..method_stack import ORTHOFINDER
from . import of_config
from .base import Adapter, AdapterError, MethodProfile, StageContext


def find_alignment_dir(native_root: Path) -> Path | None:
    """Directory that holds SpeciesIDs.txt and the Blast*.txt tables."""
    root = Path(native_root)
    if (root / "WorkingDirectory" / "SpeciesIDs.txt").is_file():
        work_dir = root / "WorkingDirectory"
    elif (root / "SpeciesIDs.txt").is_file():
        work_dir = root
    else:
        found = next(root.rglob("SpeciesIDs.txt"), None)
        if found is None:
            return None
        work_dir = found.parent
    if not any(work_dir.glob("Blast*_*.txt*")):
        return None
    return work_dir


def _split_genes(cell: str) -> list[str]:
    if not cell:
        return []
    return [token.strip() for token in cell.replace(";", ",").split(",") if token.strip()]


@dataclass
class Methods:
    """OrthoFinder search, alignment and tree methods for one run."""

    search: str = ORTHOFINDER.search
    mode: str = "msa"
    aligner: str = ORTHOFINDER.alignment
    tree: str = ORTHOFINDER.tree

    @property
    def tree_program(self) -> str:
        if self.mode != "msa":
            return "dendroblast (distance trees from the similarity graph)"
        return self.tree or "iqtree"

    def required_binaries(self) -> tuple[str, ...]:
        needed = list(of_config.binaries_for_search(self.search or ORTHOFINDER.search))
        if self.mode != "msa":
            return tuple(needed)
        needed.append(self.aligner or ORTHOFINDER.alignment)
        tree = self.tree or ORTHOFINDER.tree
        # OrthoFinder spells FastTree lowercase in its config; the shared bin
        # carries the upstream capitalisation.
        needed.append("FastTree" if tree.lower() == "fasttree" else tree)
        return tuple(needed)


def parse_methods(options: list[str] | None) -> Methods:
    """Read OrthoFinder's -S/-M/-A/-T flags out of pass-through options."""
    methods = Methods()
    tokens = list(options or [])
    flags = {"-S": "search", "-M": "mode", "-A": "aligner", "-T": "tree"}
    for index, token in enumerate(tokens[:-1]):
        attribute = flags.get(token)
        if attribute:
            setattr(methods, attribute, tokens[index + 1])
    return methods


class OrthoFinder(Adapter):
    name = "orthofinder"
    label = "F"
    repo = "https://github.com/davidemms/OrthoFinder.git"
    default_ref = "v2.5.5"
    entrypoint_candidates = (
        "orthofinder.py",
        "orthofinder/orthofinder.py",
        "orthofinder/__main__.py",
        "scripts_of/orthofinder.py",
        "*/orthofinder.py",
    )
    binaries = ("mcl",)
    # scripts_of/__main__.py prepends its own bin/ to PATH, so its bundled
    # diamond, mcl and fastme win over the shared ones unless they are linked.
    shadows_path = True
    install_notes = (
        "GPL-3.0; fetched from upstream and invoked as a subprocess. "
        "Default stack is the --alignment search (DIAMOND or BLAST), MCL, "
        "MAFFT and IQ-TREE (gene-tree reconciliation). Bundled helper binaries "
        "are linked to the shared bin at install time."
    )

    search_program = ORTHOFINDER.search
    tree_method = "msa"

    def post_install(self, root: Path, *, interpreter: Path | None = None, log=print) -> None:
        """Install numpy and scipy in the predictor virtualenv."""
        self.build_venv(root, interpreter, ("numpy", "scipy"), log=log)

    def _interpreter(self, root: Path) -> str:
        return str(self.venv_python(root) or self.python())

    def check(self, options: list[str] | None = None) -> tuple[bool, str]:
        ready, detail = super().check(options)
        if not ready:
            return ready, detail
        root = self.installed_root()
        probe = shell.run(
            [self._interpreter(root), "-c", "import numpy, scipy"],
            check=False,
        )
        if probe.returncode != 0:
            return False, "numpy/scipy unavailable; re-run suppOrth install --tools orthofinder"
        return True, detail

    def extra_binaries(self, options: list[str] | None = None) -> tuple[str, ...]:
        return parse_methods(options).required_binaries()

    def method_profile(self, context: StageContext) -> MethodProfile:
        methods = parse_methods(context.options)
        program = methods.search or self.search_program
        jobs = of_config.search_job_count(len(context.proteomes), double_blast="-1" not in context.options)
        return MethodProfile(
            search=f"{program} (via OrthoFinder)",
            clustering="MCL on the similarity graph",
            alignment=(methods.aligner or ORTHOFINDER.alignment) if methods.mode == "msa" else "none",
            tree=methods.tree_program,
            orthology=ORTHOFINDER.orthology if methods.mode == "msa" else "dendroblast graph clustering",
            extra={
                "mode": methods.mode,
                "search_jobs": str(jobs),
                "threads_per_search": str(
                    of_config.threads_per_job(context.threads, jobs)
                ),
            },
        )

    def run(self, context: StageContext) -> Path:
        root = self.installed_root()
        entrypoint = self.locate_entrypoint(root)

        previous = None
        output = context.work_dir / "results"
        if context.previous_root is not None:
            previous = find_alignment_dir(context.previous_root)
            if previous is None:
                raise AdapterError(f"no reusable OrthoFinder searches under {context.previous_root}")
            if not (previous.parent / "Log.txt").is_file():
                raise AdapterError(f"OrthoFinder update requires {previous.parent / 'Log.txt'}")
            if context.added_input_dir is None:
                raise AdapterError("OrthoFinder update requires the additional proteomes")
            number = 1
            output = previous.parent.parent / f"Results_update_{number}"
            while output.exists():
                number += 1
                output = previous.parent.parent / f"Results_update_{number}"
        elif output.exists():
            shutil.rmtree(output)
        context.work_dir.mkdir(parents=True, exist_ok=True)

        methods = parse_methods(context.options)
        requested = methods.search or self.search_program

        jobs = of_config.search_job_count(len(context.proteomes), double_blast="-1" not in context.options)
        per_job = of_config.threads_per_job(context.threads, jobs)

        with of_config.config_lock(root):
            search_program = self._configure_search(root, requested, per_job)

            argv = [
                self._interpreter(root),
                str(entrypoint),
                "-f",
                str(context.added_input_dir if previous else context.input_dir),
                "-t",
                str(context.threads),
                "-a",
                str(max(1, context.threads // 2)),
                "-S",
                search_program,
            ]
            if previous:
                # Version 2.5.5 rejects -o together with -b. Updated results
                # must be siblings of the previous Results_* directory.
                argv += ["-b", str(previous.parent), "-n", f"update_{number}"]
            else:
                argv += ["-o", str(output)]
            passthrough = of_config.strip_option(context.options, "-S")
            if previous:
                passthrough = of_config.strip_option(passthrough, "-n")
            # Defaults are only supplied when the run did not override them, so
            # a pass-through "-M dendroblast" does not collide with MSA mode.
            if "-M" not in set(passthrough):
                argv += ["-M", self.tree_method]
            if methods.mode == "msa":
                if "-A" not in set(passthrough):
                    argv += ["-A", methods.aligner or ORTHOFINDER.alignment]
                if "-T" not in set(passthrough):
                    argv += ["-T", methods.tree or ORTHOFINDER.tree]
            argv += passthrough

            # OrthoFinder runs its search jobs concurrently, so the OpenMP limit
            # is the per-job share rather than the whole machine.
            result = shell.run(
                argv,
                cwd=context.work_dir,
                log=context.log,
                timeout=None,
                threads=per_job,
                check=False,
            )

        if result.returncode != 0:
            raise shell.CommandFailed(argv, result.returncode, context.log)
        return output

    def _configure_search(self, root: Path, program: str, threads_per_search: int) -> str:
        """Register a thread-aware copy of the search method and return its name."""
        config = of_config.find_config(root)
        if config is None:
            return program
        return of_config.register_threaded_search(config, program, threads_per_search)

    def parse(self, native_root: Path, proteomes: Sequence[Proteome] | Proteome, *rest: Proteome) -> PredictorResult:
        loaded = (proteomes, *rest) if isinstance(proteomes, Proteome) else tuple(proteomes)
        pairwise = self._pairwise_tables(native_root)
        note = ""
        if pairwise:
            raw, sources = self._read_pairwise(pairwise, loaded)
        else:
            table = self._orthogroups_table(native_root)
            if table is None:
                raise AdapterError(
                    f"OrthoFinder produced no parseable output under {native_root}. "
                    f"Expected Orthologues/*__v__*.tsv or Orthogroups/Orthogroups.tsv; "
                    f"check the stage log."
                )
            raw, sources = self._read_orthogroups(table, loaded)
            note = (
                "orthogroup fallback: pairs expanded from clusters, not from "
                "OrthoFinder's pairwise orthologue tables"
            )

        result = build_result(self.name, self.label, raw, loaded, source_files=sources)
        result.note = note
        return result

    # ---- output discovery ---------------------------------------------

    def _pairwise_tables(self, native_root: Path) -> list[Path]:
        found = [p for p in native_root.rglob("*__v__*.tsv") if p.is_file()]
        preferred = [p for p in found if "Orthologues" in p.parts]
        return preferred or found

    def _orthogroups_table(self, native_root: Path) -> Path | None:
        for candidate in native_root.rglob("Orthogroups.tsv"):
            if candidate.is_file():
                return candidate
        for candidate in native_root.rglob("Orthogroups.csv"):
            if candidate.is_file():
                return candidate
        return None

    # ---- readers -------------------------------------------------------

    def _read_pairwise(
        self, tables: list[Path], proteomes: Sequence[Proteome]
    ) -> tuple[list[tuple[str, str]], list[Path]]:
        raw: list[tuple[str, str]] = []
        used: list[Path] = []
        wanted = {proteome.label for proteome in proteomes}

        for table in tables:
            stem = table.stem
            if "__v__" in stem:
                left_species, right_species = stem.split("__v__", 1)
                if wanted and not {left_species, right_species} <= wanted:
                    continue
            used.append(table)
            with table.open(newline="") as handle:
                reader = csv.reader(handle, delimiter="\t")
                header = next(reader, None)
                if header is None:
                    continue
                for row in reader:
                    if len(row) < 3:
                        continue
                    left_genes = _split_genes(row[1])
                    right_genes = _split_genes(row[2])
                    for left in left_genes:
                        for right in right_genes:
                            raw.append((left, right))
        return raw, used

    def _read_orthogroups(
        self, table: Path, proteomes: Sequence[Proteome]
    ) -> tuple[list[tuple[str, str]], list[Path]]:
        raw: list[tuple[str, str]] = []
        delimiter = "\t" if table.suffix == ".tsv" else ","
        with table.open(newline="") as handle:
            reader = csv.reader(handle, delimiter=delimiter)
            header = next(reader, None)
            if header is None:
                return raw, [table]
            for row in reader:
                genes = [gene for cell in row[1:] for gene in _split_genes(cell)]
                groups = [[g for g in genes if g in proteome.ids] for proteome in proteomes]
                for first, second in combinations(groups, 2):
                    for left in first:
                        for right in second:
                            raw.append((left, right))
        return raw, [table]

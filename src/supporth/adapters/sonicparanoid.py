"""Build SonicParanoid from source in a private package prefix. A conda-style interpreter lets upstream resolve helper binaries from PATH."""

from __future__ import annotations

import os
import shutil
import sys
import tempfile
from collections.abc import Sequence
from pathlib import Path

from ..canonical import PredictorResult, build_result
from ..proteomes import Proteome
from .. import interpreters, shell
from ..method_stack import SONICPARANOID
from .base import Adapter, AdapterError, MethodProfile, StageContext
from .tabular import id_sets, pairs_from_table


class SonicParanoid(Adapter):
    name = "sonicparanoid"
    label = "S"
    repo = "https://gitlab.com/salvo981/sonicparanoid2.git"
    # Pin the tested source release.
    default_ref = "2.0.9"
    entrypoint_candidates = (
        "site-packages/bin/sonicparanoid",
        ".venv/bin/sonicparanoid",
        "bin/sonicparanoid",
    )
    # blastp and mcl are checked at start-up even when MMseqs2 does the search.
    binaries = ("mmseqs", "diamond", "blastp", "mcl")
    # 2.0.9 rejects anything outside this range, and macOS still ships 3.9.
    python_requires = ">=3.10,<3.13"
    prefers_conda_interpreter = True
    # Pandas 3 exposes read-only arrays that the 2.0.9 Cython code cannot accept.
    pip_constraints = ("pandas<3",)
    install_notes = (
        "installed with pip --target under the tool directory, using a borrowed "
        "interpreter so its helper binaries are taken from the shared bin; "
        "compiled extensions require a C++ toolchain at install time"
    )

    sensitivity = "default"

    # Written once the profile database is usable, so an interrupted index is
    # never mistaken for a finished one.
    PFAM_MARKER = ".supporth-pfam-ready"

    PFAM_BASENAME = "pfama.mmseqs"
    # Every file SonicParanoid demands before it will use a profile database.
    # It checks these itself and deletes the whole directory when one is
    # absent; validate the complete set before adoption.
    PFAM_DATABASE_FILES = (
        "",
        ".dbtype",
        ".version",
        ".index",
        "_h",
        "_h.dbtype",
        "_h.index",
    )
    PFAM_INDEX_FILES = (".idx", ".idx.index", ".idx.dbtype")

    def post_install(
        self,
        root: Path,
        *,
        interpreter: Path | None = None,
        pfam_db: Path | None = None,
        log=print,
    ) -> None:
        chosen = interpreter or Path(self.python())
        if not interpreters.is_conda_style(chosen):
            log(
                "  note: this interpreter is not conda-style, so SonicParanoid will "
                "require exact binary versions in its own package directory"
            )
        self.install_to_prefix(root, chosen, log=log)
        if pfam_db is not None:
            self.adopt_pfam_profiles(root, pfam_db, log=log)
        else:
            self.prepare_pfam_profiles(root, chosen, log=log)

    def pfam_files(self) -> tuple[str, ...]:
        return tuple(
            f"{self.PFAM_BASENAME}{suffix}"
            for suffix in (*self.PFAM_DATABASE_FILES, *self.PFAM_INDEX_FILES)
        )

    def missing_pfam_files(self, source: Path) -> list[str]:
        """Which required files a candidate profile database lacks."""
        return [name for name in self.pfam_files() if not (source / name).is_file()]

    def resolve_pfam_source(self, source: Path) -> Path:
        """Locate a complete profile database, rejecting missing files."""
        source = source.expanduser().resolve()
        if source.name != "profile_db" and (source / "profile_db").is_dir():
            source = source / "profile_db"
        if not source.is_dir():
            raise AdapterError(f"no such profile database directory: {source}")

        missing = self.missing_pfam_files(source)
        if missing:
            partial = sorted(p.name for p in source.glob(f"{self.PFAM_BASENAME}.idx.*"))
            hint = (
                f"\n  found instead: {', '.join(partial)} -- these are shards of an "
                "unfinished index, which SonicParanoid discards"
                if partial
                else ""
            )
            raise AdapterError(
                f"{source} is not a complete profile database; missing "
                f"{len(missing)} file(s): {', '.join(missing)}{hint}"
            )
        return source

    def pfam_source_problem(self, source: Path) -> str | None:
        """Return a validation error for the profile database, or None."""
        try:
            self.resolve_pfam_source(source)
        except AdapterError as error:
            return str(error)
        return None

    def adopt_pfam_profiles(self, root: Path, source: Path, *, log=print) -> None:
        """Hardlink or copy a complete indexed database. Index compatibility is checked by the subsequent MMseqs2 search."""
        source = self.resolve_pfam_source(source)
        destination = self.profile_db(root)
        if source == destination.resolve():
            (destination / self.PFAM_MARKER).write_text(f"adopted from {source}\n")
            log(f"  using the existing PfamA profile DB at {source}")
            return
        if destination.exists():
            shutil.rmtree(destination)
        destination.mkdir(parents=True)

        linked = 0
        for name in self.pfam_files():
            target = destination / name
            try:
                os.link(source / name, target)
                linked += 1
            except OSError:
                # Cross-filesystem hardlinks require a copy.
                shutil.copyfile(source / name, target)

        how = "hardlinked" if linked == len(self.pfam_files()) else "copied"
        log(f"  {how} an existing PfamA profile DB from {source}")
        (destination / self.PFAM_MARKER).write_text(f"adopted from {source}\n")

    def prepare_pfam_profiles(
        self,
        root: Path,
        interpreter: Path,
        *,
        threads: int | None = None,
        log=print,
    ) -> None:
        """Extract and index the bundled Pfam profiles before the first analysis."""
        profile_db = self.profile_db(root)
        if self.pfam_ready(root):
            log("  PfamA profile DB already prepared")
            return

        # No need to clear a partial database first: SonicParanoid discards one
        # itself, and doing it here would destroy an index being built by
        # another process.
        threads = threads or shell.cpu_count()
        log(f"  indexing the bundled PfamA profile DB with {threads} threads (one-time)")
        # Run from a neutral directory: the source clone at the tool root holds
        # a sonicparanoid/ package without its compiled extensions, and it wins
        # the import over the installed one whenever it is the working directory.
        with tempfile.TemporaryDirectory(prefix="supporth-pfam-") as scratch:
            shell.run(
                [
                    str(interpreter),
                    "-c",
                    "from sonicparanoid.archortho import check_pfama_profile_db as check; "
                    f"print(check(threads={threads}))",
                ],
                cwd=Path(scratch),
                env=self.stage_env(root, threads),
                log=root / "pfam-index.log",
                timeout=None,
            )
        missing = self.missing_pfam_files(profile_db)
        if missing:
            raise AdapterError(
                f"indexing the PfamA profile DB failed; missing {', '.join(missing)} "
                f"(see {root / 'pfam-index.log'} and "
                f"{profile_db / 'log.obtain_precomputed_profiles.txt'})"
            )
        (profile_db / self.PFAM_MARKER).write_text("prepared by suppOrth\n")

    def check(self, options: list[str] | None = None) -> tuple[bool, str]:
        ready, detail = super().check(options)
        if not ready:
            return ready, detail
        root = self.installed_root()
        if not self.pfam_ready(root):
            return False, (
                "PfamA profile DB not indexed; domain-aware orthology cannot run "
                "(suppOrth install --tools sonicparanoid, or --pfam-db PATH to "
                "adopt an existing indexed database)"
            )
        return ready, detail

    def pfam_ready(self, root: Path) -> bool:
        """Require the completion marker and every database and index file."""
        profile_db = self.profile_db(root)
        if not (profile_db / self.PFAM_MARKER).exists():
            return False
        return not self.missing_pfam_files(profile_db)

    def profile_db(self, root: Path) -> Path:
        return root / "site-packages" / "sonicparanoid" / "pfam_files" / "profile_db"

    def leftovers(self, root: Path) -> list[tuple[Path, str]]:
        if self.pfam_ready(root):
            return []
        profile_db = self.profile_db(root)
        if not profile_db.is_dir():
            return []
        # MMseqs2 writes the index as numbered shards and merges them into
        # '.idx' only at the end. Without that merge SonicParanoid throws the
        # whole database away and starts again, so the shards of an interrupted
        # index are never picked back up.
        return [
            (path, "SonicParanoid partial Pfam index, discarded on the next attempt")
            for path in sorted(profile_db.glob("pfama.mmseqs.idx*"))
        ]

    def stage_env(self, root: Path, threads: int | None = None) -> dict[str, str]:
        """Environment exposing the private packages to the borrowed interpreter."""
        target = root / "site-packages"
        env = shell.tool_env(threads=threads)
        existing = env.get("PYTHONPATH", "")
        env["PYTHONPATH"] = f"{target}:{existing}" if existing else str(target)
        # macOS Python 3.8+ starts multiprocessing with 'spawn'. SonicParanoid's
        # Cython workers pickle a local wrapper that spawn cannot serialise.
        # fork avoids that; the ObjC flag stops the runtime aborting the fork.
        if sys.platform == "darwin":
            env["OBJC_DISABLE_INITIALIZE_FORK_SAFETY"] = "YES"
        return env

    # Run *before* sonicparanoid is imported so set_start_method still works.
    _DARWIN_FORK_LAUNCHER = (
        "import multiprocessing, sys\n"
        "multiprocessing.set_start_method('fork', force=True)\n"
        "from sonicparanoid.sonic_paranoid import main\n"
        "sys.argv[0] = 'sonicparanoid'\n"
        "raise SystemExit(main())\n"
    )

    def method_profile(self, context: StageContext) -> MethodProfile:
        return MethodProfile(
            search="MMseqs2 (via SonicParanoid)",
            clustering="MCL / InParanoid-style graph with adaptive scoring",
            alignment=SONICPARANOID.alignment,
            tree=SONICPARANOID.tree,
            orthology=SONICPARANOID.orthology,
            extra={"sensitivity": self.sensitivity},
        )

    def run(self, context: StageContext) -> Path:
        root = self.installed_root()
        entrypoint = self.locate_entrypoint(root)

        output = context.work_dir / "results"
        if context.previous_root is not None and not (output / "snapshot.tsv").is_file():
            raise AdapterError(f"SonicParanoid update requires {output / 'snapshot.tsv'}")
        output.mkdir(parents=True, exist_ok=True)

        python = str(self.prefix_interpreter(root) or self.python())
        flags = [
            "-i",
            str(context.input_dir),
            "-o",
            str(output),
            "-t",
            str(context.threads),
            *context.options,
        ]
        # Export the final pairwise predictions, rather than expanding groups.
        if not any(o in context.options for o in ("-op", "--output-pairs")):
            flags.append("-op")
        if sys.platform == "darwin":
            argv = [python, "-c", self._DARWIN_FORK_LAUNCHER, *flags]
        else:
            argv = [python, str(entrypoint), *flags]
        shell.run(
            argv,
            cwd=context.work_dir,
            log=context.log,
            env=self.stage_env(root, context.threads),
            timeout=None,
            threads=context.threads,
        )
        info = output / "last_run_info.txt"
        if info.is_file():
            for line in info.read_text().splitlines():
                if line.startswith("Run directory:"):
                    run_dir = Path(line.split("\t", 1)[1].strip())
                    if run_dir.is_dir() and run_dir.resolve().parent == (output / "runs").resolve():
                        return run_dir
        raise AdapterError(f"SonicParanoid did not record its completed run directory in {info}")

    def parse(self, native_root: Path, proteomes: Sequence[Proteome] | Proteome, *rest: Proteome) -> PredictorResult:
        loaded = (proteomes, *rest) if isinstance(proteomes, Proteome) else tuple(proteomes)
        tables = self._result_tables(native_root)
        if not tables:
            raise AdapterError(
                f"SonicParanoid produced no ortholog table under {native_root}. "
                f"Expected something under runs/*/ortholog_relations/; check the stage log."
            )
        raw: list[tuple[str, str]] = []
        for table in tables:
            raw.extend(pairs_from_table(table, *id_sets(loaded)))
        return build_result(self.name, self.label, raw, loaded, source_files=tables)

    # SonicParanoid writes its pairwise predictions to
    # 'pairwise_orthologs.<run id>.tsv'. Anything matched later is a fallback
    # for other versions: the group table is derived from the pairs and is not
    # equivalent to them, having lost some pairs by the time it is written.
    TABLE_PATTERNS = (
        "pairwise_orthologs*",
        "*ortholog*pairs*",
        "*pairwise*ortholog*",
    )
    GROUP_PATTERNS = ("*ortholog*groups*",)
    # Genes left out of every group, and a flattened restatement of the groups:
    # neither adds a prediction, and both would be counted twice.
    TABLE_EXCLUDES = ("not_assigned", "flat", "stats")

    def _result_tables(self, native_root: Path) -> list[Path]:
        for pattern in self.TABLE_PATTERNS:
            found = self._matching(native_root, pattern)
            if found:
                return found
        for pattern in self.GROUP_PATTERNS:
            found = self._matching(native_root, pattern)
            if found:
                return found
        return []

    def _matching(self, native_root: Path, pattern: str) -> list[Path]:
        return sorted(
            p
            for p in native_root.rglob(pattern)
            if p.is_file()
            and p.suffix in {".tsv", ".txt", ""}
            and not any(word in p.name for word in self.TABLE_EXCLUDES)
        )

"""Register OrthoFinder search commands with a thread allocation based on the number of proteomes."""

from __future__ import annotations

import fcntl
import json
import shutil
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator

PREFIX = "supporth_"

# Mirrors of the stock entries, differing only in the thread count. Keeping the
# rest of each command identical means results stay comparable to a stock run.
SEARCH_TEMPLATES: dict[str, dict[str, str]] = {
    "diamond": {
        "db_cmd": "diamond makedb --ignore-warnings --in INPUT -d OUTPUT",
        "search_cmd": (
            "diamond blastp --ignore-warnings -d DATABASE -q INPUT -o OUTPUT "
            "--more-sensitive -p {threads} --quiet -e 0.001 --compress 1"
        ),
    },
    "diamond_ultra_sens": {
        "db_cmd": "diamond makedb --ignore-warnings --in INPUT -d OUTPUT",
        "search_cmd": (
            "diamond blastp --ignore-warnings -d DATABASE -q INPUT -o OUTPUT "
            "--ultra-sensitive -p {threads} --quiet -e 0.001 --compress 1"
        ),
    },
    "blast": {
        "db_cmd": "makeblastdb -dbtype prot -in INPUT -out OUTPUT",
        "search_cmd": (
            "blastp -outfmt 6 -evalue 0.001 -num_threads {threads} "
            "-query INPUT -db DATABASE -out OUTPUT"
        ),
    },
    "mmseqs": {
        "db_cmd": "mmseqs createdb INPUT OUTPUT.fa ; mmseqs createindex OUTPUT.fa /tmp",
        "search_cmd": (
            "mmseqs search PATH/mmseqsDBBASENAME DATABASE.fa OUTPUT.db /tmp/tmpBASEOUTNAME "
            "--threads {threads} ; mmseqs convertalis --threads {threads} "
            "PATH/mmseqsDBBASENAME DATABASE.fa OUTPUT.db OUTPUT"
        ),
    },
}

SEARCH_BINARIES: dict[str, tuple[str, ...]] = {
    "diamond": ("diamond",),
    "diamond_ultra_sens": ("diamond",),
    "blast": ("blastp", "makeblastdb"),
    "mmseqs": ("mmseqs",),
}


def search_job_count(n_species: int = 2, *, double_blast: bool = True) -> int:
    """Number of all-versus-all search jobs OrthoFinder will create.

    With ``double_blast`` (the default) both directions of every species pair
    are searched, including self-comparisons, giving ``n**2`` jobs. The ``-1``
    option searches one direction only, giving ``n * (n + 1) / 2``.
    """
    if double_blast:
        return n_species * n_species
    return n_species * (n_species + 1) // 2


def threads_per_job(total_threads: int, jobs: int) -> int:
    """Allocate threads among search jobs, with at least one thread per job."""
    return max(1, total_threads // max(1, jobs))


def find_config(root: Path) -> Path | None:
    """Locate config.json inside an OrthoFinder source tree."""
    preferred = root / "scripts_of" / "config.json"
    if preferred.is_file():
        return preferred
    for candidate in sorted(root.rglob("config.json")):
        if candidate.is_file():
            return candidate
    return None


def variant_name(program: str) -> str:
    return f"{PREFIX}{program}"


def binaries_for_search(program: str) -> tuple[str, ...]:
    base = program[len(PREFIX):] if program.startswith(PREFIX) else program
    return SEARCH_BINARIES.get(base, ())


@contextmanager
def config_lock(root: Path) -> Iterator[None]:
    """Serialise access to one OrthoFinder install.

    The thread count is baked into config.json, so two SuppOrth runs sharing an
    install must not configure it at the same time.
    """
    lock_path = root / ".supporth.lock"
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("w") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def register_threaded_search(config_path: Path, program: str, threads: int) -> str:
    """Register a threaded variant, or preserve unrecognised search methods."""
    base = program[len(PREFIX):] if program.startswith(PREFIX) else program
    template = SEARCH_TEMPLATES.get(base)
    if template is None:
        return program

    name = variant_name(base)
    backup = config_path.with_suffix(".json.supporth-orig")
    if not backup.exists():
        shutil.copy2(config_path, backup)

    with config_path.open() as handle:
        config = json.load(handle)

    config[name] = {
        "program_type": "search",
        "db_cmd": template["db_cmd"],
        "search_cmd": template["search_cmd"].format(threads=threads),
    }

    temporary = config_path.with_suffix(".json.supporth-tmp")
    with temporary.open("w") as handle:
        json.dump(config, handle, indent=4)
        handle.write("\n")
    temporary.replace(config_path)
    return name


def strip_option(options: list[str], flag: str) -> list[str]:
    """Remove a flag and its value from pass-through options."""
    out: list[str] = []
    skip = False
    for token in options:
        if skip:
            skip = False
            continue
        if token == flag:
            skip = True
            continue
        out.append(token)
    return out

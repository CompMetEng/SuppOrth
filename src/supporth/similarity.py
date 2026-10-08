"""Read alignment statistics from OrthoFinder or run a separate BLAST/DIAMOND search."""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from collections.abc import Sequence
from pathlib import Path

from . import shell
from .proteomes import Proteome, as_proteomes

OUTFMT = ("qseqid", "sseqid", "pident", "evalue", "bitscore")


@dataclass
class Hit:
    identity: float
    evalue: float
    bitscore: float


class SimilarityIndex:
    """Best alignment per unordered protein pair."""

    def __init__(self) -> None:
        self._hits: dict[tuple[str, str], Hit] = {}

    @staticmethod
    def _key(a: str, b: str) -> tuple[str, str]:
        return (a, b) if a <= b else (b, a)

    def add(self, query: str, subject: str, hit: Hit) -> None:
        key = self._key(query, subject)
        existing = self._hits.get(key)
        if existing is None or hit.bitscore > existing.bitscore:
            self._hits[key] = hit

    def get(self, a: str, b: str) -> Hit | None:
        return self._hits.get(self._key(a, b))

    def __len__(self) -> int:
        return len(self._hits)

    def load_tabular(self, path: Path) -> None:
        with Path(path).open(errors="replace") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 5:
                    continue
                query, subject = fields[0], fields[1]
                try:
                    hit = Hit(float(fields[2]), float(fields[3]), float(fields[4]))
                except ValueError:
                    continue
                self.add(query, subject, hit)


def index_from_orthofinder(
    work_dir: Path,
    proteomes: Sequence[Proteome] | Proteome,
    *rest: Proteome,
    log=print,
) -> SimilarityIndex:
    """Read cross-species search tables and translate OrthoFinder sequence identifiers."""
    loaded = as_proteomes(proteomes, *rest)
    work_dir = Path(work_dir).expanduser()
    if (work_dir / "WorkingDirectory").is_dir():
        work_dir = work_dir / "WorkingDirectory"

    species = _orthofinder_species(work_dir / "SpeciesIDs.txt")
    wanted = {}
    for index, filename in species.items():
        for proteome in loaded:
            if filename == proteome.fasta.name or filename == f"{proteome.label}.fasta":
                wanted[index] = proteome
    if len(wanted) != len(loaded):
        names = ", ".join(p.label for p in loaded)
        raise ValueError(
            f"{work_dir / 'SpeciesIDs.txt'} lists {sorted(species.values())}, which does "
            f"not cover {names}"
        )

    names = _orthofinder_sequence_names(work_dir / "SequenceIDs.txt")
    index = SimilarityIndex()
    for first in sorted(wanted):
        for second in sorted(wanted):
            if first == second:
                continue
            table = _orthofinder_blast_table(work_dir, first, second)
            if table is None:
                continue
            log(f"  reading {table.name}")
            _load_orthofinder_table(table, names, index)
    return index


def _orthofinder_blast_table(work_dir: Path, first: str, second: str) -> Path | None:
    for name in (f"Blast{first}_{second}.txt", f"Blast{first}_{second}.txt.gz"):
        path = work_dir / name
        if path.is_file():
            return path
    return None


def _orthofinder_species(path: Path) -> dict[str, str]:
    species: dict[str, str] = {}
    with Path(path).open(errors="replace") as handle:
        for line in handle:
            if ":" not in line:
                continue
            key, filename = line.split(":", 1)
            species[key.strip()] = filename.strip()
    return species


def _orthofinder_sequence_names(path: Path) -> dict[str, str]:
    """Internal ``species_sequence`` identifiers mapped to real gene names."""
    names: dict[str, str] = {}
    with Path(path).open(errors="replace") as handle:
        for line in handle:
            if ":" not in line:
                continue
            key, description = line.split(":", 1)
            fields = description.split()
            if fields:
                names[key.strip()] = fields[0]
    return names


def _load_orthofinder_table(
    path: Path,
    names: dict[str, str],
    index: SimilarityIndex,
) -> None:
    with _open_table(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            query = names.get(fields[0])
            subject = names.get(fields[1])
            if query is None or subject is None:
                continue
            # Standard BLAST format 6: identity, e-value and bit score.
            try:
                hit = Hit(float(fields[2]), float(fields[10]), float(fields[11]))
            except ValueError:
                continue
            index.add(query, subject, hit)


def _open_table(path: Path):
    path = Path(path)
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", errors="replace")
    return path.open(errors="replace")


def build_index(
    proteomes: Sequence[Proteome] | Proteome,
    work_dir: Path,
    *rest: Proteome,
    threads: int = 1,
    program: str = "diamond",
    evalue: float = 1e-3,
    max_targets: int = 10,
    log=print,
) -> SimilarityIndex:
    loaded = as_proteomes(proteomes, *rest)
    work_dir.mkdir(parents=True, exist_ok=True)
    index = SimilarityIndex()

    if program == "none":
        return index
    if program not in {"diamond", "blast"}:
        raise ValueError(f"unsupported alignment program: {program}")

    runner = _diamond if program == "diamond" else _blastp
    for subject in loaded:
        database = _make_database(subject, work_dir, threads, program)
        for query in loaded:
            if query is subject:
                continue
            tag = f"{query.label}__{subject.label}"
            output = work_dir / f"{tag}.tsv"
            log(f"  {program} {query.label} vs {subject.label}")
            runner(query, database, work_dir, output, threads, evalue, max_targets, tag)
            index.load_tabular(output)
    return index


def _make_database(subject: Proteome, work_dir: Path, threads: int, program: str) -> Path:
    """Build one database per subject proteome for this alignment run."""
    database = work_dir / subject.label
    if program == "diamond":
        binary = shell.which("diamond")
        if binary is None:
            raise RuntimeError("diamond not available in the shared bin")
        argv = [binary, "makedb", "--in", str(subject.fasta), "-d", str(database),
                "--threads", str(threads), "--quiet"]
    else:
        binary = shell.which("makeblastdb")
        if binary is None:
            raise RuntimeError("makeblastdb not available in the shared bin")
        argv = [binary, "-in", str(subject.fasta), "-dbtype", "prot", "-out", str(database)]
    shell.run(argv, log=work_dir / f"{subject.label}.makedb.log", threads=threads)
    return database


def _diamond(
    query: Proteome,
    database: Path,
    work_dir: Path,
    output: Path,
    threads: int,
    evalue: float,
    max_targets: int,
    tag: str,
) -> None:
    binary = shell.which("diamond")
    if binary is None:
        raise RuntimeError("diamond not available in the shared bin")

    shell.run(
        [
            binary,
            "blastp",
            "-q",
            str(query.fasta),
            "-d",
            str(database),
            "-o",
            str(output),
            "--outfmt",
            "6",
            *OUTFMT,
            "--evalue",
            str(evalue),
            "--max-target-seqs",
            str(max_targets),
            "--threads",
            str(threads),
            "--quiet",
        ],
        log=work_dir / f"{tag}.search.log",
        threads=threads,
    )


def _blastp(
    query: Proteome,
    database: Path,
    work_dir: Path,
    output: Path,
    threads: int,
    evalue: float,
    max_targets: int,
    tag: str,
) -> None:
    blastp = shell.which("blastp")
    if blastp is None:
        raise RuntimeError("blastp not available in the shared bin")
    shell.run(
        [
            blastp,
            "-query",
            str(query.fasta),
            "-db",
            str(database),
            "-out",
            str(output),
            "-outfmt",
            "6 " + " ".join(OUTFMT),
            "-evalue",
            str(evalue),
            "-max_target_seqs",
            str(max_targets),
            "-num_threads",
            str(threads),
        ],
        log=work_dir / f"{tag}.search.log",
        threads=threads,
    )

"""Run predictors sequentially, import OMA output and write the consensus."""

from __future__ import annotations

import hashlib
import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import pandas as pd

from . import canonical, consensus, figures, manifest, paths, similarity, shell
from .adapters import Adapter, StageContext, resolve
from .adapters.oma import OMA
from .adapters.orthofinder import find_alignment_dir
from .proteomes import Proteome, load_proteomes, stage_input_dir


class PipelineError(RuntimeError):
    pass


@dataclass
class RunConfig:
    outdir: Path
    fastas: tuple[Path, ...] = ()
    fasta1: Path | None = None
    fasta2: Path | None = None
    tools: str = "all"
    labels: tuple[str, ...] | None = None
    threads: int = 0
    alignment: str = "diamond"
    resume: bool = False
    support_lv: int = 3
    identity_filter: float = 80.0
    collapse: tuple[Path, Path] | None = None
    tool_options: dict[str, list[str]] = field(default_factory=dict)
    oma: Path | None = None

    def __post_init__(self) -> None:
        if self.fastas:
            self.fastas = tuple(Path(path) for path in self.fastas)
        elif self.fasta1 is not None and self.fasta2 is not None:
            self.fastas = (Path(self.fasta1), Path(self.fasta2))
        if len(self.fastas) < 2:
            raise ValueError("need at least two FASTA files")
        if self.alignment in ("blast", "diamond"):
            options = list(self.tool_options.get("orthofinder", []))
            if "-S" not in options:
                self.tool_options["orthofinder"] = ["-S", self.alignment, *options]


@dataclass
class StageReport:
    tool: str
    label: str
    status: str
    seconds: float
    n_pairs: int
    method: dict
    commit: str
    dropped_unknown: int = 0
    dropped_same_species: int = 0
    note: str = ""
    error: str = ""


def warn_support_threshold(level: int, available: int, log: Callable[[str], None]) -> None:
    if level > available:
        log(
            f"WARNING: --supportLv {level} exceeds the {available} available predictors; "
            "only pairs passing the identity exception can be retained."
        )


def run_pipeline(config: RunConfig, log: Callable[[str], None] = print) -> dict:
    started = time.monotonic()
    paths.ensure_layout()

    outdir = Path(config.outdir).resolve()
    (outdir / "pairs").mkdir(parents=True, exist_ok=True)
    work_root = outdir / "work"
    work_root.mkdir(parents=True, exist_ok=True)

    proteomes = load_proteomes(config.fastas, config.labels)
    for index, proteome in enumerate(proteomes, start=1):
        log(f"Species {index}: {proteome.label} ({proteome.size} proteins)")

    adapters = resolve(config.tools)
    threads = shell.cpu_count(config.threads or None)

    _preflight(adapters, config.tool_options, log)
    selected_count = len(adapters) + int(config.oma is not None)
    warn_support_threshold(config.support_lv, selected_count, log)

    input_dir = stage_input_dir(work_root / "input", proteomes)

    reports: list[StageReport] = []
    results: list[canonical.PredictorResult] = []
    alignment_dir = find_alignment_dir(work_root / "orthofinder")

    for index, adapter in enumerate(adapters, start=1):
        log(f"\n[{index}/{len(adapters)}] {adapter.name}")
        pairs_json = outdir / "pairs" / f"{adapter.name}.json"

        if config.resume and pairs_json.exists():
            result = canonical.read_json(pairs_json)
            results.append(result)
            reports.append(
                StageReport(
                    tool=adapter.name,
                    label=adapter.label,
                    status="reused",
                    seconds=0.0,
                    n_pairs=len(result.pairs),
                    method={},
                    commit=_commit(adapter),
                )
            )
            log(f"  reusing {pairs_json.name} ({len(result.pairs)} pairs)")
            if alignment_dir is None:
                alignment_dir = find_alignment_dir(work_root / adapter.name)
            continue

        context = StageContext(
            input_dir=input_dir,
            work_dir=work_root / adapter.name,
            proteomes=proteomes,
            threads=threads,
            options=config.tool_options.get(adapter.name, []),
            alignment_dir=alignment_dir,
        )
        context.work_dir.mkdir(parents=True, exist_ok=True)

        stage_started = time.monotonic()
        try:
            native_root = adapter.run(context)
            result = adapter.parse(native_root, proteomes)
        except Exception as error:  # noqa: BLE001 - retain completed stages after a predictor failure
            elapsed = time.monotonic() - stage_started
            log(f"  FAILED after {elapsed / 60:.1f} min: {error}")
            log(f"  log: {context.log}")
            reports.append(
                StageReport(
                    tool=adapter.name,
                    label=adapter.label,
                    status="failed",
                    seconds=elapsed,
                    n_pairs=0,
                    method=adapter.method_profile(context).as_dict(),
                    commit=_commit(adapter),
                    error=str(error),
                )
            )
            continue

        elapsed = time.monotonic() - stage_started
        canonical.write_json(result, pairs_json)
        canonical.write_legacy_json(result, outdir / "pairs" / f"{adapter.name}_result_dict.json")
        results.append(result)
        reports.append(
            StageReport(
                tool=adapter.name,
                label=adapter.label,
                status="ok",
                seconds=elapsed,
                n_pairs=len(result.pairs),
                method=adapter.method_profile(context).as_dict(),
                commit=_commit(adapter),
                dropped_unknown=result.dropped_unknown,
                dropped_same_species=result.dropped_same_species,
                note=result.note,
            )
        )
        log(f"  {result.summary()} in {elapsed / 60:.1f} min")
        found = find_alignment_dir(native_root)
        if found is not None:
            alignment_dir = found

    if config.oma:
        oma_result, oma_report = _ingest_oma(config, proteomes, outdir, log)
        reports.append(oma_report)
        if oma_result is not None:
            results.append(oma_result)

    if not results:
        raise PipelineError(
            "no predictor produced a usable result; see the stage logs under "
            f"{work_root}"
        )

    successful_count = len({result.tool for result in results})
    if successful_count < selected_count:
        warn_support_threshold(config.support_lv, successful_count, log)

    log("\nAlignment statistics")
    if config.alignment == "none":
        index = similarity.SimilarityIndex()
    elif alignment_dir is not None:
        log("  reading OrthoFinder alignment tables")
        index = similarity.index_from_orthofinder(alignment_dir, proteomes, log=log)
    else:
        log(f"  OrthoFinder did not write alignment tables; running {config.alignment}")
        index = similarity.build_index(
            proteomes,
            work_root / "alignment",
            threads=threads,
            program=config.alignment,
            log=log,
        )
    log(f"  {len(index)} aligned pairs")

    table = consensus.build_table(results, index, proteomes)
    unfiltered = len(table)
    table = consensus.filter_table(
        table,
        support_lv=config.support_lv,
        identity_filter=config.identity_filter,
    )
    log(
        f"  kept {len(table)} of {unfiltered} pairs "
        f"(Lv_support ≥ {config.support_lv} or Identity > {config.identity_filter})"
    )

    output_tsv = outdir / "suppOrth.tsv"
    table.to_csv(output_tsv, sep="\t", index=False)

    summary = consensus.agreement_summary(table)
    summary.to_csv(outdir / "agreement_summary.tsv", sep="\t", index=False)

    figure_path = figures.write_support_figure(
        results, outdir / "suppOrth.png", log=log
    )

    collapsed_path = None
    if config.collapse:
        if len(proteomes) != 2:
            raise PipelineError("--collapse is only defined for exactly two proteomes")
        collapsed = consensus.collapse_to_genes(table, *config.collapse)
        collapsed_path = outdir / "suppOrth.genes.tsv"
        collapsed.to_csv(collapsed_path, sep="\t", index=False)
        log(f"Collapsed to {len(collapsed)} gene pairs")

    provenance = _provenance(config, proteomes, reports, table, threads)
    with (outdir / "provenance.json").open("w") as handle:
        json.dump(provenance, handle, indent=2)
        handle.write("\n")

    log(f"\nTotal runtime {(time.monotonic() - started) / 60:.1f} min")
    log(f"Consensus table: {output_tsv}")
    _log_summary(summary, log)

    return {
        "table": table,
        "summary": summary,
        "reports": reports,
        "output": output_tsv,
        "collapsed": collapsed_path,
        "figure": figure_path,
        "provenance": provenance,
    }


def _ingest_oma(
    config: RunConfig,
    proteomes: tuple[Proteome, ...],
    outdir: Path,
    log: Callable[[str], None],
) -> tuple[canonical.PredictorResult | None, StageReport]:
    """Parse a user-supplied OMA standalone run; never launch OMA."""
    adapter = OMA()
    source = Path(config.oma).expanduser().resolve()
    pairs_json = outdir / "pairs" / "oma.json"
    log(f"\n[oma] reading user-supplied run at {source}")

    if config.resume and pairs_json.exists():
        result = canonical.read_json(pairs_json)
        log(f"  reusing {pairs_json.name} ({len(result.pairs)} pairs)")
        return result, StageReport(
            tool=adapter.name,
            label=adapter.label,
            status="reused",
            seconds=0.0,
            n_pairs=len(result.pairs),
            method={},
            commit="user-supplied",
        )

    started = time.monotonic()
    try:
        result = adapter.parse(source, proteomes)
    except Exception as error:  # noqa: BLE001 - retain other predictions if the import fails
        elapsed = time.monotonic() - started
        log(f"  FAILED: {error}")
        return None, StageReport(
            tool=adapter.name,
            label=adapter.label,
            status="failed",
            seconds=elapsed,
            n_pairs=0,
            method=adapter.method_profile(
                StageContext(input_dir=source, work_dir=source, proteomes=proteomes)
            ).as_dict(),
            commit="user-supplied",
            error=str(error),
        )

    canonical.write_json(result, pairs_json)
    canonical.write_legacy_json(result, outdir / "pairs" / "oma_result_dict.json")
    elapsed = time.monotonic() - started
    log(f"  {result.summary()} in {elapsed:.1f} s")
    return result, StageReport(
        tool=adapter.name,
        label=adapter.label,
        status="ok",
        seconds=elapsed,
        n_pairs=len(result.pairs),
        method=adapter.method_profile(
            StageContext(input_dir=source, work_dir=source, proteomes=proteomes)
        ).as_dict(),
        commit="user-supplied",
        dropped_unknown=result.dropped_unknown,
        dropped_same_species=result.dropped_same_species,
        note=result.note,
    )


def _preflight(
    adapters: list[Adapter],
    tool_options: dict[str, list[str]],
    log: Callable[[str], None],
) -> None:
    """Check predictor dependencies before starting any stages."""
    log("\nChecking predictors")
    problems = []
    for adapter in adapters:
        ready, detail = adapter.check(tool_options.get(adapter.name, []))
        log(f"  {adapter.name:<15} {'ok' if ready else 'NOT READY'}  {detail}")
        if not ready:
            problems.append(f"{adapter.name}: {detail}")
    if problems:
        raise PipelineError(
            "predictors are not ready:\n  "
            + "\n  ".join(problems)
            + "\n\nRun: suppOrth install"
        )


def _commit(adapter: Adapter) -> str:
    record = manifest.predictor(adapter.name)
    return (record or {}).get("commit", "unknown")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _provenance(
    config: RunConfig,
    proteomes: tuple[Proteome, ...],
    reports: list[StageReport],
    table: pd.DataFrame,
    threads: int,
) -> dict:
    installed = manifest.load()
    return {
        "environment": manifest.environment(),
        "inputs": {
            proteome.label: {
                "path": str(proteome.fasta),
                "n_proteins": proteome.size,
                "sha256": _sha256(proteome.fasta),
            }
            for proteome in proteomes
        },
        "settings": {
            "tools": config.tools,
            "threads": threads,
            "alignment": config.alignment,
            "tool_options": config.tool_options,
            "support_lv": config.support_lv,
            "identity_filter": config.identity_filter,
            "oma": str(config.oma) if config.oma else None,
        },
        "predictors": [
            {
                "tool": report.tool,
                "label": report.label,
                "status": report.status,
                "commit": report.commit,
                "ref": (installed["predictors"].get(report.tool) or {}).get("ref", ""),
                "minutes": round(report.seconds / 60, 2),
                "n_pairs": report.n_pairs,
                "method": report.method,
                "dropped_unknown": report.dropped_unknown,
                "dropped_same_species": report.dropped_same_species,
                "note": report.note,
                "error": report.error,
            }
            for report in reports
        ],
        "shared_binaries": installed["binaries"],
        "result": {
            "n_pairs": int(len(table)),
            "by_support": {
                str(level): int(count)
                for level, count in table["Lv_support"].value_counts().sort_index().items()
            },
        },
    }


def _log_summary(summary: pd.DataFrame, log: Callable[[str], None]) -> None:
    log("\nAgreement between predictors:")
    for _, row in summary.iterrows():
        log(f"  {str(row['SupportedBy']):<40} {row['n_pairs']:>8} pairs")

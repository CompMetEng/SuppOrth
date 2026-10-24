"""Run predictors sequentially, import OMA output and write the consensus."""

from __future__ import annotations

import json
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import pandas as pd

from . import canonical, consensus, figures, manifest, paths, similarity, shell, updates
from .adapters import Adapter, StageContext, resolve
from .adapters.oma import OMA
from .adapters.orthofinder import find_alignment_dir
from .proteomes import Proteome, load_proteomes, stage_input_dir
from .updates import file_digest as _sha256


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
    update: bool = False
    support_lv: int = 3
    identity_filter: float = 80.0
    collapse: tuple[Path, Path] | None = None
    tool_options: dict[str, list[str]] = field(default_factory=dict)
    oma: Path | None = None

    def __post_init__(self) -> None:
        if self.resume and self.update:
            raise ValueError("--resume and --update cannot be combined")
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
    native_root: str = ""


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
    try:
        updates.recover_publication(outdir)
    except (updates.UpdateError, OSError, ValueError) as error:
        raise PipelineError(f"cannot recover interrupted update: {error}") from error
    (outdir / "pairs").mkdir(parents=True, exist_ok=True)
    work_root = outdir / "work"
    work_root.mkdir(parents=True, exist_ok=True)

    proteomes = load_proteomes(config.fastas, config.labels)
    for index, proteome in enumerate(proteomes, start=1):
        log(f"Species {index}: {proteome.label} ({proteome.size} proteins)")

    adapters = resolve(config.tools)
    previous = None
    added = ()
    config.outdir = outdir
    provenance_file = outdir / "provenance.json"
    saved = json.loads(provenance_file.read_text()) if provenance_file.is_file() else {}
    if config.resume and saved:
        current = {p.label: _sha256(p.fasta) for p in proteomes}
        recorded = {label: r.get("sha256") for label, r in saved.get("inputs", {}).items()}
        if current != recorded:
            raise PipelineError("inputs changed; use --update to add proteomes, or start a new run")
    if config.update:
        try:
            previous, added = updates.validate_update(config, proteomes, adapters)
        except updates.UpdateError as error:
            raise PipelineError(str(error)) from error
        # Retain the original species order when orienting unordered pairs.
        by_label = {p.label: p for p in proteomes}
        proteomes = tuple(by_label[label] for label in previous["inputs"]) + added
        log(f"Updating {len(previous['inputs'])} proteomes with {len(added)} additional proteome(s)")
    threads = shell.cpu_count(config.threads or None)

    _preflight(adapters, config.tool_options, log)
    selected_count = len(adapters) + int(config.oma is not None)
    warn_support_threshold(config.support_lv, selected_count, log)

    input_dir = stage_input_dir(work_root / "input", proteomes)
    added_input_dir = stage_input_dir(work_root / "added_input", added) if added else None

    reports: list[StageReport] = []
    results: list[canonical.PredictorResult] = []
    alignment_dir = find_alignment_dir(work_root / "orthofinder")
    saved_reports = {r["tool"]: r for r in saved.get("predictors", [])}
    if "orthofinder" in saved_reports:
        of_root = Path(saved_reports["orthofinder"].get("native_root") or work_root / "orthofinder" / "results")
        alignment_dir = find_alignment_dir(of_root)

    for index, adapter in enumerate(adapters, start=1):
        log(f"\n[{index}/{len(adapters)}] {adapter.name}")
        pairs_json = outdir / "pairs" / f"{adapter.name}.json"

        if config.resume and pairs_json.exists():
            saved_report = saved_reports.get(adapter.name, {})
            result = canonical.read_json(pairs_json)
            results.append(result)
            reports.append(
                StageReport(
                    tool=adapter.name,
                    label=adapter.label,
                    status="reused",
                    seconds=0.0,
                    n_pairs=len(result.pairs),
                    method=saved_report.get("method", {}),
                    commit=_commit(adapter),
                    native_root=saved_report.get("native_root", ""),
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
            previous_root=updates.native_root(previous, adapter.name, work_root / adapter.name)
            if previous else None,
            added_input_dir=added_input_dir,
        )
        context.work_dir.mkdir(parents=True, exist_ok=True)

        stage_started = time.monotonic()
        if config.update:
            log("  rerunning all four Broccoli steps on the expanded proteomes" if adapter.name == "broccoli"
                else "  reusing previous searches; computing missing comparisons")
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
            if config.update:
                raise PipelineError(
                    f"{adapter.name} update failed; consensus outputs were not replaced. "
                    f"See {context.log} and retry with --update."
                ) from error
            continue

        elapsed = time.monotonic() - stage_started
        if not config.update:
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
                native_root=str(native_root.resolve()),
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
            reuse_labels=frozenset(previous["inputs"]) if previous else frozenset(),
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

    output_dir = Path(tempfile.mkdtemp(prefix=".update-", dir=outdir)) if config.update else outdir
    output_tsv = output_dir / "suppOrth.tsv"
    if config.update:
        if any(report.status == "failed" for report in reports):
            raise PipelineError("update failed; consensus outputs were not replaced")
        for result in results:
            canonical.write_json(result, output_dir / "pairs" / f"{result.tool}.json")
            canonical.write_legacy_json(result, output_dir / "pairs" / f"{result.tool}_result_dict.json")
    table.to_csv(output_tsv, sep="\t", index=False)

    summary = consensus.agreement_summary(table)
    summary.to_csv(output_dir / "agreement_summary.tsv", sep="\t", index=False)

    figure_path = figures.write_support_figure(
        results, output_dir / "suppOrth.png", log=log
    )

    collapsed_path = None
    if config.collapse:
        if len(proteomes) != 2:
            raise PipelineError("--collapse is only defined for exactly two proteomes")
        collapsed = consensus.collapse_to_genes(table, *config.collapse)
        collapsed_path = output_dir / "suppOrth.genes.tsv"
        collapsed.to_csv(collapsed_path, sep="\t", index=False)
        log(f"Collapsed to {len(collapsed)} gene pairs")

    provenance = _provenance(config, proteomes, reports, table, threads)
    if previous:
        history = output_dir / "history"
        history.mkdir(exist_ok=True)
        number = 1
        snapshot = history / f"provenance_{number}.json"
        while (outdir / "history" / snapshot.name).exists():
            number += 1
            snapshot = history / f"provenance_{number}.json"
        snapshot.write_text(json.dumps(previous, indent=2) + "\n")
        provenance["update"] = {
            "added_proteomes": [p.label for p in added],
            "previous_provenance": str(outdir / "history" / snapshot.name),
            "broccoli": "full rerun" if any(a.name == "broccoli" for a in adapters) else "not selected",
        }
    with (output_dir / "provenance.json").open("w") as handle:
        json.dump(provenance, handle, indent=2)
        handle.write("\n")
    if config.update:
        try:
            updates.publish_outputs(output_dir, outdir)
        except OSError as error:
            raise PipelineError(f"could not publish updated outputs: {error}") from error
        output_tsv = outdir / "suppOrth.tsv"
        if figure_path:
            figure_path = outdir / figure_path.name

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

    if not config.update:
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
            "update": config.update,
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
                "native_root": report.native_root,
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

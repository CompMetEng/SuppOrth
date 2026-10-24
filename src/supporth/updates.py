"""Validate an expanded input collection before reusing predictor searches."""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

from . import manifest
from .proteomes import Proteome


class UpdateError(RuntimeError):
    pass


def file_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_update(config, proteomes: tuple[Proteome, ...], adapters) -> tuple[dict, tuple[Proteome, ...]]:
    """Require unchanged original inputs, methods and installed tool versions."""
    path = config.outdir / "provenance.json"
    try:
        previous = json.loads(path.read_text())
    except (OSError, ValueError) as error:
        raise UpdateError(f"--update requires a readable {path}") from error
    inputs = previous.get("inputs", {})
    if not inputs:
        raise UpdateError("the previous run has no recorded inputs")
    current = {p.label: p for p in proteomes}
    for label, record in inputs.items():
        if label not in current:
            raise UpdateError(f"previous proteome {label!r} is missing; supply all original and new FASTAs")
        if file_digest(current[label].fasta) != record.get("sha256"):
            raise UpdateError(f"previous proteome {label!r} has changed; start a new run")
    added = tuple(p for p in proteomes if p.label not in inputs)
    if not added:
        raise UpdateError("--update needs at least one additional proteome")

    records = {r["tool"]: r for r in previous.get("predictors", [])}
    selected = {a.name for a in adapters} | ({"oma"} if config.oma else set())
    if selected != set(records):
        raise UpdateError("--update requires the same predictors as the previous run")
    from .adapters.orthofinder import find_alignment_dir
    from .similarity import orthofinder_search_dirs

    for tool in selected - {"oma"}:
        work_dir = config.outdir / "work" / tool
        root = native_root(previous, tool, work_dir)
        if tool == "orthofinder":
            directory = find_alignment_dir(root)
            if directory is None:
                raise UpdateError(f"no reusable OrthoFinder searches under {root}")
            if not (directory.parent / "Log.txt").is_file():
                raise UpdateError(f"OrthoFinder update requires {directory.parent / 'Log.txt'}")
            try:
                orthofinder_search_dirs(directory)
            except ValueError as error:
                raise UpdateError(str(error)) from error
        elif tool == "sonicparanoid":
            for name in ("snapshot.tsv", "last_run_info.txt"):
                if not (work_dir / "results" / name).is_file():
                    raise UpdateError(f"SonicParanoid update requires {work_dir / 'results' / name}")
    installed = manifest.load()
    for tool, record in records.items():
        if record.get("status") not in {"ok", "reused"}:
            raise UpdateError(f"the previous {tool} stage was incomplete; finish that run first")
        if tool != "oma":
            commit = (installed["predictors"].get(tool) or {}).get("commit", "unknown")
            if commit != record.get("commit"):
                raise UpdateError(f"{tool} has changed version; restore it or start a new run")
    settings = previous.get("settings", {})
    if settings.get("alignment") != config.alignment:
        raise UpdateError("--update requires the previous --alignment setting")
    if _method_options(settings.get("tool_options", {})) != _method_options(config.tool_options):
        raise UpdateError("--update requires the previous --opt settings")
    for name, old in previous.get("shared_binaries", {}).items():
        now = installed["binaries"].get(name, {})
        if old.get("version") != now.get("version"):
            raise UpdateError(f"shared binary {name} has changed version; restore it or start a new run")

    reserved = {
        "orthofinder": {"-b", "-f", "-o", "-fg", "-ft", "-op", "-og", "-os", "-oa", "-ot"},
        "sonicparanoid": {"-i", "--input-directory", "-o", "--output-directory",
                          "-ow", "--overwrite", "-rs", "--remove-old-species"},
        "broccoli": {"-dir", "-steps"},
    }
    for tool, options in config.tool_options.items():
        conflicting = [o for o in options if o.split("=", 1)[0] in reserved.get(tool, set())]
        if conflicting:
            raise UpdateError(f"--update controls {tool} input and execution stages; remove {conflicting}")
    return previous, added


def _method_options(tool_options: dict[str, list[str]]) -> dict[str, list[str]]:
    options = {tool: list(values) for tool, values in tool_options.items()}
    # A SonicParanoid project names one execution, not its inference method.
    values = iter(options.get("sonicparanoid", []))
    cleaned = []
    for value in values:
        if value in {"-p", "--project-id"}:
            next(values, None)
        elif not value.startswith("--project-id="):
            cleaned.append(value)
    if cleaned:
        options["sonicparanoid"] = cleaned
    else:
        options.pop("sonicparanoid", None)
    return options


def native_root(previous: dict, tool: str, work_dir: Path) -> Path:
    record = next(r for r in previous["predictors"] if r["tool"] == tool)
    return Path(record.get("native_root") or work_dir / "results")


def recover_publication(outdir: Path) -> None:
    """Restore outputs if a process stopped while publishing an update."""
    journal = outdir / ".update-publication.json"
    if not journal.is_file():
        return
    state = json.loads(journal.read_text())
    staging = Path(state["staging"])
    if staging.resolve().parent != outdir.resolve():
        raise UpdateError(f"invalid update publication directory in {journal}")
    for name in state["files"]:
        relative = Path(name)
        if relative.is_absolute() or ".." in relative.parts:
            raise UpdateError(f"invalid output path in {journal}")
        destination = outdir / relative
        if name in state["existing"]:
            temporary = destination.with_name(destination.name + ".supporth-restore")
            shutil.copy2(staging / ".previous" / relative, temporary)
            temporary.replace(destination)
        else:
            destination.unlink(missing_ok=True)
    journal.unlink()
    shutil.rmtree(staging)


def publish_outputs(staging: Path, outdir: Path) -> None:
    """Publish matching pairs, tables, figures and provenance with rollback."""
    files = sorted(p.relative_to(staging) for p in staging.rglob("*") if p.is_file())
    existing = []
    for relative in files:
        destination = outdir / relative
        if destination.is_file():
            backup = staging / ".previous" / relative
            backup.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(destination, backup)
            existing.append(str(relative))
    journal = outdir / ".update-publication.json"
    temporary = journal.with_suffix(".tmp")
    temporary.write_text(json.dumps({"staging": str(staging), "files": [str(p) for p in files],
                                     "existing": existing}))
    temporary.replace(journal)
    try:
        for relative in files:
            destination = outdir / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            (staging / relative).replace(destination)
    except BaseException:
        recover_publication(outdir)
        raise
    journal.unlink()
    shutil.rmtree(staging)

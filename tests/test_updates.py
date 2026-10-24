"""Search reuse, expanded predictions and validation of update inputs."""

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from supporth import canonical, cli, pipeline, shell, similarity, updates
from supporth.adapters.base import Adapter, MethodProfile, StageContext
from supporth.adapters.orthofinder import OrthoFinder
from supporth.adapters.sonicparanoid import SonicParanoid
from supporth.proteomes import load_proteomes


@pytest.fixture
def fastas(tmp_path):
    files = []
    for label in ("a", "b", "c", "d"):
        path = tmp_path / f"{label}.faa"
        path.write_text(f">{label}1\nMAAAK\n>{label}2\nMCCCCK\n")
        files.append(path)
    return tuple(files)


class UpdatingPredictor(Adapter):
    name = "test_predictor"
    label = "X"

    def check(self, options=None):
        return True, "ok"

    def method_profile(self, context):
        return MethodProfile()

    def run(self, context):
        output = context.work_dir / f"run_{len(context.proteomes)}"
        output.mkdir(parents=True, exist_ok=True)
        pairs = [("a1", "b1")] if context.previous_root is None else [("a2", "c1")]
        canonical.write_json(canonical.PredictorResult(self.name, self.label, set(pairs)), output / "pairs.json")
        return output

    def parse(self, native_root, proteomes):
        return canonical.read_json(native_root / "pairs.json")


@pytest.fixture
def completed_run(tmp_path, fastas, monkeypatch):
    monkeypatch.setenv("SUPPORTH_HOME", str(tmp_path / "home"))
    adapter = UpdatingPredictor()
    monkeypatch.setattr(pipeline, "resolve", lambda _: [adapter])
    monkeypatch.setattr(pipeline.figures, "write_support_figure", lambda *args, **kwargs: None)
    config = pipeline.RunConfig(fastas=fastas[:2], outdir=tmp_path / "out", alignment="none", support_lv=1)
    pipeline.run_pipeline(config, log=lambda *_: None)
    return config, adapter


def test_update_replaces_predictions_and_can_be_extended_again(completed_run, fastas):
    config, _ = completed_run
    config.fastas = (fastas[2], fastas[1], fastas[0])
    config.update = True
    outcome = pipeline.run_pipeline(config, log=lambda *_: None)
    assert set(zip(outcome["table"].Protein1, outcome["table"].Protein2)) == {("a2", "c1")}
    assert list(outcome["provenance"]["inputs"]) == ["a", "b", "c"]
    assert outcome["provenance"]["update"]["added_proteomes"] == ["c"]
    assert Path(outcome["provenance"]["update"]["previous_provenance"]).is_file()
    config.fastas = fastas
    pipeline.run_pipeline(config, log=lambda *_: None)
    assert len(list((config.outdir / "history").glob("*.json"))) == 2


@pytest.mark.parametrize("change, message", [
    ("changed", "has changed"), ("missing", "is missing"),
    ("none_added", "additional proteome"), ("settings", "previous --opt"),
    ("version", "changed version"),
])
def test_invalid_updates_leave_outputs_untouched(completed_run, fastas, change, message, monkeypatch):
    config, adapter = completed_run
    config.update = True
    config.fastas = fastas[:3]
    before = (config.outdir / "suppOrth.tsv").read_bytes()
    if change == "changed":
        fastas[0].write_text(">a1\nMNEWK\n")
    elif change == "missing":
        config.fastas = fastas[1:3]
    elif change == "none_added":
        config.fastas = fastas[:2]
    elif change == "settings":
        config.tool_options[adapter.name] = ["--different"]
    elif change == "version":
        monkeypatch.setattr(pipeline.manifest, "load", lambda: {
            "predictors": {adapter.name: {"commit": "different"}}, "binaries": {},
        })
    monkeypatch.setattr(adapter, "run", lambda _: pytest.fail("invalid update launched a predictor"))
    with pytest.raises(pipeline.PipelineError, match=message):
        pipeline.run_pipeline(config, log=lambda *_: None)
    assert (config.outdir / "suppOrth.tsv").read_bytes() == before


def test_failed_update_preserves_consensus_and_saved_pairs(completed_run, fastas, monkeypatch):
    config, adapter = completed_run
    config.update = True
    config.fastas = fastas[:3]
    protected = [config.outdir / "suppOrth.tsv", config.outdir / "provenance.json",
                 config.outdir / "pairs" / f"{adapter.name}.json"]
    before = [p.read_bytes() for p in protected]

    def fail(context):
        raise RuntimeError("native predictor failed")

    monkeypatch.setattr(adapter, "run", fail)
    with pytest.raises(pipeline.PipelineError, match="consensus outputs were not replaced"):
        pipeline.run_pipeline(config, log=lambda *_: None)
    assert [p.read_bytes() for p in protected] == before


def test_resume_cannot_silently_ignore_added_proteomes(completed_run, fastas):
    config, _ = completed_run
    config.resume = True
    config.fastas = fastas[:3]
    with pytest.raises(pipeline.PipelineError, match="use --update"):
        pipeline.run_pipeline(config, log=lambda *_: None)


def test_resume_retains_native_root_for_a_later_update(completed_run, fastas):
    config, _ = completed_run
    config.resume = True
    resumed = pipeline.run_pipeline(config, log=lambda *_: None)
    assert resumed["reports"][0].native_root.endswith("run_2")
    config.resume = False
    config.update = True
    config.fastas = fastas[:3]
    pipeline.run_pipeline(config, log=lambda *_: None)


def test_cli_update_is_distinct_from_resume():
    parser = cli.build_parser()
    args = parser.parse_args(["run", "-f", "expanded", "-o", "results", "--update"])
    assert args.update and not args.resume
    with pytest.raises(SystemExit):
        parser.parse_args(["run", "-f", "expanded", "-o", "results", "--update", "--resume"])


def test_orthofinder_update_uses_only_new_fastas_and_keeps_previous_searches(tmp_path, fastas, monkeypatch):
    previous = tmp_path / "old" / "WorkingDirectory"
    previous.mkdir(parents=True)
    (previous / "SpeciesIDs.txt").write_text("0: a.fasta\n1: b.fasta\n")
    (previous.parent / "Log.txt").write_text("previous native run metadata\n")
    search = previous / "Blast0_1.txt"
    search.write_text("previous search\n")
    adapter = OrthoFinder()
    monkeypatch.setattr(adapter, "installed_root", lambda: tmp_path)
    monkeypatch.setattr(adapter, "locate_entrypoint", lambda _: tmp_path / "orthofinder.py")
    monkeypatch.setattr(adapter, "_interpreter", lambda _: "python")
    monkeypatch.setattr(adapter, "_configure_search", lambda _, program, threads: program)
    calls = []
    monkeypatch.setattr(shell, "run", lambda argv, **kw: calls.append(argv) or SimpleNamespace(returncode=0))
    context = StageContext(input_dir=tmp_path / "all", work_dir=tmp_path / "work",
                           proteomes=load_proteomes(fastas[:3]), previous_root=previous.parent,
                           added_input_dir=tmp_path / "new")
    output = adapter.run(context)
    argv = calls[0]
    assert argv[argv.index("-b") + 1] == str(previous.parent)
    assert argv[argv.index("-f") + 1] == str(context.added_input_dir)
    assert "-o" not in argv
    assert output.name == "Results_" + argv[argv.index("-n") + 1]
    assert output != previous.parent and search.read_text() == "previous search\n"


def test_sonicparanoid_update_parses_only_the_new_run(tmp_path, fastas, monkeypatch):
    adapter = SonicParanoid()
    monkeypatch.setattr(adapter, "installed_root", lambda: tmp_path)
    monkeypatch.setattr(adapter, "locate_entrypoint", lambda _: tmp_path / "sonicparanoid")
    monkeypatch.setattr(adapter, "prefix_interpreter", lambda _: None)
    monkeypatch.setattr(adapter, "stage_env", lambda *_: {})
    output = tmp_path / "work" / "results"
    old = output / "runs" / "old"
    old.mkdir(parents=True)
    (output / "snapshot.tsv").write_text("old snapshot\n")
    (old / "pairwise_orthologs.old.tsv").write_text("a1\tb1\n")
    calls = []

    def run(argv, **kwargs):
        calls.append(argv)
        latest = output / "runs" / "new"
        latest.mkdir()
        (latest / "pairwise_orthologs.new.tsv").write_text("a2\tc1\n")
        (output / "last_run_info.txt").write_text(f"Run directory:\t{latest}\n")

    monkeypatch.setattr(shell, "run", run)
    proteomes = load_proteomes(fastas[:3])
    context = StageContext(input_dir=tmp_path / "all", work_dir=tmp_path / "work",
                           proteomes=proteomes, previous_root=old)
    latest = adapter.run(context)
    assert adapter.parse(latest, proteomes).pairs == {("a2", "c1")}
    assert "-op" in calls[0] and "-ow" not in calls[0]
    assert (old / "pairwise_orthologs.old.tsv").read_text() == "a1\tb1\n"


def test_alignment_statistics_follow_previous_orthofinder_directories(tmp_path, fastas):
    old, new = tmp_path / "old", tmp_path / "new"
    old.mkdir()
    new.mkdir()
    (new / "SpeciesIDs.txt").write_text("0: a.fasta\n1: b.fasta\n2: c.fasta\n")
    (new / "SequenceIDs.txt").write_text("0_0: a1\n1_0: b1\n2_0: c1\n")
    (new / "previous_wd.txt").write_text(str(old) + "\n")
    line = "{}\t{}\t{}\t10\t0\t0\t1\t10\t1\t10\t1e-20\t100\n"
    (old / "Blast0_1.txt").write_text(line.format("0_0", "1_0", 90))
    (new / "Blast0_2.txt").write_text(line.format("0_0", "2_0", 85))
    index = similarity.index_from_orthofinder(new, load_proteomes(fastas[:3]), log=lambda *_: None)
    assert index.get("a1", "b1").identity == 90
    assert index.get("a1", "c1").identity == 85
    (old / "previous_wd.txt").write_text(str(new) + "\n")
    with pytest.raises(ValueError, match="cycle"):
        similarity.orthofinder_search_dirs(new)


def test_fallback_alignment_update_does_not_repeat_old_comparisons(tmp_path, fastas, monkeypatch):
    calls = []
    monkeypatch.setattr(shell, "which", lambda name: tmp_path / name)

    def run(argv, **kwargs):
        if argv[1] == "blastp":
            calls.append((Path(argv[argv.index("-q") + 1]).stem,
                          Path(argv[argv.index("-d") + 1]).name))
            Path(argv[argv.index("-o") + 1]).write_text("a1\tb1\t90\t1e-20\t100\n")

    monkeypatch.setattr(shell, "run", run)
    work = tmp_path / "searches"
    similarity.build_index(load_proteomes(fastas[:2]), work, log=lambda *_: None)
    originals = {p.name: p.read_bytes() for p in work.glob("*.tsv")}
    calls.clear()
    similarity.build_index(load_proteomes(fastas[:3]), work, reuse_labels=frozenset({"a", "b"}),
                           log=lambda *_: None)
    assert set(calls) == {("a", "c"), ("b", "c"), ("c", "a"), ("c", "b")}
    assert all((work / name).read_bytes() == content for name, content in originals.items())


def test_failed_output_generation_does_not_publish_any_updated_pairs(completed_run, fastas, monkeypatch):
    config, _ = completed_run
    config.fastas = fastas[:3]
    config.update = True
    paths = [config.outdir / "pairs" / "test_predictor.json", config.outdir / "provenance.json"]
    before = [p.read_bytes() for p in paths]

    def fail(*args, **kwargs):
        raise OSError("disk full")

    monkeypatch.setattr(canonical, "write_legacy_json", fail)
    with pytest.raises(OSError, match="disk full"):
        pipeline.run_pipeline(config, log=lambda *_: None)
    assert [p.read_bytes() for p in paths] == before


def test_publication_failure_restores_matching_previous_outputs(tmp_path, monkeypatch):
    staging = tmp_path / ".update-test"
    staging.mkdir()
    for name in ("pairs.json", "provenance.json"):
        (tmp_path / name).write_text("original")
        (staging / name).write_text("updated")
    original = Path.replace

    def fail_second_file(self, target):
        if self == staging / "provenance.json":
            raise OSError("write failed")
        return original(self, target)

    monkeypatch.setattr(Path, "replace", fail_second_file)
    with pytest.raises(OSError, match="write failed"):
        updates.publish_outputs(staging, tmp_path)
    assert (tmp_path / "pairs.json").read_text() == "original"
    assert (tmp_path / "provenance.json").read_text() == "original"
    assert not (tmp_path / ".update-publication.json").exists()


def test_recovery_after_process_stops_mid_publication(tmp_path):
    staging = tmp_path / ".update-test"
    backups = staging / ".previous"
    backups.mkdir(parents=True)
    (backups / "pairs.json").write_text("original pairs")
    (tmp_path / "pairs.json").write_text("expanded pairs")
    (tmp_path / "new-file.json").write_text("partial publication")
    (tmp_path / ".update-publication.json").write_text(json.dumps({
        "staging": str(staging), "files": ["pairs.json", "new-file.json"], "existing": ["pairs.json"],
    }))
    updates.recover_publication(tmp_path)
    assert (tmp_path / "pairs.json").read_text() == "original pairs"
    assert not (tmp_path / "new-file.json").exists()
    assert not staging.exists()


def test_project_names_can_change_without_invalidating_sonicparanoid_searches():
    assert updates._method_options({"sonicparanoid": ["-p", "first", "-m", "sensitive"]}) == \
        updates._method_options({"sonicparanoid": ["--project-id=second", "-m", "sensitive"]})

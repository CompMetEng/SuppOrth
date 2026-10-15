import gzip
from types import SimpleNamespace

import pytest

from supporth import cli, deps, pipeline, shell, similarity
from supporth.adapters.oma import OMA
from supporth.adapters.tabular import pairs_from_table
from supporth.proteomes import Proteome, ProteomeError, load_proteomes, stage_input_dir


def test_gzip_fastas_are_validated_and_staged_uncompressed(tmp_path):
    files = []
    for label in ('one', 'two'):
        path = tmp_path / f'{label}.FAA.gz'
        with gzip.open(path, 'wt') as handle:
            handle.write(f'>{label}_protein\nMAAAK\n')
        files.append(path)
    proteomes = load_proteomes(files)
    staged = stage_input_dir(tmp_path / 'input', proteomes)
    assert [p.label for p in proteomes] == ['one', 'two']
    assert (staged / 'one.fasta').read_text() == '>one_protein\nMAAAK\n'


def test_duplicate_fasta_identifiers_are_rejected(tmp_path):
    first, second = tmp_path / 'a.fa', tmp_path / 'b.fa'
    first.write_text('>a\nMAAK\n>a\nMAAK\n')
    second.write_text('>b\nMAAK\n')
    with pytest.raises(ProteomeError, match='duplicate identifier'):
        load_proteomes((first, second))
    first.write_text('>a\nMAAK\n')
    with pytest.raises(ProteomeError, match='path separators'):
        load_proteomes((first, second), labels=('../a', 'b'))


def test_large_rows_are_not_silently_dropped(tmp_path):
    table = tmp_path / 'groups.tsv'
    left = frozenset(f'a{i}' for i in range(317))
    right = frozenset(f'b{i}' for i in range(317))
    table.write_text('\t'.join(sorted(left | right)) + '\n')
    assert sum(1 for _ in pairs_from_table(table, left, right)) == 317 * 317
    with pytest.raises(ValueError, match='exceeding the requested limit'):
        list(pairs_from_table(table, left, right, max_row_pairs=10))


def test_oma_direct_folder_and_mixed_filenames_keep_all_pairs(tmp_path):
    folder = tmp_path / 'PairwiseOrthologs'
    folder.mkdir()
    (folder / 'a-b.txt').write_text('1\t1\ta1\tb1\t1:1\n')
    (folder / 'different_names.txt').write_text('1\t1\ta1\tc1\t1:1\n')
    proteomes = tuple(Proteome(label, tmp_path / f'{label}.fa', frozenset({f'{label}1'}))
                      for label in ('a', 'b', 'c'))
    assert OMA().parse(folder, proteomes).pairs == {('a1', 'b1'), ('a1', 'c1')}


def test_install_failure_returns_nonzero(tmp_path, monkeypatch):
    class Broken:
        name = 'broken'

        def install(self, **kwargs):
            raise RuntimeError('build failed')

    monkeypatch.setenv('SUPPORTH_HOME', str(tmp_path / 'home'))
    monkeypatch.setattr(cli.adapters, 'resolve', lambda _: [Broken()])
    monkeypatch.setattr(cli.deps, 'ensure', lambda *args, **kwargs: {})
    assert cli.main(['install']) == 1


def test_quoted_options_and_failed_stage_exit_status(tmp_path, monkeypatch):
    captured = []

    def run(config):
        captured.append(config)
        return {'reports': [SimpleNamespace(status='failed')]}

    monkeypatch.setattr(cli, 'run_pipeline', run)
    args = ['run', 'a.fa', 'b.fa', '-o', str(tmp_path / 'out'),
            '--opt', 'orthofinder=-X "value with spaces"']
    assert cli.main(args) == 1
    assert captured[0].tool_options['orthofinder'] == ['-S', 'diamond', '-X', 'value with spaces']


@pytest.mark.parametrize('value', ['nan', 'inf', '-1', '101'])
def test_identity_threshold_rejects_nonfinite_and_out_of_range_values(value):
    with pytest.raises(SystemExit):
        cli.build_parser().parse_args(['run', 'a.fa', 'b.fa', '-o', 'out', f'--identity_filter={value}'])


def test_shared_binary_directory_always_precedes_system_paths(tmp_path, monkeypatch):
    monkeypatch.setenv('SUPPORTH_HOME', str(tmp_path))
    monkeypatch.setenv('PATH', f'/system/bin:{tmp_path / "bin"}:/other/bin')
    assert shell.tool_env()['PATH'].split(':') == [str(tmp_path / 'bin'), '/system/bin', '/other/bin']


def test_fallback_search_builds_one_database_per_proteome(tmp_path, monkeypatch):
    proteomes = tuple(Proteome(label, tmp_path / f'{label}.fa', frozenset({f'{label}1'}))
                      for label in ('a', 'b', 'c'))
    calls = []
    monkeypatch.setattr(similarity.shell, 'which', lambda name: tmp_path / name)

    def run(argv, **kwargs):
        calls.append(argv)
        if argv[1] == 'blastp':
            output = argv[argv.index('-o') + 1]
            from pathlib import Path
            Path(output).write_text('a1\tb1\t90\t1e-20\t100\n')

    monkeypatch.setattr(similarity.shell, 'run', run)
    similarity.build_index(proteomes, tmp_path / 'alignment', log=lambda *_: None)
    assert sum(argv[1] == 'makedb' for argv in calls) == 3
    assert sum(argv[1] == 'blastp' for argv in calls) == 6


def test_interrupted_download_does_not_become_a_cached_archive(tmp_path, monkeypatch):
    class Interrupted:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def read(self, size):
            raise OSError('connection interrupted')

    monkeypatch.setattr(deps.urllib.request, 'urlopen', lambda *args, **kwargs: Interrupted())
    destination = tmp_path / 'archive.gz'
    with pytest.raises(OSError):
        deps._download('https://example.invalid/archive.gz', destination, lambda *_: None)
    assert not destination.exists()
    assert not destination.with_name('archive.gz.part').exists()


def test_failed_orthofinder_run_does_not_count_partial_groups_as_success(tmp_path, monkeypatch):
    from supporth.adapters.base import StageContext
    from supporth.adapters.orthofinder import OrthoFinder

    root = tmp_path / 'tool'
    root.mkdir()
    entrypoint = root / 'orthofinder.py'
    entrypoint.touch()
    adapter = OrthoFinder()
    monkeypatch.setattr(adapter, 'installed_root', lambda: root)
    monkeypatch.setattr(adapter, '_interpreter', lambda _: 'python')
    proteomes = tuple(Proteome(label, tmp_path / f'{label}.fa', frozenset({f'{label}1'}))
                      for label in ('a', 'b'))
    context = StageContext(input_dir=tmp_path / 'input', work_dir=tmp_path / 'work', proteomes=proteomes)

    def fail(argv, **kwargs):
        groups = context.work_dir / 'results' / 'Orthogroups'
        groups.mkdir(parents=True)
        (groups / 'Orthogroups.tsv').write_text('Orthogroup\ta\tb\nOG1\ta1\tb1\n')
        return SimpleNamespace(returncode=1)

    monkeypatch.setattr(shell, 'run', fail)
    with pytest.raises(shell.CommandFailed):
        adapter.run(context)


def test_adopting_the_installed_pfam_directory_preserves_it(tmp_path):
    from supporth.adapters.sonicparanoid import SonicParanoid

    adapter = SonicParanoid()
    database = adapter.profile_db(tmp_path)
    database.mkdir(parents=True)
    for name in adapter.pfam_files():
        (database / name).write_bytes(b'profile data')
    adapter.adopt_pfam_profiles(tmp_path, database, log=lambda *_: None)
    assert adapter.pfam_ready(tmp_path)
    assert all((database / name).read_bytes() == b'profile data' for name in adapter.pfam_files())


def test_cpu_default_respects_scheduler_affinity(monkeypatch):
    import os

    monkeypatch.setattr(os, 'cpu_count', lambda: 64)
    monkeypatch.setattr(os, 'sched_getaffinity', lambda _: {0, 1, 2, 3}, raising=False)
    assert shell.cpu_count() == 3
    assert shell.cpu_count(2) == 2


def test_cpu_default_without_affinity(monkeypatch):
    import os

    monkeypatch.delattr(os, 'sched_getaffinity', raising=False)
    monkeypatch.setattr(os, 'cpu_count', lambda: 4)
    assert shell.cpu_count() == 3


def test_sonicparanoid_install_passes_pandas_constraint(tmp_path, monkeypatch):
    from pathlib import Path
    from supporth.adapters.sonicparanoid import SonicParanoid

    calls = []
    monkeypatch.setattr(shell, 'run', lambda command, **kwargs: calls.append(command))
    SonicParanoid().install_to_prefix(tmp_path, Path('/python'), log=lambda *_: None)
    command = calls[0]
    constraints = Path(command[command.index('--constraint') + 1])
    assert constraints.read_text() == 'pandas<3\n'
    assert command[-1] == '.'


def test_failed_pfam_index_never_gets_ready_marker(tmp_path, monkeypatch):
    from pathlib import Path
    from supporth.adapters.base import AdapterError
    from supporth.adapters.sonicparanoid import SonicParanoid

    adapter = SonicParanoid()
    monkeypatch.setattr(shell, 'run', lambda *args, **kwargs: None)
    with pytest.raises(AdapterError, match='pfam-index.log'):
        adapter.prepare_pfam_profiles(tmp_path, Path('/python'), threads=1, log=lambda *_: None)
    assert not (adapter.profile_db(tmp_path) / adapter.PFAM_MARKER).exists()


def test_support_threshold_warning_mentions_identity_exception():
    messages = []
    pipeline.warn_support_threshold(3, 2, messages.append)
    assert 'identity exception' in messages[0]
    messages.clear()
    pipeline.warn_support_threshold(3, 3, messages.append)
    assert messages == []


def test_conda_detection_uses_interpreter_prefix(monkeypatch):
    from supporth import interpreters

    monkeypatch.setattr(shell, 'run', lambda *args, **kwargs:
                        SimpleNamespace(returncode=0, stdout='/opt/micromamba/envs/orthology\n'))
    assert interpreters.is_conda_style('/unrelated/python')
    monkeypatch.setattr(shell, 'run', lambda *args, **kwargs:
                        SimpleNamespace(returncode=0, stdout='/tmp/venv\n'))
    assert not interpreters.is_conda_style('/opt/anaconda3/bin/python')


@pytest.mark.parametrize('version, warns', [
    ('13.45111', False),
    ('45111b641859ed0ddd875b94d6fd1aef1a675b7e', False),
    ('15.6f452', True),
])
def test_mmseqs_adoption_version_warning(tmp_path, monkeypatch, version, warns):
    messages = []
    monkeypatch.setattr(shell, 'run', lambda *args, **kwargs:
                        SimpleNamespace(returncode=0, stdout=version, stderr=''))
    deps._warn_version(deps.SPECS['mmseqs'], tmp_path / 'mmseqs', messages.append)
    assert bool(messages) == warns


def test_oma_numeric_entry_ids_and_unknown_headers_are_not_predictions(tmp_path):
    table = tmp_path / 'pairs.txt'
    table.write_text('1\t2\ta1\tb1 mentions 2\t1:1\n3\t4\tunknown\tb1\t1:1\n')
    proteomes = (Proteome('a', tmp_path / 'a.fa', frozenset({'a1', '1'})),
                 Proteome('b', tmp_path / 'b.fa', frozenset({'b1', '2'})))
    result = OMA().parse(table, proteomes)
    assert result.pairs == {('a1', 'b1')}
    assert result.dropped_unknown == 1


def test_oma_group_fallback_retains_its_separate_format(tmp_path):
    table = tmp_path / 'OrthologousGroups.txt'
    table.write_text('OMA1\ta1\tb1\n')
    proteomes = (Proteome('a', tmp_path / 'a.fa', frozenset({'a1'})),
                 Proteome('b', tmp_path / 'b.fa', frozenset({'b1'})))
    assert OMA().parse(table, proteomes).pairs == {('a1', 'b1')}


@pytest.mark.parametrize('options', [
    [],
    ['-path_fasttree', '/custom/FastTree', '-path_diamond', '/custom/diamond'],
    ['-path_fasttree=/custom/FastTree', '-path_diamond=/custom/diamond'],
])
def test_broccoli_uses_shared_binary_paths_and_preserves_overrides(tmp_path, monkeypatch, options):
    from supporth.adapters.base import StageContext
    from supporth.adapters.broccoli import Broccoli

    monkeypatch.setenv('SUPPORTH_HOME', str(tmp_path / 'home'))
    shared = tmp_path / 'home' / 'bin'
    shared.mkdir(parents=True)
    for name in ('FastTree', 'diamond'):
        binary = shared / name
        binary.write_text('#!/bin/sh\nexit 0\n')
        binary.chmod(0o755)
    assert 'fasttree' not in {path.name for path in shared.iterdir()}
    adapter = Broccoli()
    monkeypatch.setattr(adapter, 'installed_root', lambda: tmp_path)
    monkeypatch.setattr(adapter, 'locate_entrypoint', lambda _: tmp_path / 'broccoli.py')
    monkeypatch.setattr(adapter, '_interpreter', lambda _: 'python')
    calls = []
    monkeypatch.setattr(shell, 'run', lambda command, **kwargs: calls.append(command))
    proteomes = tuple(Proteome(label, tmp_path / f'{label}.fa', frozenset({label}))
                      for label in ('a', 'b'))
    context = StageContext(input_dir=tmp_path / 'input', work_dir=tmp_path / 'work',
                           proteomes=proteomes, options=options)
    adapter.run(context)
    command = calls[0]
    if options:
        assert command[-len(options):] == options
        assert not any(str(shared) in value for value in command)
    else:
        assert command[command.index('-path_fasttree') + 1] == str(shared / 'FastTree')
        assert command[command.index('-path_diamond') + 1] == str(shared / 'diamond')

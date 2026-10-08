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
    (folder / 'a-b.txt').write_text('a1\tb1\n')
    (folder / 'different_names.txt').write_text('a1\tc1\n')
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

from pathlib import Path
import shutil
import sys
from urllib.parse import unquote

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import build_dashboard as dashboard
from dashboard_fixture import make_demo, put
from graph_evidence import source_graph_files


def test_unmapped_source_remains_downloadable_without_claiming_selected_path(tmp_path):
    root = make_demo(tmp_path/'results', samples=1)
    graph = root/'W30877/MitoLooper/getorganelle_assembly_graph.gfa'
    graph.write_text('H\tVN:Z:1.0\nS\tother\t'+'A'*100+'\n')
    summary = root/'W30877/final/summary.json'
    before = summary.read_bytes()
    mito = dashboard.build_data(root)['samples'][0]['mito']
    assert not mito['graphs']
    assert len(mito['source_graphs']) == 1
    assert (root/unquote(mito['source_graphs'][0]['href'])) == graph
    assert summary.read_bytes() == before
    assert mito['length'] == 15740  # Sequence recovery is separate from graph mapping.


def test_missing_source_is_distinct_from_retained_unmapped_source(tmp_path):
    root = make_demo(tmp_path/'results', samples=1)
    (root/'W30877/MitoLooper/getorganelle_assembly_graph.gfa').unlink()
    mito = dashboard.build_data(root)['samples'][0]['mito']
    assert not mito['graphs'] and not mito['source_graphs']


def test_source_downloads_scope_known_methods_but_keep_unattributed_graphs(tmp_path):
    root = make_demo(tmp_path/'results', samples=1)
    sample = root/'W30877'
    put(sample/'MitoLooper/mitofinder_graph.fastg', '>EDGE_1_length_4_cov_1;\nACGT\n')
    put(sample/'MitoLooper/unattributed.gfa2', 'H\tVN:Z:2.0\n')
    files = source_graph_files(sample, root, 'MitoLooper', 'getorganelle')
    assert {item['name'] for item in files} == {'getorganelle_assembly_graph.gfa','unattributed.gfa2'}
    assert all('source' not in item for item in files)


def test_manifest_restores_method_for_renamed_graph_and_links_survive_move(tmp_path):
    root = make_demo(tmp_path/'results', samples=1)
    sample = root/'W30877'
    graph = sample/'MitoLooper/getorganelle_assembly_graph.gfa'
    renamed = graph.with_name('opaque_source.gfa')
    graph.rename(renamed)
    put(sample/'final/outputs_manifest.tsv', 'source\tretained_file\n'
        '/old/location/W30877/012_getorganelle/output/assembly_graph.gfa\t'
        '/old/location/W30877/MitoLooper/opaque_source.gfa\n')
    files = source_graph_files(sample, root, 'MitoLooper', 'getorganelle')
    assert [item['name'] for item in files] == ['opaque_source.gfa']
    moved = tmp_path/'moved'
    shutil.move(root, moved)
    files = source_graph_files(moved/'W30877', moved, 'MitoLooper', 'getorganelle')
    assert (moved/unquote(files[0]['href'])).is_file()


def test_dashboard_source_links_do_not_expose_external_graph_symlinks(tmp_path):
    root = make_demo(tmp_path/'results', samples=1)
    external = put(tmp_path/'external.gfa', 'S\tx\tACGT\n')
    (root/'W30877/MitoLooper/external.gfa').symlink_to(external)
    files = source_graph_files(root/'W30877', root, 'MitoLooper', 'getorganelle')
    assert not any(item['name'] == 'external.gfa' for item in files)

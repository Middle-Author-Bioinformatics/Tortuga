import fcntl
import json
from pathlib import Path
import re
import shutil
import sys
from urllib.parse import unquote

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import build_dashboard as d
from dashboard_fixture import make_demo, put
from compact_outputs import atomic_json


def payload(path):
    text = path.read_text()
    return json.loads(re.search(r'<script type="application/json" id="report-data">(.*?)</script>', text, re.S).group(1))


def test_report_is_embedded_offline_and_links_survive_directory_move(tmp_path):
    root = make_demo(tmp_path/'original', samples=2)
    before = {path.relative_to(root): path.read_bytes() for path in root.rglob('*') if path.is_file()}
    path = d.write_dashboard(root, 'Illustrative report')
    text = path.read_text()
    assert not re.search(r'<(?:script|link)[^>]+(?:src|href)=', text)
    assert 'fetch(' not in text and 'XMLHttpRequest' not in text
    assert str(root) not in text
    data = payload(path)
    assert data['title'] == 'Illustrative report'
    assert len(data['samples']) == 2 and data['samples'][0]['its']['count'] == 2
    assert data['samples'][0]['mito']['accepted'] is True
    assert data['samples'][0]['mito']['length'] == 15740
    assert data['samples'][0]['mito']['hits'][0]['taxon'] == 'Rathbunixa occidentalis'
    assert data['samples'][0]['its']['loci'][0]['its_length_bp'] == 990
    assert data['samples'][0]['its']['loci'][0]['regions'][0]['name'] == 'SSU'
    for relative, content in before.items():
        assert (root/relative).read_bytes() == content
    moved = tmp_path/'moved'
    shutil.move(root, moved)
    data = payload(moved/'dashboard.html')
    for sample in data['samples']:
        for item in [sample['summary_file'], *sample['mito']['files'], *sample['its']['files']]:
            if item:
                assert (moved/unquote(item['href'])).is_file()


def test_gfa_graph_keeps_real_oriented_links_and_self_loops(tmp_path):
    path = put(tmp_path/'graph.gfa', 'H\tVN:Z:1.0\nS\tA\tACGT\tDP:f:22\nS\tB\t*\tLN:i:120\nL\tA\t+\tB\t-\t3M\nL\tB\t+\tB\t+\t0M\n')
    graph = d.parse_graph(path, tmp_path)
    assert graph['total_nodes'] == 2
    assert graph['nodes'][0]['length'] == 4 and graph['nodes'][0]['coverage'] == 22
    assert graph['nodes'][1]['length'] == 120
    assert graph['edges'][0]['target_strand'] == '-'
    assert graph['edges'][1]['source'] == graph['edges'][1]['target'] == 'B'
    assert all('x' in node and 'y' in node for node in graph['nodes'])


def test_fastg_reverse_complements_are_not_duplicate_segments(tmp_path):
    path = put(tmp_path/'graph.fastg', ">EDGE_1_length_60_cov_9:EDGE_2_length_70_cov_10';\nACGT\n>EDGE_1_length_60_cov_9';\nACGT\n>EDGE_2_length_70_cov_10;\nAAAA\n")
    graph = d.parse_graph(path, tmp_path)
    assert len(graph['nodes']) == 2 and graph['total_nodes'] == 2
    assert graph['edges'][0]['target_strand'] == '-'
    assert graph['nodes'][0]['coverage'] == 9


def test_graph_cap_and_oversized_graph_are_explicit(tmp_path, monkeypatch):
    path = put(tmp_path/'graph.gfa', 'S\tA\tAAAA\nS\tB\tAAAA\nL\tA\t+\tB\t+\t0M\n')
    graph = d.parse_graph(path, tmp_path, node_limit=1)
    assert graph['truncated'] and len(graph['nodes']) == 1 and not graph['edges']
    monkeypatch.setattr(d, 'GRAPH_BYTES', 1)
    graph = d.parse_graph(path, tmp_path)
    assert not graph['nodes'] and graph['file']['href'] == 'graph.gfa'
    assert '100 MiB' in graph['warnings'][0]


def test_gff_features_use_coordinates_strand_and_warn_on_multiple_sequences(tmp_path):
    path = put(tmp_path/'result.gff', '##gff-version 3\nmt\tMITOS2\tCDS\t1\t20\t.\t-\t0\tName=cox1;product=cytochrome%20oxidase\nmt\tMITOS2\ttRNA\t25\t40\t.\t+\t.\tName=trnA\nmt\tMITOS2\tCDS\t90\t110\t.\t+\t.\tName=outside\nother\tMITOS2\tCDS\t1\t10\t.\t+\t.\tName=other\n##FASTA\n>mt\nACGT\n')
    annotation = d.annotations(path, 100, tmp_path)
    assert len(annotation['features']) == 3
    assert annotation['features'][0]['strand'] == '-'
    assert annotation['features'][0]['product'] == 'cytochrome oxidase'
    assert annotation['multiple_sequences']
    assert any('exceed' in warning for warning in annotation['warnings'])


def test_missing_graph_or_annotations_are_not_fabricated(tmp_path):
    root = tmp_path/'results'
    sample = root/'SAMPLE'
    atomic_json(sample/'final/summary.json', {'sample': 'SAMPLE', 'pipeline': 'both', 'status': 'no_candidate',
                                           'itsme_ok': False, 'target_taxon': 'Expected species'})
    data = payload(d.write_dashboard(root))['samples'][0]
    assert data['mito']['sequence'] is None and not data['mito']['accepted']
    assert not data['mito']['graphs'] and not data['mito']['annotation']['features']
    assert not data['its']['graphs'] and data['its']['count'] is None
    assert data['expected_taxon'] == 'Expected species' and not data['mito']['hits']


def test_script_like_names_are_json_escaped_and_not_html(tmp_path):
    root = make_demo(tmp_path/'results', samples=1)
    title = '</script><script>alert("unsafe")</script>'
    path = d.write_dashboard(root, title)
    assert title not in path.read_text()
    assert payload(path)['title'] == title
    assert path.read_text().count('<script') == 2


def test_external_symlinks_are_not_downloadable_and_dashboard_links_refused(tmp_path):
    root = make_demo(tmp_path/'results', samples=1)
    outside = put(tmp_path/'outside.fasta', '>private\nACGT\n')
    (root/'W30877/MitoLooper/external.fasta').symlink_to(outside)
    assert not any(item['name'] == 'external.fasta' for item in d.build_data(root)['samples'][0]['mito']['files'])
    (root/'dashboard.html').symlink_to(outside)
    with pytest.raises(ValueError, match='linked'):
        d.write_dashboard(root)
    assert outside.read_text() == '>private\nACGT\n'


def test_active_controller_blocks_dashboard_without_touching_it(tmp_path):
    root = make_demo(tmp_path/'results', samples=1)
    path = d.write_dashboard(root)
    before = path.read_bytes()
    with (root/'.pipeline.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with pytest.raises(SystemExit):
            d.main([str(root)])
    assert path.read_bytes() == before


def test_reporting_refresh_also_builds_dashboard_without_analysis(tmp_path):
    from itsme_reporting import main
    root = make_demo(tmp_path/'results', samples=1)
    assert main([str(root)]) == 0
    assert (root/'cohort_its_loci.csv').is_file()
    assert payload(root/'dashboard.html')['samples'][0]['its']['count'] == 2


def test_empty_fasta_is_not_claimed_as_accepted_assembly(tmp_path):
    root = make_demo(tmp_path/'results', samples=1)
    (root/'W30877/MitoLooper/accepted_candidate.fasta').write_text('')
    sample = d.build_data(root)['samples'][0]
    assert sample['mito']['accepted'] is False
    assert sample['mito']['sequence']['name'] == 'best_sequence.fasta'


def test_gfa2_is_downloadable_but_not_misread_as_gfa1(tmp_path):
    path = put(tmp_path/'graph.gfa', 'H\tVN:Z:2.0\nS\tA\t1000\t*\n')
    result = d.parse_graph(path, tmp_path)
    assert not result['nodes']
    assert any('GFA2' in message for message in result['warnings'])
    assert result['file']['href'] == 'graph.gfa'


def test_stale_summary_does_not_show_old_recovery_as_current(tmp_path):
    import os
    root = make_demo(tmp_path/'results', samples=1)
    summary = root/'W30877/final/summary.json'
    newer = summary.stat().st_mtime_ns+1000000
    os.utime(root/'W30877/run_config.json', ns=(newer, newer))
    sample = d.build_data(root)['samples'][0]
    assert sample['status'] == 'unfinished'
    assert not sample['mito']['accepted'] and sample['mito']['length'] is None
    assert sample['its']['count'] is None and not sample['its']['loci']
    assert sample['warnings'] and sample['mito']['files']


def test_explicit_nonacceptance_and_empty_its_do_not_get_primary_accepted_downloads(tmp_path):
    root = make_demo(tmp_path/'results', samples=1)
    sample = root/'W30877'
    path = sample/'final/summary.json'
    summary = json.loads(path.read_text())
    summary['accepted_by_policy'] = False
    atomic_json(path, summary)
    (sample/'ITSME/reported_loci.fasta').write_text('')
    data = d.build_data(root)['samples'][0]
    assert not data['mito']['accepted']
    assert data['mito']['sequence']['name'] == 'best_sequence.fasta'
    assert data['its']['sequence'] is None


def test_stage_log_links_follow_recorded_paths_and_survive_relocation(tmp_path):
    from compact_outputs import table
    root = make_demo(tmp_path/'original', samples=1)
    sample = root/'W30877'
    put(sample/'logs/native error.log', 'native failure\n')
    stages = [
        {'stage': 'mitos2', 'status': 'failed', 'log': str(sample/'logs/native error.log')},
        {'stage': 'wgs', 'status': 'timeout', 'log': 'logs/wgs.log'},
        {'stage': 'missing', 'status': 'failed'},
    ]
    table(sample/'stages.tsv', stages, ['stage', 'status', 'log'])
    moved = tmp_path/'moved'
    shutil.move(root, moved)
    data = payload(d.write_dashboard(moved))['samples'][0]
    assert data['stages'][0]['log']['href'] == 'W30877/logs/native%20error.log'
    assert data['stages'][1]['log']['href'] == 'W30877/logs/wgs.log'
    assert data['stages'][2]['log'] is None
    for stage in data['stages'][:2]:
        assert (moved/unquote(stage['log']['href'])).is_file()


def test_stage_logs_do_not_link_outside_sample(tmp_path):
    root = make_demo(tmp_path/'results', samples=2)
    sample = root/'W30877'
    outside = put(tmp_path/'private.log', 'outside\n')
    (sample/'logs/linked.log').symlink_to(outside)
    for recorded in [str(outside), '../W30882/logs/wgs.log', 'logs/linked.log']:
        assert d.stage_log({'stage': 'missing', 'log': recorded}, sample, root) is None


def test_matching_reference_name_does_not_override_unresolved_target_check(tmp_path):
    root = make_demo(tmp_path/'results', samples=1)
    sample = root/'W30877'
    summary = json.loads((sample/'final/summary.json').read_text())
    summary.update(target_status='unresolved', accepted_by_policy=False)
    atomic_json(sample/'final/summary.json', summary)
    evidence_path = sample/'candidates/011_novoplasty_001/bold_evidence.json'
    evidence = json.loads(evidence_path.read_text())
    evidence['reason'] = 'Near-best BOLD records conflict at requested taxon'
    lower = dict(evidence['hits'][0], bitscore='900', metadata={'taxon': 'Other species'})
    tied = dict(evidence['hits'][0], metadata={'taxon': 'Tied species'})
    evidence['hits'] = [lower, evidence['hits'][0], tied]
    atomic_json(evidence_path, evidence)
    data = d.build_data(root)['samples'][0]['mito']
    assert data['reference_taxa'] == ['Rathbunixa occidentalis', 'Tied species']
    assert data['hits'][0]['taxon'] == 'Rathbunixa occidentalis'
    assert data['target_status'] == 'unresolved' and not data['accepted']
    assert data['target_reason'] == evidence['reason']
    assert json.loads((sample/'final/summary.json').read_text()) == summary
    # Names in older records stay visible without claiming a score ranking.
    for hit in evidence['hits']:
        hit.pop('bitscore', None)
    atomic_json(evidence_path, evidence)
    older = d.build_data(root)['samples'][0]['mito']
    assert older['reference_taxa'] == ['Other species', 'Rathbunixa occidentalis', 'Tied species']
    assert older['reference_ranked'] is False


@pytest.mark.parametrize('override', [False, True])
def test_raw_bold_matches_restore_taxa_when_confirmation_retained_no_hits(tmp_path, override):
    from compact_outputs import table
    root = make_demo(tmp_path/'results', samples=2)
    metadata = tmp_path/'metadata.tsv'
    table(metadata, [
        {'sequence_id': 'COI_A', 'processid': 'A', 'taxon': 'Observed species', 'lineage': 'Animalia;Observed species'},
        {'sequence_id': 'COI_B', 'processid': 'B', 'taxon': 'Another species', 'lineage': 'Animalia;Another species'},
    ], ['sequence_id', 'processid', 'taxon', 'lineage'])
    for name in ['W30877', 'W30882']:
        sample = root/name
        summary = json.loads((sample/'final/summary.json').read_text())
        summary.update(target_status='unresolved', accepted_by_policy=False, target_taxon='Expected species')
        atomic_json(sample/'final/summary.json', summary)
        config = {'options': {'bold_metadata': str(tmp_path/'unavailable.tsv' if override else metadata)}}
        # Preserve the checkpoint age: changing its timestamp means a new run attempt.
        config_path = sample/'run_config.json'
        import os
        old_time = config_path.stat().st_mtime_ns
        atomic_json(config_path, config)
        os.utime(config_path, ns=(old_time, old_time))
        folder = sample/'candidates/011_novoplasty_001'
        atomic_json(folder/'bold_evidence.json', {'status': 'unresolved', 'hits': [],
                                                'reason': 'No qualifying BOLD COI alignment'})
        rows = [dict(zip(d.BOLD_COLUMNS, ['COI', subject, '95', '500', '650', '650', '1', '500', '1', '500', '1e-80', score]))
                for subject, score in [('COI_B', '800'), ('COI_A', '900')]]
        put(folder/'bold_hits.tsv', '\n'.join('\t'.join(row[col] for col in d.BOLD_COLUMNS) for row in rows)+'\n')
    # One streaming pass services both samples; no full 20M-record metadata load.
    original_open = Path.open
    metadata_reads = []
    def counted_open(path, *args, **kwargs):
        if path == metadata:
            metadata_reads.append(path)
        return original_open(path, *args, **kwargs)
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(Path, 'open', counted_open)
        report = payload(d.write_dashboard(root, bold_metadata=metadata if override else None))
    assert len(metadata_reads) == 1
    for sample in report['samples']:
        mito = sample['mito']
        assert mito['reference_taxa'] == ['Observed species']
        assert mito['hits'][0]['record'] == 'A' and mito['hits'][0]['identity'] == 95
        assert mito['target_status'] == 'unresolved' and not mito['accepted']
        assert 'Expected species' not in mito['reference_taxa']
    assert str(metadata) not in (root/'dashboard.html').read_text()


def test_unknown_reference_metadata_does_not_fall_back_to_expected_identity(tmp_path):
    root = make_demo(tmp_path/'results', samples=1)
    sample = root/'W30877'
    summary = json.loads((sample/'final/summary.json').read_text())
    summary.update(target_status='unresolved')
    atomic_json(sample/'final/summary.json', summary)
    folder = sample/'candidates/011_novoplasty_001'
    atomic_json(folder/'bold_evidence.json', {'status': 'unresolved', 'hits': [
        {'sseqid': 'unknown', 'bitscore': '900', 'pident': '95', 'length': '500'}]})
    mito = d.build_data(root)['samples'][0]['mito']
    assert not mito['reference_taxa'] and mito['hits'][0]['subject'] == 'unknown'


def test_repeated_annotation_bounds_messages_are_deduplicated_in_saved_data(tmp_path):
    path = put(tmp_path/'result.gff', 'mt\tMITOS2\tCDS\t90\t110\t.\t+\t.\tName=one\n'
                                    'mt\tMITOS2\tCDS\t95\t120\t.\t+\t.\tName=two\n')
    result = d.annotations(path, 100, tmp_path)
    assert not result['features'] and len(result['warnings']) == 1

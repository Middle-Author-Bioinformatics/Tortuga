import csv
import fcntl
import json
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import itsme_reporting as r
from compact_outputs import atomic_json, table


def sample(tmp_path, rows=None):
    root = tmp_path/'SAMPLE'
    (root/'ITSME').mkdir(parents=True)
    atomic_json(root/'run_config.json', {})
    atomic_json(root/'final/summary.json', {'sample': 'SAMPLE', 'pipeline': 'both',
                'status': 'completed', 'itsme_ok': True, 'best_length_bp': 16000})
    if rows is not None:
        table(root/'ITSME/master_summary.csv', rows,
              ['contig', 'locus_type', 'length_bp', 'ITS1_coordinates', '5.8S_coordinates',
               'ITS2_coordinates', 'consensus_taxonomy', 'taxonomy_status'], ',')
    return root


def locus(name='A', **extra):
    return {'contig': name, 'locus_type': 'COMPLETE_DUAL_ANCHOR', 'length_bp': 2000,
            'ITS1_coordinates': '501-700', '5.8S_coordinates': '701-850',
            'ITS2_coordinates': '851-1100', 'consensus_taxonomy': 'Metazoa; Example; Example species',
            'taxonomy_status': 'PASS_taxonomically_coherent', **extra}


def summary(root):
    return json.loads((root/'final/summary.json').read_text())


def test_lengths_are_its_cassette_not_full_locus_and_taxonomy_uses_consensus(tmp_path):
    root = sample(tmp_path, [locus()])
    table(root/'ITSME/taxonomy_validation.tsv', [{'contig': 'A', 'consensus_species': 'Example species'}],
          ['contig', 'consensus_species'])
    result, rows = r.summarize_itsme(root, summary(root))
    assert result['itsme_locus_lengths_bp'] == [2000]
    assert result['itsme_its_lengths_bp'] == [600]
    assert result['itsme_taxa'] == ['Example species']
    assert result['itsme_taxonomy_ranks'] == ['species']
    assert [rows[0][field] for field in ['its1_length_bp', 'five_8s_length_bp', 'its2_length_bp']] == [200, 150, 250]


def test_multiple_loci_preserve_length_taxonomy_pairing_and_graph_source_mapping(tmp_path):
    root = sample(tmp_path, [locus('A'), locus('reconstructed_B', length_bp=2500,
                     ITS2_coordinates='851-1200', consensus_taxonomy='Metazoa; Otheridae')])
    table(root/'ITSME/locus_source_map.tsv', [{'locus': 'reconstructed_B', 'source_path': 'path_B'}],
          ['locus', 'source_path'])
    table(root/'ITSME/taxonomy_validation.tsv', [{'contig': 'path_B', 'consensus_family': 'Otheridae',
          'consensus_taxonomy': 'Metazoa; Otheridae'}], ['contig', 'consensus_family', 'consensus_taxonomy'])
    result, rows = r.summarize_itsme(root, summary(root))
    assert result['itsme_complete_locus_count'] == 2
    assert result['itsme_its_lengths_bp'] == [600, 700]
    assert result['itsme_taxa'] == ['Example species', 'Otheridae']
    assert result['itsme_taxonomy_ranks'] == ['', 'family']
    assert rows[1]['locus'] == 'reconstructed_B'


def test_missing_taxonomy_and_coordinates_stay_unknown_not_expected_target(tmp_path):
    root = sample(tmp_path, [locus(consensus_taxonomy='NA', ITS1_coordinates='NA')])
    result, rows = r.summarize_itsme(root, {**summary(root), 'target_taxon': 'Expected species'})
    assert result['itsme_taxa'] == [''] and result['itsme_taxonomy_ranks'] == ['']
    assert result['itsme_its_lengths_bp'] == [None]
    assert rows[0]['its1_length_bp'] is None


@pytest.mark.parametrize('coords', [('700-501', '850-701', '1100-851'),
                                   ('1100-901', '900-751', '750-501')])
def test_reverse_intervals_are_inclusive(tmp_path, coords):
    root = sample(tmp_path, [locus(**dict(zip(['ITS1_coordinates', '5.8S_coordinates', 'ITS2_coordinates'], coords)))])
    assert r.summarize_itsme(root, summary(root))[0]['itsme_its_lengths_bp'] == [600]


def test_extracted_fasta_is_length_fallback_and_legacy_results_are_read(tmp_path):
    root = sample(tmp_path)
    table(root/'final/itsme_master_summary.tsv', [locus(ITS1_coordinates='NA')], list(locus()))
    (root/'final/itsme_complete_ITS.fasta').write_text('>A\nACGT\nACGT\n')
    result, rows = r.summarize_itsme(root, summary(root))
    assert result['itsme_its_lengths_bp'] == [8]
    assert rows[0]['source_report'].endswith('itsme_master_summary.tsv')


def test_no_loci_is_distinct_from_missing_report_or_unrequested_branch(tmp_path):
    root = sample(tmp_path)
    assert r.summarize_itsme(root, summary(root))[0]['itsme_complete_locus_count'] is None
    assert r.summarize_itsme(root, {'pipeline': 'mitolooper'})[0]['itsme_reporting_status'] == 'not_requested'
    table(root/'ITSME/master_summary.csv', [], ['contig'], ',')
    cells, rows = r.summarize_itsme(root, summary(root))
    assert cells['itsme_complete_locus_count'] == 0 and cells['itsme_reporting_status'] == 'no_complete_loci'


def test_chimeras_and_partial_rows_are_not_counted(tmp_path):
    root = sample(tmp_path, [locus('A'), locus('partial', locus_type='PARTIAL_ITS'),
                            locus('chimera', taxonomy_status='CHIMERA')])
    assert r.summarize_itsme(root, summary(root))[0]['itsme_locus_ids'] == ['A']


def test_reporting_refresh_preserves_mito_columns_and_never_cleans_files(tmp_path):
    root = sample(tmp_path, [locus()])
    heavy = root/'02_itsme/heavy.fastq.gz'
    heavy.parent.mkdir()
    heavy.write_bytes(b'preserve intermediates')
    table(tmp_path/'cohort_summary.csv', [{'sample': 'SAMPLE', 'best_length_bp': 16000,
          'mito_custom': 'preserve'}], ['sample', 'best_length_bp', 'mito_custom'], ',')
    assert r.main([str(tmp_path)]) == 0
    with (tmp_path/'cohort_summary.csv').open() as handle:
        row = next(csv.DictReader(handle))
    assert row['mito_custom'] == 'preserve' and row['best_length_bp'] == '16000'
    assert json.loads(row['itsme_its_lengths_bp']) == [600]
    with (tmp_path/'cohort_its_loci.csv').open() as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 1 and rows[0]['its_length_bp'] == '600'
    assert heavy.read_bytes() == b'preserve intermediates'
    assert summary(root)['itsme_complete_locus_count'] == 1
    assert not (root/'cleanup.json').exists()


def test_active_pipeline_blocks_reporting_changes(tmp_path):
    root = sample(tmp_path, [locus()])
    before = (root/'final/summary.json').read_bytes()
    with (tmp_path/'.pipeline.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with pytest.raises(SystemExit):
            r.main([str(tmp_path)])
    assert (root/'final/summary.json').read_bytes() == before
    assert not (tmp_path/'cohort_its_loci.csv').exists()


def test_incomplete_or_stale_sample_is_not_refreshed(tmp_path):
    import os
    root = sample(tmp_path, [locus()])
    before = (root/'final/summary.json').read_bytes()
    newer = (root/'final/summary.json').stat().st_mtime_ns + 1000000
    os.utime(root/'run_config.json', ns=(newer, newer))
    assert r.main([str(tmp_path)]) == 0
    assert (root/'final/summary.json').read_bytes() == before

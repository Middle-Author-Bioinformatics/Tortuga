import csv
import json
from pathlib import Path
import subprocess
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bold_taxonomy import normalize_taxonomy
import build_dashboard as dashboard
import mito_pipeline as pipeline
from report_evidence import coi_support
from dashboard_fixture import make_demo


@pytest.mark.parametrize('missing', [None, '', 'None', 'NULL', 'NaN', 'N/A', '-'])
def test_deepest_observed_taxon_recovered_without_placeholder(missing):
    assert normalize_taxonomy(missing, 'Animalia;Porifera;None;Tethya;Tethya actinia;None') == (
        'Tethya actinia', 'Animalia;Porifera;Tethya;Tethya actinia')
    assert normalize_taxonomy(missing, 'None;null;N/A') == ('', '')


def test_legacy_metadata_recovers_name_and_preserves_target_comparison(tmp_path):
    path = tmp_path/'metadata.tsv'
    pipeline.write_table(path, [{'sequence_id':'A', 'processid':'A', 'taxon':'None',
        'lineage':'Animalia;Porifera;Tethya;Tethya actinia;None', 'source':'test'}])
    metadata = pipeline.load_bold_metadata(path)
    assert metadata['A']['taxon'] == 'Tethya actinia'
    row = dict(qseqid='cox1', sseqid='A', pident='96.598', length='1558', qlen='1563',
               slen='1560', bitscore='2571')
    strict = pipeline.bold_decision([row], metadata, 'Tethya', 97, .8, 150, .02, 500)
    relaxed = pipeline.bold_decision([row], metadata, 'Tethya', 95, .8, 150, .02, 500)
    assert strict['status'] == 'unresolved' and relaxed['status'] == 'confirmed_target'
    # A lower identity gate does not substitute the expected species for the hit.
    assert pipeline.bold_decision([row], metadata, 'Tethya californiana', 95, .8, 150, .02, 500)['status'] == 'non_target'


@pytest.mark.parametrize('saved', [True, False])
def test_dashboard_recovers_existing_names_without_blast_or_metadata_rewrite(tmp_path, saved):
    root = make_demo(tmp_path/'results', samples=1)
    sample = root/'W30877'
    evidence_path = sample/'candidates/011_novoplasty_001/bold_evidence.json'
    evidence = json.loads(evidence_path.read_text())
    old_meta = {'processid':'A', 'taxon':'None',
                'lineage':'Animalia;Arthropoda;Rathbunixa;Rathbunixa occidentalis;None'}
    if saved:
        evidence['hits'][0]['metadata'] = old_meta
    else:
        metadata = tmp_path/'metadata.tsv'
        pipeline.write_table(metadata, [dict(old_meta, sequence_id='COI_A')])
        evidence['hits'] = []
        row = ['cox1','COI_A','96.598','1558','1563','1560','1','1558','1','1558','0','2571']
        (evidence_path.parent/'bold_hits.tsv').write_text('\t'.join(row)+'\n')
    evidence_path.write_text(json.dumps(evidence))
    original = evidence_path.read_bytes()
    result = dashboard.build_data(root, bold_metadata=metadata if not saved else None)['samples'][0]['mito']
    assert result['reference_taxa'] == ['Rathbunixa occidentalis']
    assert result['hits'][0]['taxon'] == 'Rathbunixa occidentalis'
    assert 'None' not in result['hits'][0]['lineage']
    assert evidence_path.read_bytes() == original


def test_confidence_uses_shared_observed_genus_with_missing_subspecies():
    hit = dict(subject='A', bitscore=2571, taxon='None',
               lineage='Animalia;Porifera;Tethya;Tethya actinia;None',
               identity=96.598, aligned_bp=1558, coverage=.997)
    other = dict(hit, subject='B', bitscore=2556,
                 lineage='Animalia;Porifera;Tethya;Tethya minuta;None')
    result = coi_support([hit, other], {})
    assert result['taxon'] == 'Tethya' and result['grade'] == 'Mixed'
    assert hit['taxon'] == 'None'  # Reporting does not mutate stored evidence.


def test_prepare_ignores_missing_terminal_ranks(tmp_path):
    source = tmp_path/'bold.tsv'
    source.write_text('processid\tmarker_code\tnuc\tkingdom\tgenus\tspecies\tsubspecies\n'
                      'A\tCOI\t'+'ACGT'*50+'\tAnimalia\tTethya\tTethya actinia\tNone\n')
    output = tmp_path/'prepared'
    result = subprocess.run([sys.executable, str(Path(pipeline.__file__).parent/'prepare_bold.py'),
                             str(source), '--source','test','-o',str(output)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    with (output/'metadata.tsv').open() as handle:
        row = next(csv.DictReader(handle, delimiter='\t'))
    assert row['taxon'] == 'Tethya actinia'
    assert row['lineage'] == 'Animalia;Tethya;Tethya actinia'

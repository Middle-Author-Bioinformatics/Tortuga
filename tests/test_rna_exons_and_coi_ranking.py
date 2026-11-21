import json
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import build_dashboard as dashboard
import mito_pipeline as pipeline
from report_evidence import read_annotation, reverse_complement, coi_support
from dashboard_fixture import put


def test_mitos_rrna_parent_envelope_uses_explicit_exons_before_and_after_parent(tmp_path):
    path = put(tmp_path/'result.gff',
        'W31021\tmitfi\texon\t1\t879\t0\t+\t.\tParent=transcript_rrnL_0;Name=rrnL_0-b\n'
        'W31021\tmitos\tgene\t1107\t1619\t.\t+\t.\tID=gene_nad6;Name=nad6\n'
        'W31021\tmitos\texon\t1107\t1619\t.\t+\t0\tParent=transcript_nad6;Name=nad6\n'
        'W31021\tmitfi\tncRNA_gene\t1\t18475\t.\t+\t.\tID=gene_rrnL_0;Name=rrnL_0\n'
        'W31021\tmitfi\trRNA\t1\t18475\t.\t+\t.\tID=transcript_rrnL_0;Name=rrnL_0;Parent=gene_rrnL_0\n'
        'W31021\tmitfi\texon\t17751\t18475\t0\t+\t.\tParent=transcript_rrnL_0;Name=rrnL_0-a\n')
    sequence = 'ACGT'*4737+'A'
    before = path.read_bytes()
    features = read_annotation(path, 18949, circular=True, sequence=sequence)['features']
    rna = next(feature for feature in features if feature['type'] == 'rRNA')
    assert rna['segments'] == [[17751,18475],[1,879]]
    assert rna['sequence'] == sequence[17750:18475]+sequence[:879]
    assert len(rna['sequence']) == 1604
    assert (rna['native_start'],rna['native_end']) == (1,18475)
    assert rna['coordinate_source'] == 'child_exons'
    assert len([feature for feature in features if feature['type'] == 'rRNA']) == 1
    assert any(feature['gene'] == 'nad6' and feature['type'] == 'CDS' for feature in features)
    assert path.read_bytes() == before


@pytest.mark.parametrize('strand', ['+','-'])
def test_split_rna_copies_only_explicit_parts_in_transcript_order(tmp_path,strand):
    path = put(tmp_path/'rna.gff',
        f'mt\tmitfi\trRNA\t1\t50\t.\t{strand}\t.\tID=rna;Name=rrnL\n'
        f'mt\tmitfi\texon\t2\t6\t.\t{strand}\t.\tParent=rna;Name=rrnL-b\n'
        f'mt\tmitfi\texon\t40\t44\t.\t{strand}\t.\tParent=rna;Name=rrnL-a\n')
    nt = 'ACGTTGCAAA'*5
    rna = read_annotation(path,50,True,nt)['features'][0]
    expected = nt[39:44]+nt[1:6] if strand == '+' else reverse_complement(nt[39:44])+reverse_complement(nt[1:6])
    assert rna['sequence'] == expected


def test_invalid_rna_exon_does_not_fall_back_to_parent_envelope(tmp_path):
    path = put(tmp_path/'rna.gff',
        'mt\tmitfi\trRNA\t1\t50\t.\t+\t.\tID=rna;Name=rrnL\n'
        'mt\tmitfi\texon\t40\t51\t.\t+\t.\tParent=rna;Name=rrnL\n')
    result = read_annotation(path,50,True,'A'*50)
    assert not result['features'] and result['warnings']


def bold_rows():
    return [
        dict(qseqid='COI',sseqid='weak',pident='80.796',length='1557',qlen='1593',slen='1593',bitscore='1461'),
        dict(qseqid='COI',sseqid='target',pident='99.696',length='658',qlen='1593',slen='658',bitscore='1178'),
        dict(qseqid='COI',sseqid='target2',pident='99.543',length='657',qlen='1593',slen='657',bitscore='1172')]


def metadata():
    return {row['sseqid']:{'taxon':'Meretrix lamarckii' if row['sseqid']=='weak' else 'Tivela stultorum',
                          'lineage':'Animalia;Mollusca;'+('Meretrix;Meretrix lamarckii' if row['sseqid']=='weak' else 'Tivela;Tivela stultorum')}
            for row in bold_rows()}


def test_long_divergent_homolog_does_not_veto_high_identity_short_barcode():
    result = pipeline.bold_decision(bold_rows(),metadata(),'Tivela stultorum',97,.8,150,.02,500)
    assert result['status'] == 'confirmed_target'
    assert {row['sseqid'] for row in result['hits']} == {'target','target2'}
    # Genuine near-best conflicts and incomplete search windows still block confirmation.
    meta = metadata();meta['target2']['taxon'] = 'Other species';meta['target2']['lineage']='Animalia;Other species'
    assert pipeline.bold_decision(bold_rows(),meta,'Tivela stultorum',97,.8,150,.02,500)['reason_code'] == 'taxon_conflict'
    assert pipeline.bold_decision(bold_rows(),metadata(),'Tivela stultorum',97,.8,150,.02,3)['reason_code'] == 'near_best_truncated'
    assert pipeline.bold_decision(bold_rows(),metadata(),'Tivela stultorum',100,.8,150,.02,500)['reason_code'] == 'no_qualifying_alignment'


@pytest.mark.parametrize('raw', [True,False])
def test_dashboard_uses_qualifying_barcodes_and_keeps_no_hit_fallback(tmp_path,raw):
    rows = bold_rows();meta=metadata()
    evidence = {'hits':[dict(row,metadata=meta[row['sseqid']]) for row in rows],
                'alignment_thresholds':{'identity_percent':95,'shorter_sequence_coverage':.8,'minimum_aligned_bp':150}}
    path = None
    if raw:
        path = put(tmp_path/'hits.tsv','\n'.join('\t'.join(str(row.get(key,1)) for key in dashboard.BOLD_COLUMNS) for row in rows)+'\n')
    original = json.dumps(evidence)
    hits,count = dashboard.coi_matches(evidence,path)
    assert count == 2 and [hit['subject'] for hit in hits] == ['target','target2']
    support = coi_support(hits,evidence)
    assert support['taxon'] == 'Tivela stultorum' and support['grade'] == 'High'
    assert json.dumps(evidence) == original
    evidence['alignment_thresholds']['identity_percent'] = 100
    hits,count = dashboard.coi_matches(evidence,path)
    assert count == 3 and hits[0]['taxon'] == 'Meretrix lamarckii'

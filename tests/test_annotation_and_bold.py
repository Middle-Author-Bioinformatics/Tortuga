import csv
import fcntl
import json
from pathlib import Path
import sys

import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import mito_pipeline as m
import annotate_existing as helper
from test_pipeline import fixture, summary


def hit(subject,score=1000,**extra):
    return {'qseqid':'COI','sseqid':subject,'pident':'99','length':'600','qlen':'650','slen':'650',
            'bitscore':str(score),**extra}


def test_full_hit_page_with_score_gap_can_confirm_target():
    rows=[hit('A'),hit('B',400),hit('C',300)]
    meta={'A':{'taxon':'Target','lineage':'Animalia;Target'}}
    decision=m.bold_decision(rows,meta,'Target',97,.8,400,.02,3)
    assert decision['status']=='confirmed_target' and decision['hit_limit_reached']
    assert not decision['near_best_truncated']


def test_tied_page_remains_unresolved_at_limit_and_duplicate_hsps_are_not_subjects():
    meta={key:{'taxon':'Target','lineage':'Animalia;Target'} for key in ['A','B']}
    assert m.bold_decision([hit('A'),hit('B')],meta,'Target',97,.8,400,.02,2)['reason_code']=='near_best_truncated'
    result=m.bold_decision([hit('A'),hit('A')],meta,'Target',97,.8,400,.02,2)
    assert result['status']=='confirmed_target' and result['returned_subjects']==1


def test_bold_expands_only_when_near_best_evidence_is_cut_off(fixture):
    code,events=fixture['run']({'novoplasty':{'circular':True,'complete':True,'target':'expand'}},
                             ['--no-assembly','--bold-max-hits','1','--bold-search-max-hits','10'])
    assert code==0 and summary(fixture)['target_status']=='confirmed_target'
    calls=[event for event in events if event['tool']=='blastn' and '-db' in event['args']]
    assert [event['args'][event['args'].index('-max_target_seqs')+1] for event in calls]==['1','10']
    path=fixture['out']/'W30877/candidates/011_novoplasty_001/bold_evidence.json'
    evidence=json.loads(path.read_text())
    assert len(evidence['search_attempts'])==2 and evidence['search_limit']==10


def test_bold_hard_limit_does_not_turn_ties_into_false_confirmation(fixture):
    code,events=fixture['run']({'novoplasty':{'circular':True,'complete':True}},
                             ['--no-assembly','--bold-max-hits','1','--bold-search-max-hits','1'])
    result=summary(fixture)
    assert code==0 and result['target_status']=='unresolved' and not result['accepted_by_policy']
    assert result['annotation_ok'] and result['annotation_scope']=='exploratory_candidate'


def test_unclosed_unresolved_candidate_is_annotated_by_default_without_promotion(fixture):
    code,events=fixture['run']({'novoplasty':{'complete':True,'target':'mixed'}},['--no-assembly'])
    result=summary(fixture)
    assert code==0 and result['annotation_ok'] and result['annotation_status']=='completed'
    assert not result['accepted_by_policy'] and result['target_status']=='unresolved'
    annotation=next(event for event in events if event['tool']=='mitos2')
    assert '--linear' in annotation['args']
    assert not (fixture['out']/'W30877/MitoLooper/accepted_candidate.fasta').exists()
    assert (fixture['out']/'W30877/MitoLooper/mitos_result.gff').is_file()


@pytest.mark.parametrize('policy',['accepted','none'])
def test_annotation_policy_can_explicitly_disable_exploratory_annotation(fixture,policy):
    code,events=fixture['run']({'novoplasty':{'complete':True}},['--no-assembly','--annotation-policy',policy])
    assert code==0 and not any(event['tool']=='mitos2' for event in events)
    assert not summary(fixture)['accepted_by_policy']


def test_existing_compacted_candidates_are_annotated_without_reassembly_or_itsme(fixture):
    fixture['argv'].remove('--keep-intermediates')
    code,events=fixture['run']({'novoplasty':{'complete':True,'target':'mixed'}},
                             ['--pipeline','both','--itsme-db',str(fixture['itsmedb']),'--no-assembly','--annotation-policy','accepted'])
    assert code==0
    before=len(events)
    seq=fixture['out']/'W30877/MitoLooper/best_sequence.fasta'
    content=seq.read_bytes()
    assert helper.main([str(fixture['out'])])==0
    trace=[json.loads(line) for line in fixture['trace'].read_text().splitlines()]
    assert [row['tool'] for row in trace[before:]]==['mitos2']
    assert seq.read_bytes()==content
    result=summary(fixture)
    assert result['annotation_ok'] and result['annotation_scope']=='exploratory_candidate'
    assert result['itsme_ok'] and not result['accepted_by_policy'] and result['target_status']=='unresolved'
    assert not list((fixture['out']/'W30877/05_mitos2').glob('attempt_*'))
    with (fixture['out']/'cohort_summary.csv').open() as handle:
        row=next(csv.DictReader(handle))
    assert row['annotation_ok']=='True' and row['annotation_scope']=='exploratory_candidate'
    count=len(trace)
    assert helper.main([str(fixture['out'])])==0
    assert len(fixture['trace'].read_text().splitlines())==count


def test_existing_candidate_can_be_revalidated_and_accepted_under_original_policy(fixture):
    fixture['argv'].remove('--keep-intermediates')
    code,events=fixture['run']({'novoplasty':{'circular':True,'complete':True,'target':'mixed'}},
                             ['--no-assembly','--annotation-policy','accepted'])
    assert code==0 and not summary(fixture)['accepted_by_policy']
    config_path=fixture['sheet'].parent/'scenario.json'
    scenario=json.loads(config_path.read_text())
    scenario['methods']['novoplasty']['target']='yes'
    config_path.write_text(json.dumps(scenario))
    # The helper queries the retained full sequence, so the double's default target is affirmative.
    assert helper.main([str(fixture['out']),'--recheck-bold'])==0
    result=summary(fixture)
    assert result['accepted_by_policy'] and result['target_status']=='confirmed_target' and result['annotation_ok']
    evidence=json.loads((fixture['out']/'W30877/candidates/011_novoplasty_001/bold_evidence.json').read_text())
    assert evidence['query_source']=='saved_protein_localization'
    assert (fixture['out']/'W30877/MitoLooper/accepted_candidate.fasta').is_file()
    assert not list((fixture['out']/'postprocess_bold_reference').glob('coi.*'))
    trace=[json.loads(line) for line in fixture['trace'].read_text().splitlines()]
    assert all(event['tool'] not in {'novoplasty','getorganelle','mitofinder','mtgrasp','metaspades','itsme'} for event in trace[len(events):])


def test_active_run_blocks_existing_annotation(fixture):
    code,events=fixture['run']({'novoplasty':{'complete':True}},['--no-assembly','--annotation-policy','accepted'])
    with (fixture['out']/'.pipeline.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        with pytest.raises(SystemExit):
            helper.main([str(fixture['out'])])
    assert len(fixture['trace'].read_text().splitlines())==len(events)


def test_disabling_annotation_on_an_accepted_genome_is_not_an_execution_failure(fixture):
    code,events=fixture['run']({'novoplasty':{'circular':True,'complete':True}},
                             ['--no-assembly','--annotation-policy','none'])
    assert code==0 and summary(fixture)['status']=='accepted_annotation_disabled'
    assert summary(fixture)['accepted_by_policy'] and not any(event['tool']=='mitos2' for event in events)


def test_missing_code_prevents_guessing_for_exploratory_annotation(fixture):
    from Bio import SeqIO
    row={'sample':'W30877','target_taxon':'Target species'}
    m.write_table(fixture['sheet'],[row],list(row))
    code,events=fixture['run']({'getorganelle':{'complete':True}},
                             ['--genbank-dir',str(fixture['sheet'].parent/'absent'),'--no-assembly'])
    assert code==0 and summary(fixture)['annotation_status']=='skipped_unknown_code'
    assert not any(event['tool']=='mitos2' for event in events)


def test_bold_recheck_failure_does_not_block_existing_candidate_annotation(fixture):
    fixture['argv'].remove('--keep-intermediates')
    code,events=fixture['run']({'novoplasty':{'complete':True,'target':'mixed'}},
                             ['--no-assembly','--annotation-policy','accepted'])
    scenario_path=fixture['sheet'].parent/'scenario.json'
    scenario=json.loads(scenario_path.read_text())
    scenario['fail']=['blastn']
    scenario_path.write_text(json.dumps(scenario))
    assert helper.main([str(fixture['out']),'--recheck-bold'])==1
    result=summary(fixture)
    assert result['annotation_ok'] and result['annotation_scope']=='exploratory_candidate'
    assert result['bold_recheck_status']=='failed' and not result['accepted_by_policy']


def test_force_repeats_annotation_even_when_prior_helper_scratch_was_retained(fixture):
    code,events=fixture['run']({'novoplasty':{'complete':True}},['--no-assembly','--annotation-policy','accepted'])
    assert helper.main([str(fixture['out']),'--keep-intermediates'])==0
    assert helper.main([str(fixture['out']),'--force','--keep-intermediates'])==0
    rows=[json.loads(line) for line in fixture['trace'].read_text().splitlines()]
    assert sum(row['tool']=='mitos2' for row in rows)==2

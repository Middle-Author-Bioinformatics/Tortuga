import json
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import mito_pipeline as pipeline
from compact_outputs import finished_sample
from test_pipeline import fixture, summary


def run_circle(fixture, methods, extra=()):
    return fixture['run'](methods, ['--recovery-policy','circle-first',*extra])


def test_defaults_stop_on_circle_and_leave_assembly_unlimited():
    args = pipeline.parser().parse_args([])
    assert args.recovery_policy == 'circle-first'
    assert args.assembly_timeout == 0 and args.assembler == 'metaspades'
    assert args.linear_ranking == 'genes'
    assert pipeline.parser().parse_args(['--no-assemble']).no_assembly


def test_first_read_circle_stops_mito_despite_taxonomy_and_gene_gates_but_runs_itsme(fixture):
    code, events = run_circle(fixture, {'novoplasty':{'circular':True,'complete':False,'target':'mixed'},
                                      'getorganelle':{'circular':True,'complete':True}},
                              ['--branch','full','--itsme-db',str(fixture['itsmedb'])])
    tools = [event['tool'] for event in events]
    assert code == 0 and 'itsme' in tools and 'mitos2' in tools
    assert all(tool not in tools for tool in ['getorganelle','mitofinder','mtgrasp','metaspades','megahit'])
    result = summary(fixture)
    assert result['circular'] and not result['accepted_by_policy']
    assert result['target_status'] == 'unresolved'
    assert result['recovery_selection'] == 'first_circular' and result['winning_method'] == 'novoplasty'
    assert result['status'] == 'circular_recovered_annotated'
    assert finished_sample(fixture['out']/'W30877')


def test_second_read_circle_stops_remaining_methods_and_wgs(fixture):
    code, events = run_circle(fixture, {'novoplasty':{'circular':False,'complete':True},
                                      'getorganelle':{'circular':True,'complete':False,'target':'no'}})
    tools = [event['tool'] for event in events]
    assert code == 0 and 'novoplasty' in tools and 'getorganelle' in tools
    assert all(tool not in tools for tool in ['mitofinder','mtgrasp','metaspades'])
    assert summary(fixture)['winning_method'] == 'getorganelle'


def test_wgs_runs_without_deadline_and_stops_on_first_assembly_circle(fixture):
    methods = {method:{'circular':False,'complete':False} for method in pipeline.READ_METHODS}
    methods['getorganelle_assembly'] = {'circular':True,'complete':False,'target':'mixed'}
    code, events = run_circle(fixture,methods)
    tools = [event['tool'] for event in events]
    assert code == 0 and 'metaspades' in tools and 'getorganelle_assembly' in tools
    assert not any(event['tool']=='mitofinder' and '-a' in event['args'] for event in events)
    assert summary(fixture)['winning_method'] == 'getorganelle_assembly'
    state = json.loads((fixture['out']/'W30877/state/wgs.json').read_text())
    assert state['timeout_seconds'] is None


@pytest.mark.parametrize('flag', ['--no-assembly','--no-assemble'])
def test_no_assembly_retains_gene_rich_linear_without_stopping_on_legacy_length_acceptance(fixture,flag):
    methods = {'novoplasty':{'circular':False,'length':11000,'genes':['cox1']},
               'getorganelle':{'circular':False,'length':10000,'complete':True}}
    code, events = run_circle(fixture,methods,[flag,'--accept-topology','any','--accept-rule','length'])
    tools = [event['tool'] for event in events]
    assert code == 0 and 'getorganelle' in tools and 'metaspades' not in tools
    result = summary(fixture)
    assert result['best_method'] == 'getorganelle' and result['best_length_bp'] == 10000
    assert not result['circular'] and result['recovery_selection'] == 'best_linear'
    assert result['status'] == 'linear_recovered_annotated'


def test_no_circle_after_wgs_uses_best_assembly_derived_linear(fixture):
    methods = {'novoplasty':{'circular':False,'length':11000,'genes':['cox1']},
               'getorganelle_assembly':{'circular':False,'length':10000,'complete':True}}
    code, events = run_circle(fixture,methods)
    assert code == 0 and any(event['tool']=='metaspades' for event in events)
    assert any(event['tool']=='mitofinder' and '-a' in event['args'] for event in events)
    assert summary(fixture)['best_method'] == 'getorganelle_assembly'
    assert not summary(fixture)['circular']


def test_explicit_length_priority_changes_linear_selection(fixture):
    methods = {'novoplasty':{'circular':False,'length':11000,'genes':['cox1']},
               'getorganelle':{'circular':False,'length':10000,'complete':True}}
    code, events = run_circle(fixture,methods,['--no-assembly','--linear-ranking','length'])
    assert code == 0 and summary(fixture)['best_method'] == 'novoplasty'


def test_linear_ranking_has_no_hidden_taxonomy_or_acceptance_override():
    args = pipeline.parser().parse_args([])
    short = dict(length_bp=9000,genes_found=['cox1'],rrnas_found=[],ambiguous_bp=0,accepted=True,target_status='confirmed_target')
    rich = dict(short,length_bp=10000,genes_found=['cox1','nad1'],accepted=False,target_status='unresolved')
    assert pipeline.selected_candidate([short,rich],args) is rich
    # Without a reference-gene screen the fallback naturally uses sequence length.
    assert pipeline.selected_candidate([dict(short,genes_found=[]),dict(rich,genes_found=[])],args)['length_bp'] == 10000

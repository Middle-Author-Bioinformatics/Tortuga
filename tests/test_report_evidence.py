import json
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import build_dashboard as dashboard
from report_evidence import read_annotation, core_gene_counts, reverse_complement, coi_support
from graph_evidence import SourceGraph
from dashboard_fixture import make_demo, put


def test_mitos_gene_features_are_cds_even_with_rna_and_aliases(tmp_path):
    path = put(tmp_path/'result.gff',
        'mt\tMITOS2\tgene\t1\t9\t.\t+\t.\tName=COI\n'
        'mt\tMITOS2\trRNA\t10\t18\t.\t+\t.\tName=rrnL\n'
        'mt\tMITOS2\ttRNA\t19\t21\t.\t+\t.\tName=trnA\n')
    result = read_annotation(path, 30, sequence='ACGT'*7+'AC')
    assert [row['type'] for row in result['features']] == ['CDS', 'rRNA', 'tRNA']
    assert result['features'][0]['sequence'] == 'ACGTACGTA'


def test_parent_child_cds_and_split_features_are_not_double_counted(tmp_path):
    path = put(tmp_path/'result.gff',
        'mt\tMITOS2\tgene\t1\t20\t.\t-\t.\tID=g1;Name=cox1\n'
        'mt\tMITOS2\tCDS\t1\t4\t.\t-\t0\tID=c1;Parent=g1\n'
        'mt\tMITOS2\tCDS\t9\t12\t.\t-\t0\tID=c1;Parent=g1\n')
    sequence = 'ACGTTGCAAACC'+'A'*8
    result = read_annotation(path, len(sequence), sequence=sequence)
    assert len(result['features']) == 1
    assert result['features'][0]['sequence'] == reverse_complement(sequence[:4]+sequence[8:12])


def test_origin_crossing_feature_extracts_correct_stranded_sequence(tmp_path):
    path = put(tmp_path/'result.gff', 'mt\tMITOS2\tgene\t17\t25\t.\t-\t.\tName=cox1\n')
    sequence = 'ACGTAGCATCGATTCGACTA'
    feature = read_annotation(path, 20, circular=True, sequence=sequence)['features'][0]
    assert feature['segments'] == [[17, 20], [1, 5]]
    assert feature['sequence'] == reverse_complement(sequence[16:]+sequence[:5])
    assert not read_annotation(path, 20, circular=False, sequence=sequence)['features']


def test_core_counts_use_unique_13_cds_and_two_rrnas_not_extra_reference_genes():
    candidate = {'genes_found': ['COI','cox1','COII','COIII','cob','ND1','ND2','ND3','ND4','ND4L','ND5','ND6','atp6','atp8','mutS'],
                 'rrnas_found':['rrnS','rrnL','12S_rRNA']}
    result = core_gene_counts(candidate)
    assert result == {'cds':13, 'rrna':2, 'extra':['muts']}
    assert core_gene_counts({})['cds'] is None


def test_selected_paths_explain_locus_length_and_shared_nodes(tmp_path):
    root = make_demo(tmp_path/'results', samples=1)
    data = dashboard.build_data(root, export_graphs=True)['samples'][0]
    graph = data['its']['graphs'][0]
    assert len(graph['nodes']) == 5
    assert len([node for node in graph['nodes'] if len(node['loci']) == 2]) == 1
    path = graph['paths']['rDNA_locus_1']
    assert path['sequence_verified']
    assert sum(step['added_bp'] for step in path['steps']) == path['locus_bp'] == 7230
    assert sum(step['overlap'] for step in path['steps']) == 100
    assert all(node['taxa'] == ['Rathbunixa occidentalis'] for node in graph['nodes'])
    assert all(node['taxonomy_scope'] == 'Locus consensus' for node in graph['nodes'])
    assert not any(node['id'] == '99' for graph in data['mito']['graphs'] for node in graph['nodes'])
    # Retained reduced graphs still work after the full source graphs disappear.
    sample = root/'W30877'
    for path in [sample/'ITSME/assembly_graph_with_scaffolds.gfa', sample/'ITSME/assembly_graph.fastg', sample/'MitoLooper/getorganelle_assembly_graph.gfa']:
        path.unlink()
    later = dashboard.build_data(root)['samples'][0]
    assert later['its']['graphs'][0]['paths'] == graph['paths']
    assert later['mito']['graphs']


def test_graph_ids_without_matching_sequence_do_not_establish_a_path(tmp_path):
    root = make_demo(tmp_path/'results', samples=1)
    path = root/'W30877/ITSME/assembly_graph_with_scaffolds.gfa'
    path.write_text(path.read_text().replace('A','T'))
    assert not dashboard.build_data(root)['samples'][0]['its']['graphs']


def test_fastg_overlap_is_inferred_only_after_exact_reconstruction(tmp_path):
    path = put(tmp_path/'graph.fastg', '>EDGE_1_length_6_cov_2:EDGE_2_length_6_cov_2;\nACGTAC\n>EDGE_2_length_6_cov_2;\nTACGGA\n')
    graph = SourceGraph(path)
    evidence = graph.path_evidence([['1+','2+']], 'ACGTACGGA')
    assert evidence['sequence_verified'] and evidence['assembled_bp'] == 9
    assert evidence['steps'][1]['overlap'] == 3
    assert not graph.path_evidence([['1+','2+']], 'TTTTTTTTT')['sequence_verified']


def test_no_observed_coi_does_not_show_expected_taxonomy_even_if_saved_confirmed(tmp_path):
    root = make_demo(tmp_path/'results', samples=1)
    put(root/'W30877/candidates/011_novoplasty_001/bold_evidence.json', '{}')
    data = dashboard.build_data(root)['samples'][0]
    assert data['expected_taxon'] == 'Rathbunixa occidentalis'
    assert not data['mito']['confidence']['taxon'] and not data['mito']['reference_taxa']
    assert data['mito']['confidence']['grade'] == 'Unavailable'


def test_confidence_accounts_for_short_barcodes_conflicts_and_truncation():
    short = {'subject':'a','bitscore':400,'identity':99,'aligned_bp':200,'coverage':1,
             'taxon':'Genus species','lineage':'Animalia;Genus;Genus species'}
    assert coi_support([short], {})['grade'] == 'Moderate'
    assert coi_support([dict(short, aligned_bp=600)], {})['grade'] == 'High'
    limited = coi_support([short], {'near_best_truncated':True})
    assert limited['grade'] == 'Limited' and not limited['taxon']
    other = dict(short, subject='b', taxon='Genus other', lineage='Animalia;Genus;Genus other')
    mixed = coi_support([short, other], {})
    assert mixed['grade'] == 'Mixed' and mixed['taxon'] == 'Genus'


def test_duplicate_fasta_aliases_hidden_without_deleting_sources(tmp_path):
    root = make_demo(tmp_path/'results', samples=1)
    sample = root/'W30877'
    original = sample/'MitoLooper/best_sequence.fasta'
    alias = put(sample/'MitoLooper/annotation_input.fasta', original.read_text().replace('candidate','alias'))
    data = dashboard.build_data(root)['samples'][0]
    assert len([file for file in data['mito']['files'] if file['name'].endswith('.fasta')]) == 1
    assert original.is_file() and alias.is_file()


def test_compaction_preserves_selected_graphs_before_native_work_is_removed(tmp_path):
    import compact_outputs as compact
    root = make_demo(tmp_path/'results', samples=1)
    sample = root/'W30877'
    for branch, filename, native in [
        ('ITSME','assembly_graph_with_scaffolds.gfa','02_itsme/attempt_001/output/assembly/spades/assembly_graph_with_scaffolds.gfa'),
        ('MitoLooper','getorganelle_assembly_graph.gfa','012_getorganelle/attempt_001/output/extended_spades/assembly_graph.gfa')]:
        target = sample/native
        target.parent.mkdir(parents=True, exist_ok=True)
        (sample/branch/filename).replace(target)
    compact.compact_sample(sample)
    compact.validate_retained(sample)
    assert not (sample/'012_getorganelle').exists()
    data = dashboard.build_data(root)['samples'][0]
    assert data['mito']['graphs'] and data['its']['graphs']
    assert data['its']['graphs'][0]['paths']['rDNA_locus_1']['sequence_verified']
    assert not list((sample/'MitoLooper').glob('*.gfa'))
    assert list((sample/'MitoLooper/selected_graphs').glob('*.gfa'))


def test_short_barcode_default_and_preparation_override(tmp_path):
    import subprocess
    import mito_pipeline as pipeline
    from Bio import SeqIO
    assert pipeline.parser().parse_args([]).coi_min_bp == 150
    source = put(tmp_path/'bold.tsv', 'processid\tmarker_code\tnuc\tspecies\nSHORT\tCOI\t'+'ACGT'*50+'\tObserved species\n')
    result = subprocess.run([sys.executable, str(Path(pipeline.__file__).parent/'prepare_bold.py'),
                             str(source), '--source', 'test', '-o', str(tmp_path/'prepared')], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert len(next(SeqIO.parse(tmp_path/'prepared/coi.fasta', 'fasta')).seq) == 200
    assert json.loads((tmp_path/'prepared/provenance.json').read_text())['minimum_length_bp'] == 150


def test_annotations_from_a_different_input_are_not_projected_onto_best_sequence(tmp_path):
    root = make_demo(tmp_path/'results', samples=1)
    put(root/'W30877/MitoLooper/annotation_input.fasta', '>different\n'+'A'*15740+'\n')
    data = dashboard.build_data(root)['samples'][0]['mito']
    assert not data['annotation']['features']
    assert data['core_genes']['source'] == 'Reference screen'


def test_report_refuses_linked_graph_output_directory(tmp_path):
    root = make_demo(tmp_path/'results', samples=1)
    outside = tmp_path/'outside'
    outside.mkdir()
    (root/'W30877/ITSME/selected_graphs').symlink_to(outside, target_is_directory=True)
    with pytest.raises(ValueError, match='linked graph'):
        dashboard.write_dashboard(root)
    assert not list(outside.iterdir())

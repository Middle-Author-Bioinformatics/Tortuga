"""Illustrative data for dashboard tests and visual checks, not biological results."""
import json
import random
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from compact_outputs import atomic_json, table


def put(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(value)
    return path


def make_demo(root, samples=8):
    names = ['W30877', 'W30882', 'W30892', 'W30897', 'W30924', 'W30936', 'W31027', 'W31051']
    taxa = ['Rathbunixa occidentalis', 'Rhizophysa filiformis', 'Theodisca dendritica', 'Balanus trigonus',
            'Patiria miniata', 'Ophiothrix spiculata', 'Geodia', 'Corambe pacifica']
    results = []
    for index in range(samples):
        name, taxon = names[index], taxa[index]
        sample = root/name
        length = 15740+index*213
        accepted = index < 5
        its_count = 2 if index == 0 else 0 if index == 6 else 1
        atomic_json(sample/'run_config.json', {'version': 'demo'})
        result = {'sample': name, 'pipeline': 'both', 'status': 'recovered_annotated' if accepted else 'partial_target_confirmed',
                  'target_taxon': taxon, 'itsme_ok': True, 'best_length_bp': length, 'circular': accepted,
                  'genes_recovered': 15 if accepted else 9, 'genes_expected': 15, 'complete_by_reference': accepted,
                  'best_method': 'getorganelle', 'winning_method': 'getorganelle' if accepted else '',
                  'target_status': 'unresolved' if index == 6 else 'confirmed_target', 'annotation_ok': accepted,
                  'taxonomy': {'kingdom':'Animalia', 'phylum':'Arthropoda', 'class':'Malacostraca', 'order':'Decapoda', 'family':'Pinnotheridae', 'genus':taxon.split()[0], 'species':taxon},
                  'best_candidate': {'candidate_id': '011_novoplasty_001',
                                     'genes_found':['atp6','atp8','cox1','cox2','cox3','cytb','nad1','nad2','nad3','nad4','nad4l','nad5','nad6']+(['mutS'] if index==1 else []),
                                     'rrnas_found':['rrns','rrnl']}}
        if index == 1:
            result['genes_expected'] = 16
        atomic_json(sample/'final/summary.json', result)
        rng = random.Random(index)
        seq = ''.join(rng.choice('ACGT') for _ in range(length))
        put(sample/'MitoLooper/best_sequence.fasta', '>candidate\n'+seq+'\n')
        if accepted:
            put(sample/'MitoLooper/accepted_candidate.fasta', '>candidate\n'+seq+'\n')
            genes = [('cox1', 'CDS', 1, 1533, '+'), ('cox2', 'CDS', 1800, 2490, '+'),
                     ('atp8', 'CDS', 2700, 2870, '+'), ('atp6', 'CDS', 2870, 3560, '+'),
                     ('cox3', 'CDS', 3800, 4560, '+'), ('nad3', 'CDS', 4800, 5150, '+'),
                     ('nad5', 'CDS', 5400, 7050, '-'), ('nad4', 'CDS', 7300, 8650, '-'),
                     ('nad4l', 'CDS', 8800, 9100, '-'), ('nad6', 'CDS', 9300, 9770, '+'),
                     ('cytb', 'CDS', 9990, 11130, '+'), ('nad1', 'CDS', 11300, 12260, '-'),
                     ('rrnL', 'rRNA', 12500, 13800, '-'), ('rrnS', 'rRNA', 14000, 14700, '-'),
                     ('nad2', 'CDS', 14780, length-30, '+')]
            lines = ['##gff-version 3']
            for i, (gene, kind, start, end, strand) in enumerate(genes):
                lines.append(f'candidate\tMITOS2\t{"gene" if kind=="CDS" else kind}\t{start}\t{end}\t.\t{strand}\t.\tID={gene};Name={gene}')
                if i % 2 == 0 and start > 80:
                    lines.append(f'candidate\tMITOS2\ttRNA\t{start-75}\t{start-10}\t.\t+\t.\tID=trn{i};Name=trn{i}')
            put(sample/'MitoLooper/mitos_result.gff', '\n'.join(lines)+'\n')
        graph = ['H\tVN:Z:1.0', f'S\t1\t{seq[:5000]}\tDP:f:22', f'S\t2\t{seq[4950:10000]}\tDP:f:23', f'S\t3\t{seq[9950:]}\tDP:f:24', 'S\t99\tACGT\tDP:f:1',
                 'L\t1\t+\t2\t+\t50M', 'L\t2\t+\t3\t+\t50M', 'L\t3\t+\t1\t+\t0M']
        put(sample/'MitoLooper/getorganelle_assembly_graph.gfa', '\n'.join(graph)+'\n')
        fastg = []
        for segment in range(1, 9):
            target = f':EDGE_{segment+1}_length_{segment*140}_cov_38' if segment < 8 else ''
            fastg += [f'>EDGE_{segment}_length_{segment*140}_cov_38{target};', 'ACGT']
        fastg += [">EDGE_3_length_420_cov_38:EDGE_6_length_840_cov_38';", 'ACGT']
        put(sample/'ITSME/assembly_graph.fastg', '\n'.join(fastg)+'\n')
        rows, path_rows, locus_sequences = [], [], []
        its_graph = ['H\tVN:Z:1.0']
        shared = ''.join(rng.choice('ACGT') for _ in range(2000))
        its_graph.append(f'S\t101\t{shared}\tDP:f:38')
        tax_rows = []
        for locus in range(its_count):
            locus_id = 'rDNA_locus_'+str(locus+1)
            shift = index*17+locus*50
            row = {'contig': locus_id, 'locus_type': 'COMPLETE_DUAL_ANCHOR', 'length_bp': 7230+shift,
                   'SSU_coordinates': '1-1800', 'ITS1_coordinates': '1801-2170', '5.8S_coordinates': '2171-2330',
                   'ITS2_coordinates': f'2331-{2790+shift}', 'LSU_coordinates': f'{2791+shift}-{7230+shift}',
                   'taxonomy_status': 'PASS_taxonomically_coherent', 'consensus_taxonomy': 'Metazoa; '+taxon}
            rows.append(row)
            nt = shared+''.join(rng.choice('ACGT') for _ in range(row['length_bp']-2000))
            a, b = str(102+locus*2), str(103+locus*2)
            its_graph += [f'S\t{a}\t{nt[1950:5000]}\tDP:f:38', f'S\t{b}\t{nt[4950:]}\tDP:f:38',
                          f'L\t101\t+\t{a}\t+\t50M', f'L\t{a}\t+\t{b}\t+\t50M']
            path_rows.append({'locus':locus_id,'source_path':f'ITSME_PATH_{locus+1:06d}', 'node_path':f'101+,{a}+,{b}+', 'path_status':'PROMOTED'})
            locus_sequences.append('>'+locus_id+'\n'+nt+'\n')
            tax_rows.append({'contig': locus_id, 'consensus_species': taxon if len(taxon.split()) > 1 else 'NA',
                             'consensus_genus': taxon.split()[0]})
        columns = ['contig','locus_type','length_bp','SSU_coordinates','ITS1_coordinates','5.8S_coordinates',
                   'ITS2_coordinates','LSU_coordinates','taxonomy_status','consensus_taxonomy']
        table(sample/'ITSME/master_summary.csv', rows, columns, ',')
        table(sample/'ITSME/taxonomy_validation.tsv', tax_rows, ['contig','consensus_genus','consensus_species'])
        put(sample/'ITSME/reported_loci.fasta', ''.join(locus_sequences))
        put(sample/'ITSME/assembly_graph_with_scaffolds.gfa', '\n'.join(its_graph)+'\n')
        table(sample/'ITSME/locus_source_map.tsv', path_rows, ['locus','source_path','node_path','path_status'])
        atomic_json(sample/'candidates/011_novoplasty_001/bold_evidence.json', {
            'status': result['target_status'],
            'reason': 'BOLD near-best score window reaches the search limit' if index == 6 else 'All near-best qualifying BOLD records match requested taxon',
            'hits': [{'bitscore': None if index == 6 else '1200', 'pident': '99.3', 'length': '648', 'qlen':'650', 'slen':'650', 'metadata': {'taxon': taxon, 'processid': 'DEMO0001', 'lineage': 'Metazoa; '+taxon}}]})
        if index == 0:
            # Reproduce multiple out-of-bounds feature warnings in the native GFF.
            with (sample/'MitoLooper/mitos_result.gff').open('a') as handle:
                for start in [length+1, length+50]:
                    handle.write(f'candidate\tMITOS2\tCDS\t{start}\t{start+20}\t.\t+\t.\tName=outside\n')
        stages = [{'stage': step, 'status': status, 'elapsed_seconds': secs, 'reason': reason} for step, status, secs, reason in [
            ('qc','ok',182,''), ('novoplasty','ok',740,''), ('getorganelle','skipped',0,'Candidate already met acceptance policy'),
            ('itsme','ok',1840,''), ('wgs','timeout',7200,'Stopped at the configured wall-clock limit'), ('mtgrasp','failed',12,''), ('mitos2','ok' if accepted else 'skipped',215,'')]]
        table(sample/'stages.tsv', stages, ['stage','status','elapsed_seconds','reason'])
        atomic_json(sample/'qc/fastp.json', {'summary': {'before_filtering': {'total_reads': 82304166, 'q30_rate': .93},
                                                        'after_filtering': {'total_reads': 79421224, 'q30_rate': .962}}})
        put(sample/'logs/wgs.log', 'Illustrative timeout log; no assembly run was executed.\n')
        if index == 0:
            put(sample/'logs/mtgrasp.log', 'Illustrative failed stage log.\n')
        results.append(result)
    table(root/'cohort_summary.csv', results, ['sample','status','target_taxon','best_length_bp'], ',')
    atomic_json(root/'last_invocation.json', {'results': results})
    return root

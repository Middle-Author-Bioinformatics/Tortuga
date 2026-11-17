"""Illustrative data for dashboard tests and visual checks, not biological results."""
import json
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
                  'best_method': 'novoplasty' if accepted else 'getorganelle', 'winning_method': 'novoplasty' if accepted else '',
                  'target_status': 'confirmed_target', 'annotation_ok': accepted,
                  'best_candidate': {'candidate_id': '011_novoplasty_001'}}
        atomic_json(sample/'final/summary.json', result)
        seq = 'ACGT'*(length//4)+'A'*(length%4)
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
                lines.append(f'candidate\tMITOS2\t{kind}\t{start}\t{end}\t.\t{strand}\t.\tID={gene};Name={gene}')
                if i % 2 == 0 and start > 80:
                    lines.append(f'candidate\tMITOS2\ttRNA\t{start-75}\t{start-10}\t.\t+\t.\tID=trn{i};Name=trn{i}')
            put(sample/'MitoLooper/mitos_result.gff', '\n'.join(lines)+'\n')
        graph = ['H\tVN:Z:1.0']
        for segment in range(1, 13):
            graph.append(f'S\t{segment}\t*\tLN:i:{segment*213}\tDP:f:{20+segment}')
        for start, end in [(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,9),(9,10),(10,11),(11,12),(12,1),(3,7)]:
            graph.append(f'L\t{start}\t+\t{end}\t+\t55M')
        put(sample/'MitoLooper/getorganelle_assembly_graph.gfa', '\n'.join(graph)+'\n')
        fastg = []
        for segment in range(1, 9):
            target = f':EDGE_{segment+1}_length_{segment*140}_cov_38' if segment < 8 else ''
            fastg += [f'>EDGE_{segment}_length_{segment*140}_cov_38{target};', 'ACGT']
        fastg += [">EDGE_3_length_420_cov_38:EDGE_6_length_840_cov_38';", 'ACGT']
        put(sample/'ITSME/assembly_graph.fastg', '\n'.join(fastg)+'\n')
        rows = []
        tax_rows = []
        for locus in range(its_count):
            locus_id = 'rDNA_locus_'+str(locus+1)
            shift = index*17+locus*50
            row = {'contig': locus_id, 'locus_type': 'COMPLETE_DUAL_ANCHOR', 'length_bp': 7230+shift,
                   'SSU_coordinates': '1-1800', 'ITS1_coordinates': '1801-2170', '5.8S_coordinates': '2171-2330',
                   'ITS2_coordinates': f'2331-{2790+shift}', 'LSU_coordinates': f'{2791+shift}-{7230+shift}',
                   'taxonomy_status': 'PASS_taxonomically_coherent', 'consensus_taxonomy': 'Metazoa; '+taxon}
            rows.append(row)
            tax_rows.append({'contig': locus_id, 'consensus_species': taxon if len(taxon.split()) > 1 else 'NA',
                             'consensus_genus': taxon.split()[0]})
        columns = ['contig','locus_type','length_bp','SSU_coordinates','ITS1_coordinates','5.8S_coordinates',
                   'ITS2_coordinates','LSU_coordinates','taxonomy_status','consensus_taxonomy']
        table(sample/'ITSME/master_summary.csv', rows, columns, ',')
        table(sample/'ITSME/taxonomy_validation.tsv', tax_rows, ['contig','consensus_genus','consensus_species'])
        put(sample/'ITSME/reported_loci.fasta', ''.join('>'+row['contig']+'\n'+'ACGT'*1800+'\n' for row in rows))
        atomic_json(sample/'candidates/011_novoplasty_001/bold_evidence.json', {'status': 'confirmed_target', 'hits': [
            {'pident': '99.3', 'length': '648', 'metadata': {'taxon': taxon, 'processid': 'DEMO0001', 'lineage': 'Metazoa; '+taxon}}]})
        stages = [{'stage': step, 'status': status, 'elapsed_seconds': secs, 'reason': reason} for step, status, secs, reason in [
            ('qc','ok',182,''), ('novoplasty','ok',740,''), ('getorganelle','skipped',0,'Candidate already met acceptance policy'),
            ('itsme','ok',1840,''), ('wgs','timeout',7200,'Stopped at the configured wall-clock limit'), ('mitos2','ok' if accepted else 'skipped',215,'')]]
        table(sample/'stages.tsv', stages, ['stage','status','elapsed_seconds','reason'])
        atomic_json(sample/'qc/fastp.json', {'summary': {'before_filtering': {'total_reads': 82304166, 'q30_rate': .93},
                                                        'after_filtering': {'total_reads': 79421224, 'q30_rate': .962}}})
        put(sample/'logs/wgs.log', 'Illustrative timeout log; no assembly run was executed.\n')
        results.append(result)
    table(root/'cohort_summary.csv', results, ['sample','status','target_taxon','best_length_bp'], ',')
    atomic_json(root/'last_invocation.json', {'results': results})
    return root

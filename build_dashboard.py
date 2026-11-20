#!/usr/bin/env python3
"""Build an offline Tortuga dashboard beside the saved cohort results."""
from __future__ import annotations

import argparse
from collections import defaultdict, deque
from datetime import datetime, timezone
import csv
import fcntl
import heapq
import json
import math
from pathlib import Path
import re
from urllib.parse import quote, unquote
from bold_taxonomy import normalize_taxonomy
from report_evidence import (read_annotation, fasta_records, core_gene_counts,
                             distinct_files, coi_support)

GRAPH_NODES = 240
GRAPH_EDGES = 720
GRAPH_BYTES = 100 * 1024**2
FEATURE_LIMIT = 3000
BOLD_COLUMNS = 'qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore'.split()
BOLD_PREVIEW = 12


def number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (ValueError, TypeError):
        return None


def inside(path, root):
    try:
        path.resolve().relative_to(root.resolve())
        return path.is_file()
    except (ValueError, OSError):
        return False


def relative_file(path, root):
    if path is None or not inside(path, root):
        return None
    return {'name': path.name, 'href': quote(path.relative_to(root).as_posix(), safe='/'),
            'bytes': path.stat().st_size}


def read_json(path, warnings):
    if not path.is_file():
        return {}
    try:
        value = json.loads(path.read_text())
        return value if isinstance(value, dict) else {}
    except (ValueError, OSError) as error:
        warnings.append(f'{path.name}: could not read report ({type(error).__name__})')
        return {}


def read_table(path):
    if not path or not path.is_file():
        return []
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter=',' if path.suffix == '.csv' else '\t'))


def first(paths, root):
    return next((path for path in paths if inside(path, root)), None)


def stage_log(stage, sample, root):
    """Resolve native logs, including absolute paths saved before a folder move."""
    recorded = Path(stage['log']) if stage.get('log') else None
    candidates = []
    if recorded:
        if '..' in recorded.parts or (recorded.is_absolute() and recorded.exists() and not inside(recorded, sample)):
            return None
        candidates.extend([recorded] if recorded.is_absolute() else [sample/recorded, root/recorded])
        candidates.append(sample/'logs'/recorded.name)
    name = stage.get('stage', '')
    if name and Path(name).name == name and name not in {'.', '..'}:
        candidates.append(sample/'logs'/(name+'.log'))
    return relative_file(next((path for path in candidates if inside(path, sample)), None), root)


def choose_product(sample, branch, name, root):
    prefix = 'itsme_' if branch == 'ITSME' else ''
    return first([sample/branch/name, sample/'final'/(prefix+name)], root)


def fastg_node(token):
    token = token.strip().rstrip(';')
    reverse = token.endswith("'")
    match = re.match(r'(?:EDGE|NODE)_(.+?)_length_(\d+)(?:_cov_([\d.eE+-]+))?', token)
    if match:
        return match.group(1), '-' if reverse else '+', int(match.group(2)), number(match.group(3))
    return token.rstrip("'"), '-' if reverse else '+', None, None


def layout_graph(nodes, edges):
    """Deterministic graph-distance layers; geometry is a layout, not genomic coordinates."""
    adjacency = {node['id']: set() for node in nodes}
    for edge in edges:
        adjacency[edge['source']].add(edge['target'])
        adjacency[edge['target']].add(edge['source'])
    remaining, components = set(adjacency), []
    while remaining:
        seed = min(remaining, key=lambda key: (len(adjacency[key]), key))
        queue, levels, seen = deque([(seed, 0)]), defaultdict(list), {seed}
        while queue:
            current, depth = queue.popleft()
            levels[depth].append(current)
            for neighbor in sorted(adjacency[current]):
                if neighbor not in seen:
                    seen.add(neighbor)
                    queue.append((neighbor, depth+1))
        remaining -= seen
        components.append(levels)
    lookup, offset = {}, 0
    for component in components:
        height = max(len(level) for level in component.values()) * 78 + 70
        for level, members in component.items():
            for index, key in enumerate(members):
                lookup[key] = (85+level*150, offset+height/2+(index-(len(members)-1)/2)*78)
        offset += height
    for node in nodes:
        node['x'], node['y'] = lookup[node['id']]
    return len(components)


def parse_graph(path, root, node_limit=GRAPH_NODES, edge_limit=GRAPH_EDGES):
    result = {'file': relative_file(path, root), 'nodes': [], 'edges': [], 'warnings': [],
              'total_nodes': 0, 'total_links': 0, 'truncated': False}
    if path.stat().st_size > GRAPH_BYTES:
        result['warnings'].append('Graph exceeds the 100 MiB preview limit; use the original graph download.')
        result['truncated'] = True
        return result
    nodes, links, seen_links = {}, [], set()
    def add_node(key, length=None, coverage=None):
        if key not in nodes and len(nodes) < node_limit:
            nodes[key] = {'id': key, 'length': length, 'coverage': coverage}
        elif key in nodes and length is not None:
            nodes[key].update(length=length, coverage=coverage)
    def add_edge(source, target, source_strand='+', target_strand='+', overlap=''):
        result['total_links'] += 1
        key = (source, target, source_strand, target_strand)
        reverse = (target, source, '+' if target_strand == '-' else '-', '+' if source_strand == '-' else '-')
        if key not in seen_links and reverse not in seen_links and len(links) < edge_limit:
            links.append({'source': source, 'target': target, 'source_strand': source_strand,
                          'target_strand': target_strand, 'overlap': overlap})
            seen_links.add(key)
    with path.open(errors='replace') as handle:
        for line in handle:
            if path.suffix.lower() == '.gfa2' or (line.startswith('H\t') and 'VN:Z:2.' in line):
                result['warnings'].append('GFA2 is not supported in the inline viewer; download the original graph.')
                return result
            if path.suffix.lower() == '.fastg':
                if not line.startswith('>'):
                    continue
                source, _, neighbors = line[1:].strip().rstrip(';').partition(':')
                key, strand, length, coverage = fastg_node(source)
                if strand == '+':
                    result['total_nodes'] += 1
                add_node(key, length, coverage)
                for neighbor in filter(None, neighbors.split(',')):
                    target, target_strand, _, _ = fastg_node(neighbor)
                    add_edge(key, target, strand, target_strand)
            else:
                fields = line.rstrip('\n').split('\t')
                if fields[0] == 'S' and len(fields) >= 3:
                    result['total_nodes'] += 1
                    tags = {item.split(':', 2)[0]: item.split(':', 2)[-1] for item in fields[3:] if ':' in item}
                    length = len(fields[2]) if fields[2] != '*' else number(tags.get('LN'))
                    add_node(fields[1], length, number(tags.get('DP', tags.get('dp'))))
                elif fields[0] == 'L' and len(fields) >= 6:
                    add_edge(fields[1], fields[3], fields[2], fields[4], fields[5])
                elif fields[0] == 'E':
                    result['warnings'].append('GFA2 edge records are not rendered; download the original graph.')
    result['nodes'] = list(nodes.values())
    result['edges'] = [edge for edge in links if edge['source'] in nodes and edge['target'] in nodes]
    result['truncated'] = result['total_nodes'] > len(nodes) or result['total_links'] > edge_limit
    if result['truncated']:
        result['warnings'].append(f'Bounded preview: up to {node_limit} segments and {edge_limit} stored links. Download the full graph for all paths.')
    if not nodes:
        result['warnings'].append('No supported GFA1 segments or FASTG headers were found.')
    result['components'] = layout_graph(result['nodes'], result['edges']) if nodes else 0
    result['warnings'] = list(dict.fromkeys(result['warnings']))
    return result


def annotations(path, length, root, circular=False, sequence=''):
    result = read_annotation(path, length, circular, sequence, FEATURE_LIMIT)
    result['file'] = relative_file(path, root)
    return result


def coi_matches(evidence, raw):
    """Observed reference matches are independent of the target acceptance filter."""
    saved = evidence.get('hits', [])
    metadata = {hit.get('sseqid'): hit.get('metadata', {}) for hit in saved if hit.get('sseqid')}
    rows = saved
    count = len(saved)
    if raw:
        with raw.open(newline='') as handle:
            # Keep a bounded preview even for saturated searches. No sequences or
            # full BOLD metadata database are loaded into the dashboard builder.
            reader = csv.DictReader(handle, fieldnames=BOLD_COLUMNS, delimiter='\t')
            count = 0
            def valid_rows():
                nonlocal count
                for row in reader:
                    if row.get('sseqid') and number(row.get('bitscore')) is not None:
                        count += 1
                        yield row
            rows = heapq.nlargest(BOLD_PREVIEW, valid_rows(), key=lambda hit: number(hit.get('bitscore')))
    else:
        rows = sorted(saved, key=lambda hit: number(hit.get('bitscore')) or 0, reverse=True)[:BOLD_PREVIEW]
    hits = []
    for hit in rows:
        meta = hit.get('metadata') or metadata.get(hit.get('sseqid'), {})
        taxon, lineage = normalize_taxonomy(meta.get('taxon'), meta.get('lineage'))
        shorter = min(number(hit.get('qlen')) or 0, number(hit.get('slen')) or 0)
        hits.append({'subject': hit.get('sseqid', ''), 'bitscore': number(hit.get('bitscore')),
                     'taxon': taxon, 'lineage': lineage,
                     'identity': number(hit.get('pident')), 'aligned_bp': number(hit.get('length')),
                     'coverage': min(1, (number(hit.get('length')) or 0)/shorter) if shorter else None,
                     'record': meta.get('processid') or hit.get('sseqid', '')})
    return hits, count


def resolve_coi_names(samples, metadata_sources):
    """Scan each metadata TSV once, retaining only the small set of needed IDs."""
    groups = defaultdict(list)
    for sample, source in zip(samples, metadata_sources):
        if source:
            groups[Path(source)].append(sample)
    for path, group in groups.items():
        pending = {hit['subject'] for sample in group for hit in sample['mito']['hits']
                   if hit['subject'] and not hit['taxon']}
        if not pending or not path.is_file():
            continue
        keys = {key: subject for subject in pending for key in [subject, subject.removeprefix('lcl|')]}
        found = {}
        try:
            with path.open(newline='') as handle:
                reader = csv.DictReader(handle, delimiter='\t')
                if not {'sequence_id', 'taxon'} <= set(reader.fieldnames or []):
                    continue
                for row in reader:
                    subject = keys.get(row['sequence_id'])
                    if subject:
                        row['taxon'], row['lineage'] = normalize_taxonomy(row.get('taxon'), row.get('lineage'))
                    if subject and row.get('taxon'):
                        found[subject] = row
                        if len(found) == len(pending):
                            break
        except (OSError, csv.Error) as error:
            for sample in group:
                sample['warnings'].append(f'BOLD metadata lookup failed ({type(error).__name__})')
        for sample in group:
            for hit in sample['mito']['hits']:
                meta = found.get(hit['subject'])
                if meta:
                    hit.update(taxon=meta['taxon'], lineage=meta.get('lineage', ''),
                               record=meta.get('processid') or hit['subject'])
    for sample in samples:
        mito = sample['mito']
        hits = mito['hits']
        top = max((hit['bitscore'] for hit in hits if hit['bitscore'] is not None), default=None)
        mito['reference_ranked'] = top is not None
        mito['reference_taxa'] = sorted({hit['taxon'] for hit in hits if hit['taxon']
                                        and (top is None or hit['bitscore'] == top)})


def collect_sample(sample, root):
    from itsme_reporting import summarize_itsme, first_file, interval
    warnings = []
    summary_path = first([sample/'final/summary.json', sample/'MitoLooper/summary.json'], root)
    summary = read_json(summary_path, warnings) if summary_path else {}
    name = summary.get('sample', sample.name)
    mito_files, its_files = [], []
    for branch, files in [('MitoLooper', mito_files), ('ITSME', its_files)]:
        directory = sample/branch
        if directory.is_dir() and not directory.is_symlink():
            for path in sorted(directory.iterdir()):
                if path.is_file() and path.suffix.lower() in {'.fasta', '.fa', '.fna', '.csv', '.tsv', '.gff', '.gff3', '.gb', '.gbk', '.gfa', '.gfa2', '.fastg', '.html'}:
                    item = relative_file(path, root)
                    if item:
                        files.append(item)
    accepted = choose_product(sample, 'MitoLooper', 'accepted_candidate.fasta', root)
    best = choose_product(sample, 'MitoLooper', 'best_sequence.fasta', root)
    if accepted and (accepted.stat().st_size == 0 or summary.get('accepted_by_policy') is False):
        accepted = None
    if best and best.stat().st_size == 0:
        best = None
    selected = accepted or best
    length = number(summary.get('best_length_bp'))
    if selected:
        # Streaming lengths only; the HTML never embeds genomic sequences.
        from itsme_reporting import fasta_lengths
        lengths = list(fasta_lengths(selected).values())
        if len(lengths) == 1 and lengths[0] > 0:
            length = lengths[0]
        else:
            warnings.append('Selected mitochondrial FASTA does not contain exactly one nonempty sequence; inspect the download.')
            accepted = None
            length = None
    annotation = first([sample/'MitoLooper/mitos_result.gff', sample/'MitoLooper/mitos_result.gff3',
                        sample/'final/mitos_result.gff', sample/'final/result.gff'], root)
    if annotation is None:
        annotation = first(sorted((sample/'MitoLooper').glob('mitos_*.gff*')), root)
    annotation_report=first([sample/'MitoLooper/annotation_summary.json',sample/'final/annotation_summary.json'],root)
    if annotation_report:
        gff=read_json(annotation_report,warnings).get('gff','')
        if gff:
            recorded=first([Path(gff),sample/'MitoLooper'/Path(gff).name],root)
            annotation=recorded or annotation
    selected_records = fasta_records(selected)
    mitochondrial_nt = next(iter(selected_records.values()), '') if len(selected_records) == 1 else ''
    gene_data = annotations(annotation, length, root, bool(summary.get('circular')), mitochondrial_nt)
    annotation_input = first([sample/'MitoLooper/annotation_input.fasta'], root)
    if annotation_input:
        annotation_sequences = list(fasta_records(annotation_input).values())
        if annotation_sequences != [mitochondrial_nt]:
            gene_data['features'] = []
            warnings.append('Annotation input differs from the selected sequence')
    cells, loci = summarize_itsme(sample, summary)
    report = first_file(sample, ['master_summary.csv', 'master_summary.tsv'])
    native_rows = {row.get('contig', row.get('locus', '')): row for row in read_table(report)}
    its_nt = fasta_records(first_file(sample, ['reported_loci.fasta', 'complete_rDNA_loci.fasta']))
    for locus in loci:
        native = native_rows.get(locus['locus'], {})
        locus['regions'] = [{'name': region, 'start': positions[0], 'end': positions[1]}
            for region in ['SSU', 'ITS1', '5.8S', 'ITS2', 'LSU']
            for positions in [interval(native.get(region+'_coordinates'), locus.get('locus_length_bp'))] if positions]
        locus.pop('source_report', None)
        locus['sequence'] = its_nt.get(locus['locus'], '')
        for region in locus['regions']:
            region['sequence'] = locus['sequence'][region['start']-1:region['end']]
    graphs = {'ITSME': [], 'MitoLooper': []}
    best_candidate = summary.get('best_candidate') or {}
    evidence = {}
    candidate_id = best_candidate.get('candidate_id', '')
    if candidate_id and '/' not in candidate_id and '\\' not in candidate_id:
        candidate_dir = sample/'candidates'/candidate_id
        if inside(candidate_dir/'bold_evidence.json', root):
            evidence = read_json(candidate_dir/'bold_evidence.json', warnings)
    raw = first([sample/'candidates'/candidate_id/'bold_hits.tsv'], root) if candidate_id and '/' not in candidate_id and '\\' not in candidate_id else None
    hits, hit_count = coi_matches(evidence, raw)
    qc_path = first([sample/'qc/fastp.json', sample/'ITSME/fastp.json'], root)
    qc = read_json(qc_path, warnings) if qc_path else {}
    qsummary = qc.get('summary', {})
    before, after = qsummary.get('before_filtering', {}), qsummary.get('after_filtering', {})
    stages = read_table(sample/'stages.tsv')
    its_sequence = first_file(sample, ['reported_loci.fasta', 'complete_rDNA_loci.fasta', 'complete_ITS.fasta'])
    its_cassette = first_file(sample, ['complete_ITS.fasta'])
    if its_sequence and its_sequence.stat().st_size == 0:
        its_sequence = None
    if its_cassette and its_cassette.stat().st_size == 0:
        its_cassette = None
    logs = [item for path in sorted((sample/'logs').glob('*')) if path.suffix in {'.log', '.err'}
            for item in [relative_file(path, root)] if item]
    config = read_json(sample/'run_config.json', warnings)
    evidence.setdefault('alignment_thresholds', {}).setdefault('near_best_fraction', config.get('options', {}).get('bold_near_fraction', .02))
    saved_sample = config.get('sample') or {}
    taxonomy = summary.get('taxonomy') or (saved_sample.get('taxonomy', {}) if isinstance(saved_sample, dict) else {})
    gene_rows = read_table(first([sample/'MitoLooper/gene_recovery.csv', sample/'final/gene_recovery.csv'], root))
    counts = core_gene_counts(best_candidate, gene_rows)
    counts['source'] = 'Reference screen'
    if gene_data['features'] and not gene_data.get('multiple_sequences') and summary.get('annotation_ok') is not False:
        counts = core_gene_counts({'genes_found': [feature['gene'] for feature in gene_data['features'] if feature['type'] == 'CDS'],
                                   'rrnas_found': [feature['gene'] for feature in gene_data['features'] if feature['type'] == 'rRNA']})
        counts['source'] = 'MITOS annotation'
    result = {'id': name, 'expected_taxon': summary.get('target_taxon', ''), 'expected_taxonomy': taxonomy,
            'status': summary.get('status', 'unknown'), 'pipeline': summary.get('pipeline', 'both'),
            'summary_file': relative_file(summary_path, root), 'warnings': warnings,
            'mito': {'length': length, 'accepted': bool(accepted), 'sequence': relative_file(selected, root),
                'circular': summary.get('circular'), 'method': summary.get('winning_method') or summary.get('best_method', ''),
                'genes_recovered': summary.get('genes_recovered'), 'genes_expected': summary.get('genes_expected'),
                'core_genes': counts, 'reference_genes': gene_rows,
                'complete_by_reference': summary.get('complete_by_reference'), 'target_status': summary.get('target_status', 'unresolved'),
                'annotation_ok': summary.get('annotation_ok'), 'annotation_status': summary.get('annotation_status',''),
                'annotation_scope':summary.get('annotation_scope',''), 'target_reason':summary.get('target_reason') or evidence.get('reason',''), 'annotation': gene_data, 'hits': hits,
                'reference_taxa': [], 'hit_count': hit_count, 'bold_evidence': evidence,
                'raw_hits_file': relative_file(raw, root),
                'evidence_file': relative_file(sample/'candidates'/candidate_id/'bold_evidence.json', root) if candidate_id else None,
                'files': distinct_files(mito_files, root), 'graphs': graphs['MitoLooper']},
            'its': {'status': cells.get('itsme_reporting_status'), 'run_status':cells.get('itsme_run_status'),
                    'partial_count':cells.get('itsme_partial_locus_count'), 'diagnostic_reason':cells.get('itsme_diagnostic_reason'), 'count': cells.get('itsme_complete_locus_count'),
                    'loci': loci, 'files': distinct_files(its_files, root), 'graphs': graphs['ITSME'],
                    'summary_file': relative_file(report, root),
                    'sequence': relative_file(its_sequence, root), 'cassette_sequence': relative_file(its_cassette, root)},
            'qc': {'reads_before': before.get('total_reads'), 'reads_after': after.get('total_reads'),
                   'q30_before': before.get('q30_rate'), 'q30_after': after.get('q30_rate'),
                   'report': relative_file(first([sample/'qc/fastp.html', sample/'ITSME/fastp.html'], root), root)},
            'stages': [dict({key: stage.get(key, '') for key in ['stage', 'status', 'reason', 'elapsed_seconds']},
                            log=stage_log(stage, sample, root)) for stage in stages],
            'logs': logs}
    checkpoint = sample/'run_config.json'
    if summary_path and checkpoint.is_file() and checkpoint.stat().st_mtime_ns > summary_path.stat().st_mtime_ns:
        result['status'] = 'unfinished'
        result['warnings'].append('The saved summary predates the latest sample attempt. Old recovery claims are hidden; retained files remain available for inspection.')
        result['mito'].update(length=None, accepted=False, sequence=None, circular=None, method='',
                              genes_recovered=None, genes_expected=None, target_status='unresolved',
                              complete_by_reference=None, annotation_ok=None, hits=[], reference_taxa=[], hit_count=0)
        result['mito']['core_genes'] = {'cds':None, 'rrna':None, 'extra':[], 'source':'Unavailable'}
        result['mito']['reference_genes'] = []
        result['mito']['annotation']['features'] = []
        result['its'].update(status='not_run', run_status='not_run', count=None, loci=[])
    return result


def build_data(root, title=None, bold_metadata=None, export_graphs=False):
    root = Path(root).absolute()
    samples, metadata_sources = [], []
    for sample in sorted(root.iterdir()):
        if sample.is_dir() and not sample.is_symlink() and ((sample/'final/summary.json').is_file() or (sample/'MitoLooper/summary.json').is_file()):
            samples.append(collect_sample(sample, root))
            config = read_json(sample/'run_config.json', samples[-1]['warnings'])
            metadata_sources.append(bold_metadata or config.get('options', {}).get('bold_metadata'))
    if not samples:
        raise ValueError('No saved Tortuga sample summaries found in this directory')
    resolve_coi_names(samples, metadata_sources)
    from graph_evidence import report_graphs
    for sample in samples:
        mito = sample['mito']
        confidence_evidence = dict(mito['bold_evidence'])
        if mito['hit_count'] > len(mito['hits']) and mito['hits']:
            high, low = mito['hits'][0]['bitscore'], mito['hits'][-1]['bitscore']
            fraction = confidence_evidence.get('alignment_thresholds', {}).get('near_best_fraction', .02)
            if high is None or low is None or low >= (1-fraction)*high:
                confidence_evidence['near_best_truncated'] = True
        mito['confidence'] = coi_support(mito['hits'], confidence_evidence)
        graphs = report_graphs(root/sample['id'], root, sample, export_graphs)
        mito['graphs'], sample['its']['graphs'] = graphs['MitoLooper'], graphs['ITSME']
        its_taxa = {locus['locus']: locus['taxon'] for locus in sample['its']['loci']}
        for graph in sample['its']['graphs']:
            for node in graph['nodes']:
                node['taxa'] = sorted({its_taxa.get(locus, '') for locus in node['loci']} - {''})
        for graph in mito['graphs']:
            for node in graph['nodes']:
                node['taxa'] = [mito['confidence']['taxon']] if mito['confidence']['taxon'] else mito['reference_taxa']
                node['taxonomy_scope'] = 'Assembly COI match' if mito['confidence']['taxon'] else 'Closest COI references'
        for branch in [mito, sample['its']]:
            branch['files'] = [item for item in branch['files'] if Path(item['name']).suffix.lower() not in {'.gfa', '.gfa2', '.fastg'}]
        mito['files'] = [item for item in mito['files'] if Path(item['name']).suffix.lower() not in {'.fasta', '.fa', '.fna'}
                         or item['name'] in {'accepted_candidate.fasta', 'best_sequence.fasta', 'circular_target_candidate.fasta'}
                         or item['name'].startswith('mitos_')]
    return {'title': title or root.name, 'built': datetime.now(timezone.utc).isoformat(timespec='seconds'),
            'samples': samples, 'files': [item for name in ['cohort_summary.csv', 'cohort_its_loci.csv', 'cohort_gene_recovery.csv', 'reference_selection.tsv']
                                        for item in [relative_file(root/name, root)] if item]}


def write_dashboard(root, title=None, bold_metadata=None):
    """Read saved results and replace only dashboard.html. Caller owns the run lock."""
    root = Path(root).absolute()
    target, temporary = root/'dashboard.html', root/'dashboard.html.tmp'
    if target.is_symlink() or temporary.is_symlink():
        raise ValueError('Refusing a linked dashboard output')
    data = build_data(root, title, bold_metadata, export_graphs=True)
    assets = Path(__file__).parent/'dashboard_assets'
    template = (assets/'dashboard.html').read_text()
    payload = json.dumps(data, ensure_ascii=True, allow_nan=False).replace('<', '\\u003c').replace('>', '\\u003e').replace('&', '\\u0026')
    html = template.replace('/*__STYLE__*/', (assets/'dashboard.css').read_text())
    html = html.replace('/*__SCRIPT__*/', (assets/'dashboard.js').read_text()).replace('__REPORT_DATA__', payload)
    temporary.write_text(html)
    temporary.replace(target)
    return target


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run_dir', type=Path, help='Completed Tortuga results directory')
    parser.add_argument('--title', help='Display title; defaults to the results directory name')
    parser.add_argument('--bold-metadata', type=Path, help='BOLD metadata TSV; defaults to the saved run configuration')
    args = parser.parse_args(argv)
    if args.bold_metadata is not None and not args.bold_metadata.is_file():
        parser.error('BOLD metadata file not found: '+str(args.bold_metadata))
    root = args.run_dir.absolute()
    if root.is_symlink() or not root.is_dir():
        parser.error('Not a regular results directory')
    with (root/'.pipeline.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            parser.error('The pipeline is running; retry after it stops. No dashboard was changed.')
        try:
            path = write_dashboard(root, args.title, args.bold_metadata)
        except (ValueError, OSError, csv.Error) as error:
            parser.error(str(error))
    print(f'Dashboard: {path}\nOpen it in a browser. Move the entire results folder to preserve download links.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

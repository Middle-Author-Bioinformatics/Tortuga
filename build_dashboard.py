#!/usr/bin/env python3
"""Build an offline Tortuga dashboard beside the saved cohort results."""
from __future__ import annotations

import argparse
from collections import defaultdict, deque
from datetime import datetime, timezone
import csv
import fcntl
import json
import math
from pathlib import Path
import re
from urllib.parse import quote, unquote

GRAPH_NODES = 240
GRAPH_EDGES = 720
GRAPH_BYTES = 100 * 1024**2
FEATURE_LIMIT = 3000


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


def annotations(path, length, root):
    result = {'file': relative_file(path, root), 'features': [], 'warnings': []}
    if path is None:
        return result
    with path.open(errors='replace') as handle:
        for line in handle:
            if line.startswith('##FASTA'):
                break
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                continue
            start, end = number(fields[3]), number(fields[4])
            if start is None or end is None or start < 1 or end < start:
                result['warnings'].append('Malformed feature coordinates were omitted.')
                continue
            if length and end > length:
                result['warnings'].append('Some feature coordinates exceed the displayed sequence length; those features were omitted.')
                continue
            attrs = {key: unquote(value) for item in fields[8].split(';') if '=' in item for key, value in [item.split('=', 1)]}
            result['features'].append({'sequence_id': fields[0], 'type': fields[2], 'start': int(start), 'end': int(end),
                'strand': fields[6], 'gene': attrs.get('Name') or attrs.get('gene') or attrs.get('gene_id') or attrs.get('ID', fields[2]),
                'product': attrs.get('product', ''), 'source': fields[1]})
            if len(result['features']) >= FEATURE_LIMIT:
                result['warnings'].append(f'Annotation preview is limited to {FEATURE_LIMIT} features; download the full GFF.')
                break
    if len({feature['sequence_id'] for feature in result['features']}) > 1:
        result['warnings'].append('Annotations contain multiple sequence IDs; the gene map is disabled to avoid mixing contigs.')
        result['multiple_sequences'] = True
    return result


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
    gene_data = annotations(annotation, length, root)
    cells, loci = summarize_itsme(sample, summary)
    report = first_file(sample, ['master_summary.csv', 'master_summary.tsv'])
    native_rows = {row.get('contig', row.get('locus', '')): row for row in read_table(report)}
    for locus in loci:
        native = native_rows.get(locus['locus'], {})
        locus['regions'] = [{'name': region, 'start': positions[0], 'end': positions[1]}
            for region in ['SSU', 'ITS1', '5.8S', 'ITS2', 'LSU']
            for positions in [interval(native.get(region+'_coordinates'), locus.get('locus_length_bp'))] if positions]
        locus.pop('source_report', None)
    graphs = {'ITSME': [], 'MitoLooper': []}
    for branch in graphs:
        paths = sorted(set((sample/branch).glob('*.gfa')) | set((sample/branch).glob('*.fastg')) | set((sample/branch).glob('*.gfa2')))
        # Compatibility with old compact files that have not yet been reorganized.
        if not paths:
            paths = [path for path in sorted((sample/'final').glob('*')) if path.suffix in {'.gfa', '.gfa2', '.fastg'}
                     and (path.name.startswith('itsme_') == (branch == 'ITSME'))]
        for path in paths:
            if not inside(path, root):
                continue
            try:
                graph = parse_graph(path, root)
                graph['scope'] = 'Whole-genome graph; not specific to the selected mitogenome' if path.name.startswith('wgs_') else 'Assembly-method graph; may include alternative or off-target paths'
                graphs[branch].append(graph)
            except (ValueError, OSError) as error:
                warnings.append(f'{path.name}: graph unavailable ({type(error).__name__})')
    best_candidate = summary.get('best_candidate') or {}
    evidence = {}
    candidate_id = best_candidate.get('candidate_id', '')
    if candidate_id and '/' not in candidate_id and '\\' not in candidate_id:
        candidate_dir = sample/'candidates'/candidate_id
        if inside(candidate_dir/'bold_evidence.json', root):
            evidence = read_json(candidate_dir/'bold_evidence.json', warnings)
    hits = []
    evidence_hits = sorted(evidence.get('hits', []),
                           key=lambda hit: number(hit.get('bitscore')) or 0, reverse=True)
    top_score = number(evidence_hits[0].get('bitscore')) if evidence_hits else None
    # Display observed reference names, never substitute the expected identity.
    # Without scores, older evidence can supply the table but no claimed top hit.
    reference_taxa = sorted({hit.get('metadata', {}).get('taxon', '') for hit in evidence_hits
                             if top_score is not None and number(hit.get('bitscore')) == top_score}
                            - {''})
    for hit in evidence_hits[:12]:
        meta = hit.get('metadata', {})
        hits.append({'taxon': meta.get('taxon', ''), 'lineage': meta.get('lineage', ''),
                     'identity': number(hit.get('pident')), 'aligned_bp': number(hit.get('length')),
                     'record': meta.get('processid', hit.get('sseqid', ''))})
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
    result = {'id': name, 'expected_taxon': summary.get('target_taxon', ''),
            'status': summary.get('status', 'unknown'), 'pipeline': summary.get('pipeline', 'both'),
            'summary_file': relative_file(summary_path, root), 'warnings': warnings,
            'mito': {'length': length, 'accepted': bool(accepted), 'sequence': relative_file(selected, root),
                'circular': summary.get('circular'), 'method': summary.get('winning_method') or summary.get('best_method', ''),
                'genes_recovered': summary.get('genes_recovered'), 'genes_expected': summary.get('genes_expected'),
                'complete_by_reference': summary.get('complete_by_reference'), 'target_status': summary.get('target_status', 'unresolved'),
                'annotation_ok': summary.get('annotation_ok'), 'annotation_status': summary.get('annotation_status',''),
                'annotation_scope':summary.get('annotation_scope',''), 'target_reason':summary.get('target_reason') or evidence.get('reason',''), 'annotation': gene_data, 'hits': hits,
                'reference_taxa': reference_taxa, 'hit_count': len(evidence_hits),
                'evidence_file': relative_file(sample/'candidates'/candidate_id/'bold_evidence.json', root) if candidate_id else None,
                'files': mito_files, 'graphs': graphs['MitoLooper']},
            'its': {'status': cells.get('itsme_reporting_status'), 'run_status':cells.get('itsme_run_status'),
                    'partial_count':cells.get('itsme_partial_locus_count'), 'diagnostic_reason':cells.get('itsme_diagnostic_reason'), 'count': cells.get('itsme_complete_locus_count'),
                    'loci': loci, 'files': its_files, 'graphs': graphs['ITSME'],
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
        result['mito']['annotation']['features'] = []
        result['its'].update(status='not_run', count=None, loci=[])
    return result


def build_data(root, title=None):
    root = Path(root).absolute()
    samples = []
    for sample in sorted(root.iterdir()):
        if sample.is_dir() and not sample.is_symlink() and ((sample/'final/summary.json').is_file() or (sample/'MitoLooper/summary.json').is_file()):
            samples.append(collect_sample(sample, root))
    if not samples:
        raise ValueError('No saved Tortuga sample summaries found in this directory')
    return {'title': title or root.name, 'built': datetime.now(timezone.utc).isoformat(timespec='seconds'),
            'samples': samples, 'files': [item for name in ['cohort_summary.csv', 'cohort_its_loci.csv', 'cohort_gene_recovery.csv', 'reference_selection.tsv']
                                        for item in [relative_file(root/name, root)] if item]}


def write_dashboard(root, title=None):
    """Read saved results and replace only dashboard.html. Caller owns the run lock."""
    root = Path(root).absolute()
    target, temporary = root/'dashboard.html', root/'dashboard.html.tmp'
    if target.is_symlink() or temporary.is_symlink():
        raise ValueError('Refusing a linked dashboard output')
    data = build_data(root, title)
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
    args = parser.parse_args(argv)
    root = args.run_dir.absolute()
    if root.is_symlink() or not root.is_dir():
        parser.error('Not a regular results directory')
    with (root/'.pipeline.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            parser.error('The pipeline is running; retry after it stops. No dashboard was changed.')
        try:
            path = write_dashboard(root, args.title)
        except (ValueError, OSError, csv.Error) as error:
            parser.error(str(error))
    print(f'Dashboard: {path}\nOpen it in a browser. Move the entire results folder to preserve download links.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

"""Link final sequences to their source graphs without assembly or classification."""
from __future__ import annotations

from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
import re
from urllib.parse import unquote

from report_evidence import fasta_records, reverse_complement


def table(path):
    if not path or not path.is_file():
        return []
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter=',' if path.suffix == '.csv' else '\t'))


def path_tokens(value):
    pieces = []
    for piece in value.strip().split(';'):
        tokens = [token.strip() for token in piece.split(',') if token.strip()]
        if not tokens or any(not re.fullmatch(r'[^\s,;]+[+-]', token) for token in tokens):
            return []
        pieces.append(tokens)
    return pieces


def contig_paths(path):
    result, current = {}, None
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith(('NODE_', 'EDGE_')):
                current = line
                result[current] = []
            elif current and line:
                result[current].extend(path_tokens(line.rstrip(';')))
    return result


def local_file(value, sample, root):
    if not value:
        return None
    path = Path(unquote(str(value)))
    candidates = [path, sample/path, root/path, sample/'ITSME'/path.name, sample/'MitoLooper'/path.name]
    for candidate in candidates:
        try:
            candidate.resolve().relative_to(sample.resolve())
            if candidate.is_file():
                return candidate
        except (ValueError, OSError):
            pass
    return None


def inventory(sample, root, branch):
    """Use original manifest paths to distinguish numeric node IDs in different bins."""
    entries = {}
    for row in table(sample/'final/outputs_manifest.tsv'):
        path = local_file(row.get('retained_file'), sample, root)
        if path and path.suffix.lower() in {'.gfa', '.gfa2', '.fastg', '.paths'}:
            if branch in path.parts or (branch == 'ITSME') == path.name.startswith('itsme_'):
                entries[path.resolve()] = {'path': path, 'source': row.get('source') or str(path)}
    for path in sorted((sample/branch).glob('*')):
        if path.is_file() and path.suffix.lower() in {'.gfa', '.gfa2', '.fastg', '.paths'}:
            entries.setdefault(path.resolve(), {'path': path, 'source': str(path)})
    # Before compact cleanup, preserve final pooled and per-bin assembly provenance.
    for stage in ['02_itsme', 'itsme'] if branch == 'ITSME' else ['011_novoplasty', '012_getorganelle', '013_mitofinder', '014_mtgrasp', '041_getorganelle_assembly', '042_mitofinder_assembly', '03_wgs']:
        for path in sorted((sample/stage).rglob('*')):
            if not path.is_file() or path.suffix.lower() not in {'.gfa', '.gfa2', '.fastg', '.paths'}:
                continue
            if any(re.fullmatch(r'K\d+', part) or part in {'tmp', 'temp'} for part in path.parts):
                continue
            if local_file(str(path), sample, root):
                entries.setdefault(path.resolve(), {'path': path, 'source': str(path)})
    return list(entries.values())


def source_graph_files(sample, root, branch, method=''):
    """Download inventory independent of whether a selected path can be mapped.

    Include sources attributed to the winning method and unattributed retained
    graphs. Never imply that an unmapped source reconstructs the final sequence.
    """
    from build_dashboard import relative_file
    files = []
    methods = ('getorganelle', 'novoplasty', 'mitofinder', 'mtgrasp', 'wgs')
    wanted = re.sub(r'[^a-z0-9]', '', method.lower())
    for entry in inventory(sample, root, branch):
        if entry['path'].suffix.lower() not in {'.gfa', '.gfa2', '.fastg'}:
            continue
        source = re.sub(r'[^a-z0-9]', '', (entry['source']+' '+entry['path'].name).lower())
        attributed = next((name for name in methods if name in source), '')
        if branch == 'MitoLooper' and wanted and attributed and attributed not in wanted:
            if not (wanted.endswith('assembly') and attributed == 'wgs'):
                continue
        item = relative_file(entry['path'], root)
        if item:
            files.append(item)
    return files


class SourceGraph:
    def __init__(self, path):
        self.path, self.nodes, self.edges = path, {}, []
        self.named_paths = {}
        self.too_large = path.stat().st_size > 100*1024**2
        if self.too_large:
            return
        def node(label):
            label = label.rstrip("';")
            match = re.match(r'(?:EDGE|NODE)_(.+?)_length_(\d+)_cov_([\d.eE+-]+)', label)
            return (match.group(1), int(match.group(2)), float(match.group(3))) if match else (label, 0, None)
        with path.open(errors='replace') as handle:
            header, sequence = None, []
            def save():
                if header is None:
                    return
                source, _, neighbors = header.rstrip(';').partition(':')
                key, length, coverage = node(source)
                strand = '-' if source.endswith("'") else '+'
                seq = ''.join(sequence).upper()
                if strand == '-':
                    seq = reverse_complement(seq)
                self.nodes.setdefault(key, {'id': key, 'length': length or len(seq), 'coverage': coverage, 'sequence': seq})
                for target in filter(None, neighbors.split(',')):
                    self.edges.append({'source': key, 'source_strand': strand, 'target': node(target)[0],
                                       'target_strand': '-' if target.endswith("'") else '+', 'overlap': None})
            for line in handle:
                if line.startswith('>'):
                    save(); header, sequence = line[1:].strip(), []
                elif header is not None:
                    sequence.append(line.strip())
                else:
                    f = line.rstrip('\n').split('\t')
                    if f[0] == 'H' and any(x.startswith('VN:Z:2') for x in f):
                        break
                    if f[0] == 'S' and len(f) >= 3:
                        tags = {x.split(':', 2)[0]: x.split(':', 2)[-1] for x in f[3:] if ':' in x}
                        seq = f[2] if f[2] != '*' else ''
                        self.nodes[f[1]] = {'id': f[1], 'length': len(seq) if seq else int(tags.get('LN', 0)),
                            'coverage': float(tags.get('DP', tags.get('dp', 0))) or None, 'sequence': seq.upper()}
                    elif f[0] == 'L' and len(f) >= 6:
                        overlap = int(f[5][:-1]) if re.fullmatch(r'\d+M', f[5]) else None
                        self.edges.append({'source': f[1], 'source_strand': f[2], 'target': f[3], 'target_strand': f[4], 'overlap': overlap})
                    elif f[0] == 'P' and len(f) >= 3:
                        self.named_paths[f[1]] = path_tokens(f[2])
            save()

    def path_evidence(self, pieces, final_sequence):
        edges = {}
        for edge in self.edges:
            a, b = edge['source']+edge['source_strand'], edge['target']+edge['target_strand']
            edges[a, b] = edge['overlap']
            flip = lambda x: x[:-1]+('+' if x[-1] == '-' else '-')
            edges[flip(b), flip(a)] = edge['overlap']
        steps, reconstructed, total, known, junctions = [], '', 0, len(pieces) == 1, []
        for piece in pieces:
            for i, token in enumerate(piece):
                record = self.nodes[token[:-1]]
                overlap = 0 if i == 0 else edges.get((piece[i-1], token))
                if overlap is None:
                    known = False
                seq = record['sequence']
                if token[-1] == '-':
                    seq = reverse_complement(seq)
                if not seq or len(seq) != record['length']:
                    known = False
                added = record['length']-overlap if overlap is not None else None
                steps.append({'node': token[:-1], 'orientation': token[-1], 'length': record['length'],
                              'overlap': overlap, 'added_bp': added, 'start': total+1 if known else None})
                if i:
                    junctions.append([piece[i-1], token])
                if overlap is not None and seq:
                    if overlap and reconstructed and reconstructed[-overlap:] != seq[:overlap]:
                        known = False
                    reconstructed += seq[overlap:]
                total += added or 0
        forward = reconstructed.find(final_sequence) if known and final_sequence else -1
        reverse = reverse_complement(reconstructed).find(final_sequence) if known and final_sequence else -1
        # FASTG lacks overlap lengths. Infer a common SPAdes k overlap only when
        # recorded links exist and it reconstructs the entire final sequence exactly.
        if not known and len(pieces) == 1 and len(steps) > 1 and final_sequence:
            difference = sum(step['length'] for step in steps)-len(final_sequence)
            count = len(steps)-1
            if difference >= 0 and difference % count == 0:
                overlap = difference//count
                linked = all((a, b) in edges for a, b in zip(pieces[0], pieces[0][1:]))
                sequences = [self.nodes[token[:-1]]['sequence'] if token[-1] == '+' else
                             reverse_complement(self.nodes[token[:-1]]['sequence']) for token in pieces[0]]
                valid = linked and all(len(nt) == step['length'] and overlap < len(nt)
                                       for nt, step in zip(sequences, steps))
                valid = valid and all(not overlap or a[-overlap:] == b[:overlap] for a, b in zip(sequences, sequences[1:]))
                candidate = sequences[0]+''.join(nt[overlap:] for nt in sequences[1:]) if valid else ''
                if candidate and final_sequence in {candidate, reverse_complement(candidate)}:
                    known, total, reconstructed = True, len(candidate), candidate
                    forward, reverse = (0, -1) if candidate == final_sequence else (-1, 0)
                    offset = 0
                    for i, step in enumerate(steps):
                        step.update(overlap=overlap if i else 0, added_bp=step['length']-(overlap if i else 0), start=offset+1)
                        offset += step['added_bp']
        verified = forward >= 0 or reverse >= 0
        return {'steps': steps, 'assembled_bp': total if known else None, 'sequence_verified': verified,
                'locus_bp': len(final_sequence) or None,
                'trimmed_bp': total-len(final_sequence) if verified else None,
                'reverse_to_final': forward < 0 and reverse >= 0, 'junctions': junctions}


def source_bin(source):
    match = re.search(r'ITSME_(BIN_[A-Za-z0-9]+)_', source, re.I)
    return match.group(1).lower() if match else 'pooled'


def report_graphs(sample, root, data, write=False):
    """Only graphs with explicit path provenance or an exact selected-sequence placement."""
    from build_dashboard import first, relative_file, layout_graph
    from itsme_reporting import first_file
    results = {'ITSME': [], 'MitoLooper': []}
    mapping = table(first_file(sample, ['locus_nodes.tsv', 'locus_source_map.tsv']))
    validation = {row.get('path_id'): row for row in table(first_file(sample, ['graph_locus_validation.tsv']))}
    sequences = fasta_records(first_file(sample, ['reported_loci.fasta', 'complete_rDNA_loci.fasta']))
    mito_sequence = next(iter(fasta_records(local_file(data['mito']['sequence']['href'], sample, root)).values()), '') if data['mito']['sequence'] else ''
    summary_path = first([sample/'final/summary.json', sample/'MitoLooper/summary.json'], root)
    candidate = (json.loads(summary_path.read_text()).get('best_candidate') or {}) if summary_path else {}
    candidate_source = candidate.get('source_path', '')
    for record in table(sample/'final/outputs_manifest.tsv'):
        if candidate_source and Path(record.get('retained_file', '')).name == Path(candidate_source).name:
            candidate_source = record.get('source') or candidate_source
            break
    for branch in results:
        if (sample/branch).is_symlink() or (sample/branch/'selected_graphs').is_symlink():
            raise ValueError('Refusing linked graph result directory')
        signature = hashlib.sha256(json.dumps(
            [[(locus['locus'], locus.get('sequence', '')) for locus in data['its']['loci']], mapping]
            if branch == 'ITSME' else [mito_sequence, candidate.get('candidate_id'), candidate_source]).encode()).hexdigest()
        cache = sample/branch/'selected_graphs/index.json'
        if cache.is_file() and not cache.is_symlink():
            cached = json.loads(cache.read_text())
            if cached.get('signature') == signature:
                cached_graphs = cached.get('graphs', [])
                if all(local_file(graph['file']['href'], sample, root) for graph in cached_graphs if graph.get('file')):
                    results[branch] = cached_graphs
                    for target in data['its']['loci'] if branch == 'ITSME' else []:
                        if target['locus'] in cached.get('selections', {}):
                            target.update(cached['selections'][target['locus']])
                    continue
        entries = inventory(sample, root, branch)
        path_tables = [entry for entry in entries if entry['path'].suffix == '.paths']
        graphs = [entry for entry in entries if entry['path'].suffix in {'.gfa', '.fastg'}]
        graphs.sort(key=lambda e: (e['path'].suffix != '.gfa', 'with_scaffolds' not in e['path'].name, e['source']))
        if branch == 'ITSME':
            targets = data['its']['loci']
        else:
            targets = [{'locus': 'mitogenome', 'taxon': '', 'sequence': mito_sequence}] if mito_sequence else []
            method = data['mito']['method']
            graphs = [e for e in graphs if method and (method in e['source'].lower() or (method.endswith('_assembly') and '03_wgs' in e['source']))]
            scoped = [e for e in graphs if candidate_source and e['source'].startswith(str(Path(candidate_source).parent)+'/')]
            if scoped:
                graphs = scoped
        grouped, used = {}, set()
        for target in targets:
            locus = target['locus']
            sequence = target.get('sequence') or sequences.get(locus, '')
            rows = [row for row in mapping if (row.get('locus') or row.get('sequence_id')) == locus
                    and row.get('path_status') != 'COLLAPSED_SAME_TAXONOMY']
            if len(rows) > 1:
                continue
            row = rows[0] if rows else {}
            source = row.get('source_path') or row.get('source') or locus
            tokens = path_tokens(row.get('node_path', ''))
            selected = None
            for entry in graphs:
                scope = source_bin(source)
                text = entry['source'].lower()
                if branch == 'ITSME':
                    if scope != 'pooled' and scope not in text:
                        continue
                    if scope == 'pooled' and 'taxon_bins/' in text:
                        continue
                graph = entry.get('parsed')
                if graph is None:
                    try:
                        graph = entry['parsed'] = SourceGraph(entry['path'])
                    except (OSError, ValueError, KeyError):
                        continue
                if not graph.nodes:
                    continue
                pieces = tokens
                if branch == 'ITSME' and not pieces:
                    original = re.sub(r'^ITSME_BIN_[A-Za-z0-9]+_', '', source, flags=re.I)
                    pieces = graph.named_paths.get(original, [])
                    for path_entry in path_tables:
                        if pieces:
                            break
                        if str(Path(path_entry['source']).parent) == str(Path(entry['source']).parent):
                            pieces = contig_paths(path_entry['path']).get(original, [])
                if branch == 'MitoLooper' or not pieces:
                    placements = []
                    search = sequence+sequence if branch == 'MitoLooper' and data['mito']['circular'] else sequence
                    for key, node in graph.nodes.items():
                        nt = node['sequence']
                        if len(nt) < 31 or len(nt) > len(sequence) or set(nt) - set('ACGT'):
                            continue
                        for strand, value in [('+', nt), ('-', reverse_complement(nt))]:
                            pos = search.find(value)
                            next_pos = search.find(value, pos+1)
                            if 0 <= pos < len(sequence) and (next_pos < 0 or next_pos >= len(sequence)):
                                placements.append((pos, key+strand)); break
                    pieces = [[token for pos, token in sorted(placements)]] if placements else []
                if not pieces or any(token[:-1] not in graph.nodes for piece in pieces for token in piece):
                    continue
                evidence = graph.path_evidence(pieces, sequence)
                if evidence['reverse_to_final'] and len(pieces) == 1:
                    pieces = [[token[:-1]+('+' if token[-1] == '-' else '-') for token in reversed(pieces[0])]]
                    evidence = graph.path_evidence(pieces, sequence)
                # Numeric node IDs recur across bins/rounds: never pick an arbitrary
                # graph when both its source and sequence reconstruction disagree.
                explicit_source = row.get('source_graph')
                reduced_file = row.get('graph_file')
                if explicit_source and str(Path(explicit_source)) != entry['source'] and not (
                        reduced_file and Path(entry['source']).name == reduced_file):
                    continue
                if sequence and evidence['assembled_bp'] is not None and not evidence['sequence_verified']:
                    continue
                selected = entry, graph, pieces, evidence
                break
            if selected is None:
                target['graph_status'] = 'unavailable'
                continue
            entry, graph, pieces, evidence = selected
            key = str(entry['path'].resolve())
            grouped.setdefault(key, {'entry': entry, 'source': graph, 'paths': {}, 'taxa': {}, 'label': source_bin(source)})
            grouped[key]['paths'][locus] = evidence
            grouped[key]['taxa'][locus] = target.get('taxon', '')
            decision = validation.get(source, {})
            target['graph_status'] = 'mapped'
            target['selection'] = {'source': source, 'supported_junctions': decision.get('supported_junctions'),
                'junction_count': decision.get('junction_count'), 'status': decision.get('path_status') or row.get('path_status', ''),
                'alternative_count': sum(1 for v in validation.values() if v.get('promoted_locus') == locus and v.get('path_status') == 'COLLAPSED_SAME_TAXONOMY')}
        for item in grouped.values():
            graph, paths = item['source'], item['paths']
            keep = {step['node'] for path in paths.values() for step in path['steps']}
            tag = hashlib.sha256((item['entry']['source']+','.join(sorted(paths))).encode()).hexdigest()[:12]
            output = sample/branch/'selected_graphs'/f'{tag}.gfa'
            edges, seen = [], set()
            for edge in graph.edges:
                key = (edge['source'], edge['source_strand'], edge['target'], edge['target_strand'])
                if edge['source'] in keep and edge['target'] in keep and key not in seen:
                    seen.add(key); edges.append(edge)
            if write:
                if output.is_symlink() or output.parent.is_symlink():
                    raise ValueError('Refusing linked graph output')
                output.parent.mkdir(exist_ok=True, parents=True)
                lines = ['H\tVN:Z:1.0']
                for key in sorted(keep):
                    node = graph.nodes[key]
                    coverage = f"\tDP:f:{node['coverage']}" if node['coverage'] is not None else ''
                    lines.append(f"S\t{key}\t{node['sequence'] or '*'}\tLN:i:{node['length']}"+coverage)
                for edge in edges:
                    overlap = str(edge['overlap'])+'M' if edge['overlap'] is not None else '*'
                    lines.append(f"L\t{edge['source']}\t{edge['source_strand']}\t{edge['target']}\t{edge['target_strand']}\t{overlap}")
                temp = output.with_suffix('.tmp')
                if temp.is_symlink():
                    raise ValueError('Refusing linked graph temporary output')
                temp.write_text('\n'.join(lines)+'\n'); temp.replace(output)
            nodes = []
            for key in sorted(keep):
                loci = [locus for locus, path in paths.items() if any(step['node'] == key for step in path['steps'])]
                node = dict(graph.nodes[key], loci=loci, taxa=sorted({item['taxa'][locus] for locus in loci if item['taxa'][locus]}),
                            taxonomy_scope='Locus consensus' if branch == 'ITSME' else 'Assembly COI match')
                nodes.append(node)
            layout_graph(nodes, edges)
            results[branch].append({'file': relative_file(output, root), 'source_file': relative_file(item['entry']['path'], root),
                'label': item['label'] if branch == 'ITSME' else data['mito']['method'],
                'nodes': nodes, 'edges': edges, 'paths': paths, 'loci': list(paths), 'warnings': [],
                'total_nodes': len(nodes), 'total_links': len(edges), 'truncated': False})
        if write and results[branch]:
            if cache.is_symlink() or cache.with_suffix('.tmp').is_symlink():
                raise ValueError('Refusing linked graph provenance output')
            content = {'signature': signature, 'graphs': results[branch],
                       'selections': {target['locus']: {key: target[key] for key in ['graph_status', 'selection'] if key in target}
                                      for target in targets}}
            temp = cache.with_suffix('.tmp')
            temp.write_text(json.dumps(content)); temp.replace(cache)
    return results

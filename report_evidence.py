"""Sequence and annotation evidence for portable result reports."""
from __future__ import annotations

from collections import defaultdict
import hashlib
from pathlib import Path
import re
from urllib.parse import unquote
from bold_taxonomy import normalize_taxonomy

CORE_CDS = {'atp6', 'atp8', 'cox1', 'cox2', 'cox3', 'cytb',
            'nad1', 'nad2', 'nad3', 'nad4', 'nad4l', 'nad5', 'nad6'}
CORE_RRNA = {'rrns', 'rrnl'}
COMPLEMENT = str.maketrans('ACGTRYMKBDHVSWNacgtrymkbdhvswn',
                          'TGCAYRKMVHDBSWNtgcayrkmvhdbswn')


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def fasta_records(path, limit=100*1024**2):
    """Bounded reader for final sequences; never load recruited reads or WGS FASTAs."""
    if path is None or not Path(path).is_file() or Path(path).stat().st_size > limit:
        return {}
    records, name, parts = {}, None, []
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith('>'):
                if name:
                    records[name] = ''.join(parts).upper()
                name, parts = line[1:].strip().split()[0], []
            elif name:
                parts.append(''.join(line.split()))
    if name:
        records[name] = ''.join(parts).upper()
    return records


def gene_name(value):
    from mito_pipeline import normalize_gene
    value = re.sub(r'^(gene|cds|rna)[:_-]', '', value, flags=re.I)
    value = re.sub(r'(?:[_.-](?:\d+|[ab]))+$', '', value, flags=re.I)
    return normalize_gene(value)


def gff_attributes(text):
    result = {}
    for part in text.strip().split(';'):
        if '=' in part:
            key, value = part.split('=', 1)
        else:
            match = re.match(r'\s*(\S+)\s+"(.*)"\s*$', part)
            if not match:
                continue
            key, value = match.groups()
        result[key.strip()] = unquote(value.strip())
    return result


def read_annotation(path, length, circular=False, sequence='', limit=3000):
    result = {'features': [], 'warnings': [], 'multiple_sequences': False}
    if path is None:
        return result
    rows, names = [], {}
    with Path(path).open(errors='replace') as handle:
        for line in handle:
            if line.startswith('##FASTA'):
                break
            fields = line.rstrip('\n').split('\t')
            if line.startswith('#') or len(fields) != 9:
                continue
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError:
                continue
            attrs = gff_attributes(fields[8])
            name = attrs.get('gene') or attrs.get('Name') or attrs.get('gene_id') or attrs.get('ID', '')
            if attrs.get('ID'):
                names[attrs['ID']] = name
            rows.append((fields, attrs, name, start, end))
            if len(rows) >= limit:
                break
    result['multiple_sequences'] = len({row[0][0] for row in rows}) > 1
    def resolved_name(fields, attrs, name):
        direct = attrs.get('gene') or attrs.get('Name') or attrs.get('gene_id')
        return direct or names.get(attrs.get('Parent', '').split(',')[0]) or name or fields[2]
    explicit = {gene_name(resolved_name(fields, attrs, name)) for fields, attrs, name, start, end in rows
                if fields[2].lower() in {'cds', 'rrna', 'trna'}}
    seen = set()
    for fields, attrs, name, start, end in rows:
        name = resolved_name(fields, attrs, name)
        canonical = gene_name(name)
        native_type = fields[2]
        kind = native_type
        if native_type.lower() in {'gene', 'mrna', 'exon'}:
            if canonical in explicit:
                continue
            if canonical in CORE_CDS | {'atp9', 'muts', 'tatc'} or attrs.get('gene_biotype') == 'protein_coding' or attrs.get('biotype') == 'protein_coding':
                kind = 'CDS'
            elif canonical in CORE_RRNA:
                kind = 'rRNA'
            elif canonical.startswith(('trn', 'trna')):
                kind = 'tRNA'
            elif native_type.lower() != 'gene':
                continue
        segments = []
        if start >= 1 and end >= start and (not length or end <= length):
            segments = [[start, end]]
        elif circular and length and 1 <= start <= length and end > length and end-start+1 <= length:
            segments = [[start, length], [1, end-length]]
        elif circular and length and 1 <= end < start <= length:
            segments = [[start, length], [1, end]]
        if not segments:
            result['warnings'].append('Some feature coordinates exceed the displayed sequence length; those features were omitted.')
            continue
        key = (fields[0], kind, canonical, start, end, fields[6])
        if key in seen:
            continue
        seen.add(key)
        nt = ''.join(sequence[a-1:b] for a, b in segments) if sequence and not result['multiple_sequences'] else ''
        if fields[6] == '-':
            nt = reverse_complement(nt)
        result['features'].append({'sequence_id': fields[0], 'type': kind, 'native_type': native_type,
            'start': start, 'end': end, 'segments': segments, 'strand': fields[6], 'gene': name,
            'canonical_gene': canonical, 'sequence': nt,
            'feature_id': attrs.get('ID') or attrs.get('Parent', ''),
            'product': attrs.get('product', ''), 'source': fields[1]})
    result['warnings'] = list(dict.fromkeys(result['warnings']))
    combined, index = [], {}
    for feature in result['features']:
        key = (feature['sequence_id'], feature['feature_id'], feature['type'], feature['gene'], feature['strand'])
        if feature['type'] == 'CDS' and feature['feature_id'] and key in index:
            previous = index[key]
            previous['segments'].extend(feature['segments'])
            previous['segments'].sort()
            previous['start'], previous['end'] = previous['segments'][0][0], previous['segments'][-1][1]
            nt = ''.join(sequence[a-1:b] for a, b in previous['segments']) if sequence and not result['multiple_sequences'] else ''
            previous['sequence'] = reverse_complement(nt) if previous['strand'] == '-' else nt
        else:
            combined.append(feature); index[key] = feature
    result['features'] = combined
    return result


def core_gene_counts(candidate, gene_rows=()):
    genes = {gene_name(name) for name in candidate.get('genes_found', [])}
    rrnas = {gene_name(name) for name in candidate.get('rrnas_found', [])}
    available = candidate.get('gene_screen_available') is not False and ('genes_found' in candidate or bool(gene_rows))
    for row in gene_rows:
        if str(row.get('screen_pass', '')).lower() == 'true':
            (rrnas if row.get('type') == 'rRNA' else genes).add(gene_name(row['gene']))
    return {'cds': len(genes & CORE_CDS) if available else None,
            'rrna': len(rrnas & CORE_RRNA) if available else None,
            'extra': sorted((genes | rrnas) - CORE_CDS - CORE_RRNA)}


def distinct_files(files, root):
    """Deduplicate report downloads by sequence content, retaining native files on disk."""
    priority = {'accepted_candidate.fasta': 0, 'best_sequence.fasta': 1,
                'reported_loci.fasta': 0, 'complete_ITS.fasta': 1}
    output, seen = [], set()
    from urllib.parse import unquote
    for item in sorted(files, key=lambda f: (priority.get(f['name'], 5), f['name'])):
        path = Path(root)/unquote(item['href'])
        signature = None
        if path.suffix.lower() in {'.fa', '.fna', '.fasta', '.ffn', '.faa'}:
            sequences = fasta_records(path)
            if sequences:
                signature = ('protein' if path.suffix.lower() == '.faa' else 'nucleotide',
                             tuple(sorted(hashlib.sha256(seq.encode()).hexdigest() for seq in sequences.values())))
        if signature and signature in seen:
            continue
        if signature:
            seen.add(signature)
        output.append(item)
    return output


def coi_support(hits, evidence):
    """An explicit evidence grade, never a calibrated species probability."""
    normalized = []
    for hit in hits:
        taxon, lineage = normalize_taxonomy(hit.get('taxon'), hit.get('lineage'))
        normalized.append(dict(hit, taxon=taxon, lineage=lineage))
    hits = normalized
    if not hits:
        return {'grade': 'Unavailable', 'taxon': '', 'agreement': None}
    scored = [hit for hit in hits if hit.get('bitscore') is not None]
    top = max((hit['bitscore'] for hit in scored), default=None)
    fraction = evidence.get('alignment_thresholds', {}).get('near_best_fraction', .02)
    near = [hit for hit in hits if top is None or (hit.get('bitscore') or 0) >= (1-fraction)*top]
    subjects = {hit.get('subject') or hit.get('record') or str(i): hit for i, hit in enumerate(near)}
    near = list(subjects.values())
    known = [hit for hit in near if hit.get('taxon')]
    names = {hit['taxon'] for hit in known}
    common = ''
    if len(known) == len(near) and len(names) == 1:
        common = next(iter(names))
    elif len(known) == len(near) and all(hit.get('lineage') for hit in near):
        lineages = [[x.strip() for x in hit['lineage'].split(';') if x.strip()] for hit in near]
        for level in zip(*lineages):
            if len({x.casefold() for x in level}) != 1:
                break
            common = level[0]
    best = max(hits, key=lambda h: h.get('bitscore') or 0)
    identity = best.get('identity') or 0
    length = best.get('aligned_bp') or 0
    coverage = best.get('coverage')
    incomplete = (len(known) != len(near) or top is None or coverage is None or evidence.get('near_best_truncated')
                  or evidence.get('reason_code') in {'near_best_truncated', 'missing_metadata', 'search_failed'})
    grade = ('Limited' if incomplete or not common else
             'High' if identity >= 99 and length >= 400 and coverage is not None and coverage >= .9 else
             'Moderate' if identity >= 97 and length >= 150 and coverage is not None and coverage >= .8 else 'Low')
    if len(names) > 1:
        grade = 'Mixed'
    return {'grade': grade, 'taxon': '' if incomplete else common, 'agreement': f'{len(known)}/{len(near)}',
            'identity': best.get('identity'), 'aligned_bp': best.get('aligned_bp'), 'coverage': coverage}

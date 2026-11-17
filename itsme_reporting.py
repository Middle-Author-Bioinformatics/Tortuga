#!/usr/bin/env python3
"""Summarize reported ITSME loci without rerunning analysis or cleaning files."""
from __future__ import annotations

import argparse
import csv
import fcntl
import json
from pathlib import Path
import re

RANKS = ('domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
ITS_FIELDS = ['itsme_reporting_status', 'itsme_complete_locus_count', 'itsme_locus_ids',
              'itsme_locus_lengths_bp', 'itsme_its_lengths_bp', 'itsme_taxa', 'itsme_taxonomy_ranks']
LOCUS_FIELDS = ['sample', 'locus', 'locus_type', 'locus_length_bp', 'its_length_bp',
                'its1_length_bp', 'five_8s_length_bp', 'its2_length_bp', 'taxon',
                'taxonomy_rank', 'consensus_taxonomy', 'taxonomy_status', 'source_report']
MISSING = {'', 'na', 'n/a', 'none', 'unclassified', 'unresolved', 'unknown'}


def present(value):
    return str(value or '').strip().casefold() not in MISSING


def read_rows(path):
    with path.open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t' if path.suffix == '.tsv' else ',')
        return list(reader), list(reader.fieldnames or [])


def first_file(root, names):
    # Prefer current compact outputs; legacy and native work directories are read-only fallbacks.
    for directory, prefix in [(root/'ITSME', ''), (root/'final', 'itsme_'),
                              (root/'02_itsme/output', ''), (root/'itsme/output', '')]:
        for name in names:
            path = directory/(prefix+name)
            if path.is_file():
                return path
            path = directory/'final'/name
            if not prefix and path.is_file():
                return path
    return None


def fasta_lengths(path):
    lengths = {}
    if path is None:
        return lengths
    current = None
    with path.open() as handle:
        for line in handle:
            if line.startswith('>'):
                current = line[1:].split()[0]
                lengths[current] = 0
            elif current:
                lengths[current] += len(''.join(line.split()))
    return lengths


def integer(value):
    try:
        number = int(value)
        return number if number >= 0 else None
    except (ValueError, TypeError):
        return None


def interval(value, limit=None):
    match = re.fullmatch(r'\s*(\d+)\s*-\s*(\d+)\s*', str(value or ''))
    if not match:
        return None
    start, end = sorted(map(int, match.groups()))
    return (start, end) if start >= 1 and (limit is None or end <= limit) else None


def summarize_itsme(root, summary):
    """Return cohort cells and one row per complete locus; never infer taxonomy from the target."""
    root = Path(root)
    cells = {field: [] for field in ITS_FIELDS[2:]}
    cells.update(itsme_complete_locus_count=None, itsme_reporting_status='missing_report')
    if summary.get('pipeline') == 'mitolooper':
        cells['itsme_reporting_status'] = 'not_requested'
        return cells, []
    report = first_file(root, ['master_summary.csv', 'master_summary.tsv'])
    if report is None:
        if summary.get('itsme_ok') is False:
            cells['itsme_reporting_status'] = 'failed_no_report'
        elif summary.get('itsme_ok') is None:
            cells['itsme_reporting_status'] = 'not_run'
        return cells, []
    rows, columns = read_rows(report)
    if not ({'contig', 'locus'} & set(columns)):
        cells['itsme_reporting_status'] = 'unsupported_report'
        return cells, []
    taxonomy_path = first_file(root, ['taxonomy_validation.tsv'])
    taxonomy = read_rows(taxonomy_path)[0] if taxonomy_path else []
    by_id = {row.get('contig', ''): row for row in taxonomy}
    mapping_path = first_file(root, ['locus_source_map.tsv'])
    mapping = read_rows(mapping_path)[0] if mapping_path else []
    its_lengths = fasta_lengths(first_file(root, ['complete_ITS.fasta']))
    locus_lengths = fasta_lengths(first_file(root, ['reported_loci.fasta', 'complete_rDNA_loci.fasta']))
    results = []
    for row in rows:
        locus = row.get('contig') or row.get('locus') or ''
        if not locus:
            continue
        kind = row.get('locus_type', '')
        status = row.get('taxonomy_status', '')
        # Older summaries may contain incomplete rows. Do not promote them to complete loci.
        if 'partial' in kind.casefold() or 'chimera' in status.casefold() or 'partial' in row.get('status', '').casefold():
            continue
        lineage = row.get('consensus_taxonomy', '')
        tax = by_id.get(locus, {})
        if not tax:
            sources = [by_id.get(item.get('source_path'), {}) for item in mapping if item.get('locus') == locus]
            tax = next((item for item in sources if item and item.get('consensus_taxonomy') == lineage), {})
        if not present(lineage):
            lineage = tax.get('consensus_taxonomy', '')
        rank, name = '', ''
        for candidate in reversed(RANKS):
            value = row.get('consensus_'+candidate) or tax.get('consensus_'+candidate)
            if present(value):
                rank, name = candidate, value.strip()
                break
        if not name and present(lineage):
            name = next((item.strip() for item in reversed(lineage.split(';')) if present(item)), '')
        length = integer(row.get('length_bp'))
        if length is None:
            length = locus_lengths.get(locus)
        regions = [interval(row.get(field), length) for field in
                   ['ITS1_coordinates', '5.8S_coordinates', 'ITS2_coordinates']]
        its_length = None
        if all(regions):
            ordered = sorted(regions)
            if all(left[1] < right[0] for left, right in zip(ordered, ordered[1:])):
                its_length = ordered[-1][1] - ordered[0][0] + 1
        if its_length is None:
            its_length = its_lengths.get(locus)
        region_lengths = [region[1]-region[0]+1 if region else None for region in regions]
        results.append(dict(zip(LOCUS_FIELDS, [summary.get('sample', root.name), locus, kind, length, its_length,
            *region_lengths,
            name, rank, lineage if present(lineage) else '', status or tax.get('taxonomy_status', ''), str(report)])))
    cells.update(itsme_complete_locus_count=len(results),
                 itsme_reporting_status=('available_from_failed_run' if summary.get('itsme_ok') is False else
                                         'available' if results else 'no_complete_loci'))
    for field, source in [('itsme_locus_ids', 'locus'), ('itsme_locus_lengths_bp', 'locus_length_bp'),
                          ('itsme_its_lengths_bp', 'its_length_bp'), ('itsme_taxa', 'taxon'),
                          ('itsme_taxonomy_ranks', 'taxonomy_rank')]:
        cells[field] = [row[source] for row in results]
    return cells, results


def refresh_itsme_reports(root):
    """Refresh only saved reporting files. Caller must hold the cohort lock."""
    from compact_outputs import atomic_json, finished_sample, read_json, table
    root = Path(root)
    updates, loci, summaries = {}, [], {}
    for path in sorted(root.glob('*/final/summary.json')):
        sample = path.parent.parent
        if sample.is_symlink() or not (sample/'run_config.json').is_file() or not finished_sample(sample):
            continue
        summary = read_json(path)
        cells, details = summarize_itsme(sample, summary)
        summary.update(cells)
        atomic_json(path, summary)
        mito = sample/'MitoLooper/summary.json'
        if mito.is_file():
            content = read_json(mito)
            content.update(cells)
            atomic_json(mito, content)
        summaries[summary.get('sample', sample.name)] = summary
        updates[summary.get('sample', sample.name)] = cells
        loci.extend(details)
    cohort = root/'cohort_summary.csv'
    if cohort.exists():
        rows, fields = read_rows(cohort)
    else:
        rows = list(summaries.values())
        fields = sorted({field for row in rows for field in row if field not in ITS_FIELDS})
    existing = {row.get('sample') for row in rows}
    for name, summary in summaries.items():
        if name not in existing:
            rows.append(summary)
            fields += [field for field in summary if field not in fields]
    fields += [field for field in ITS_FIELDS if field not in fields]
    for row in rows:
        row.update(updates.get(row.get('sample'), {}))
    serialized = [{key: json.dumps(value) if isinstance(value, (list, dict)) else value
                   for key, value in row.items()} for row in rows]
    table(cohort, serialized, fields, ',')
    table(root/'cohort_its_loci.csv', loci, LOCUS_FIELDS, ',')
    invocation = root/'last_invocation.json'
    if invocation.is_file():
        content = read_json(invocation)
        for row in content.get('results', []):
            row.update(updates.get(row.get('sample'), {}))
        atomic_json(invocation, content)
    # The dashboard only reads saved evidence; report generation must never rerun analysis.
    if updates:
        from build_dashboard import write_dashboard
        try:
            write_dashboard(root)
        except (OSError, ValueError, TypeError, csv.Error) as error:
            import sys
            print(f"WARNING: Dashboard could not be generated: {error}", file=sys.stderr)
    return len(updates), len(loci)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run_dir', type=Path, help='Existing Tortuga cohort directory; reporting only')
    args = parser.parse_args(argv)
    root = args.run_dir.absolute()
    if root.is_symlink() or not root.is_dir() or not any(root.glob('*/final/summary.json')):
        parser.error('Not a Tortuga results directory with saved sample summaries')
    with (root/'.pipeline.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            parser.error('The pipeline is still running; reports were not changed. Retry after it stops.')
        samples, loci = refresh_itsme_reports(root)
    print(f'Updated ITS reporting for {samples} completed samples ({loci} complete loci); no analyses or cleanup ran.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

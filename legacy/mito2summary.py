#!/usr/bin/env python3

"""
Summarize a single GetOrganelle sample from three inputs:
  1) Hits CSV (tab- or comma-delimited) with 'details' column listing genes and coordinates
  2) Assembled FASTA
  3) Log file (to read total cost time in seconds)

Outputs:
  <out_prefix>_run_summary.csv   (one row: sample, total_mito_bp, num_fragments, total_cost_seconds)
  <out_prefix>_gene_lengths_long.csv
      columns: sample, gene, present, max_length_bp, pct_expected_length

Usage:
  python getorganelle_triplet2csv.py \
      -c B10156-COI-SKQ23_getOrganelleReads_NOseed.csv \
      -f B10156-COI-SKQ23_getOrganelleReads_NOseed.fasta \
      -l B10156-COI-SKQ23_getOrganelleReads_NOseed.log \
      -o B10156-COI-SKQ23 \
      --expected-csv mitochondrial_gene_average_lengths.csv   # optional
"""

import re
import csv
import argparse
from pathlib import Path
from typing import Optional, Dict, Tuple, List
from collections import defaultdict
import pandas as pd

# ---------- FASTA + LOG ----------
def fasta_total_length_and_count(fasta_path: Path) -> Tuple[int, int]:
    total = 0
    count = 0
    seq_len = 0
    with fasta_path.open("r", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_len > 0:
                    total += seq_len
                    seq_len = 0
                count += 1
            else:
                seq_len += len(line)
        if seq_len > 0:
            total += seq_len
    return total, count

def parse_total_cost_seconds(log_text: str) -> Optional[float]:
    # Matches e.g., "Total cost 2810.65 s" (or min/hr variants)
    m = re.search(r"Total\s+cost\s+([\d.]+)\s*([A-Za-z]+)", log_text)
    if not m:
        return None
    val = float(m.group(1))
    unit = m.group(2).lower()
    to_seconds = {
        "s": 1, "sec": 1, "secs": 1, "second": 1, "seconds": 1,
        "m": 60, "min": 60, "mins": 60, "minute": 60, "minutes": 60,
        "h": 3600, "hr": 3600, "hrs": 3600, "hour": 3600, "hours": 3600
    }
    if unit not in to_seconds:
        return None
    return val * to_seconds[unit]

# ---------- HITS CSV → longest per gene ----------
GENE_COORD_RE = re.compile(r'^([A-Za-z0-9_.+-]+)\((\d+)-(\d+),')

def parse_details_field(details_str: str) -> List[Tuple[str, int, int]]:
    """
    Parse a 'details' string like:
      'cox1(2020-3545,animal_mt)>>cox2(3665-4347,animal_mt)>>...'
    into a list of (gene, start, end).
    """
    if not isinstance(details_str, str) or not details_str.strip():
        return []
    parts = [p.strip() for p in details_str.split('>>') if p.strip()]
    hits = []
    for p in parts:
        m = GENE_COORD_RE.match(p)
        if m:
            gene = m.group(1)
            start = int(m.group(2))
            end = int(m.group(3))
            hits.append((gene, start, end))
    return hits

def autodetect_delimiter(first_line: str) -> str:
    if '\t' in first_line and ',' not in first_line:
        return '\t'
    if '\t' in first_line:
        return '\t'
    return ','

def longest_gene_lengths_from_hits(hits_csv: Path) -> Dict[str, int]:
    """
    Return {gene_lower: max_len} within this sample from the hits CSV.
    """
    with hits_csv.open("r", encoding="utf-8", errors="replace") as f:
        first = f.readline()
    delim = autodetect_delimiter(first)

    longest: Dict[str, int] = defaultdict(lambda: -1)
    with hits_csv.open("r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter=delim)
        for row in reader:
            details = row.get("details", "")
            for gene, start, end in parse_details_field(details):
                gl = gene.lower()
                length = abs(end - start) + 1  # inclusive
                if length > longest[gl]:
                    longest[gl] = length
    return dict(longest)

# ---------- Expected gene lengths ----------
def load_expected_lengths(expected_csv: Optional[Path]) -> Dict[str, int]:
    """
    Load expected lengths from CSV if provided, else use defaults.
    CSV columns required: Gene, Average_Length_bp
    """
    defaults: Dict[str, int] = {
        # Common animal mitochondrial genes (approx.)
        "cox1": 1545, "cox2": 684, "cox3": 786,
        "cytb": 1140,
        "atp6": 681, "atp8": 207,
        "nad1": 957, "nad2": 1041, "nad3": 354, "nad4": 1377, "nad4l": 297, "nad5": 1812, "nad6": 525,
        "rrnl": 1200,  # 16S rRNA
        "rrns": 950,   # 12S rRNA
        # Any trn* (tRNAs) will be treated as ~70 bp below
    }
    if expected_csv is None:
        return defaults

    try:
        df = pd.read_csv(expected_csv)
        for _, row in df.iterrows():
            g = str(row["Gene"]).strip().lower()
            L = float(row["Average_Length_bp"])
            if L > 0:
                defaults[g] = int(round(L))
    except Exception:
        pass
    return defaults

def canonical_label(g_lower: str) -> str:
    # Prefer mixed-case labels for some genes
    special = {"rrnl": "rrnL", "rrns": "rrnS", "nad4l": "nad4L"}
    return special.get(g_lower, g_lower)

def expected_length_for_gene(g_lower: str, expected_map: Dict[str, int]) -> Optional[int]:
    if g_lower.startswith("trn"):
        return 70  # generic tRNA length
    return expected_map.get(g_lower)

# ---------- Main CLI ----------
def main():
    ap = argparse.ArgumentParser(description="Summarize a GetOrganelle sample from hits CSV, FASTA, and LOG.")
    ap.add_argument("-c", "--csv", required=True, help="GetOrganelle hits CSV/TSV (with 'details' column).")
    ap.add_argument("-f", "--fasta", required=True, help="Assembled organelle FASTA.")
    ap.add_argument("-l", "--log", required=True, help="GetOrganelle run log file.")
    ap.add_argument("-o", "--output-prefix", required=True, help="Output prefix for CSVs.")
    ap.add_argument("--expected-csv", help="Optional CSV with columns: Gene,Average_Length_bp (overrides defaults).")
    args = ap.parse_args()

    csv_path = Path(args.csv)
    fasta_path = Path(args.fasta)
    log_path = Path(args.log)

    # Sample name
    sample = fasta_path.name.split("_getOrganelleReads_NOseed")[0]

    # FASTA + LOG
    total_bp, n_frags = fasta_total_length_and_count(fasta_path)
    log_text = log_path.read_text(encoding="utf-8", errors="ignore")
    total_secs = parse_total_cost_seconds(log_text)

    # Gene lengths
    expected_map = load_expected_lengths(Path(args.expected_csv) if args.expected_csv else None)
    longest = longest_gene_lengths_from_hits(csv_path)

    # Build long-format gene table (presence & % expected)
    gene_rows = []
    # Expected genes = keys in expected_map + any trn* that appear (handle trn if present)
    expected_genes = set(expected_map.keys()) | {g for g in longest.keys() if g.startswith("trn")}
    for g in sorted(expected_genes):
        max_len = int(longest.get(g, 0))
        present = max_len > 0
        exp = expected_length_for_gene(g, expected_map)
        pct = round((max_len / exp) * 100.0, 1) if exp and exp > 0 else None
        gene_rows.append({
            "sample": sample,
            "gene": canonical_label(g),
            "present": present,
            "max_length_bp": max_len,
            "pct_expected_length": pct
        })
    gene_df = pd.DataFrame(gene_rows).sort_values(["sample", "gene"])

    # Run summary table
    run_df = pd.DataFrame([{
        "sample": sample,
        "total_mito_bp": int(total_bp),
        "num_fragments": int(n_frags),
        "total_cost_seconds": total_secs
    }])

    # Write
    out_run = f"{args.output_prefix}_run_summary.csv"
    out_gene = f"{args.output_prefix}_gene_lengths_long.csv"
    run_df.to_csv(out_run, index=False)
    gene_df.to_csv(out_gene, index=False)
    print(f"Wrote {out_run}")
    print(f"Wrote {out_gene}")

if __name__ == "__main__":
    main()

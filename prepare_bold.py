#!/usr/bin/env python3
"""Normalize a BOLD COI TSV export into a FASTA and auditable taxonomy table."""
import argparse
import csv
# Public BOLD exports can contain fields exceeding csv's 128 KiB default.
# This changes the parser limit; it does not preallocate a buffer of this size.
csv.field_size_limit(2**31 - 1)
import gzip
import json
from pathlib import Path
import re
import sys
from mito_pipeline import digest_file, write_table, fasta_write


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("export", type=Path)
    p.add_argument("--source", required=True, help="BOLD dataset URL/version/date; saved with every record")
    p.add_argument("-o", "--outdir", required=True, type=Path)
    p.add_argument("--sequence-column")
    p.add_argument("--processid-column")
    p.add_argument("--marker-column")
    args = p.parse_args()
    if args.outdir.exists() and any(args.outdir.iterdir()):
        raise ValueError("Output directory must be new or empty")
    args.outdir.mkdir(parents=True, exist_ok=True)
    op = gzip.open if str(args.export).endswith(".gz") else open
    rows, ids, skipped = [], set(), 0
    with op(args.export, "rt", encoding="utf-8-sig", newline="") as handle, (args.outdir/"coi.fasta").open("w") as fasta:
        reader = csv.DictReader(handle, delimiter="\t")
        def choose(explicit, names):
            column = explicit or next((x for x in names if x in (reader.fieldnames or [])), None)
            if column not in (reader.fieldnames or []):
                raise ValueError(f"Missing column; specify an override. Tried: {names}")
            return column
        scol = choose(args.sequence_column,["nuc","nucleotides","nucleotide_sequence","sequence"])
        pcol = choose(args.processid_column,["processid","process_id"])
        mcol = choose(args.marker_column,["marker_code","marker","markercode"])
        ranks = ["kingdom","phylum","class","order","family","subfamily","tribe","genus","species","subspecies"]
        for row in reader:
            marker = row[mcol].strip().upper()
            if marker not in {"COI-5P","COI","CO1","COX1"}:
                skipped += 1
                continue
            seq = re.sub(r"[\s.\-]","",row[scol]).upper().replace("U","T")
            pid = row[pcol].strip()
            taxa = [row.get(rank,"").strip() for rank in ranks if row.get(rank,"").strip()]
            if not pid or not re.fullmatch(r"[A-Za-z0-9_.:-]+",pid) or not taxa or len(seq)<400 or set(seq)-set("ACGTRYSWKMBDHVN"):
                skipped += 1
                continue
            sid = pid+"_"+re.sub(r"[^A-Za-z0-9]","",marker)
            if sid in ids:
                raise ValueError(f"Duplicate barcode {sid}; curate duplicate export records")
            ids.add(sid)
            fasta.write(f">{sid}\n{seq}\n")
            rows.append({"sequence_id":sid,"processid":pid,"taxon":taxa[-1],
                         "lineage":";".join(taxa),"source":args.source})
    if not rows:
        raise ValueError("No eligible COI records; check export columns and marker labels")
    write_table(args.outdir/"metadata.tsv",rows,["sequence_id","processid","taxon","lineage","source"])
    (args.outdir/"provenance.json").write_text(json.dumps({"source":args.source,"input":str(args.export.resolve()),
        "input_sha256":digest_file(args.export),"retained_records":len(rows),"skipped_records":skipped},indent=2)+"\n")
    print(f"Prepared {len(rows)} COI references; skipped {skipped} records")


if __name__ == "__main__":
    try:
        main()
    except (ValueError,OSError,KeyError,csv.Error) as e:
        print(f"ERROR: {e}",file=sys.stderr)
        sys.exit(1)

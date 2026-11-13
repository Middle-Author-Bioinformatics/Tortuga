#!/usr/bin/env python3
"""Sequential SCCWRP mitogenome workflow. See README.md for evidence rules."""
from __future__ import annotations

import argparse
import csv
import dataclasses
import gzip
import hashlib
import json
import os
from pathlib import Path
import random
import re
import shlex
import shutil
import signal
import subprocess
import sys
import time

VERSION = "2.0.8"
READ_METHODS = ("novoplasty", "getorganelle", "mitofinder", "mtgrasp")
DEFAULT_COMMANDS = {
    "novoplasty": ["perl", "NOVOPlasty4.3.5.pl"],
    "getorganelle": ["get_organelle_from_reads.py"],
    "getorganelle_assembly": ["get_organelle_from_assembly.py"],
    "mitofinder": ["mitofinder"], "mtgrasp": ["mtgrasp.py"],
    "metaspades": ["spades.py"], "megahit": ["megahit"],
    "mitos2": ["runmitos.py"], "itsme": ["bash", str(Path.home()/"bin/ITSME/itsme.v2.sh")],
    "blastn": ["blastn"], "tblastn": ["tblastn"],
    "makeblastdb": ["makeblastdb"], "fastp": ["fastp"],
    "minimap2": ["minimap2"],
}


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, default=str) + "\n")
    tmp.replace(path)


def digest_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def stamp(path):
    p = Path(path).resolve()
    s = p.stat()
    return {"path": str(p), "size": s.st_size, "mtime_ns": s.st_mtime_ns}


def output_exists(path):
    p = Path(path)
    return p.exists() or (p.suffix == ".nsq" and any(p.parent.glob(p.stem+".*.nsq")))


def write_table(path, rows, fields=None, delimiter="\t"):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or sorted({k for r in rows for k in r})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter=delimiter,
                                extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v) if isinstance(v, (list, dict)) else v
                             for k, v in row.items()})


def fasta_write(path, name, sequence):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write(f">{name}\n")
        for i in range(0, len(sequence), 80):
            f.write(str(sequence[i:i+80]) + "\n")


def normalize_gene(name):
    s = re.sub(r"[^a-z0-9]", "", name.lower())
    s = re.sub(r"^mt", "", s)
    aliases = {"coi": "cox1", "co1": "cox1", "coii": "cox2", "co2": "cox2",
               "coiii": "cox3", "co3": "cox3", "cob": "cytb", "cib": "cytb",
               "12srrna": "rrns", "rrns": "rrns", "srrna": "rrns", "rrn12": "rrns",
               "16srrna": "rrnl", "rrnl": "rrnl", "lrrna": "rrnl", "rrn16": "rrnl"}
    if s.startswith("nd") and not s.startswith("nad"):
        s = "nad" + s[2:]
    return aliases.get(s, s)


def attempt_dir(root, runner, stage, base):
    """Reuse only a successful attempt; failed tool directories are never overwritten."""
    state = runner.root/"state"/(stage+".json")
    if runner.resume and state.exists():
        old = json.loads(state.read_text())
        previous = Path(old.get("cwd", ""))
        if (old.get("status") == "ok" and previous.is_dir()
            and all(output_exists(p) for p in old.get("required_outputs",[])+old.get("produced_paths",[]))):
            return previous
    base = Path(base)
    base.mkdir(parents=True, exist_ok=True)
    i = 1
    while (base/f"attempt_{i:03d}").exists():
        i += 1
    path = base/f"attempt_{i:03d}"
    path.mkdir()
    return path


def fastq_records(path):
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as f:
        while True:
            a = f.readline()
            if not a:
                return
            b, c, d = f.readline(), f.readline(), f.readline()
            if not a.startswith("@") or not c.startswith("+") or len(b.rstrip()) != len(d.rstrip()):
                raise ValueError(f"Malformed four-line FASTQ: {path}")
            yield (a, b, c, d)


def template_id(record):
    return re.sub(r"/[12]$", "", record[0].split()[0][1:])


def subsample_pairs(r1, r2, out1, out2, fraction, seed):
    """One paired pass with reproducible Bernoulli sampling and mate validation."""
    from itertools import zip_longest
    rng, total, kept = random.Random(seed), 0, 0
    with gzip.open(out1, "wt") as a, gzip.open(out2, "wt") as b:
        for left, right in zip_longest(fastq_records(r1), fastq_records(r2)):
            if left is None or right is None or template_id(left) != template_id(right):
                raise ValueError("FASTQ mates differ in count or template order")
            total += 1
            if rng.random() < fraction:
                kept += 1
                a.writelines(left)
                b.writelines(right)
    if not kept:
        raise ValueError("Assembly read subset is empty; increase --assembly-percent")
    return {"input_templates": total, "retained_templates": kept,
            "requested_percent": fraction*100, "sampling_seed": seed}


@dataclasses.dataclass
class Sample:
    sample: str
    r1: Path
    r2: Path
    target_taxon: str
    mt_code: int
    reference_gb: Path
    seed_fasta: Path | None = None
    taxonomy: dict = dataclasses.field(default_factory=dict)


def discover_pairs(folder):
    patterns = [("_R1_001.fastq.gz", "_R2_001.fastq.gz"),
                ("_R1.fastq.gz", "_R2.fastq.gz"), ("_1.fq.gz", "_2.fq.gz"),
                ("_1.fastq.gz", "_2.fastq.gz"), ("_R1.fq.gz", "_R2.fq.gz"),
                ("_1.fq", "_2.fq"), ("_R1.fastq", "_R2.fastq")]
    found = {}
    for p in sorted(Path(folder).iterdir()):
        for left, right in patterns:
            if p.name.endswith(left):
                name = p.name[:-len(left)]
                mate = p.with_name(name + right)
                if not mate.is_file():
                    raise ValueError(f"Missing mate for {p}: {mate}")
                if name in found:
                    raise ValueError(f"Multiple read pairs for sample {name}; merge lanes first")
                found[name] = (p.resolve(), mate.resolve())
                break
    if not found:
        raise ValueError(f"No supported FASTQ pairs found in {folder}")
    return found


class Runner:
    def __init__(self, root, commands, resume=False, dry_run=False):
        self.root, self.commands = Path(root), commands
        self.resume, self.dry_run = resume, dry_run
        self.events = []

    def skip(self, stage, reason):
        row = {"stage": stage, "status": "skipped", "reason": reason}
        self.events.append(row)
        print(f"  {stage}: SKIPPED ({reason})", flush=True)

    def run(self, stage, tool, args, cwd, timeout, outputs=()):
        cwd = Path(cwd)
        cwd.mkdir(parents=True, exist_ok=True)
        if not self.commands.get(tool):
            reason=f"No executable found for {tool} in the configured Conda environments"
            row={"stage":stage,"status":"planned_missing" if self.dry_run else "failed","reason":reason,
                 "returncode":127,"command":[],"cwd":str(cwd)}
            self.events.append(row)
            if not self.dry_run:
                atomic_json(self.root/"state"/(stage+".json"),row)
            print(f"  {stage}: {reason}",flush=True)
            return False
        cmd = self.commands[tool] + list(map(str, args))
        state = self.root / "state" / (stage + ".json")
        payload = {"command": cmd, "cwd": str(cwd), "timeout_seconds": timeout,
                   "required_outputs":[str(x) for x in outputs]}
        files = {}
        for i, flag in enumerate(args[:-1]):
            if str(flag) in {"-query","-subject","-i","-I","-1","-2","-g","-a","-r","-r1","-r2","-c","-in","-s","--genes"}:
                p = Path(str(args[i+1]))
                if p.is_file():
                    files[str(p)] = (stamp(p) if re.search(r"\.(fastq|fq)(\.gz)?$",p.name)
                                     else {"path":str(p.resolve()),"sha256":digest_file(p)})
        payload["input_fingerprints"] = files
        signature = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
        if state.exists() and self.resume:
            old = json.loads(state.read_text())
            if (old.get("signature") == signature and old.get("status") == "ok"
                and all(output_exists(x) for x in list(outputs)+old.get("produced_paths",[]))):
                self.events.append(dict(old, status="cached"))
                print(f"  {stage}: CACHED", flush=True)
                return True
        if self.dry_run:
            print(f"  {stage}: {shlex.join(cmd)}", flush=True)
            self.events.append(dict(payload, stage=stage, status="planned"))
            return True
        log = self.root / "logs" / (stage + ".log")
        log.parent.mkdir(parents=True, exist_ok=True)
        print(f"  {stage}: RUNNING", flush=True)
        row = dict(payload, stage=stage, signature=signature, started=time.time(), log=str(log))
        try:
            with log.open("w") as f:
                f.write(shlex.join(cmd) + "\n")
                f.flush()
                proc = subprocess.Popen(cmd, cwd=cwd, stdout=f, stderr=subprocess.STDOUT,
                                        start_new_session=True)
                try:
                    code = proc.wait(timeout=timeout or None)
                except (subprocess.TimeoutExpired, KeyboardInterrupt):
                    os.killpg(proc.pid, signal.SIGTERM)
                    try:
                        proc.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        os.killpg(proc.pid, signal.SIGKILL)
                        proc.wait()
                    if sys.exc_info()[0] is KeyboardInterrupt:
                        raise
                    row["status"] = "timeout"
                    code = -signal.SIGTERM
                row.setdefault("status", "ok" if code == 0 else "failed")
                row["returncode"] = code
                if row["status"] == "ok" and not all(output_exists(x) for x in outputs):
                    row.update(status="failed", reason="Expected output missing")
        except OSError as e:
            row.update(status="failed", reason=str(e))
        row["elapsed_seconds"] = round(time.time() - row["started"], 3)
        atomic_json(state, row)
        self.events.append(row)
        print(f"  {stage}: {row['status'].upper()}", flush=True)
        return row["status"] == "ok"


def reference_gene_name(feature):
    gene = feature.qualifiers.get("gene", [""])[0]
    if not gene and feature.type == "rRNA":
        product = feature.qualifiers.get("product", [""])[0]
        compact = re.sub(r"[^a-z0-9]", "", product.lower())
        if compact in {"srrna","ssurrna","smallsubunitribosomalrna"} or "12s" in compact:
            gene = "rrns"
        elif compact in {"lrrna","lsurrna","largesubunitribosomalrna"} or "16s" in compact:
            gene = "rrnl"
    return normalize_gene(gene) if gene else ""


def reference_boundary_warnings(record):
    from Bio.SeqFeature import ExactPosition
    return [{"gene":reference_gene_name(feature),"location":str(feature.location),
             "reason":"Uncertain annotation boundary; usable as a homology query, not a full-length boundary standard"}
            for feature in record.features if feature.type in {"CDS","rRNA"} and feature.location is not None
            and any(not isinstance(pos,ExactPosition) for part in feature.location.parts for pos in (part.start,part.end))]


def reference_data(record, mt_code):
    """Validate and extract reference genes without writing files."""
    from Bio.Seq import Seq
    from Bio.SeqFeature import ExactPosition, BeforePosition, AfterPosition
    if not record.seq or set(str(record.seq).upper()) - set("ACGTRYSWKMBDHVN"):
        raise ValueError("Invalid nucleotide reference")
    proteins, rnas, labels, names = [], [], [], set()
    seed = None
    for feature in record.features:
        if feature.type not in {"CDS", "rRNA"}:
            continue
        gene = reference_gene_name(feature)
        if not gene:
            raise ValueError(f"Reference {record.id} has {feature.type} without a gene name")
        if gene in names:
            raise ValueError(f"Duplicate reference gene {gene}; use a single curated mitogenome")
        names.add(gene)
        nt = str(feature.extract(record.seq))
        if set(nt.upper())-set("ACGT"):
            raise ValueError(f"Reference gene {gene} contains ambiguous DNA bases; curate reference")
        positions = [pos for part in feature.location.parts for pos in (part.start,part.end)]
        if any(not isinstance(pos,(ExactPosition,BeforePosition,AfterPosition)) for pos in positions):
            raise ValueError(f"Reference gene {gene} has unsupported uncertain coordinates; curate reference")
        if any(int(part.start)<0 or int(part.end)>len(record.seq) for part in feature.location.parts) or not nt:
            raise ValueError(f"Reference gene {gene} has invalid coordinates")
        uncertain = any(not isinstance(pos,ExactPosition) for pos in positions)
        if uncertain and feature.type=="CDS" and not feature.qualifiers.get("translation"):
            raise ValueError(f"Reference gene {gene} has an uncertain boundary without an annotated translation")
        labels.append((gene, nt))
        if feature.type == "CDS":
            aa = feature.qualifiers.get("translation", [""])[0]
            if not aa:
                offset = int(feature.qualifiers.get("codon_start", [1])[0])-1
                aa = str(Seq(nt[offset:]).translate(table=mt_code)).rstrip("*")
            if "*" in aa.rstrip("*") or not aa:
                raise ValueError(f"Invalid protein translation for reference {gene}")
            proteins.append((gene, aa.rstrip("*")))
            if gene == "cox1":
                seed = nt
        else:
            rnas.append((gene, nt))
    if seed is None or len(proteins) < 3 or not rnas:
        raise ValueError("Reference must annotate COX1, at least 3 CDS genes, and rRNA genes")
    return proteins, rnas, labels, seed


def reference_files(sample, folder):
    """Derive seeds and per-gene completeness queries from an annotated reference."""
    from Bio import SeqIO
    records = list(SeqIO.parse(sample.reference_gb, "genbank"))
    if len(records) != 1:
        raise ValueError(f"{sample.sample}: reference_gb must contain ONE mitogenome")
    record = records[0]
    proteins, rnas, labels, seed = reference_data(record, sample.mt_code)
    folder.mkdir(parents=True, exist_ok=True)
    fasta_write(folder/"reference.fasta", record.id, str(record.seq))
    SeqIO.write(record,folder/"reference.gb","genbank")
    warnings = reference_boundary_warnings(record)
    atomic_json(folder/"reference_quality.json",{"accession":record.id,"annotation_warnings":warnings,
                "interpretation":"Homology reference; uncertain gene boundaries cannot establish full-length recovery"})
    for filename, rows in [("proteins.faa", proteins), ("rrna.fasta", rnas), ("labels.fasta", labels)]:
        with (folder/filename).open("w") as f:
            for name, seq in rows:
                f.write(f">{name}\n{seq}\n")
    fasta_write(folder/"seed.fasta", "reference_cox1", seed)
    return {"seed": sample.seed_fasta or folder/"seed.fasta", "reference": folder/"reference.fasta",
            "annotation_warnings":warnings,
            "genbank": folder/"reference.gb",
            "proteins": folder/"proteins.faa", "rrna": folder/"rrna.fasta", "labels": folder/"labels.fasta",
            "required_proteins": [x[0] for x in proteins], "required_rrnas": [x[0] for x in rnas]}


BLAST_FIELDS = "qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore"


def blast_rows(path):
    rows = []
    if not Path(path).exists():
        return rows
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        v = line.split("\t")
        if len(v) != 12:
            raise ValueError(f"Unexpected BLAST columns in {path}")
        rows.append(dict(zip(BLAST_FIELDS.split(), v)))
    return rows


def qualifying_genes(rows, min_identity, min_coverage):
    best = {}
    for row in rows:
        coverage = (int(row["qend"])-int(row["qstart"])+1)/int(row["qlen"])
        if float(row["pident"]) >= min_identity and coverage >= min_coverage:
            name = row["qseqid"]
            if name not in best or float(row["bitscore"]) > float(best[name]["bitscore"]):
                best[name] = row
    return best


def bold_decision(rows, metadata, target, identity, coverage, min_bp, near_fraction, max_hits):
    good = [r for r in rows if float(r["pident"]) >= identity and int(r["length"]) >= min_bp
            and int(r["length"])/min(int(r["qlen"]), int(r["slen"])) >= coverage]
    if not good:
        return {"status": "unresolved", "reason": "No qualifying BOLD COI alignment", "hits": []}
    best = max(float(r["bitscore"]) for r in good)
    near = [r for r in good if float(r["bitscore"]) >= best*(1-near_fraction)]
    subjects = set(r["sseqid"] for r in rows)
    if any(float(r["bitscore"]) > best and r not in good for r in rows):
        status, reason = "unresolved", "A higher-scoring BOLD alignment fails quality gates"
    elif len(subjects) >= max_hits:
        status, reason = "unresolved", "BOLD hit ceiling reached; increase --bold-max-hits"
    elif any(r["sseqid"] not in metadata for r in near):
        status, reason = "unresolved", "Missing BOLD record metadata"
    else:
        matches = []
        for r in near:
            m = metadata[r["sseqid"]]
            tokens = {t.strip().casefold() for t in (m["lineage"] + ";" + m["taxon"]).split(";")}
            matches.append(target.casefold() in tokens)
        if all(matches):
            status, reason = "confirmed_target", "All near-best qualifying BOLD records match requested taxon"
        elif any(matches):
            status, reason = "unresolved", "Near-best BOLD records conflict at requested taxon"
        else:
            status, reason = "non_target", "No near-best qualifying BOLD record matches requested taxon"
    hits = [dict(r, metadata=metadata.get(r["sseqid"], {})) for r in near]
    return {"status": status, "reason": reason, "requested_taxon": target, "hits": hits}


def terminal_overlap(seq, minimum=100, maximum=2000):
    maximum = min(maximum, len(seq)//3)
    for size in range(maximum, minimum-1, -1):
        tip = seq[:size]
        if tip == seq[-size:] and len(set(tip)) == 4:
            # Reject low-complexity repeat ends: at least half of 15-mers must be unique.
            words = [tip[i:i+15] for i in range(len(tip)-14)]
            if len(set(words))/len(words) >= .5:
                return size
    return 0


def candidate_files(method, folder):
    if method == "novoplasty":
        patterns = ["Circularized_assembly*.fasta", "Contigs*.fasta", "Option*.fasta"]
    elif method.startswith("getorganelle"):
        patterns = ["*.path_sequence.fasta"]
    elif method.startswith("mitofinder"):
        patterns = ["*_mtDNA_contig.fasta"]
    else:
        patterns = ["*.final-mtgrasp-assembly.fa", "*.postsealer.postpilon.fasta"]
    return sorted({p for pattern in patterns for p in folder.rglob(pattern) if p.is_file()})


def native_circular(method, filename, description):
    if method.startswith("getorganelle") and re.search(r"\.complete\.graph", filename):
        return "getorganelle_closed_graph"
    if method == "novoplasty" and filename.startswith("Circularized_assembly"):
        return "novoplasty_circularized"
    if (method == "mtgrasp" and re.search(r"(?<![A-Za-z])circular(?![A-Za-z])", description, re.I)
        and not re.search(r"(?<![A-Za-z])linear(?![A-Za-z])|circular[=:_ -]+(?:false|no|0)(?![A-Za-z])",description,re.I)):
        return "mtgrasp_header_circular"
    return "not_established"


def junction_support(sam_path, junction, flank=20, min_mapq=20):
    """Count distinct templates with a continuous alignment across the artificial join."""
    templates, starts = set(), set()
    with open(sam_path) as f:
        for line in f:
            if line.startswith("@"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields)<11:
                continue
            flag = int(fields[1])
            if fields[2] != "junction" or flag & (4|256|2048|512|1024) or int(fields[4]) < min_mapq:
                continue
            pos = int(fields[3])-1
            origin = pos
            for length,op in re.findall(r"(\d+)([MIDNSHP=X])",fields[5]):
                length = int(length)
                if op in "M=X":
                    if pos <= junction-flank and pos+length >= junction+flank:
                        templates.add(re.sub(r"/[12]$","",fields[0]))
                        starts.add((origin,bool(flag&16)))
                    pos += length
                elif op in "DN":
                    pos += length
    return {"spanning_templates":len(templates),"distinct_alignment_starts":len(starts),
            "minimum_mapq":min_mapq,"required_flank_bp":flank}


def summarize_annotation(gff):
    from urllib.parse import unquote
    features, genes = {}, set()
    with open(gff) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            v = line.rstrip().split("\t")
            if len(v) != 9:
                continue
            features[v[2]] = features.get(v[2],0)+1
            attributes = dict(x.split("=",1) for x in v[8].split(";") if "=" in x)
            gene = attributes.get("Name") or attributes.get("gene") or attributes.get("gene_id")
            if gene:
                genes.add(unquote(gene))
    return {"feature_counts":features,"gene_labels":sorted(genes),"gff":str(gff),
            "review_required":"Inspect MITOS2 tRNA/rRNA content and coding boundaries; reference-gene recovery is a screening criterion, not a finished annotation guarantee."}


def acceptance_policy(args):
    """Stopping thresholds are independent of reference-database length filters."""
    return {"topology":args.accept_topology,"rule":args.accept_rule,
            "min_detected_genes":args.min_detected_genes,
            "min_length_bp":args.min_accepted_bp if args.min_accepted_bp is not None else args.min_genome_bp,
            "max_length_bp":args.max_genome_bp,
            "candidate_floor_bp":args.min_candidate_bp,
            "gene_count_definition":"Unique reference CDS plus rRNA genes passing homology identity/coverage; excludes tRNAs"}


def validate_acceptance_options(args):
    policy=acceptance_policy(args)
    if not 0 < args.min_candidate_bp <= args.min_genome_bp <= args.max_genome_bp:
        raise ValueError("Require 0 < --min-candidate-bp <= --min-genome-bp <= --max-genome-bp")
    if not args.min_candidate_bp <= policy["min_length_bp"] <= args.max_genome_bp:
        raise ValueError("--min-accepted-bp must be between --min-candidate-bp and --max-genome-bp")
    count=policy["min_detected_genes"]
    if count is not None and count<1:
        raise ValueError("--min-detected-genes must be positive")
    if policy["rule"] in {"genes","either","both"} and count is None:
        raise ValueError("The selected --accept-rule requires --min-detected-genes")
    if policy["rule"] in {"all-genes","length"} and count is not None:
        raise ValueError("Use --accept-rule genes, either or both with --min-detected-genes")
    if policy["rule"]=="genes" and args.min_accepted_bp is not None:
        raise ValueError("Use --accept-rule both or either to combine gene count with --min-accepted-bp")


def apply_acceptance(result,args):
    """An explicit policy can stop recovery without claiming complete mitogenomes."""
    policy=acceptance_policy(args)
    count=len(set(result.get("genes_found",[])) | set(result.get("rrnas_found",[])))
    long_enough=result["length_bp"]>=policy["min_length_bp"]
    enough_genes=(count>=policy["min_detected_genes"] if policy["min_detected_genes"] is not None else False)
    criteria={"all-genes":bool(result["complete_by_reference"]) and long_enough,
              "genes":enough_genes,"length":long_enough,
              "either":enough_genes or long_enough,"both":enough_genes and long_enough}
    checks={"validation_succeeded":not bool(result.get("validation_error")),
            "target_confirmed":result["target_status"]=="confirmed_target",
            "unambiguous_sequence":result["ambiguous_bp"]==0,
            "candidate_length_bounds":policy["candidate_floor_bp"]<=result["length_bp"]<=policy["max_length_bp"],
            "topology":policy["topology"]=="any" or result["circular"],
            "selected_rule":criteria[policy["rule"]]}
    rule_description={"all-genes":f"all reference CDS/rRNA genes and at least {policy['min_length_bp']} bp",
                      "genes":f"at least {policy['min_detected_genes']} CDS/rRNA genes",
                      "length":f"at least {policy['min_length_bp']} bp",
                      "either":f"at least {policy['min_detected_genes']} CDS/rRNA genes OR {policy['min_length_bp']} bp",
                      "both":f"at least {policy['min_detected_genes']} CDS/rRNA genes AND {policy['min_length_bp']} bp"}[policy["rule"]]
    reasons={"validation_succeeded":"validation command failed",
             "target_confirmed":"target identity is "+result["target_status"],
             "unambiguous_sequence":f"{result['ambiguous_bp']} ambiguous bases",
             "candidate_length_bounds":f"length must be {policy['candidate_floor_bp']}-{policy['max_length_bp']} bp",
             "topology":"circularity required but not established",
             "selected_rule":"requires "+rule_description}
    result.update(detected_gene_count=count,detected_cds_count=len(set(result.get("genes_found",[]))),
                  detected_rrna_count=len(set(result.get("rrnas_found",[]))),
                  acceptance_policy=policy,acceptance_checks=checks,
                  acceptance_failures=[name for name,ok in checks.items() if not ok],
                  acceptance_reason=("met "+rule_description if all(checks.values()) else
                                     "; ".join(reasons[name] for name,ok in checks.items() if not ok)),
                  gene_count_threshold_met=enough_genes,length_threshold_met=long_enough,
                  accepted=all(checks.values()))
    result["strict_complete_candidate"]=(checks["validation_succeeded"] and checks["target_confirmed"]
        and checks["unambiguous_sequence"] and result["circular"] and result["complete_by_reference"]
        and args.min_genome_bp<=result["length_bp"]<=args.max_genome_bp)
    return result


def assess_candidate(sequence, description, filename, method, sample, refs, runner, folder,
                     args, bold_db, metadata, mitofinder_circle=False):
    from Bio.Seq import Seq
    folder.mkdir(parents=True, exist_ok=True)
    seq = str(sequence).upper()
    overlap = terminal_overlap(seq)
    evidence = native_circular(method, filename, description)
    result = {"candidate_id": folder.name, "method": method, "source_file": filename,
              "source_record": description, "length_bp": len(seq), "ambiguous_bp": sum(x not in "ACGT" for x in seq),
              "terminal_overlap_bp": overlap, "circular_evidence": evidence,
              "circular": evidence != "not_established", "fasta": str(folder/"sequence.fasta"),
              "target_status": "unresolved", "complete_by_reference": False, "accepted": False}
    # Establish topology before homology screening. Unclosed contigs retain their
    # ends and are screened once, so linear acceptance cannot count artificial
    # end-to-start gene alignments or silently remove terminal repeats.
    if not result["circular"] and (overlap or mitofinder_circle):
        junction_seq=seq[:-overlap] if overlap else seq
        flank_size=min(500,len(junction_seq)//3)
        reference=folder/"junction_reference.fasta"
        reference.write_text(">linear_candidate\n"+junction_seq+"\n>junction\n"+
                             junction_seq[-flank_size:]+junction_seq[:flank_size]+"\n")
        sam=folder/"junction.sam"
        if runner.run(folder.name+"_junction","minimap2",
                      ["-ax","sr","--secondary=no","--sam-hit-only","-t",args.threads,
                       "-o",sam,reference,sample.r1,sample.r2],folder,
                      args.validation_timeout*60,[sam]):
            support=junction_support(sam,flank_size)
            result["junction_support"]=support
            if support["spanning_templates"]>=args.min_junction_templates and support["distinct_alignment_starts"]>=2:
                result["circular"]=True
                result["circular_evidence"]=("mitofinder_log_and_raw_junction_reads" if mitofinder_circle
                                                else "terminal_overlap_and_raw_junction_reads")
            atomic_json(folder/"junction_support.json",support)
    removed=overlap if result["circular"] else 0
    if removed:
        seq=seq[:-removed]
    result.update(length_bp=len(seq),ambiguous_bp=sum(x not in "ACGT" for x in seq),
                  terminal_overlap_removed_bp=removed,
                  homology_sequence_mode="circular_doubled" if result["circular"] else "linear")
    fasta_write(folder/"sequence.fasta",sample.sample,seq)
    fasta_write(folder/"screen.fasta","candidate",seq+seq if result["circular"] else seq)
    # Overlap alone is retained as a review flag; it is not proof of circularity.
    ok = runner.run(folder.name+"_proteins", "tblastn",
                    ["-query", refs["proteins"], "-subject", folder/"screen.fasta",
                     "-db_gencode", sample.mt_code, "-evalue", "1e-10", "-num_threads", args.threads,
                     "-outfmt", "6 "+BLAST_FIELDS, "-out", folder/"proteins.tsv"],
                    folder, args.validation_timeout*60, [folder/"proteins.tsv"])
    ok2 = runner.run(folder.name+"_rrna", "blastn",
                     ["-query", refs["rrna"], "-subject", folder/"screen.fasta", "-task", "blastn",
                      "-evalue", "1e-10", "-num_threads", args.threads,
                      "-outfmt", "6 "+BLAST_FIELDS, "-out", folder/"rrna.tsv"],
                     folder, args.validation_timeout*60, [folder/"rrna.tsv"])
    if not ok or not ok2:
        result["validation_error"] = "Gene screening failed; candidate cannot pass"
        apply_acceptance(result,args)
        atomic_json(folder/"assessment.json", result)
        return result
    prows = blast_rows(folder/"proteins.tsv")
    genes = qualifying_genes(prows, args.protein_identity, args.gene_coverage)
    rnas = qualifying_genes(blast_rows(folder/"rrna.tsv"), args.rrna_identity, args.gene_coverage)
    result["genes_found"] = sorted(genes)
    result["rrnas_found"] = sorted(rnas)
    result["missing_genes"] = sorted(set(refs["required_proteins"])-set(genes))
    result["missing_rrnas"] = sorted(set(refs["required_rrnas"])-set(rnas))
    result["complete_by_reference"] = not result["missing_genes"] and not result["missing_rrnas"]
    coi_hits = [r for r in prows if r["qseqid"] == "cox1"]
    if coi_hits:
        hit = max(coi_hits, key=lambda r: float(r["bitscore"]))
        start, end = int(hit["sstart"]), int(hit["send"])
        # Protein HSP coordinates delimit the recovered marker; BOLD aligns to its barcode region.
        coi = (seq+seq if result["circular"] else seq)[min(start,end)-1:max(start,end)]
        if start > end:
            coi = str(Seq(coi).reverse_complement())
        fasta_write(folder/"coi.fasta", "COI", coi)
        if runner.run(folder.name+"_bold", "blastn",
                      ["-query", folder/"coi.fasta", "-db", bold_db, "-task", "blastn",
                       "-evalue", "1e-20", "-max_target_seqs", args.bold_max_hits,
                       "-num_threads", args.threads, "-outfmt", "6 "+BLAST_FIELDS,
                       "-out", folder/"bold_hits.tsv"], folder, args.validation_timeout*60,
                      [folder/"bold_hits.tsv"]):
            decision = bold_decision(blast_rows(folder/"bold_hits.tsv"), metadata,
                                     sample.target_taxon, args.coi_identity, args.coi_coverage,
                                     args.coi_min_bp, args.bold_near_fraction, args.bold_max_hits)
            result["target_status"] = decision["status"]
            atomic_json(folder/"bold_evidence.json", decision)
    apply_acceptance(result,args)
    atomic_json(folder/"assessment.json", result)
    return result


def score_candidate(c):
    return (c["accepted"], c["target_status"] == "confirmed_target",
            len(c.get("genes_found", [])) + len(c.get("rrnas_found", [])),
            -c["ambiguous_bp"], c["length_bp"])


GENE_FIELDS = ["sample","method","gene","type","reference_length","aligned_reference_length",
               "query_coverage","identity_percent","screen_pass"]


def gene_report(sample, refs, candidate, args):
    from Bio import SeqIO
    rows = []
    for kind, filename, queries, min_identity in [
        ("CDS","proteins.tsv",refs["proteins"],args.protein_identity),
        ("rRNA","rrna.tsv",refs["rrna"],args.rrna_identity)]:
        hits = blast_rows(Path(candidate["fasta"]).parent/filename) if candidate else []
        best = {}
        for hit in hits:
            if hit["qseqid"] not in best or float(hit["bitscore"])>float(best[hit["qseqid"]]["bitscore"]):
                best[hit["qseqid"]] = hit
        best.update(qualifying_genes(hits,min_identity,args.gene_coverage))
        for gene in SeqIO.parse(queries,"fasta"):
            hit = best.get(gene.id)
            length = int(hit["qend"])-int(hit["qstart"])+1 if hit else 0
            cov = length/len(gene.seq)
            identity = float(hit["pident"]) if hit else 0
            rows.append({"sample":sample.sample,"method":candidate["method"] if candidate else "",
                         "gene":gene.id,"type":kind,"reference_length":len(gene.seq),
                         "aligned_reference_length":length,"query_coverage":round(cov,4),
                         "identity_percent":identity,"screen_pass":cov>=args.gene_coverage and identity>=min_identity})
    return rows


def best_rescue(candidates):
    eligible = [c for c in candidates if c["target_status"] == "confirmed_target" and not c["circular"]]
    return max(eligible, key=score_candidate) if eligible else None


def novo_config(sample, refs, r1, r2, folder, args):
    first = next(fastq_records(sample.r1 if args.dry_run else r1), None)
    if first is None:
        raise ValueError("Empty FASTQ")
    fields = [("Project name", sample.sample), ("Type", "mito"),
              ("Genome Range", f"{args.min_genome_bp}-{args.max_genome_bp}"),
              ("K-mer", args.novo_k), ("Max memory", args.memory), ("Extended log", 0),
              ("Save assembled reads", "no"), ("Seed Input", refs["seed"]),
              ("Extend seed directly", "no"), ("Reference sequence", ""),
              ("Variance detection", ""), ("Chloroplast sequence", "")]
    text = "Project:\n-----------------------\n" + "".join(f"{k:<22}= {v}\n" for k,v in fields)
    text += "Dataset 1:\n-----------------------\n"
    text += "".join(f"{k:<22}= {v}\n" for k,v in [
        ("Read Length", len(first[1].strip())), ("Insert size", args.insert_size),
        ("Platform", "illumina"), ("Single/Paired", "PE"), ("Combined reads", ""),
        ("Forward reads", r1), ("Reverse reads", r2), ("Store Hash", "")])
    text += ("Heteroplasmy:\n-----------------------\nMAF =\nHP exclude list =\nPCR-free =\n"
             "Optional:\n-----------------------\nInsert size auto = yes\nUse Quality Scores = yes\n"
             "Reduce ambigious N's = no\nOutput path = " + str(folder) + "/\n")
    path = folder/"config.txt"
    path.write_text(text)
    return path


def method_args(method, sample, refs, r1, r2, folder, args, graph=None, contigs=None, rescue=None):
    seed = Path(rescue["fasta"]) if rescue else refs["seed"]
    if method == "novoplasty":
        return ["-c", novo_config(sample, refs, r1, r2, folder, args)]
    if method == "getorganelle":
        return ["-1", r1, "-2", r2, "-s", seed, "-F", "animal_mt", "-o", folder/"output",
                "-t", args.threads, "-R", args.getorg_rounds, "--disentangle-time-limit", args.disentangle_seconds]
    if method == "getorganelle_assembly":
        # Assembly mode has no read seed (-s). A verified rescue sequence can label/slim the graph.
        labels = refs["labels"]
        if rescue:
            labels = folder/"rescue_labels.fasta"
            labels.write_text(refs["labels"].read_text() + Path(rescue["fasta"]).read_text())
        return ["-g", graph, "--genes", labels, "-F", "animal_mt", "-o", folder/"output",
                "-t", args.threads, "--disentangle-time-limit", args.disentangle_seconds,
                "--expected-min-size", args.min_genome_bp, "--expected-max-size", args.max_genome_bp]
    if method.startswith("mitofinder"):
        options = ["-j", sample.sample, "-r", refs["genbank"], "-o", sample.mt_code,
                   "-p", args.threads, "-m", args.memory]
        if contigs:
            return options + ["-a", contigs]
        return options + ["-1", r1, "-2", r2, "--megahit"]
    if method == "mtgrasp":
        return ["-r1", r1, "-r2", r2, "-r", seed if rescue else refs["reference"],
                "-m", sample.mt_code, "-o", folder/"output", "-t", args.threads]
    raise ValueError(method)


def load_bold_metadata(path):
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not {"sequence_id", "processid", "taxon", "lineage", "source"} <= set(reader.fieldnames or []):
            raise ValueError("BOLD metadata requires sequence_id, processid, taxon, lineage, source columns")
        rows = {}
        for row in reader:
            if row["sequence_id"] in rows or not all(row[k].strip() for k in ["sequence_id", "processid", "taxon", "source"]):
                raise ValueError("Duplicate or incomplete BOLD metadata record")
            rows[row["sequence_id"]] = row
    if not rows:
        raise ValueError("BOLD metadata is empty")
    return rows


def process_sample(sample, args, commands, metadata, bold_db):
    from Bio import SeqIO
    root = args.outdir/sample.sample
    root.mkdir(parents=True, exist_ok=True)
    runner = Runner(root, commands, args.resume, args.dry_run)
    checkpoint = root/"run_config.json"
    inputs = {"sample": dataclasses.asdict(sample), "reads": [stamp(sample.r1), stamp(sample.r2)],
              "reference_sha256": digest_file(sample.reference_gb),
              "seed_sha256": digest_file(sample.seed_fasta) if sample.seed_fasta else None,
              "bold_fasta": stamp(args.bold_fasta), "bold_metadata_sha256": digest_file(args.bold_metadata),
              "commands": commands, "version": VERSION,
              "options": {k: v for k,v in vars(args).items() if k not in {"resume", "dry_run", "make_sample_sheet", "samples", "input_dir", "only_sample","keep_intermediates"}}}
    serialized = json.loads(json.dumps(inputs, default=str))
    if checkpoint.exists():
        if not args.resume:
            raise ValueError(f"{root} already has a run; use --resume or a new output directory")
        if json.loads(checkpoint.read_text()) != serialized:
            raise ValueError(f"{sample.sample}: inputs/options changed; use a new output directory")
        compacted=root/"cleanup.json"
        if compacted.exists() and json.loads(compacted.read_text()).get("compacted"):
            from compact_outputs import validate_retained
            validate_retained(root)
            print("  outputs: CACHED finalized compact results",flush=True)
            return json.loads((root/"final/summary.json").read_text())
    elif any(root.iterdir()):
        raise ValueError(f"{root} is nonempty without a matching run checkpoint")
    atomic_json(checkpoint, serialized)
    refs = reference_files(sample, root/"references")
    r1, r2 = sample.r1, sample.r2
    if not args.no_qc:
        qc = root/"qc"
        a, b = qc/"R1.fq.gz", qc/"R2.fq.gz"
        opts = ["-i", r1, "-I", r2, "-o", a, "-O", b,
                "--thread", min(args.threads, 16), "--detect_adapter_for_pe",
                "--length_required", 50, "--json", qc/"fastp.json", "--html", qc/"fastp.html"]
        if not runner.run("qc", "fastp", opts, qc, args.method_timeout*60, [a,b]):
            raise RuntimeError("Read QC failed")
        r1, r2 = a, b
    candidates = []
    winner = None
    final = root/"final"
    final.mkdir(exist_ok=True)
    if args.resume and not args.dry_run:
        for name in ["best_sequence.fasta","accepted_candidate.fasta","circular_target_candidate.fasta","annotation_summary.json"]:
            (final/name).unlink(missing_ok=True)
    annotation_ok = None
    mitos = None

    def annotate_winner():
        nonlocal annotation_ok, mitos
        accepted_path=final/"accepted_candidate.fasta"
        shutil.copyfile(winner["fasta"],accepted_path)
        if winner["strict_complete_candidate"]:
            shutil.copyfile(winner["fasta"],final/"circular_target_candidate.fasta")
        mitos = attempt_dir(root,runner,"mitos2",root/"05_mitos2")
        annotation_args=["-i",accepted_path,"-c",sample.mt_code,
                         "--refdir",args.mitos_refdir,"-r",args.mitos_refseqver,
                         "-o",mitos/"output"]
        if not winner["circular"]:
            annotation_args.append("--linear")
        annotation_ok = runner.run("mitos2","mitos2",
                          annotation_args,mitos,args.annotation_timeout*60,
                          [mitos/"output/result.gff"])
        if annotation_ok:
            report = summarize_annotation(mitos/"output/result.gff")
            atomic_json(final/"annotation_summary.json",report)
            if not report["feature_counts"]:
                annotation_ok = False
                runner.events.append({"stage":"annotation_validation","status":"failed",
                                      "reason":"MITOS2 GFF has no features"})

    def recover(method, number, graph=None, contigs=None):
        nonlocal winner
        tool = "mitofinder" if method.startswith("mitofinder") else method
        if tool in args.skip_method:
            runner.skip(method, "disabled by --skip-method")
            return
        folder = attempt_dir(root, runner, method, root/f"{number}_{method}")
        rescue = best_rescue(candidates) if method in {"getorganelle","getorganelle_assembly","mtgrasp"} else None
        if rescue:
            atomic_json(folder/"rescue_seed.json", {"candidate_id": rescue["candidate_id"],
                        "fasta": rescue["fasta"], "target_status": rescue["target_status"]})
        params = method_args(method, sample, refs, r1, r2, folder, args, graph, contigs, rescue)
        success = runner.run(method, tool, params, folder, args.method_timeout*60)
        if args.dry_run or not success:
            return
        files = candidate_files(method,folder)
        state_path = root/"state"/(method+".json")
        stage_state = json.loads(state_path.read_text())
        stage_state["produced_paths"] = list(map(str,files))
        atomic_json(state_path,stage_state)
        record_count = sum(1 for f in files for _ in SeqIO.parse(f,"fasta"))
        method_log = root/"logs"/(method+".log")
        logtext = method_log.read_text(errors="replace") if method_log.exists() else ""
        mitofinder_circle = (method.startswith("mitofinder") and record_count == 1
                            and "Evidences of circularization were found!" in logtext
                            and "Evidences of circularization could not be found" not in logtext)
        count = 0
        seen = set()
        for file in files:
            for rec in SeqIO.parse(file, "fasta"):
                seq = str(rec.seq).upper()
                if len(seq) < args.min_candidate_bp or seq in seen:
                    continue
                seen.add(seq)
                count += 1
                cdir = root/"candidates"/f"{number}_{method}_{count:03d}"
                candidate = assess_candidate(seq, rec.description, file.name, method, sample, refs,
                                             runner, cdir, args, bold_db, metadata,mitofinder_circle)
                candidate["source_path"] = str(file)
                candidates.append(candidate)
                print(f"  Candidate {candidate['candidate_id']}: {candidate['length_bp']} bp; "
                      f"{candidate['detected_gene_count']} genes; circular={candidate['circular']}; "
                      f"target={candidate['target_status']}; "
                      +("ACCEPTED: " if candidate["accepted"] else "continue: ")
                      +candidate["acceptance_reason"],flush=True)
                if count >= args.max_candidates:
                    runner.skip(method+"_candidate_limit", f"Limit {args.max_candidates}; review extra paths manually")
                    break
            if count >= args.max_candidates:
                break
        eligible = [c for c in candidates if c["accepted"]]
        if eligible:
            winner = max(eligible, key=score_candidate)

    for i, method in enumerate(READ_METHODS, 1):
        if winner:
            runner.skip(method, "Candidate already meets configured acceptance policy")
        else:
            recover(method, f"01{i}")
    if winner:
        annotate_winner()

    # Nuclear rDNA remains a distinct optional branch and never gates mito recovery.
    itsme_ok = None
    if args.branch == "full":
        idir = attempt_dir(root, runner, "itsme", root/"02_itsme")
        itsme_ok = runner.run("itsme", "itsme", ["-1", sample.r1, "-2", sample.r2,
                              "--db-dir", args.itsme_db, "-o", idir/"output",
                              "-t", args.threads, "-m", args.memory, "--sensitivity", args.itsme_sensitivity],
                              idir, args.itsme_timeout*60,
                              [idir/"output"/"run_summary.txt"])
    else:
        runner.skip("itsme", "--branch mito")

    wgs_ok = None
    if args.no_assembly or (winner and args.skip_wgs_if_accepted):
        runner.skip("wgs", "--no-assembly" if args.no_assembly else "Candidate meets acceptance policy; --skip-wgs-if-accepted")
    else:
        wgs = attempt_dir(root, runner, "wgs", root/"03_wgs")
        a, b = r1, r2
        if args.assembly_percent < 100:
            a, b = wgs/"subset_R1.fq.gz", wgs/"subset_R2.fq.gz"
            if not args.dry_run:
                counts = subsample_pairs(r1,r2,a,b,args.assembly_percent/100,args.sampling_seed)
                atomic_json(wgs/"sampling.json", counts)
        if args.assembler == "metaspades":
            opts = ["--meta", "--only-assembler", "-1", a, "-2", b, "-o", wgs/"output",
                    "-t", args.threads, "-m", args.memory]
            graph, contigs = wgs/"output/assembly_graph.fastg", wgs/"output/contigs.fasta"
        else:
            opts = ["-1", a, "-2", b, "-o", wgs/"output", "-t", args.threads,
                    "-m", args.memory*1024**3]
            graph, contigs = None, wgs/"output/final.contigs.fa"
        wgs_ok = runner.run("wgs", args.assembler, opts, wgs, args.assembly_timeout*60,
                            [contigs] + ([graph] if graph else []))
        if not winner and wgs_ok:
            if graph:
                recover("getorganelle_assembly", "041", graph=graph)
            else:
                runner.skip("getorganelle_assembly", "Selected assembler supplies no supported FASTG/GFA graph")
            if not winner:
                recover("mitofinder_assembly", "042", contigs=contigs)
            else:
                runner.skip("mitofinder_assembly", "Candidate meets configured acceptance policy")
        elif not winner:
            runner.skip("assembly_recovery", "WGS failed or timed out")
        else:
            runner.skip("assembly_recovery", "Read-based recovery already succeeded")

    best = winner or (max(candidates, key=score_candidate) if candidates else None)
    if best:
        shutil.copyfile(best["fasta"], final/"best_sequence.fasta")
    if winner and annotation_ok is None:
        annotate_winner()
    elif not winner:
        runner.skip("mitos2", "No candidate met configured acceptance policy and target gates")
    strict_winner=bool(winner and winner["strict_complete_candidate"])
    status = ("recovered_annotated" if strict_winner and annotation_ok else
              "recovered_annotation_failed" if strict_winner else
              "accepted_by_policy_annotated" if winner and annotation_ok else
              "accepted_by_policy_annotation_failed" if winner else
              "partial_target_confirmed" if best and best["target_status"] == "confirmed_target" else
              "unresolved_candidate" if best else "no_candidate")
    summary = {"sample": sample.sample, "pipeline":args.pipeline or ("both" if args.branch=="full" else "mitolooper"),
               "status": "planned" if args.dry_run else status,
               "target_taxon": sample.target_taxon, "taxonomy": sample.taxonomy, "mt_code": sample.mt_code,
               "reference_gb":str(sample.reference_gb),"reference_annotation_warnings":refs["annotation_warnings"],
               "winning_method": winner["method"] if winner else "", "best_candidate": best,
               "accepted_by_policy":bool(winner),"strict_complete_candidate":strict_winner,
               "acceptance_policy":acceptance_policy(args),
               "complete_by_reference":bool(best and best["complete_by_reference"]),
               "candidate_count": len(candidates), "itsme_ok": itsme_ok, "wgs_ok": wgs_ok,
               "annotation_ok": annotation_ok, "annotation_dir": str(mitos/"output") if winner else "",
               "best_method": best["method"] if best else "", "best_length_bp": best["length_bp"] if best else 0,
               "circular": best["circular"] if best else False,
               "target_status": best["target_status"] if best else "unresolved",
               "genes_recovered": len(best.get("genes_found",[]))+len(best.get("rrnas_found",[])) if best else 0,
               "genes_expected": len(refs["required_proteins"])+len(refs["required_rrnas"]),
               "failed_stages": [e["stage"] for e in runner.events if e["status"] in {"failed", "timeout"}],
               "completeness_definition": "All reference CDS/rRNA genes pass configured homology coverage; no ambiguous bases; expected genome length. Policy acceptance can stop recovery without meeting this screen or establishing circularity. Neither screen proves full tRNA content or absence of NUMTs/chimeras."}
    atomic_json(final/"summary.json", summary)
    write_table(final/"gene_recovery.csv",gene_report(sample,refs,best,args),GENE_FIELDS,",")
    atomic_json(root/"stages.json", runner.events)
    write_table(root/"candidates.tsv", candidates,
                ["candidate_id","method","length_bp","circular","circular_evidence","terminal_overlap_bp",
                 "target_status","complete_by_reference","genes_found","missing_genes","missing_rrnas",
                 "ambiguous_bp","detected_gene_count","detected_cds_count","detected_rrna_count",
                 "accepted","strict_complete_candidate","acceptance_policy","acceptance_checks",
                 "acceptance_failures","acceptance_reason","homology_sequence_mode","terminal_overlap_removed_bp","source_path","fasta"])
    write_table(root/"stages.tsv", runner.events,
                ["stage","status","reason","returncode","elapsed_seconds","log","command"])
    return summary


def parser(require_pipeline=False):
    p = argparse.ArgumentParser(description="Run ITSME, MitoLooper, or both with automatic Conda environment switching." if require_pipeline else __doc__,
                                prog="pipeline_controller.sh" if require_pipeline else None)
    p.add_argument("input_dir", type=Path, nargs="?", help="Paired FASTQ directory (including W*_1.fq.gz)")
    p.add_argument("-i","--input-dir",type=Path,dest="input_dir_option",help="Alternative to the positional FASTQ directory")
    p.add_argument("--pipeline",choices=["itsme","mitolooper","both"],required=require_pipeline,
                   help="Run nuclear rDNA only, mitochondrial recovery only, or both in workflow order")
    p.add_argument("--version", action="version", version=VERSION)
    p.add_argument("--samples", type=Path, help="TSV with sample,target_taxon; optional mt_code,reference_gb,seed_fasta overrides")
    p.add_argument("--genbank-dir", type=Path, default=Path.home()/"NOAA/mitos_genbank",
                   help="Local annotated GenBank collection for automatic reference selection")
    p.add_argument("--reference-match", choices=["exact","genus","hierarchy"], default="hierarchy",
                   help="Prefer target, then genus and supplied higher ranks; default hierarchy stops at family")
    p.add_argument("--reference-min-rank", choices=["family","order","class","phylum"], default="family",
                   help="Broadest allowed hierarchy fallback; does not change the expected target for BOLD")
    p.add_argument("--reference-index", type=Path, help="Optional SQLite index path; default: genbank-dir/tortuga_reference_index.sqlite")
    p.add_argument("--make-sample-sheet", type=Path, help="Write a blank metadata TSV and exit")
    p.add_argument("--only-sample", help="Run just one sample from the sheet")
    p.add_argument("-o", "--outdir", type=Path,help="Output directory; defaults to itsme_results, mito_results, or both_results for the selected pipeline")
    p.add_argument("-t", "--threads", type=int, default=24)
    p.add_argument("-m", "--memory", type=int, default=128, help="Assembler memory GB; not a global process limit")
    p.add_argument("--config", type=Path, help="Optional explicit command mapping; otherwise environments are selected automatically")
    p.add_argument("--env-root",type=Path,default=Path(os.environ.get("MITO_ENV_ROOT","/home/ark/miniconda3/envs")),
                   help="Root containing existing Conda environments")
    p.add_argument("--conda",type=Path,help="Optional absolute Conda executable override")
    p.add_argument("--itsme-script",type=Path,default=Path.home()/"bin/ITSME/itsme.v2.sh")
    p.add_argument("--branch", choices=["mito","full"], default=None,help="Legacy mode flag; prefer --pipeline on the controller")
    p.add_argument("--itsme-db", type=Path)
    p.add_argument("--sensitivity","--itsme-sensitivity",type=int,choices=[1,2,3],default=2,
                   dest="itsme_sensitivity",help="ITSME recruitment preset (default 2)")
    p.add_argument("--bold-fasta", type=Path, help="Local BOLD-derived COI reference FASTA, including alternative taxa")
    p.add_argument("--bold-metadata", type=Path, help="Matching BOLD provenance/taxonomy TSV")
    p.add_argument("--mitos-refdir", type=Path)
    p.add_argument("--mitos-refseqver", help="MITOS2 reference release subdirectory")
    p.add_argument("--no-assembly", action="store_true")
    p.add_argument("--skip-wgs-if-accepted",action="store_true",
                   help="Skip separate WGS when read recovery already meets acceptance policy; ITSME still runs")
    p.add_argument("--assembler", choices=["metaspades","megahit"], default="metaspades")
    p.add_argument("--assembly-percent", type=float, default=100)
    p.add_argument("--sampling-seed", type=int, default=17)
    p.add_argument("--no-qc", action="store_true", help="Use supplied reads without fastp QC")
    p.add_argument("--skip-method", choices=["novoplasty","getorganelle","getorganelle_assembly","mitofinder","mtgrasp"],
                   action="append", default=[])
    p.add_argument("--resume", "--skip-existing", action="store_true", dest="resume")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--keep-intermediates",action="store_true",
                   help="Keep large stage work files for debugging/resume; default exports essentials and deletes intermediates")
    p.add_argument("--method-timeout", type=float, default=60, help="Minutes per mitochondrial recovery method")
    p.add_argument("--assembly-timeout", type=float, default=120, help="Minutes for WGS process group")
    p.add_argument("--annotation-timeout", type=float, default=120)
    p.add_argument("--validation-timeout", type=float, default=20)
    p.add_argument("--itsme-timeout", type=float, default=180)
    p.add_argument("--getorg-rounds", type=int, default=15)
    p.add_argument("--disentangle-seconds", type=int, default=300)
    p.add_argument("--novo-k", type=int, default=33)
    p.add_argument("--insert-size", type=int, default=300)
    p.add_argument("--min-genome-bp", type=int, default=10000)
    p.add_argument("--max-genome-bp", type=int, default=30000)
    p.add_argument("--min-candidate-bp", type=int, default=500)
    p.add_argument("--max-candidates", type=int, default=100)
    p.add_argument("--accept-topology",choices=["circular","any"],default="circular",
                   help="Require circularity, or allow candidates without established circularity")
    p.add_argument("--accept-rule",choices=["all-genes","genes","length","either","both"],default="all-genes",
                   help="Stop mito recovery on all reference genes, gene count, length, either threshold or both thresholds")
    p.add_argument("--min-detected-genes",type=int,
                   help="Unique qualifying reference CDS+rRNA count for genes/either/both; excludes tRNAs")
    p.add_argument("--min-accepted-bp",type=int,
                   help="Contig length threshold for all-genes/length/either/both (default --min-genome-bp); independent of reference filtering")
    p.add_argument("--protein-identity", type=float, default=30)
    p.add_argument("--rrna-identity", type=float, default=70)
    p.add_argument("--gene-coverage", type=float, default=.8)
    p.add_argument("--coi-identity", type=float, default=97)
    p.add_argument("--coi-coverage", type=float, default=.8, help="Aligned fraction of shorter COI/barcode sequence")
    p.add_argument("--coi-min-bp", type=int, default=400)
    p.add_argument("--bold-near-fraction", type=float, default=.02)
    p.add_argument("--bold-max-hits", type=int, default=500)
    p.add_argument("--min-junction-templates", type=int, default=3)
    return p


def resolve_commands(args, tools=None):
    tools=list(DEFAULT_COMMANDS if tools is None else tools)
    commands={tool:list(DEFAULT_COMMANDS[tool]) for tool in tools}
    rows=[]
    if args.config:
        config=json.loads(args.config.read_text())
        for key,value in config.get("commands",{}).items():
            if key not in DEFAULT_COMMANDS or not isinstance(value,list) or not value or not all(isinstance(x,str) for x in value):
                raise ValueError(f"Invalid command prefix: {key}")
            if key in commands:
                commands[key]=[os.path.expanduser(os.path.expandvars(x)) for x in value]
    else:
        from configure_belly import runtime_commands
        commands,rows=runtime_commands(args.env_root,args.conda,args.itsme_script,tools,args.branch)
    return commands,rows


def finalize_outputs(root,args,result):
    if args.dry_run:
        return result
    from compact_outputs import compact_sample
    try:
        compact_sample(root,remove=not args.keep_intermediates)
        return json.loads((root/"final/summary.json").read_text())
    except (ValueError,OSError) as error:
        print(f"ERROR [{root.name}, cleanup]: {error}; inspect retained stage files",file=sys.stderr,flush=True)
        result["cleanup_error"]=str(error)
        atomic_json(root/"final/summary.json",result)
        return result


def run_itsme_only(args,pairs):
    """Standalone rDNA batches need no mitochondrial metadata or databases."""
    if not args.itsme_db or not args.itsme_db.expanduser().is_dir():
        raise ValueError("ITSME requires --itsme-db pointing to the existing rRNA database root")
    if args.threads<1 or args.memory<1 or args.itsme_timeout<=0:
        raise ValueError("Threads, memory and ITSME timeout must be positive")
    args.itsme_db=args.itsme_db.expanduser().resolve()
    selected=dict(pairs)
    if args.samples:
        sheet=args.samples.expanduser().resolve()
        with sheet.open(newline="") as f:
            reader=csv.DictReader(f,delimiter="\t")
            if "sample" not in (reader.fieldnames or []):
                raise ValueError("Optional ITSME sample sheet needs a sample column")
            names=[row["sample"].strip() for row in reader]
        if len(names)!=len(set(names)) or any(name not in pairs for name in names):
            raise ValueError("Duplicate sample name or missing read pair in sample sheet")
        selected={name:pairs[name] for name in names}
    if args.only_sample:
        if args.only_sample not in selected:
            raise ValueError(f"Sample not found: {args.only_sample}")
        selected={args.only_sample:selected[args.only_sample]}
    if not selected:
        raise ValueError("No samples selected")
    if any(not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*",name) for name in selected):
        raise ValueError("Sample names must begin with a letter/digit and contain only letters, digits, _, . or -")
    commands,environment_rows=resolve_commands(args,["itsme"])
    import fcntl
    args.outdir.mkdir(parents=True,exist_ok=True)
    with (args.outdir/".pipeline.lock").open("w") as lock:
        try:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            raise ValueError("Another pipeline is writing this output directory")
        if environment_rows:
            atomic_json(args.outdir/"environment_selection.json",environment_rows)
            write_table(args.outdir/"environment_selection.tsv",environment_rows,
                        ["tool","status","environment","entrypoint","alternatives","note"])
        results=[]
        for name,(r1,r2) in selected.items():
            print(f"\n=== {name}: ITSME ===",flush=True)
            root=args.outdir/name
            root.mkdir(parents=True,exist_ok=True)
            try:
                payload={"version":VERSION,"pipeline":"itsme","sample":name,
                         "reads":[stamp(r1),stamp(r2)],"itsme_db":str(args.itsme_db),
                         "commands":commands,"threads":args.threads,"memory":args.memory,
                         "timeout":args.itsme_timeout,"sensitivity":args.itsme_sensitivity}
                checkpoint=root/"run_config.json"
                if checkpoint.exists():
                    if not args.resume:
                        raise ValueError(f"{root} already has a run; use --resume or a new output directory")
                    if json.loads(checkpoint.read_text())!=payload:
                        raise ValueError("Inputs/settings/pipeline changed; use a new output directory")
                    compacted=root/"cleanup.json"
                    if compacted.exists() and json.loads(compacted.read_text()).get("compacted"):
                        from compact_outputs import validate_retained
                        validate_retained(root)
                        results.append(json.loads((root/"final/summary.json").read_text()))
                        print("  outputs: CACHED finalized compact results",flush=True)
                        continue
                elif any(root.iterdir()):
                    raise ValueError(f"{root} is nonempty without a matching checkpoint")
                atomic_json(checkpoint,payload)
                runner=Runner(root,commands,args.resume,args.dry_run)
                work=attempt_dir(root,runner,"itsme",root/"itsme")
                output=work/"output"
                ok=runner.run("itsme","itsme",
                              ["-1",r1,"-2",r2,"--db-dir",args.itsme_db,"-o",output,
                               "-t",args.threads,"-m",args.memory,"--sensitivity",args.itsme_sensitivity],
                              work,args.itsme_timeout*60,[output/"run_summary.txt"])
                result={"sample":name,"pipeline":"itsme","status":"planned" if args.dry_run else "completed" if ok else "failed",
                        "itsme_ok":ok,"itsme_output":str(output),"failed_stages":[] if ok else ["itsme"]}
                atomic_json(root/"stages.json",runner.events)
                write_table(root/"stages.tsv",runner.events,["stage","status","reason","returncode","elapsed_seconds","log","command"])
                atomic_json(root/"final/summary.json",result)
                result=finalize_outputs(root,args,result)
            except (ValueError,RuntimeError,OSError) as error:
                print(f"ERROR [{name}]: {error}",file=sys.stderr)
                result={"sample":name,"pipeline":"itsme","status":"failed","error":str(error)}
            results.append(result)
        previous={p.parent.parent.name:json.loads(p.read_text()) for p in args.outdir.glob("*/final/summary.json")}
        previous.update({r["sample"]:r for r in results})
        write_table(args.outdir/"cohort_summary.csv",list(previous.values()),
                    ["sample","pipeline","status","itsme_ok","itsme_output","final_directory","failed_stages","error","cleanup_error"],",")
        atomic_json(args.outdir/"last_invocation.json",{"version":VERSION,"pipeline":"itsme","results":results})
        return 1 if any(r["status"]=="failed" or r.get("cleanup_error") for r in results) else 0


def resolve_sample_reference(row, args, pair, catalog):
    """Resolve one sample independently so reference errors do not abort a cohort."""
    from Bio import SeqIO
    from Bio.Data import CodonTable
    from reference_selection import select_reference, taxonomy_from_row
    name, target = row["sample"].strip(), row["target_taxon"].strip()
    lineage = taxonomy_from_row(row)
    def refpath(value):
        path = Path(os.path.expandvars(value)).expanduser()
        path = (args.samples.parent/path).resolve() if not path.is_absolute() else path.resolve()
        if not path.is_file():
            raise ValueError(f"Missing reference: {path}")
        return path
    code_text = (row.get("mt_code") or "").strip()
    code = int(code_text) if code_text else None
    if code is not None and code not in CodonTable.unambiguous_dna_by_id:
        raise ValueError(f"Invalid genetic code for {name}")
    reference_text = (row.get("reference_gb") or "").strip()
    if reference_text:
        reference = refpath(reference_text)
        record = SeqIO.read(reference,"genbank")
        if code is None:
            codes = {int(x) for feature in record.features if feature.type=="CDS"
                     for x in feature.qualifiers.get("transl_table", [])}
            if len(codes)!=1:
                raise ValueError(f"{name}: reference does not specify a unique genetic code; fill mt_code")
            code = next(iter(codes))
        reference_data(record,code)
        choice = {"accession":record.id,"organism":record.annotations.get("organism", ""),
                  "match_level":"manual_override","match_rank":"manual","matched_taxon":"",
                  "eligible_candidates":1,"annotation_warnings":reference_boundary_warnings(record)}
    else:
        choice = select_reference(catalog,target,code,args.reference_match,
                                  args.min_genome_bp,args.max_genome_bp,lineage,args.reference_min_rank)
        reference = Path(choice["path"])
        if code is None:
            code = choice["mt_code"]
    if code not in CodonTable.unambiguous_dna_by_id:
        raise ValueError(f"Invalid inferred genetic code for {name}")
    seed = refpath(row["seed_fasta"]) if row.get("seed_fasta") else None
    sample = Sample(name,*pair,target,code,reference,seed,lineage)
    report = {"sample":name,"target_taxon":target,"reference_gb":str(reference),
              "reference_accession":choice["accession"],"reference_taxon":choice["organism"],
              "match_level":choice["match_level"],"match_rank":choice["match_rank"],
              "matched_taxon":choice["matched_taxon"],"mt_code":code,
              "code_source":"sample_sheet" if code_text else "GenBank_CDS",
              "eligible_candidates":choice["eligible_candidates"],
              "annotation_warnings":choice.get("annotation_warnings",[]),
              "selection_status":"available","reference_error":""}
    resolved = {**row,"sample":name,"target_taxon":target,"mt_code":code,
                "reference_gb":str(reference),"seed_fasta":str(seed) if seed else ""}
    return sample, report, resolved


def main(argv=None,require_pipeline=False):
    cli=parser(require_pipeline=require_pipeline)
    args = cli.parse_args(argv)
    if require_pipeline and not args.pipeline:
        cli.error("--pipeline is required: choose itsme, mitolooper, or both")
    if args.input_dir and args.input_dir_option:
        cli.error("Provide the FASTQ directory either positionally or with --input-dir, not both")
    args.input_dir=args.input_dir or args.input_dir_option
    del args.input_dir_option
    if not args.input_dir:
        cli.error("A FASTQ directory is required (positional or --input-dir)")
    if args.pipeline:
        branch="full" if args.pipeline in {"itsme","both"} else "mito"
        if args.branch and (args.pipeline=="itsme" or args.branch!=branch):
            cli.error("--pipeline conflicts with --branch; use --pipeline alone")
        args.branch=branch
        if args.pipeline in {"itsme","both"} and args.itsme_db is None:
            args.itsme_db=Path(os.environ.get("ITSME_DB","/home/ark/databases/itsme_db"))
    else:
        args.branch=args.branch or "mito"
    args.input_dir = args.input_dir.expanduser().resolve()
    default_output={"itsme":"itsme_results","mitolooper":"mito_results","both":"both_results"}.get(args.pipeline,"mito_results")
    args.outdir = (args.outdir or Path(default_output)).expanduser().resolve()
    if args.dry_run and args.outdir.exists() and any(args.outdir.iterdir()):
        raise ValueError("Dry runs require a new/empty output directory to protect existing results")
    pairs = discover_pairs(args.input_dir)
    if args.make_sample_sheet:
        if args.make_sample_sheet.exists():
            raise ValueError("Sample sheet already exists; refusing to overwrite")
        write_table(args.make_sample_sheet, [{"sample": s} for s in pairs],
                    ["sample","target_taxon","mt_code","reference_gb","seed_fasta",
                     "kingdom","phylum","class","order","family","genus","species"])
        print(f"Wrote {args.make_sample_sheet}; fill metadata before running")
        return 0
    if args.pipeline=="itsme":
        return run_itsme_only(args,pairs)
    for key in ["samples","bold_fasta","bold_metadata","mitos_refdir","mitos_refseqver"]:
        if not getattr(args,key):
            raise ValueError(f"--{key.replace('_','-')} is required")
    for key in ["samples","bold_fasta","bold_metadata","mitos_refdir"]:
        setattr(args,key,getattr(args,key).expanduser().resolve())
        if not getattr(args,key).exists():
            raise ValueError(f"Missing {key}: {getattr(args,key)}")
    if not (args.mitos_refdir/args.mitos_refseqver).is_dir():
        raise ValueError("MITOS2 reference release not found under --mitos-refdir")
    if args.branch == "full":
        if not args.itsme_db or not args.itsme_db.expanduser().is_dir():
            raise ValueError("--branch full requires --itsme-db")
        args.itsme_db = args.itsme_db.expanduser().resolve()
    if args.threads < 1 or args.memory < 1 or not 0 < args.assembly_percent <= 100:
        raise ValueError("Invalid threads, memory or assembly percentage")
    validate_acceptance_options(args)
    if any(getattr(args,k) <= 0 for k in ["method_timeout","assembly_timeout","annotation_timeout","validation_timeout","itsme_timeout"]):
        raise ValueError("Timeouts must be positive")
    if not 0 < args.gene_coverage <= 1 or not 0 < args.coi_coverage <= 1 or not 0 <= args.bold_near_fraction < 1:
        raise ValueError("Invalid coverage or near-best fraction")
    commands,environment_rows=resolve_commands(args)
    samples, reference_choices, resolved_samples, reference_failures = [], [], [], []
    catalog = None
    with args.samples.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        input_columns = list(reader.fieldnames or [])
        if not {"sample","target_taxon"}.issubset(input_columns):
            raise ValueError("Sample sheet needs sample and target_taxon columns")
        rows = list(reader)
    seen = set()
    selected = []
    for row in rows:
        name = (row.get("sample") or "").strip()
        if name in seen or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", name):
            raise ValueError(f"Duplicate or unsafe sample name: {name}")
        seen.add(name)
        if args.only_sample and name != args.only_sample:
            continue
        if name not in pairs or not (row.get("target_taxon") or "").strip():
            raise ValueError(f"Missing FASTQs or target_taxon for {name}")
        selected.append(row)
    if not selected:
        raise ValueError("No samples selected")
    if any(not (row.get("reference_gb") or "").strip() for row in selected):
        from reference_selection import load_catalog
        catalog = load_catalog(args.genbank_dir, reference_data, args.reference_index)
    for row in selected:
        name = row["sample"].strip()
        try:
            sample, report, resolved = resolve_sample_reference(row,args,pairs[name],catalog)
        except (ValueError,OSError) as error:
            print(f"ERROR [{name}, reference]: {error}; continuing other samples",file=sys.stderr,flush=True)
            failure = {"sample":name,"target_taxon":row["target_taxon"].strip(),
                       "selection_status":"unavailable","reference_error":str(error)}
            reference_failures.append(failure)
            reference_choices.append(failure)
            continue
        samples.append(sample)
        reference_choices.append(report)
        resolved_samples.append(resolved)
        print(f"Reference [{name}]: {report['reference_accession']} ({report['reference_taxon']}); "
              f"{report['match_level']} [{report['matched_taxon']}]; code {report['mt_code']}",flush=True)
        for warning in report["annotation_warnings"]:
            print(f"  Reference warning [{name}]: {warning['gene']} {warning['location']}: {warning['reason']}",flush=True)
    # A scoped output lock prevents two workers writing the shared BOLD index or cohort files.
    import fcntl
    args.outdir.mkdir(parents=True,exist_ok=True)
    with (args.outdir/".pipeline.lock").open("w") as lock:
        try:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            raise ValueError("Another pipeline is writing this output directory")
        write_table(args.outdir/"reference_selection.tsv",reference_choices,
                    ["sample","target_taxon","reference_gb","reference_accession","reference_taxon",
                     "match_level","match_rank","matched_taxon","mt_code","code_source","eligible_candidates",
                     "selection_status","annotation_warnings","reference_error"])
        write_table(args.outdir/"reference_failures.tsv",reference_failures,
                    ["sample","target_taxon","selection_status","reference_error"])
        atomic_json(args.outdir/"acceptance_policy.json",acceptance_policy(args))
        core_columns = ["sample","target_taxon","mt_code","reference_gb","seed_fasta"]
        write_table(args.outdir/"resolved_samples.tsv",resolved_samples,
                    core_columns + [column for column in input_columns if column not in core_columns])
        if environment_rows:
            atomic_json(args.outdir/"environment_selection.json",environment_rows)
            write_table(args.outdir/"environment_selection.tsv",environment_rows,
                        ["tool","status","environment","entrypoint","alternatives","note"])
        results = []
        # ITSME does not depend on a mitochondrial reference. Keep its checkpoint
        # tree separate so reference failures can be repaired without overwriting it.
        if reference_failures and args.branch == "full":
            fallback_args = argparse.Namespace(**vars(args))
            fallback_args.pipeline = "itsme"
            fallback_args.samples = None
            fallback_args.only_sample = None
            fallback_args.outdir = args.outdir/"unresolved_references_itsme"
            run_itsme_only(fallback_args,{r["sample"]:pairs[r["sample"]] for r in reference_failures})
            independent = json.loads((fallback_args.outdir/"last_invocation.json").read_text())
            itsme_results = {r["sample"]:r for r in independent["results"]}
        else:
            itsme_results = {}
        for failure in reference_failures:
            result = {"sample":failure["sample"],"target_taxon":failure["target_taxon"],
                      "pipeline":args.pipeline or ("both" if args.branch=="full" else "mitolooper"),
                      "status":"reference_unavailable","error":failure["reference_error"],
                      "failed_stages":["reference_selection"]}
            independent = itsme_results.get(failure["sample"])
            if independent:
                result.update(itsme_ok=independent.get("itsme_ok",False),
                              itsme_output=independent.get("itsme_output",""))
                if independent["status"] == "failed":
                    result["failed_stages"].append("itsme")
                    result["error"] += "; ITSME: " + independent.get("error","execution failed")
                if not args.dry_run:
                    from compact_outputs import promote_independent_final
                    result=promote_independent_final(args.outdir/"unresolved_references_itsme"/failure["sample"],
                                                     args.outdir/failure["sample"],result)
            results.append(result)
        pending=[]
        for sample in samples:
            compacted=args.outdir/sample.sample/"cleanup.json"
            if args.resume and compacted.exists() and json.loads(compacted.read_text()).get("compacted"):
                try:
                    # Finalized samples need neither the copied BOLD index nor stage work files.
                    results.append(process_sample(sample,args,commands,None,None))
                except (ValueError,RuntimeError,OSError) as error:
                    print(f"ERROR [{sample.sample}]: {error}",file=sys.stderr)
                    results.append({"sample":sample.sample,"status":"failed","error":str(error)})
            else:
                pending.append(sample)
        samples=pending
        if samples:
            metadata = load_bold_metadata(args.bold_metadata)
            from Bio import SeqIO
            ids = set()
            for rec in SeqIO.parse(args.bold_fasta,"fasta"):
                if rec.id in ids or rec.id not in metadata:
                    raise ValueError(f"Duplicate BOLD sequence ID or missing metadata: {rec.id}")
                ids.add(rec.id)
            if not ids:
                raise ValueError("Empty BOLD reference FASTA")
            dbdir = args.outdir/"bold_reference"
            dbdir.mkdir(exist_ok=True)
            bdb = dbdir/"coi"
            brunner = Runner(args.outdir,commands,args.resume,args.dry_run)
            bold_stamp = {"fasta":stamp(args.bold_fasta),"metadata_sha256":digest_file(args.bold_metadata)}
            if (dbdir/"source.json").exists() and json.loads((dbdir/"source.json").read_text()) != bold_stamp:
                raise ValueError("BOLD references changed; choose a new output directory")
            if not brunner.run("bold_index","makeblastdb",["-in",args.bold_fasta,"-dbtype","nucl","-parse_seqids","-out",bdb],
                               dbdir,args.validation_timeout*60,[dbdir/"coi.nsq"]):
                raise RuntimeError("Could not build BOLD BLAST database")
            atomic_json(dbdir/"source.json",bold_stamp)
            for sample in samples:
                print(f"\n=== {sample.sample} ===",flush=True)
                try:
                    result=process_sample(sample,args,commands,metadata,bdb)
                    results.append(finalize_outputs(args.outdir/sample.sample,args,result))
                except (ValueError,RuntimeError,OSError) as e:
                    print(f"ERROR [{sample.sample}]: {e}",file=sys.stderr)
                    results.append({"sample":sample.sample,"status":"failed","error":str(e)})
        # Rebuild cohort table from saved per-sample summaries, preserving --only-sample runs.
        previous = {p.parent.parent.name:json.loads(p.read_text()) for p in args.outdir.glob("*/final/summary.json")}
        previous.update({r["sample"]:r for r in results})
        write_table(args.outdir/"cohort_summary.csv",list(previous.values()),
                    ["sample","pipeline","status","target_taxon","mt_code","winning_method","best_method","best_length_bp",
                     "circular","target_status","genes_recovered","genes_expected","candidate_count",
                     "accepted_by_policy","strict_complete_candidate","complete_by_reference","acceptance_policy",
                     "itsme_ok","itsme_output","final_directory","wgs_ok","annotation_ok","failed_stages","error","cleanup_error"],",")
        gene_rows = []
        for name, row in previous.items():
            path = args.outdir/name/"final/gene_recovery.csv"
            if row["status"] not in {"failed","reference_unavailable"} and path.exists():
                with path.open(newline="") as f:
                    gene_rows.extend(csv.DictReader(f))
        write_table(args.outdir/"cohort_gene_recovery.csv",gene_rows,GENE_FIELDS,",")
        atomic_json(args.outdir/"last_invocation.json",{"version":VERSION,
                    "pipeline":args.pipeline or ("both" if args.branch=="full" else "mitolooper"),"results":results})
        if not args.keep_intermediates and not args.dry_run:
            from compact_outputs import cleanup_bold_index
            removed=cleanup_bold_index(args.outdir)
            atomic_json(args.outdir/"cleanup_summary.json",{"bold_index_bytes_removed":removed,
                        "sample_cleanup_reports":[r.get("cleanup_report",str(args.outdir/r["sample"]/"cleanup.json")) for r in results]})
            print(f"Removed {removed/1024**3:.2f} GiB of copied BOLD index",flush=True)
        return 1 if any(r["status"] in {"failed","reference_unavailable","recovered_annotation_failed","accepted_by_policy_annotation_failed"} or r.get("wgs_ok") is False
                       or r.get("itsme_ok") is False or r.get("cleanup_error") or (not r.get("winning_method") and r.get("failed_stages"))
                       for r in results) else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError,RuntimeError,OSError,KeyError) as error:
        print(f"ERROR: {error}",file=sys.stderr)
        sys.exit(1)

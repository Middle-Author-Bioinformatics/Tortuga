#!/usr/bin/env python3
"""Export concise Tortuga deliverables and remove owned, disposable stage files."""
from __future__ import annotations

import argparse
import csv
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import sys

STAGE_RE=re.compile(r"^(?:01[1-4]_(?:novoplasty|getorganelle|mitofinder|mtgrasp)|02_itsme|03_wgs|04[12]_(?:getorganelle_assembly|mitofinder_assembly)|05_mitos2|itsme)$")
SEQUENCES={".fa",".fasta",".fna",".faa",".ffn",".fas"}
GRAPHS={".gfa",".fastg",".gfa2"}
ANNOTATIONS={".csv",".tsv",".gff",".gff3",".gtf",".gb",".gbk",".gbff",".bed"}
LOGS={".log",".err"}
WGS_PRODUCTS={"contigs.fasta","scaffolds.fasta","final.contigs.fa","assembly_graph.fastg",
              "assembly_graph.gfa","assembly_graph_with_scaffolds.gfa","scaffolds.paths","contigs.paths"}


def read_json(path):
    return json.loads(Path(path).read_text())


def atomic_json(path,value):
    path=Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_name(path.name+".tmp")
    if path.is_symlink() or temporary.is_symlink():
        raise ValueError(f"Refusing linked metadata file: {path}")
    temporary.write_text(json.dumps(value,indent=2)+"\n")
    temporary.replace(path)


def table(path,rows,fields,delimiter="\t"):
    path=Path(path)
    temporary=path.with_name(path.name+".tmp")
    if path.is_symlink() or temporary.is_symlink():
        raise ValueError(f"Refusing linked metadata file: {path}")
    with temporary.open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter=delimiter,extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def owned_files(root):
    """Walk links themselves, never external targets or linked directories."""
    for base,dirs,names in os.walk(root,followlinks=False):
        for name in list(dirs):
            path=Path(base)/name
            if path.is_symlink():
                dirs.remove(name)
                yield path
        for name in names:
            yield Path(base)/name


def protected_inputs(root):
    config=read_json(root/"run_config.json")
    protected={Path(row["path"]).absolute() for row in config.get("reads",[]) if isinstance(row,dict) and row.get("path")}
    sample=config.get("sample")
    if isinstance(sample,dict):
        protected.update(Path(sample[key]).absolute() for key in ["r1","r2","reference_gb","seed_fasta"] if sample.get(key))
    return protected


def kind(path):
    name=path.name.lower()
    if name=="run_summary.txt":
        return "annotation"
    if name in {"fastp.json","fastp.html","sampling.json"}:
        return "qc_or_report"
    suffix=Path(name[:-3] if name.endswith(".gz") else name).suffix
    if suffix in SEQUENCES:
        return "sequence"
    if suffix in GRAPHS:
        return "graph"
    if suffix in ANNOTATIONS:
        return "annotation"
    if path.suffix.lower() in LOGS or name in {"run_summary.txt","params.txt","warnings.log","spades.log","config.txt"}:
        return "log"
    if path.suffix.lower() in {".html",".htm",".svg",".png",".pdf"}:
        return "qc_or_report"
    return ""


def products(stage):
    """Select final products, excluding k-mer/round/read-recruitment intermediates."""
    label="itsme" if stage.name in {"itsme","02_itsme"} else re.sub(r"^\d+_","",stage.name)
    for path in owned_files(stage):
        if not path.is_file():
            continue
        if path.is_symlink():
            try:
                path.resolve().relative_to(stage.resolve())
            except ValueError:
                continue  # Do not materialize links to raw inputs or external databases.
        rel=path.relative_to(stage)
        category=kind(path)
        if not category:
            continue
        if category in {"sequence","graph"} and path.stat().st_size==0:
            continue
        # Text logs are useful even after a failed assembler or annotation job.
        if category=="log":
            yield path,label,category
            continue
        if label=="itsme":
            # Native final products, root summaries, and final SPAdes graphs.
            # Do not export seed/recruitment/round FASTAs or intermediate BLASTs.
            in_final="final" in rel.parts
            output_root=path.parent.name=="output"
            final_graph=(path.name in WGS_PRODUCTS and path.parent.name=="spades")
            qc_report=category=="qc_or_report" and ("qc" in rel.parts or path.name.startswith(("fastp.","multiqc")))
            if in_final or output_root or final_graph or qc_report:
                yield path,label,category
        elif label=="wgs":
            if path.name in WGS_PRODUCTS and path.parent.name=="output":
                yield path,label,category or "graph"
            elif path.name=="sampling.json":
                yield path,label,category
        elif label=="mitos2":
            if "output" in rel.parts:
                yield path,"mitos",category
        else:
            # Terminal mitochondrial FASTAs/graphs; exclude nested assembler K* files.
            if any(re.fullmatch(r"K\d+",part) or part.lower() in {"tmp","temp","intermediate_contigs"} for part in rel.parts[:-1]):
                continue
            if category=="graph" or (category=="sequence" and
                (path.name.endswith(".path_sequence.fasta") or
                 path.name.startswith(("Circularized_assembly","Contigs","Option")) or
                 path.name.endswith(("_mtDNA_contig.fasta",".final-mtgrasp-assembly.fa",".postsealer.postpilon.fasta")))):
                yield path,label,category


def same_bytes(a,b):
    if a.stat().st_size!=b.stat().st_size:
        return False
    with a.open("rb") as left,b.open("rb") as right:
        while True:
            block=left.read(1024*1024)
            other=right.read(1024*1024)
            if block!=other:
                return False
            if not block:
                return True


def export_file(source,directory,label,dry_run=False):
    """Install a real retained file before deleting any stage data; never overwrite."""
    safe=re.sub(r"[^A-Za-z0-9_.-]","_",source.name)
    prefix=label+"_" if label else ""
    target=directory/(prefix+safe)
    if target.is_symlink():
        raise ValueError(f"Refusing deliverable symlink: {target}")
    if target.exists() and not same_bytes(source,target):
        tag=hashlib.sha256(str(source).encode()).hexdigest()[:12]
        target=directory/(prefix+tag+"_"+safe)
        if target.is_symlink() or (target.exists() and not same_bytes(source,target)):
            raise ValueError(f"Deliverable collision: {target}")
    if not dry_run and not target.exists():
        directory.mkdir(parents=True,exist_ok=True)
        temp=target.with_name(target.name+".exporting")
        if temp.exists() or temp.is_symlink():
            raise ValueError(f"Unfinished export exists: {temp}")
        try:
            try:
                os.link(source.resolve(),temp)
            except OSError:
                shutil.copy2(source,temp)
                if not same_bytes(source,temp):
                    raise OSError(f"Export verification failed: {source}")
            temp.replace(target)
        finally:
            temp.unlink(missing_ok=True)
    return target


def rewrite(value,mapping):
    if isinstance(value,str):
        return mapping.get(value,value)
    if isinstance(value,list):
        return [rewrite(item,mapping) for item in value]
    if isinstance(value,dict):
        return {key:rewrite(item,mapping) for key,item in value.items()}
    return value


def annotation_table(source):
    """A readable native GFF feature table, without changing annotation calls."""
    from urllib.parse import unquote
    target=source.with_name(source.stem+"_features.tsv")
    if target.exists() or target.is_symlink():
        return
    rows=[]
    fields=["sequence_id","source","feature_type","start","end","score","strand","phase","gene","attributes"]
    with source.open() as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            parts=line.rstrip("\n").split("\t")
            if len(parts)!=9:
                continue
            attrs=dict(item.split("=",1) for item in parts[8].split(";") if "=" in item)
            gene=unquote(attrs.get("Name") or attrs.get("gene") or attrs.get("gene_id") or attrs.get("ID", ""))
            rows.append(dict(zip(fields,parts[:8]+[gene,parts[8]])))
    table(target,rows,fields)


def rewrite_reports(root,mapping,summary):
    candidates=root/"candidates"
    assessments=[p for p in owned_files(candidates) if p.name=="assessment.json"] if candidates.is_dir() else []
    for path in [root/"final/annotation_summary.json",*assessments]:
        if path.exists() and not path.is_symlink():
            atomic_json(path,rewrite(read_json(path),mapping))
    path=root/"candidates.tsv"
    if path.exists() and not path.is_symlink():
        with path.open(newline="") as handle:
            reader=csv.DictReader(handle,delimiter="\t")
            fields=reader.fieldnames
            rows=[{key:mapping.get(value,value) for key,value in row.items()} for row in reader]
        table(path,rows,fields)
    atomic_json(root/"final/summary.json",summary)


def compact_sample(root,remove=True,dry_run=False):
    root=Path(root).absolute()
    if root.is_symlink() or not (root/"run_config.json").is_file() or not (root/"final/summary.json").is_file():
        raise ValueError(f"Not a completed Tortuga sample: {root}")
    final=root/"final"
    for directory in [final,root/"logs",root/"qc",root/"candidates"]:
        if directory.is_symlink():
            raise ValueError(f"Refusing linked output directory: {directory}")
    previous=read_json(root/"cleanup.json") if (root/"cleanup.json").exists() else {}
    if previous.get("compacted"):
        validate_retained(root)
        return previous
    protected=protected_inputs(root)
    stages=sorted(p for p in root.iterdir() if STAGE_RE.fullmatch(p.name) and p.is_dir() and not p.is_symlink())
    selected=[]
    for stage in stages:
        selected.extend(products(stage))
    # Keep candidate sequences but remove doubled subjects, junction mappings and temporary queries.
    candidates=root/"candidates"
    if candidates.is_dir():
        selected.extend((path,"mito_"+path.parent.name,"sequence") for path in owned_files(candidates) if path.name=="sequence.fasta" and not path.is_symlink())
    old_manifest=final/"outputs_manifest.tsv"
    manifest=[]
    if old_manifest.exists():
        with old_manifest.open(newline="") as handle:
            manifest=list(csv.DictReader(handle,delimiter="\t"))
    mapping={}
    for source,label,category in sorted(selected,key=lambda row:str(row[0])):
        destination=root/"logs" if category=="log" else final
        target=export_file(source,destination,label,dry_run)
        mapping[str(source)]=str(target)
        manifest.append({"source":str(source),"retained_file":str(target),"kind":category,"bytes":source.stat().st_size})
        if not dry_run and label=="mitos" and target.suffix.lower() in {".gff",".gff3"}:
            annotation_table(target)
    # Cleanup is limited to generated stage trees and explicitly known scratch files.
    disposable=[]
    for stage in stages:
        disposable.extend(owned_files(stage))
    for path in owned_files(root/"qc") if (root/"qc").is_dir() else []:
        if re.search(r"\.(?:fastq|fq)(?:\.gz)?$",path.name,re.I):
            disposable.append(path)
    for path in owned_files(candidates) if candidates.is_dir() else []:
        if path.name in {"sequence.fasta","screen.fasta","doubled.fasta","junction_reference.fasta","coi.fasta","junction.sam"}:
            disposable.append(path)
    disposable=sorted(set(p for p in disposable if p.absolute() not in protected),key=str)
    # Never erase protected inputs, even if a user placed one inside a stage tree.
    report={"sample":root.name,"compacted":bool(remove),"dry_run":dry_run,
            "files_removed":len(disposable) if remove else 0,
            "bytes_removed":sum(p.lstat().st_size for p in disposable if str(p) not in mapping) if remove else 0,
            "exported_files":len({row["retained_file"] for row in manifest}),
            "final_directory":str(final),"resume_note":"Finalized compact results are reused as a whole; removed stage intermediates cannot be resumed. Retry analyses in a new output directory."}
    if dry_run:
        return report
    # Export manifests and repaired report paths first. Any export failure leaves stage files intact.
    unique={row["source"]:row for row in manifest}
    table(old_manifest,list(unique.values()),["source","retained_file","kind","bytes"])
    summary=rewrite(read_json(final/"summary.json"),mapping)
    summary["final_directory"]=str(final)
    summary["output_retention"]="compact" if remove else "all"
    summary["cleanup_report"]=str(root/"cleanup.json")
    itsme_files=[row["retained_file"] for row in manifest if Path(row["retained_file"]).name.startswith("itsme_") and row["kind"]!="log"]
    if itsme_files:
        summary["itsme_work_output"]=summary.get("itsme_work_output",summary.get("itsme_output",""))
        summary["itsme_output"]=str(final)
        summary["itsme_files"]=sorted(set(itsme_files))
    if summary.get("annotation_dir"):
        summary.setdefault("annotation_work_dir",summary["annotation_dir"])
        summary["annotation_dir"]=str(final)
    rewrite_reports(root,mapping,summary)
    # Include the wrapper's own final sequences/tables in the retained-output check.
    retained={row["retained_file"] for row in unique.values()}
    for path in final.iterdir():
        if path.is_file() and not path.is_symlink() and path.name not in {"summary.json","outputs_manifest.tsv"} and str(path) not in retained:
            unique[str(path)]={"source":str(path),"retained_file":str(path),"kind":kind(path) or "metadata","bytes":path.stat().st_size}
    table(old_manifest,list(unique.values()),["source","retained_file","kind","bytes"])
    deletion_rows=[]
    if remove:
        for path in disposable:
            deletion_rows.append({"path":str(path),"bytes":path.lstat().st_size,"reason":"exported product or disposable stage intermediate"})
            path.unlink()
        for stage in stages:
            for base,dirs,names in os.walk(stage,topdown=False,followlinks=False):
                path=Path(base)
                if not any(path.iterdir()):
                    path.rmdir()
    table(root/"cleanup_removed.tsv",deletion_rows,["path","bytes","reason"])
    atomic_json(root/"cleanup.json",report)
    print(f"  outputs: {final}; removed {report['bytes_removed']/1024**3:.2f} GiB of stage files",flush=True)
    return report


def validate_retained(root):
    manifest=Path(root)/"final/outputs_manifest.tsv"
    if not manifest.is_file():
        raise ValueError(f"Missing retained-output manifest: {manifest}")
    with manifest.open(newline="") as handle:
        for row in csv.DictReader(handle,delimiter="\t"):
            path=Path(row["retained_file"])
            if not path.is_file() or path.is_symlink() or path.stat().st_size!=int(row["bytes"]):
                raise ValueError(f"Retained deliverable missing or changed: {path}")


def cleanup_bold_index(root,dry_run=False):
    directory=Path(root)/"bold_reference"
    if not directory.is_dir() or directory.is_symlink():
        return 0
    removed=0
    extensions={"nhr","nin","nsq","ndb","not","ntf","nto","nal","nsd","nsi","nog"}
    for path in directory.iterdir():
        if path.name.startswith("coi.") and path.suffix[1:] in extensions and (path.is_file() or path.is_symlink()):
            removed+=path.lstat().st_size
            if not dry_run:
                path.unlink()
    return removed


def promote_independent_final(source,destination,result):
    """Place reference-unavailable ITSME products in the ordinary parent sample final."""
    source,destination=Path(source),Path(destination)
    final=destination/"final"
    if destination.is_symlink() or final.is_symlink():
        raise ValueError(f"Refusing linked sample final directory: {destination}")
    final.mkdir(parents=True,exist_ok=True)
    mapping={}
    manifest=[]
    for path in (source/"final").iterdir():
        if path.is_file() and not path.is_symlink() and path.name not in {"summary.json","outputs_manifest.tsv"}:
            target=export_file(path,final,"")
            mapping[str(path)]=str(target)
            manifest.append({"source":str(path),"retained_file":str(target),"kind":kind(path) or "metadata","bytes":target.stat().st_size})
    independent=read_json(source/"final/summary.json")
    result={**result,"final_directory":str(final),"itsme_output":str(final),
            "itsme_files":[mapping.get(path,path) for path in independent.get("itsme_files",[])],
            "cleanup_report":independent.get("cleanup_report",""),
            "output_retention":independent.get("output_retention","")}
    atomic_json(final/"summary.json",result)
    table(final/"outputs_manifest.tsv",manifest,["source","retained_file","kind","bytes"])
    return result


def refresh_cohort(root):
    path=root/"cohort_summary.csv"
    summaries={p.parent.parent.name:read_json(p) for p in root.glob("*/final/summary.json")}
    if path.exists():
        with path.open(newline="") as handle:
            reader=csv.DictReader(handle)
            fields=list(reader.fieldnames or [])
            rows=list(reader)
        for field in ["itsme_output","final_directory"]:
            if field not in fields:
                fields.append(field)
        for row in rows:
            summary=summaries.get(row["sample"],{})
            for field in ["itsme_output","final_directory"]:
                if field in summary:
                    row[field]=summary[field]
        table(path,rows,fields,",")
    path=root/"last_invocation.json"
    if path.exists():
        invocation=read_json(path)
        invocation["results"]=[{**row,**summaries.get(row["sample"],{})} for row in invocation.get("results",[])]
        atomic_json(path,invocation)


def finished_sample(sample):
    """A fresh terminal sample report, not merely an existing final directory."""
    config=sample/"run_config.json"
    summary=sample/"final/summary.json"
    if not summary.is_file() or summary.is_symlink():
        return False
    result=read_json(summary)
    if result.get("status") not in {"completed","failed","recovered_annotated","recovered_annotation_failed",
            "accepted_by_policy_annotated","accepted_by_policy_annotation_failed","partial_target_confirmed",
            "unresolved_candidate","no_candidate","reference_unavailable"}:
        return False
    # A newer checkpoint can indicate an interrupted retry with a stale old summary.
    return summary.stat().st_mtime_ns>=config.stat().st_mtime_ns


def sample_directories(root):
    return sorted(path for path in root.iterdir() if path.is_dir() and not path.is_symlink()
                  and (path/"run_config.json").is_file() and not (path/"run_config.json").is_symlink())


def compact_run(root,dry_run=False):
    """Clean finished samples in one inactive cohort, preserving unfinished work."""
    root=Path(root).expanduser().absolute()
    if not root.exists():
        raise ValueError(f"Output directory does not exist: {root}")
    if root.is_symlink() or not root.is_dir():
        raise ValueError(f"Not a regular Tortuga output directory: {root}")
    marker=root/"last_invocation.json"
    fallback=root/"unresolved_references_itsme"
    children=sample_directories(root)
    fallback_children=sample_directories(fallback) if fallback.is_dir() and not fallback.is_symlink() else []
    if not marker.is_file() and not children and not fallback_children:
        raise ValueError(f"Not a recognized Tortuga output directory: {root}")
    reports=[]
    with (root/".pipeline.lock").open("a") as lock:
        try:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            raise ValueError("The pipeline is still running in this directory; nothing was cleaned. Run cleanup again after it stops")
        finished=[sample for sample in children if finished_sample(sample)]
        unfinished=[sample.name for sample in children if sample not in finished]
        expected=set()
        sheet=root/"resolved_samples.tsv"
        if sheet.is_file():
            with sheet.open(newline="") as handle:
                expected={row.get("sample","").strip() for row in csv.DictReader(handle,delimiter="\t")}
            expected.discard("")
        unstarted=sorted(expected-{sample.name for sample in children})
        if not marker.is_file():
            print("No batch completion record: inspecting individual completed samples; shared BOLD index will be retained",flush=True)
        for sample in finished:
            reports.append(compact_sample(sample,dry_run=dry_run))
        for name in unfinished:
            print(f"SKIPPED {name}: no fresh terminal sample summary; work files retained",flush=True)
        for name in unstarted:
            print(f"SKIPPED {name}: sample has not started; shared work retained",flush=True)
        fallback_pending=any(not finished_sample(sample) for sample in fallback_children)
        if fallback_children:
            reports.extend(compact_run(fallback,dry_run))
            if not dry_run:
                parent_results={row["sample"]:row for row in read_json(marker).get("results",[])} if marker.is_file() else {}
                for sample in fallback_children:
                    if not finished_sample(sample):
                        continue
                    destination=root/sample.name
                    if (destination/"run_config.json").exists():
                        continue
                    independent=read_json(sample/"final/summary.json")
                    summary={**independent,**parent_results.get(sample.name,{})}
                    if sample.name not in parent_results:
                        summary.update(pipeline="both",status="reference_unavailable",
                                       error="Mitochondrial reference unavailable; batch completion record absent",
                                       failed_stages=sorted(set(independent.get("failed_stages",[])+["reference_selection"])))
                    promote_independent_final(sample,destination,summary)
        batch_finished=marker.is_file() and not unfinished and not unstarted and not fallback_pending
        removed=cleanup_bold_index(root,dry_run) if batch_finished else 0
        if not dry_run:
            refresh_cohort(root)
            atomic_json(root/"cleanup_summary.json",{"samples":reports,"bold_index_bytes_removed":removed,
                        "batch_completion_record_present":marker.is_file(),"batch_finished":batch_finished,
                        "unfinished_samples":unfinished,"not_started_samples":unstarted,
                        "shared_bold_index_retained":not batch_finished})
        if batch_finished:
            print(f"{'Would remove' if dry_run else 'Removed'} {removed/1024**3:.2f} GiB of copied BOLD index",flush=True)
        else:
            print("Shared BOLD index retained: full batch completion is not established",flush=True)
        print(f"{'Would compact' if dry_run else 'Compacted'} {len(reports)} completed sample(s); unfinished samples were not cleaned",flush=True)
    return reports


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir",type=Path,help="Tortuga cohort output directory; finished samples only, while the controller is inactive")
    parser.add_argument("--dry-run",action="store_true",help="Show estimated cleanup without exporting or deleting")
    args=parser.parse_args(argv)
    compact_run(args.run_dir,args.dry_run)
    return 0


if __name__=="__main__":
    try:
        sys.exit(main())
    except (ValueError,OSError,KeyError) as error:
        print(f"ERROR: {error}",file=sys.stderr)
        sys.exit(1)

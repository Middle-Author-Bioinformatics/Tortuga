#!/usr/bin/env python3
"""Find existing mitochondrial tools and write a prefix-based pipeline config.

Read-only with respect to Conda environments and references: never installs,
deletes, activates, moves, or merges them. Writes only a config and audit reports.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import re
import shutil
import sys
import time

# Runtime choice order: the FIRST environment containing the executable wins.
# Other entries are relevant alternative Conda installations, kept here as notes.
# GetOrganelle: use getorganelle first; mito-getorganelle is another relevant package environment.
# mtGrasp: use mito-mtgrasp first; mtgrasp is another relevant package environment.
# SPAdes: use itsme first; mito-spades is another relevant package environment.
# BLAST+: use itsme if present; mito-blast is another relevant package environment.
# This order follows the first suitable choices in the user's supplied inventory.
SPECS = {
    "fastp": (["itsme"], ["fastp"]),
    "minimap2": (["itsme"], ["minimap2"]),
    "metaspades": (["itsme", "mito-spades"], ["spades.py", "spades"]),
    "blastn": (["itsme", "mito-blast"], ["blastn"]),
    "tblastn": (["itsme", "mito-blast"], ["tblastn"]),
    "makeblastdb": (["itsme", "mito-blast"], ["makeblastdb"]),
    "getorganelle": (["getorganelle", "mito-getorganelle"], ["get_organelle_from_reads.py"]),
    "getorganelle_assembly": (["getorganelle", "mito-getorganelle"], ["get_organelle_from_assembly.py"]),
    "mitofinder": (["mito-mitofinder"], ["mitofinder", "MitoFinder", "mitofinder.py"]),
    "mtgrasp": (["mito-mtgrasp", "mtgrasp"], ["mtgrasp.py", "mtgrasp"]),
    "mitos2": (["mito-mitos2", "mitos2"], ["runmitos.py", "runmitos"]),
    "megahit": (["mito-megahit", "itsme"], ["megahit"]),
}
MISSING_NOTE = "Executable not found; inspect environment/package paths before installing"


def natural_key(value):
    return [int(x) if x.isdigit() else x.lower() for x in re.split(r"(\d+)", str(value))]


def packages(prefix):
    result = []
    for file in sorted((prefix/"conda-meta").glob("*.json")):
        try:
            result.append(json.loads(file.read_text()))
        except (ValueError, OSError):
            pass
    return result


def executable(prefix, names):
    for name in names:
        path = prefix/"bin"/name
        if path.is_file() and os.access(path, os.X_OK):
            return [str(path)]
        if path.is_file() and path.suffix == ".py" and (prefix/"bin/python").is_file():
            return [str(prefix/"bin/python"), str(path)]
    return None


def novo_executable(prefix):
    # Prefer explicit Perl entry points, not similarly named config/example files.
    paths = list((prefix/"bin").glob("*NOVOPlasty*.pl")) + list((prefix/"bin").glob("*novoplasty*.pl"))
    for package in packages(prefix):
        if "novoplasty" in package.get("name", "").lower():
            for filename in package.get("files", []):
                path = prefix/filename
                if re.fullmatch(r"novoplasty[\d.]*\.pl", path.name, re.I):
                    paths.append(path)
    for path in sorted(set(paths), key=natural_key, reverse=True):
        if path.is_file():
            perl = prefix/"bin/perl"
            return [str(perl) if perl.is_file() else "perl", str(path)]
    return executable(prefix, ["NOVOPlasty", "novoplasty", "NOVOPlasty.pl", "novoplasty.pl"])


def locate(env_root, names, tool_names, special=None):
    alternatives = []
    selected = None
    for name in names:
        prefix = env_root/name
        command = special(prefix) if special else executable(prefix, tool_names)
        if command:
            alternatives.append(str(prefix))
            if selected is None:
                selected = (prefix, command)
    return selected, alternatives


def reference_inventory(roots, seconds=20, max_dirs=5000):
    """Bounded filename inventory only; no sequences or arbitrary config contents read."""
    started, visited = time.monotonic(), set()
    hits = []
    pruned = {"envs", "conda_envs", ".snakemake", ".git", "pkgs", "raw_reads", "node_modules"}
    truncated = False
    for root in roots:
        root = root.expanduser().resolve()
        if not root.is_dir():
            continue
        for base, dirs, files in os.walk(root, followlinks=False):
            base = Path(base)
            if base in visited:
                dirs[:] = []
                continue
            visited.add(base)
            depth = len(base.relative_to(root).parts)
            dirs[:] = sorted(d for d in dirs if d not in pruned and not d.startswith(".")) if depth < 6 else []
            if "auxinfo.json" in files:
                hits.append({"type":"MITOS_reference_candidate", "path":str(base),
                             "hint":f"--mitos-refdir {base.parent} --mitos-refseqver {base.name}"})
            for name in files:
                lower = name.lower()
                if re.search(r"\.(fasta|fa|fna|tsv|csv|json|zip)(\.gz)?$",lower) and "bold" in lower:
                    hits.append({"type":"BOLD_candidate", "path":str(base/name), "hint":"Check COI sequences, matching taxonomy metadata and BOLD provenance"})
                if lower.startswith("animal_mt") and lower.endswith((".fasta", ".fa", ".fasta.gz")):
                    hits.append({"type":"GetOrganelle_animal_mt", "path":str(base/name), "hint":"GetOrganelle seed/label resource; not a BOLD or MITOS database"})
                if lower.endswith((".gb", ".gbk", ".gbff", ".genbank")) and re.search(r"mito|mtref|reference|refseq", str(base/name), re.I):
                    hits.append({"type":"GenBank_reference_candidate", "path":str(base/name), "hint":"Check that this is one annotated mitochondrial genome appropriate for the sample"})
            if len(visited) >= max_dirs or time.monotonic()-started >= seconds:
                truncated = True
                break
        if truncated:
            break
    return {"roots":[str(x.expanduser().resolve()) for x in roots], "directories_visited":len(visited),
            "truncated":truncated,"candidates":hits,
            "note":"Filename candidates are not validated databases; no match does not prove absence."}


def build_config(env_root, conda, itsme_script, branch="mito"):
    commands, rows = {}, []
    for tool, (env_names, tool_names) in SPECS.items():
        selected, alternatives = locate(env_root, env_names, tool_names)
        if selected:
            prefix, native = selected
            commands[tool] = [str(conda), "run", "--no-capture-output", "-p", str(prefix)] + native
            rows.append({"tool":tool,"status":"found", "environment":str(prefix),
                         "entrypoint":" ".join(native),"alternatives":";".join(alternatives),
                         "note":"Executable found; run check_tools.py to probe launchability"})
        else:
            rows.append({"tool":tool,"status":"missing", "environment":";".join(str(env_root/n) for n in env_names),
                         "entrypoint":";".join(tool_names),"alternatives":"", "note":MISSING_NOTE})
    selected, alternatives = locate(env_root,["mito-novoplasty"],[],novo_executable)
    if selected:
        prefix,native=selected
        commands["novoplasty"]=[str(conda),"run","--no-capture-output","-p",str(prefix)]+native
        rows.append({"tool":"novoplasty","status":"found","environment":str(prefix),
                     "entrypoint":" ".join(native),"alternatives":";".join(alternatives),
                     "note":"Detected actual entry point; check version/config compatibility with help probe"})
    else:
        rows.append({"tool":"novoplasty","status":"missing","environment":str(env_root/"mito-novoplasty"),
                     "entrypoint":"NOVOPlasty*.pl or NOVOPlasty wrapper","alternatives":"","note":MISSING_NOTE})
    if itsme_script.is_file():
        commands["itsme"]=[str(conda),"run","--no-capture-output","-p",str(env_root/"itsme"),"bash",str(itsme_script)]
        rows.append({"tool":"itsme","status":"found","environment":str(env_root/"itsme"),
                     "entrypoint":str(itsme_script),"alternatives":"","note":"Only used with --branch full"})
    elif branch=="full":
        rows.append({"tool":"itsme","status":"missing","environment":str(env_root/"itsme"),
                     "entrypoint":str(itsme_script),"alternatives":"","note":"Provide --itsme-script"})
    return {"commands":commands},rows


def runtime_commands(env_root, conda, itsme_script, tools, branch="mito"):
    """Select stage environments in-process; no configuration generation required."""
    env_root=Path(env_root).expanduser().resolve()
    preferred=env_root.parent/"bin/conda"
    conda=conda or os.environ.get("MITO_CONDA_EXE") or (
        str(preferred) if preferred.is_file() else os.environ.get("CONDA_EXE") or shutil.which("conda"))
    if not conda or not Path(conda).expanduser().is_file():
        raise ValueError("Conda executable not found; set MITO_CONDA_EXE or use --conda /path/to/conda")
    conda=Path(conda).expanduser().resolve()
    config,rows=build_config(env_root,conda,Path(itsme_script).expanduser().resolve(),branch)
    commands={tool:config["commands"].get(tool,[]) for tool in tools}
    rows=[row for row in rows if row["tool"] in commands]
    for row in rows:
        if row["status"]=="found":
            alternatives=[p for p in row["alternatives"].split(";") if p and p!=row["environment"]]
            row["alternatives"]=";".join(alternatives)
            row["note"]="First available executable selected at startup"
            if alternatives:
                row["note"]+="; other relevant Conda installation(s): "+", ".join(alternatives)
            print(f"[environment] {row['tool']}: {row['environment']}"+
                  ("; alternative(s): "+", ".join(alternatives) if alternatives else ""),flush=True)
        elif row["tool"]!="itsme" or branch=="full":
            print(f"[environment] {row['tool']}: not found in checked environments",flush=True)
    return commands,rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--env-root",type=Path,default=Path("/home/ark/miniconda3/envs"))
    p.add_argument("--conda",type=Path)
    p.add_argument("--itsme-script",type=Path,default=Path.home()/"bin/ITSME/itsme.v2.sh")
    p.add_argument("--branch",choices=["mito","full"],default="mito")
    p.add_argument("--no-assembly",action="store_true")
    p.add_argument("--skip-method",action="append",default=[],choices=["novoplasty","getorganelle","getorganelle_assembly","mitofinder","mtgrasp"])
    p.add_argument("--config-out",type=Path,default=Path("config.belly.json"))
    p.add_argument("--report-prefix",type=Path,default=Path("belly_audit"))
    p.add_argument("--reference-root",type=Path,action="append")
    p.add_argument("--no-reference-scan",action="store_true")
    args=p.parse_args()
    args.env_root=args.env_root.expanduser().resolve()
    args.itsme_script=args.itsme_script.expanduser().resolve()
    preferred=args.env_root.parent/"bin/conda"
    conda=args.conda or (preferred if preferred.is_file() else Path(os.environ.get("CONDA_EXE") or shutil.which("conda") or ""))
    conda=conda.expanduser().resolve()
    if not conda.is_file():
        raise ValueError("Conda executable not found; specify --conda /absolute/path/to/conda")
    if args.config_out.exists():
        raise ValueError(f"{args.config_out} already exists; choose another --config-out to preserve it")
    config,rows=build_config(args.env_root,conda,args.itsme_script,args.branch)
    optional={"megahit"}
    if args.no_assembly:
        optional.update(["metaspades","getorganelle_assembly"])
    if args.branch=="mito":
        optional.add("itsme")
    optional.update(args.skip_method)
    if "getorganelle" in args.skip_method:
        # Read and assembly modes can still be independently enabled by the caller.
        optional.add("getorganelle")
    missing=[r["tool"] for r in rows if r["status"]=="missing" and r["tool"] not in optional]
    args.report_prefix.parent.mkdir(parents=True,exist_ok=True)
    with Path(str(args.report_prefix)+".tsv").open("w",newline="") as f:
        fields=["tool","status","environment","entrypoint","alternatives","note"]
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t")
        w.writeheader();w.writerows(rows)
    info={"conda":str(conda),"environment_root":str(args.env_root),"tools":rows,"missing_required":missing}
    if not args.no_reference_scan:
        roots=args.reference_root or [Path("/home/ark/databases"),Path.home()/".GetOrganelle",
                                      Path.home()/"bin/mitogenome_pipeline"]
        info["references"]=reference_inventory(roots)
    Path(str(args.report_prefix)+".json").write_text(json.dumps(info,indent=2)+"\n")
    for row in rows:
        print(f"{row['tool']:24s} {row['status']:8s} {row['environment']}")
    if missing:
        print("\nNo runnable config written. Resolve missing tools: "+", ".join(missing))
        print("Audit reports were saved; missing here means not found in the checked locations, not necessarily absent everywhere.")
        return 1
    args.config_out.parent.mkdir(parents=True,exist_ok=True)
    args.config_out.write_text(json.dumps(config,indent=2)+"\n")
    print(f"\nWrote {args.config_out}. No environments or databases were changed.")
    print("Next: python check_tools.py --config "+str(args.config_out)+(" --branch full" if args.branch=="full" else ""))
    print("Reference filename candidates, if any, are in "+str(args.report_prefix)+".json")
    return 0


if __name__=="__main__":
    try:
        sys.exit(main())
    except (OSError,ValueError) as error:
        print(f"ERROR: {error}",file=sys.stderr)
        sys.exit(1)

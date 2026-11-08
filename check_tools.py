#!/usr/bin/env python3
"""Probe configured command prefixes with help flags; does not run assemblies."""
import argparse
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
from mito_pipeline import DEFAULT_COMMANDS


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config",type=Path)
    p.add_argument("--env-root",type=Path,default=Path(os.environ.get("MITO_ENV_ROOT","/home/ark/miniconda3/envs")))
    p.add_argument("--conda",type=Path)
    p.add_argument("--itsme-script",type=Path,default=Path.home()/"bin/ITSME/itsme.v2.sh")
    p.add_argument("--branch",choices=["mito","full"],default="mito")
    p.add_argument("--no-assembly",action="store_true")
    p.add_argument("--no-qc",action="store_true")
    p.add_argument("-o","--out",type=Path,default=Path("tool_checks.json"))
    args=p.parse_args()
    commands=dict(DEFAULT_COMMANDS)
    if args.config:
        commands.update(json.loads(args.config.read_text()).get("commands",{}))
    else:
        from configure_belly import runtime_commands
        commands,_=runtime_commands(args.env_root,args.conda,args.itsme_script,commands,args.branch)
    skip=set()
    if args.branch=="mito": skip.add("itsme")
    if args.no_assembly: skip.update(["metaspades","megahit","getorganelle_assembly"])
    else: skip.add("megahit")
    if args.no_qc: skip.add("fastp")
    results=[]
    for tool,cmd in commands.items():
        if tool in skip: continue
        if not cmd:
            results.append({"tool":tool,"ok":False,"error":"Executable not found in checked Conda environments"})
            print(f"{tool}: NOT FOUND")
            continue
        cmd=[os.path.expanduser(os.path.expandvars(x)) for x in cmd]+["-h"]
        try:
            result=subprocess.run(cmd,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=60)
            # Some tools print usage with status 1. Preserve it for inspection.
            ok=result.returncode==0
            row={"tool":tool,"command":cmd,"returncode":result.returncode,"help_output":result.stdout,"ok":ok}
        except (OSError,subprocess.TimeoutExpired) as error:
            row={"tool":tool,"command":cmd,"ok":False,"error":str(error)}
        results.append(row)
        print(f"{tool}: {'OK' if row['ok'] else 'CHECK FAILED'}")
    args.out.write_text(json.dumps(results,indent=2)+"\n")
    print(f"Saved detailed help/version output to {args.out}")
    return 0 if all(x["ok"] for x in results) else 1


if __name__=="__main__":
    sys.exit(main())

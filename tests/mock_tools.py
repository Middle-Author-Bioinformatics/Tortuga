#!/usr/bin/env python3
"""CLI doubles for integration tests; no real biological inference."""
import json
from pathlib import Path
import sys

scenario_path, tool, *args = sys.argv[1:]
scenario = json.loads(Path(scenario_path).read_text())
with Path(scenario["trace"]).open("a") as f:
    f.write(json.dumps({"tool":tool,"args":args,"cwd":str(Path.cwd())})+"\n")


def val(flag):
    return args[args.index(flag)+1]


def put(path, value):
    p=Path(path)
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(value)


if tool in scenario.get("fail",[]):
    sys.exit(2)
if tool == "makeblastdb":
    put(val("-out")+".nsq","mock index")
elif tool == "fastp":
    import shutil
    shutil.copyfile(val("-i"),val("-o"))
    shutil.copyfile(val("-I"),val("-O"))
elif tool in {"novoplasty","getorganelle","getorganelle_assembly","mitofinder","mtgrasp"}:
    key = tool
    if tool == "mitofinder" and "-a" in args:
        key = "mitofinder_assembly"
    state = scenario.get("methods",{}).get(key)
    if state:
        seq=Path(scenario["sequence"]).read_text().strip()
        if state.get("overlap"):
            seq += seq[:150]
        if tool == "novoplasty":
            folder=Path.cwd()
            filename="Circularized_assembly_test.fasta" if state.get("circular") else "Contigs_test.fasta"
        elif tool.startswith("getorganelle"):
            folder=Path(val("-o"))
            filename="animal_mt.complete.graph1.1.path_sequence.fasta" if state.get("circular") else "animal_mt.incomplete.graph1.1.path_sequence.fasta"
        elif tool == "mitofinder":
            folder=Path.cwd()/"job_MitoFinder_Final_Results"
            filename="job_mtDNA_contig.fasta"
            if state.get("circular"):
                print("Evidences of circularization were found!")
        else:
            folder=Path(val("-o"))/"final_output"
            filename="test.final-mtgrasp-assembly.fa"
        put(folder/filename,">candidate "+("Circular" if state.get("circular") and tool=="mtgrasp" else "Linear")+"\n"+seq+"\n")
elif tool in {"tblastn","blastn"}:
    output=Path(val("-out"))
    key = next((x for x in sorted(scenario.get("methods",{}),key=len,reverse=True) if ("_"+x+"_") in str(output)), None)
    state=scenario.get("methods",{}).get(key,{})
    if tool == "tblastn":
        genes = ["cox1","nad1","cytb"] if state.get("complete") else ["cox1"]
        rows=[]
        for i,g in enumerate(genes):
            rows.append(f"{g}\tcandidate\t90\t200\t200\t12000\t1\t200\t{1+i*1000}\t{600+i*1000}\t1e-80\t500")
        put(output,"\n".join(rows)+"\n")
    elif "-db" in args:
        target = state.get("target","yes")
        hits = ["barcode1"] if target=="yes" else ["barcode2"] if target=="no" else ["barcode1","barcode2"]
        put(output,"\n".join(f"COI\t{x}\t99\t600\t600\t650\t1\t600\t1\t600\t1e-100\t1000" for x in hits)+"\n")
    else:
        put(output,"rrns\tcandidate\t90\t1000\t1000\t24000\t1\t1000\t8000\t8999\t1e-80\t500\n" if state.get("complete") else "")
elif tool in {"metaspades","megahit"}:
    folder=Path(val("-o"))
    put(folder/("contigs.fasta" if tool=="metaspades" else "final.contigs.fa"),">wgs\nACGT\n")
    if tool == "metaspades":
        put(folder/"assembly_graph.fastg",">NODE_1;\nACGT\n")
elif tool == "itsme":
    put(Path(val("-o"))/"run_summary.txt","mock ITS complete\n")
elif tool == "mitos2":
    put(Path(val("-o"))/"result.gff","##gff-version 3\ncandidate\tMITOS2\tgene\t1\t600\t.\t+\t.\tName=cox1\n")
elif tool == "minimap2":
    put(val("-o"),"@SQ\tSN:junction\tLN:1000\n"+"".join(
        f"read{i}\t0\tjunction\t{450+i}\t40\t150M\t*\t0\t0\tACGT\tIIII\n" for i in range(3)))
else:
    sys.exit("Unknown mock tool: "+tool)

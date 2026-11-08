import csv
import gzip
import importlib.util
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import time

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import mito_pipeline as m


@pytest.fixture
def fixture(tmp_path):
    reads=tmp_path/"reads with spaces"
    reads.mkdir()
    for mate in [1,2]:
        with gzip.open(reads/f"W30877_{mate}.fq.gz","wt") as f:
            for i in range(20):
                f.write(f"@read{i}/{mate}\n"+"ACGT"*40+"\n+\n"+"I"*160+"\n")
    rng=random.Random(12)
    sequence="".join(rng.choices("ACGT",k=12000))
    seqpath=tmp_path/"sequence.txt"
    seqpath.write_text(sequence)
    rec=SeqRecord(Seq(sequence),id="TEST",description="test reference")
    rec.annotations["molecule_type"]="DNA"
    for i,gene in enumerate(["COI","ND1","CYTB"]):
        rec.features.append(SeqFeature(FeatureLocation(i*1000,i*1000+600),type="CDS",
                                      qualifiers={"gene":[gene],"translation":["M"+"A"*199]}))
    rec.features.append(SeqFeature(FeatureLocation(8000,9000),type="rRNA",qualifiers={"gene":["rrnS"]}))
    gb=tmp_path/"reference.gb"
    SeqIO.write(rec,gb,"genbank")
    sheet=tmp_path/"samples.tsv"
    m.write_table(sheet,[{"sample":"W30877","target_taxon":"Target species","mt_code":5,"reference_gb":str(gb)}],
                  ["sample","target_taxon","mt_code","reference_gb"])
    bold=tmp_path/"bold.fasta"
    bold.write_text(">barcode1\n"+sequence[:650]+"\n>barcode2\n"+sequence[1000:1650]+"\n")
    meta=tmp_path/"bold.tsv"
    m.write_table(meta,[{"sequence_id":"barcode1","processid":"BOLD1","taxon":"Target species",
                        "lineage":"Animalia;Target genus;Target species","source":"BOLD test fixture"},
                       {"sequence_id":"barcode2","processid":"BOLD2","taxon":"Other species",
                        "lineage":"Animalia;Other genus;Other species","source":"BOLD test fixture"}],
                  ["sequence_id","processid","taxon","lineage","source"])
    refdir=tmp_path/"mitos_refs"
    (refdir/"release").mkdir(parents=True)
    itsmedb=tmp_path/"itsme_db"
    itsmedb.mkdir()
    trace=tmp_path/"trace.jsonl"
    scenario=tmp_path/"scenario.json"
    config=tmp_path/"config.json"
    commands={k:[sys.executable,str(ROOT/"tests/mock_tools.py"),str(scenario),k] for k in m.DEFAULT_COMMANDS}
    config.write_text(json.dumps({"commands":commands}))
    out=tmp_path/"results"
    argv=[str(reads),"--samples",str(sheet),"--bold-fasta",str(bold),"--bold-metadata",str(meta),
          "--mitos-refdir",str(refdir),"--mitos-refseqver","release","--config",str(config),"-o",str(out)]

    def run(methods,extra=(),fail=()):
        scenario.write_text(json.dumps({"trace":str(trace),"sequence":str(seqpath),"methods":methods,"fail":list(fail)}))
        code=m.main(argv+list(extra))
        events=[json.loads(x) for x in trace.read_text().splitlines()] if trace.exists() else []
        return code,events
    return {"run":run,"out":out,"reads":reads,"sheet":sheet,"trace":trace,"itsmedb":itsmedb,"argv":argv}


def summary(f):
    return json.loads((f["out"]/"W30877/final/summary.json").read_text())


def test_early_stop_keeps_full_branch_and_wgs(fixture):
    code,events=fixture["run"]({"novoplasty":{"circular":True,"complete":True}},
                              ["--branch","full","--itsme-db",str(fixture["itsmedb"])])
    assert code==0
    tools=[x["tool"] for x in events]
    assert "getorganelle" not in tools and "mitofinder" not in tools and "mtgrasp" not in tools
    assert tools.index("novoplasty")<tools.index("mitos2")<tools.index("itsme")<tools.index("metaspades")
    assert summary(fixture)["status"]=="recovered_annotated"
    assert summary(fixture)["genes_recovered"]==4
    with (fixture["out"]/"cohort_gene_recovery.csv").open() as f:
        rows=list(csv.DictReader(f))
    assert len(rows)==4 and all(r["screen_pass"]=="True" for r in rows)


def test_all_fallbacks_then_graph_recovery_with_confirmed_rescue(fixture):
    methods={"novoplasty":{"circular":False,"complete":False},
             "getorganelle":{"circular":False,"complete":True,"target":"no"},
             "mitofinder":{"complete":False},"mtgrasp":{"circular":False,"complete":False},
             "getorganelle_assembly":{"circular":True,"complete":True}}
    code,events=fixture["run"](methods)
    assert code==0
    names=[e["tool"] for e in events if e["tool"] in m.READ_METHODS]
    assert names==list(m.READ_METHODS)
    graph=next(e for e in events if e["tool"]=="getorganelle_assembly")
    assert "-g" in graph["args"] and "-s" not in graph["args"]
    assert "mitofinder_assembly" not in [e["tool"] for e in events]
    assert summary(fixture)["winning_method"]=="getorganelle_assembly"
    seeds=list((fixture["out"]/"W30877").rglob("rescue_seed.json"))
    assert seeds
    assert all(json.loads(p.read_text())["target_status"]=="confirmed_target" for p in seeds)
    assert all("getorganelle_" not in json.loads(p.read_text())["candidate_id"] for p in seeds)


def test_no_assembly_exports_longest_target_partial(fixture):
    code,events=fixture["run"]({"novoplasty":{"circular":False,"complete":True}},["--no-assembly"])
    assert code==0
    assert "metaspades" not in [e["tool"] for e in events]
    assert "mitos2" not in [e["tool"] for e in events]
    assert summary(fixture)["status"]=="partial_target_confirmed"
    assert (fixture["out"]/"W30877/final/best_sequence.fasta").exists()


def test_conflicting_bold_blocks_acceptance_and_rescue(fixture):
    code,events=fixture["run"]({"novoplasty":{"circular":True,"complete":True,"target":"mixed"}},["--no-assembly"])
    assert code==0
    assert summary(fixture)["status"]=="unresolved_candidate"
    assert not summary(fixture)["winning_method"]
    assert not list((fixture["out"]/"W30877").rglob("rescue_seed.json"))


def test_annotation_failure_reported(fixture):
    code,_=fixture["run"]({"novoplasty":{"circular":True,"complete":True}},["--no-assembly"],fail=["mitos2"])
    assert code==1
    assert summary(fixture)["status"]=="recovered_annotation_failed"


def test_resume_caches_success_and_rejects_changed_options(fixture):
    methods={"novoplasty":{"circular":True,"complete":True}}
    fixture["run"](methods,["--no-assembly"])
    count=len(fixture["trace"].read_text().splitlines())
    code,_=fixture["run"](methods,["--no-assembly","--resume"])
    assert code==0
    assert len(fixture["trace"].read_text().splitlines())==count
    code,_=fixture["run"](methods,["--no-assembly","--resume","--coi-identity","98"])
    assert code==1


def test_failed_method_gets_new_workdir_on_resume(fixture):
    methods={"novoplasty":{"circular":True,"complete":True}}
    fixture["run"](methods,["--no-assembly"],fail=["novoplasty"])
    code,events=fixture["run"](methods,["--no-assembly","--resume"])
    assert code==0
    calls=[x for x in events if x["tool"]=="novoplasty"]
    assert len(calls)==2 and calls[0]["cwd"]!=calls[1]["cwd"]


def test_megahit_skips_graph_extraction_but_uses_mitofinder(fixture):
    code,events=fixture["run"]({"mitofinder_assembly":{"complete":True}},["--assembler","megahit"])
    assert code==0
    assert "getorganelle_assembly" not in [e["tool"] for e in events]
    mf=[e for e in events if e["tool"]=="mitofinder"]
    assert any("-a" in e["args"] for e in mf)


def test_dry_run_uses_no_external_commands(fixture):
    code,events=fixture["run"]({},["--dry-run","--assembly-percent","10"])
    assert code==0 and not events
    assert summary(fixture)["status"]=="planned"


def test_mitofinder_circularity_requires_raw_junction_support(fixture):
    code,events=fixture["run"]({"mitofinder":{"circular":True,"complete":True}},["--no-assembly"])
    assert code==0
    assert "minimap2" in [e["tool"] for e in events]
    result=summary(fixture)
    assert result["winning_method"]=="mitofinder"
    assert result["best_candidate"]["circular_evidence"]=="mitofinder_log_and_raw_junction_reads"


def test_terminal_overlap_does_not_pass_when_junction_verification_fails(fixture):
    code,_=fixture["run"]({"novoplasty":{"complete":True,"overlap":True}},["--no-assembly"],fail=["minimap2"])
    assert code==1
    assert not summary(fixture)["winning_method"]
    assert summary(fixture)["best_candidate"]["terminal_overlap_bp"]==150
    assert not summary(fixture)["best_candidate"]["circular"]


def test_prepare_bold_preserves_provenance_and_excludes_other_markers(tmp_path):
    export=tmp_path/"export.tsv"
    m.write_table(export,[{"processid":"PID1","marker_code":"COI-5P","nuc":"ACGT"*150,
                          "kingdom":"Animalia","genus":"Target","species":"Target species"},
                         {"processid":"PID2","marker_code":"ITS","nuc":"ACGT"*150,
                          "kingdom":"Fungi","genus":"Other","species":"Other species"}],
                  ["processid","marker_code","nuc","kingdom","genus","species"])
    out=tmp_path/"bold_prepared"
    p=subprocess.run([sys.executable,str(ROOT/"prepare_bold.py"),str(export),"--source","BOLD test export","-o",str(out)],capture_output=True,text=True)
    assert p.returncode==0,p.stderr
    records=list(SeqIO.parse(out/"coi.fasta","fasta"))
    assert [r.id for r in records]==["PID1_COI5P"]
    metadata=m.load_bold_metadata(out/"metadata.tsv")
    assert metadata["PID1_COI5P"]["source"]=="BOLD test export"


def test_pair_subsampling_is_reproducible_and_checks_mates(fixture,tmp_path):
    paths=[tmp_path/x for x in ["a.fq.gz","b.fq.gz","c.fq.gz","d.fq.gz"]]
    a,b,c,d=paths
    r1,r2=m.discover_pairs(fixture["reads"])["W30877"]
    first=m.subsample_pairs(r1,r2,a,b,.5,17)
    second=m.subsample_pairs(r1,r2,c,d,.5,17)
    assert first==second
    assert list(m.fastq_records(a))==list(m.fastq_records(c))
    assert [m.template_id(r) for r in m.fastq_records(a)]==[m.template_id(r) for r in m.fastq_records(b)]
    with gzip.open(r2,"wt") as f:
        f.write("@wrong/2\nACGT\n+\nIIII\n")
    with pytest.raises(ValueError,match="mates differ"):
        m.subsample_pairs(r1,r2,c,d,.5,17)


def test_native_circularity_never_uses_incomplete_or_length():
    assert m.native_circular("getorganelle","x.incomplete.graph1.path_sequence.fasta","")=="not_established"
    assert m.native_circular("getorganelle","x.complete.graph1.path_sequence.fasta","")!="not_established"
    assert m.native_circular("mitofinder","x_mtDNA_contig.fasta","")=="not_established"
    assert m.native_circular("mtgrasp","x.final-mtgrasp-assembly.fa","Linear")=="not_established"
    assert m.native_circular("mtgrasp","x.final-mtgrasp-assembly.fa","sample_Circular_Standardized")!="not_established"
    assert m.native_circular("mtgrasp","x.final-mtgrasp-assembly.fa","circular=false")=="not_established"


def test_bold_hit_limit_and_missing_metadata_remain_unresolved():
    row={"qseqid":"COI","sseqid":"id","pident":"99","length":"600","qlen":"600","slen":"650","bitscore":"1000"}
    params=("Target",97,.8,400,.02)
    assert m.bold_decision([row],{},*params,500)["status"]=="unresolved"
    meta={"id":{"taxon":"Target","lineage":"Animalia;Target"}}
    assert m.bold_decision([row],meta,*params,1)["status"]=="unresolved"
    assert m.bold_decision([row],meta,*params,500)["status"]=="confirmed_target"


def test_junction_spanning_requires_continuous_supported_alignment(tmp_path):
    sam=tmp_path/"junction.sam"
    rows=["@SQ\tSN:junction\tLN:1000"]
    for name,flag,pos,mapq,cigar in [("a/1",0,450,30,"150M"),("a/2",0,450,30,"150M"),
                                    ("b",0,460,30,"150M"),("c",0,450,0,"150M"),
                                    ("d",256,450,30,"150M"),("e",0,450,30,"40M40D60M")]:
        rows.append(f"{name}\t{flag}\tjunction\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\tACGT\tIIII")
    sam.write_text("\n".join(rows)+"\n")
    result=m.junction_support(sam,500)
    assert result["spanning_templates"]==2
    assert result["distinct_alignment_starts"]==2


def test_timeout_kills_owned_child_process(tmp_path):
    marker=tmp_path/"child_finished"
    script=tmp_path/"spawn.py"
    script.write_text("import subprocess,sys,time\nsubprocess.Popen([sys.executable,'-c',\"import time;from pathlib import Path;time.sleep(.6);Path("+repr(str(marker))+").write_text('bad')\"])\ntime.sleep(20)\n")
    runner=m.Runner(tmp_path,{"spawn":[sys.executable,str(script)]})
    assert not runner.run("timeout","spawn",[],tmp_path,.15)
    time.sleep(.8)
    assert not marker.exists()
    assert runner.events[-1]["status"]=="timeout"


def test_main_selects_environments_without_a_config_file(fixture,monkeypatch):
    import configure_belly
    argv=fixture["argv"]
    index=argv.index("--config")
    commands=json.loads(Path(argv[index+1]).read_text())["commands"]
    del argv[index:index+2]
    calls=[]
    def resolve(*args):
        calls.append(args)
        return commands,[{"tool":"novoplasty","status":"found","environment":"/chosen/env",
                          "entrypoint":"mock","alternatives":"/other/env","note":"Alternative relevant package"}]
    monkeypatch.setattr(configure_belly,"runtime_commands",resolve)
    code,_=fixture["run"]({"novoplasty":{"circular":True,"complete":True}},["--no-assembly"])
    assert code==0 and len(calls)==1
    notes=json.loads((fixture["out"]/"environment_selection.json").read_text())
    assert notes[0]["environment"]=="/chosen/env"
    assert notes[0]["alternatives"]=="/other/env"


def test_missing_auto_selected_tool_is_reported_without_using_current_path(tmp_path,monkeypatch):
    runner=m.Runner(tmp_path,{"novoplasty":[]})
    def forbidden(*args,**kwargs):
        raise AssertionError("Missing tool must not launch an unintended PATH executable")
    monkeypatch.setattr(subprocess,"Popen",forbidden)
    assert runner.run("novoplasty","novoplasty",[],tmp_path,10) is False
    assert runner.events[-1]["returncode"]==127

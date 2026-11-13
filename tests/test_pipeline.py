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
          "--mitos-refdir",str(refdir),"--mitos-refseqver","release","--config",str(config),"-o",str(out),"--keep-intermediates"]

    def run(methods,extra=(),fail=()):
        scenario.write_text(json.dumps({"trace":str(trace),"sequence":str(seqpath),"methods":methods,"fail":list(fail)}))
        code=m.main(argv+list(extra))
        events=[json.loads(x) for x in trace.read_text().splitlines()] if trace.exists() else []
        return code,events
    return {"run":run,"out":out,"reads":reads,"sheet":sheet,"trace":trace,"itsmedb":itsmedb,"argv":argv}


def summary(f):
    return json.loads((f["out"]/"W30877/final/summary.json").read_text())


def test_automatic_reference_and_code_from_taxon_only_sheet(fixture,tmp_path):
    catalog=tmp_path/"genbank"
    catalog.mkdir()
    source=tmp_path/"reference.gb"
    record=SeqIO.read(source,"genbank")
    record.annotations.update(organism="Target species",taxonomy=["Eukaryota","Metazoa","Target"],topology="circular")
    record.description="Target species mitochondrion, complete genome"
    for feature in record.features:
        if feature.type=="CDS":
            feature.qualifiers["transl_table"]=["5"]
    SeqIO.write(record,catalog/"selected.gb","genbank")
    m.write_table(fixture["sheet"],[{"sample":"W30877","target_taxon":"Target species"}],
                  ["sample","target_taxon"])
    code,events=fixture["run"]({"novoplasty":{"circular":True,"complete":True}},
                              ["--genbank-dir",str(catalog)])
    assert code==0
    with (fixture["out"]/"reference_selection.tsv").open() as handle:
        rows=list(csv.DictReader(handle,delimiter="\t"))
    assert rows[0]["reference_gb"]==str(catalog/"selected.gb")
    assert rows[0]["code_source"]=="GenBank_CDS"
    assert rows[0]["match_level"]=="exact_taxon"
    assert summary(fixture)["target_taxon"]=="Target species"
    with (fixture["out"]/"resolved_samples.tsv").open() as handle:
        resolved=list(csv.DictReader(handle,delimiter="\t"))
    assert resolved[0]["mt_code"]=="5"


def test_family_reference_preserves_target_taxonomy_and_extra_metadata(fixture,tmp_path):
    catalog=tmp_path/"genbank"
    catalog.mkdir()
    record=SeqIO.read(tmp_path/"reference.gb","genbank")
    record.annotations.update(organism="Relative species",
                              taxonomy=["Eukaryota","Metazoa","Annelida","Targetida","Targetidae","Relative"])
    record.description="Relative species mitochondrion, complete genome"
    for feature in record.features:
        if feature.type=="CDS":
            feature.qualifiers["transl_table"]=["5"]
    SeqIO.write(record,catalog/"family.gb","genbank")
    input_row={"sample":"W30877","target_taxon":"Target species","Genus":"Target",
               "Species":"species","Family":"Targetidae","Order":"Targetida",
               "Phylum":"Annelida","AphiaID":"1234","original_identification":"Target (Subgenus) species"}
    m.write_table(fixture["sheet"],[input_row],list(input_row))
    code,events=fixture["run"]({"novoplasty":{"circular":True,"complete":True}},
                              ["--genbank-dir",str(catalog)])
    assert code==0
    with (fixture["out"]/"reference_selection.tsv").open() as handle:
        selected=list(csv.DictReader(handle,delimiter="\t"))[0]
    assert selected["match_rank"]=="family" and selected["matched_taxon"]=="Targetidae"
    with (fixture["out"]/"resolved_samples.tsv").open() as handle:
        resolved=list(csv.DictReader(handle,delimiter="\t"))[0]
    assert all(resolved[key]==value for key,value in input_row.items())
    result=summary(fixture)
    assert result["target_taxon"]=="Target species"
    assert result["taxonomy"]["family"]=="Targetidae"
    assert result["best_candidate"]["target_status"]=="confirmed_target"


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


def test_prepare_bold_accepts_large_fields(tmp_path):
    export = tmp_path/"large.tsv"
    m.write_table(export, [{"processid":"PID1", "marker_code":"COI-5P",
                           "nuc":"ACGT"*150, "species":"Target species",
                           "notes":"x"*200_000}],
                  ["processid","marker_code","nuc","species","notes"])
    out = tmp_path/"prepared"
    p = subprocess.run([sys.executable,str(ROOT/"prepare_bold.py"),str(export),
                        "--source","large field test","-o",str(out)],
                       capture_output=True,text=True)
    assert p.returncode == 0, p.stderr
    records = list(SeqIO.parse(out/"coi.fasta","fasta"))
    assert len(records) == 1 and str(records[0].seq) == "ACGT"*150
    metadata = m.load_bold_metadata(out/"metadata.tsv")
    assert metadata[records[0].id]["taxon"] == "Target species"


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


def itsme_argv(fixture):
    argv=fixture["argv"]
    config=argv[argv.index("--config")+1]
    argv[:]=["--pipeline","itsme","-i",str(fixture["reads"]),"--itsme-db",str(fixture["itsmedb"]),
             "--config",config,"-o",str(fixture["out"]),"--keep-intermediates"]


def test_controller_itsme_only_needs_no_mito_references_or_sample_sheet(fixture):
    itsme_argv(fixture)
    code,events=fixture["run"]({})
    assert code==0
    assert [e["tool"] for e in events]==["itsme"]
    assert summary(fixture)["pipeline"]=="itsme"
    assert summary(fixture)["status"]=="completed"
    assert not (fixture["out"]/"bold_reference").exists()


def test_controller_both_runs_itsme_once_in_order(fixture):
    code,events=fixture["run"]({"novoplasty":{"circular":True,"complete":True}},
                              ["--pipeline","both","--itsme-db",str(fixture["itsmedb"]),"--sensitivity","3"])
    assert code==0
    names=[e["tool"] for e in events]
    assert names.count("itsme")==1
    assert names.index("novoplasty")<names.index("itsme")<names.index("metaspades")
    call=next(e for e in events if e["tool"]=="itsme")
    assert call["args"][call["args"].index("--sensitivity")+1]=="3"
    assert summary(fixture)["pipeline"]=="both"


def test_controller_mitolooper_never_calls_itsme(fixture):
    code,events=fixture["run"]({"novoplasty":{"circular":True,"complete":True}},["--pipeline","mitolooper","--no-assembly"])
    assert code==0
    assert "itsme" not in [e["tool"] for e in events]
    assert summary(fixture)["pipeline"]=="mitolooper"


def test_controller_itsme_resume_and_retry(fixture):
    itsme_argv(fixture)
    code,_=fixture["run"]({},fail=["itsme"])
    assert code==1
    code,events=fixture["run"]({},["--resume"])
    assert code==0
    assert len(events)==2 and events[0]["cwd"]!=events[1]["cwd"]
    count=len(events)
    code,events=fixture["run"]({},["--resume"])
    assert code==0 and len(events)==count


def test_controller_itsme_continues_after_failures_and_accepts_name_only_sheet(fixture):
    import shutil
    for mate in [1,2]:
        shutil.copyfile(fixture["reads"]/f"W30877_{mate}.fq.gz",fixture["reads"]/f"W30882_{mate}.fq.gz")
    itsme_argv(fixture)
    m.write_table(fixture["sheet"],[{"sample":"W30877"},{"sample":"W30882"}],["sample"])
    code,events=fixture["run"]({},["--samples",str(fixture["sheet"])],fail=["itsme"])
    assert code==1 and len(events)==2
    with (fixture["out"]/"cohort_summary.csv").open() as f:
        rows=list(csv.DictReader(f))
    assert {r["sample"] for r in rows}=={"W30877","W30882"}
    assert all(r["status"]=="failed" for r in rows)


def test_controller_requires_valid_selector_and_rejects_conflicting_branch(fixture):
    with pytest.raises(SystemExit) as error:
        m.main([str(fixture["reads"])],require_pipeline=True)
    assert error.value.code==2
    with pytest.raises(SystemExit) as error:
        m.main([str(fixture["reads"]),"--pipeline","bogus"],require_pipeline=True)
    assert error.value.code==2
    with pytest.raises(SystemExit) as error:
        m.main([str(fixture["reads"]),"--pipeline","mitolooper","--branch","full"],require_pipeline=True)
    assert error.value.code==2


def test_controller_shell_uses_shared_launcher(tmp_path):
    env=dict(os.environ,MITO_PYTHON=sys.executable)
    p=subprocess.run(["bash",str(ROOT/"pipeline_controller.sh"),"--version"],env=env,capture_output=True,text=True)
    assert p.returncode==0 and p.stdout.strip()==m.VERSION


@pytest.mark.parametrize("combined",[False,True])
def test_reference_failure_does_not_block_other_samples_or_independent_itsme(fixture,tmp_path,combined):
    import shutil
    for mate in [1,2]:
        shutil.copyfile(fixture["reads"]/f"W30877_{mate}.fq.gz",fixture["reads"]/f"NOREF_{mate}.fq.gz")
    m.write_table(fixture["sheet"],[
        {"sample":"NOREF","target_taxon":"Unknown species","mt_code":5,"reference_gb":str(tmp_path/"missing.gb")},
        {"sample":"W30877","target_taxon":"Target species","mt_code":5,"reference_gb":str(tmp_path/"reference.gb")}
    ],["sample","target_taxon","mt_code","reference_gb"])
    extra=["--pipeline","both" if combined else "mitolooper"]
    if combined:
        extra += ["--itsme-db",str(fixture["itsmedb"])]
    code,events=fixture["run"]({"novoplasty":{"circular":True,"complete":True}},extra)
    assert code==1
    assert summary(fixture)["status"]=="recovered_annotated"
    invocation=json.loads((fixture["out"]/"last_invocation.json").read_text())
    unavailable=next(row for row in invocation["results"] if row["sample"]=="NOREF")
    assert unavailable["status"]=="reference_unavailable"
    assert "reference_selection" in unavailable["failed_stages"]
    with (fixture["out"]/"reference_selection.tsv").open() as f:
        choices={r["sample"]:r for r in csv.DictReader(f,delimiter="\t")}
    assert choices["NOREF"]["selection_status"]=="unavailable"
    assert choices["W30877"]["selection_status"]=="available"
    assert sum(event["tool"]=="itsme" for event in events)==(2 if combined else 0)
    if combined:
        assert unavailable["itsme_ok"] is True
        assert Path(unavailable["itsme_output"]).joinpath("itsme_run_summary.txt").exists()
        assert unavailable["itsme_output"]==str(fixture["out"]/"NOREF/final")
    assert not (fixture["out"]/"NOREF/01_read_recovery").exists()


def test_no_usable_references_still_runs_itsme_without_building_bold(fixture,tmp_path):
    m.write_table(fixture["sheet"],[{"sample":"W30877","target_taxon":"Target species",
                                      "mt_code":5,"reference_gb":str(tmp_path/"missing.gb")}],
                  ["sample","target_taxon","mt_code","reference_gb"])
    code,events=fixture["run"]({},["--pipeline","both","--itsme-db",str(fixture["itsmedb"])])
    assert code==1
    assert [event["tool"] for event in events]==["itsme"]
    assert not (fixture["out"]/"bold_reference").exists()
    result=json.loads((fixture["out"]/"last_invocation.json").read_text())["results"][0]
    assert result["status"]=="reference_unavailable" and result["itsme_ok"] is True
    code,events=fixture["run"]({},["--pipeline","both","--itsme-db",str(fixture["itsmedb"]),"--resume"])
    assert code==1 and len(events)==1


def test_linear_all_gene_candidate_stops_and_annotates_with_linear_mode(fixture):
    code,events=fixture["run"]({"novoplasty":{"complete":True}},
                              ["--accept-topology","any","--no-assembly"])
    assert code==0
    tools=[e["tool"] for e in events]
    assert "getorganelle" not in tools and "mtgrasp" not in tools
    annotation=next(e for e in events if e["tool"]=="mitos2")
    assert "--linear" in annotation["args"]
    result=summary(fixture)
    assert result["status"]=="accepted_by_policy_annotated"
    assert result["accepted_by_policy"] and not result["strict_complete_candidate"]
    assert result["complete_by_reference"] and not result["circular"]
    assert result["best_candidate"]["homology_sequence_mode"]=="linear"
    final=fixture["out"]/"W30877/final"
    assert (final/"accepted_candidate.fasta").exists()
    assert not (final/"circular_target_candidate.fasta").exists()
    candidate_dir=fixture["out"]/"W30877/candidates"/result["best_candidate"]["candidate_id"]
    assert len(SeqIO.read(candidate_dir/"screen.fasta","fasta").seq)==12000


def test_gene_threshold_stops_on_partial_and_skips_wgs_without_skipping_itsme(fixture):
    code,events=fixture["run"]({"novoplasty":{"genes":["cox1","nad1"]}},
        ["--pipeline","both","--itsme-db",str(fixture["itsmedb"]),
         "--accept-topology","any","--accept-rule","genes","--min-detected-genes","2",
         "--skip-wgs-if-accepted"])
    assert code==0
    tools=[e["tool"] for e in events]
    assert tools.count("itsme")==1 and "metaspades" not in tools
    assert not any(tool in tools for tool in ["getorganelle","mitofinder","mtgrasp"])
    result=summary(fixture)
    assert result["status"]=="accepted_by_policy_annotated"
    assert result["genes_recovered"]==2 and result["genes_expected"]==4
    assert not result["complete_by_reference"] and result["itsme_ok"]
    assert result["wgs_ok"] is None
    assert result["best_candidate"]["missing_genes"]==["cytb"]
    assert result["best_candidate"]["missing_rrnas"]==["rrns"]
    with (fixture["out"]/"cohort_summary.csv").open() as handle:
        cohort=list(csv.DictReader(handle))[0]
    assert cohort["accepted_by_policy"]=="True" and cohort["strict_complete_candidate"]=="False"


def test_insufficient_genes_continues_to_later_method(fixture):
    code,events=fixture["run"]({"novoplasty":{"genes":["cox1"]},
                              "getorganelle":{"genes":["cox1","nad1"]}},
        ["--accept-topology","any","--accept-rule","genes","--min-detected-genes","2","--no-assembly"])
    assert code==0 and summary(fixture)["winning_method"]=="getorganelle"
    assert "mtgrasp" not in [e["tool"] for e in events]
    assessments=sorted((fixture["out"]/"W30877/candidates").glob("*/assessment.json"))
    first=json.loads(assessments[0].read_text())
    assert not first["accepted"] and first["acceptance_failures"]==["selected_rule"]


@pytest.mark.parametrize("rule,min_genes,min_bp,accepted",[
    ("length",None,12000,True),
    ("length",None,12001,False),
    ("either",2,12000,True),
    ("either",1,12001,True),
    ("either",2,12001,False),
    ("both",1,12000,True),
    ("both",2,12000,False),
    ("both",1,12001,False),
])
def test_length_and_combined_rules_use_one_contig(fixture,rule,min_genes,min_bp,accepted):
    extra=["--accept-topology","any","--accept-rule",rule,"--min-accepted-bp",str(min_bp),"--no-assembly"]
    if min_genes is not None:
        extra += ["--min-detected-genes",str(min_genes)]
    code,events=fixture["run"]({"novoplasty":{"genes":["cox1"]}},extra)
    result=summary(fixture)
    assert code==0 and result["accepted_by_policy"] is accepted
    assert ("mitos2" in [e["tool"] for e in events]) is accepted
    assert result["best_candidate"]["length_threshold_met"] is (12000>=min_bp)
    assert not result["complete_by_reference"]


def test_gene_count_does_not_sum_fragments_from_different_assemblers(fixture):
    code,events=fixture["run"]({"novoplasty":{"genes":["cox1","nad1"]},
                              "getorganelle":{"genes":["cox1","cytb"]}},
        ["--accept-topology","any","--accept-rule","genes","--min-detected-genes","3","--no-assembly"])
    assert code==0 and not summary(fixture)["accepted_by_policy"]
    assert summary(fixture)["genes_recovered"]==2


@pytest.mark.parametrize("target",["no","mixed"])
def test_relaxed_length_policy_still_rejects_non_target_and_conflicting_bold(fixture,target):
    code,events=fixture["run"]({"novoplasty":{"complete":True,"target":target}},
        ["--accept-topology","any","--accept-rule","length","--min-accepted-bp","10000","--no-assembly"])
    result=summary(fixture)
    assert code==0 and not result["accepted_by_policy"]
    assert "target_confirmed" in result["best_candidate"]["acceptance_failures"]
    assert "mitos2" not in [e["tool"] for e in events]


def test_relaxed_policy_does_not_override_topology_requirement(fixture):
    code,events=fixture["run"]({"novoplasty":{"complete":True}},
        ["--accept-rule","length","--min-accepted-bp","10000","--no-assembly"])
    result=summary(fixture)
    assert code==0 and not result["accepted_by_policy"]
    assert result["best_candidate"]["acceptance_failures"]==["topology"]


def test_annotation_failure_under_relaxed_policy_returns_nonzero(fixture):
    code,events=fixture["run"]({"novoplasty":{"genes":["cox1"]}},
        ["--accept-topology","any","--accept-rule","length","--no-assembly"],fail=["mitos2"])
    assert code==1 and summary(fixture)["status"]=="accepted_by_policy_annotation_failed"


def test_linear_terminal_repeat_is_preserved_without_junction_support(fixture):
    code,events=fixture["run"]({"novoplasty":{"complete":True,"overlap":True}},
        ["--accept-topology","any","--no-assembly"],fail=["minimap2"])
    result=summary(fixture)
    assert code==0 and result["accepted_by_policy"] and not result["circular"]
    candidate=result["best_candidate"]
    assert candidate["terminal_overlap_bp"]==150 and candidate["terminal_overlap_removed_bp"]==0
    sequence=SeqIO.read(candidate["fasta"],"fasta").seq
    screen=SeqIO.read(fixture["out"]/"W30877/candidates"/candidate["candidate_id"]/"screen.fasta","fasta").seq
    assert len(sequence)==12150 and sequence==screen


def test_strict_circular_candidate_uses_doubled_screen_and_circular_annotation(fixture):
    code,events=fixture["run"]({"novoplasty":{"circular":True,"complete":True}},["--no-assembly"])
    result=summary(fixture)
    assert code==0 and result["status"]=="recovered_annotated"
    assert result["strict_complete_candidate"]
    annotation=next(e for e in events if e["tool"]=="mitos2")
    assert "--linear" not in annotation["args"]
    candidate_dir=fixture["out"]/"W30877/candidates"/result["best_candidate"]["candidate_id"]
    assert len(SeqIO.read(candidate_dir/"screen.fasta","fasta").seq)==24000
    assert (fixture["out"]/"W30877/final/circular_target_candidate.fasta").exists()


def test_acceptance_options_validate_before_launching_tools(fixture):
    for options in [
        ["--accept-rule","genes"],
        ["--accept-rule","both","--min-detected-genes","0"],
        ["--accept-rule","length","--min-detected-genes","3"],
        ["--accept-rule","genes","--min-detected-genes","2","--min-accepted-bp","10000"],
        ["--min-accepted-bp","400"],
        ["--min-accepted-bp","31000"],
    ]:
        with pytest.raises(ValueError):
            fixture["run"]({},options)
    assert not fixture["trace"].exists()


def test_policy_safety_guards_and_unique_gene_counts():
    args=m.parser().parse_args(["--accept-topology","any","--accept-rule","genes","--min-detected-genes","2"])
    base={"length_bp":12000,"ambiguous_bp":0,"circular":False,"complete_by_reference":False,
          "target_status":"confirmed_target","genes_found":["cox1","cox1"],"rrnas_found":[]}
    assert not m.apply_acceptance(dict(base),args)["accepted"]
    base["rrnas_found"]=["rrns"]
    assert m.apply_acceptance(dict(base),args)["accepted"]
    for override in [{"ambiguous_bp":1},{"length_bp":31000},{"validation_error":"failed"},
                     {"target_status":"unresolved"},{"length_bp":499}]:
        assert not m.apply_acceptance({**base,**override},args)["accepted"]


def test_relaxed_policy_also_stops_post_wgs_fallback(fixture):
    code,events=fixture["run"]({"getorganelle_assembly":{"genes":["cox1"]}},
        ["--accept-topology","any","--accept-rule","length","--skip-wgs-if-accepted"])
    assert code==0 and summary(fixture)["winning_method"]=="getorganelle_assembly"
    assert summary(fixture)["status"]=="accepted_by_policy_annotated"
    assert "mitofinder" in [e["tool"] for e in events]  # Read-based attempt still ran.
    assert not any(e["tool"]=="mitofinder" and "-a" in e["args"] for e in events)


def test_default_compaction_preserves_deliverables_and_resumes_without_reassembly(fixture):
    fixture["argv"].remove("--keep-intermediates")
    raw_before={p:p.read_bytes() for p in fixture["reads"].iterdir()}
    methods={"novoplasty":{"circular":True,"complete":True}}
    code,events=fixture["run"](methods,["--pipeline","both","--itsme-db",str(fixture["itsmedb"])])
    assert code==0
    root=fixture["out"]/"W30877"
    final=root/"final"
    assert (final/"itsme_complete_ITS.fasta").exists()
    assert (final/"itsme_master_summary.csv").exists()
    assert (final/"itsme_run_summary.txt").exists()
    assert (final/"mitos_result.gff").exists()
    assert (final/"wgs_assembly_graph.fastg").exists()
    assert (final/"wgs_contigs.fasta").exists()
    assert not (root/"03_wgs").exists() and not (root/"02_itsme").exists()
    assert not list(root.rglob("*.fq.gz")) and not list(root.rglob("*.fastq.gz"))
    assert not (fixture["out"]/"bold_reference/coi.nsq").exists()
    assert (root/"qc/fastp.html").exists()
    result=summary(fixture)
    assert result["itsme_output"]==str(final) and Path(result["best_candidate"]["fasta"]).exists()
    assert result["annotation_dir"]==str(final)
    assert {p:p.read_bytes() for p in fixture["reads"].iterdir()}==raw_before
    before=len(events)
    code,events=fixture["run"](methods,["--pipeline","both","--itsme-db",str(fixture["itsmedb"]),"--resume"])
    assert code==0 and len(events)==before


def test_default_compaction_applies_to_standalone_itsme_and_caches_final_results(fixture):
    itsme_argv(fixture)
    fixture["argv"].remove("--keep-intermediates")
    code,events=fixture["run"]({})
    assert code==0
    root=fixture["out"]/"W30877"
    assert (root/"final/itsme_complete_ITS.fasta").exists()
    assert not (root/"itsme").exists()
    code,events=fixture["run"]({},["--resume"])
    assert code==0 and len(events)==1


def test_default_compaction_preserves_failed_stage_status_logs_and_final_itsme(fixture):
    fixture["argv"].remove("--keep-intermediates")
    code,events=fixture["run"]({"novoplasty":{"complete":True}},
        ["--pipeline","both","--itsme-db",str(fixture["itsmedb"])],fail=["metaspades","mtgrasp"])
    assert code==1
    result=summary(fixture)
    assert result["wgs_ok"] is False and result["itsme_ok"] is True
    assert "wgs" in result["failed_stages"] and "mtgrasp" in result["failed_stages"]
    root=fixture["out"]/"W30877"
    assert (root/"logs/wgs.log").exists() and (root/"logs/mtgrasp.log").exists()
    assert (root/"final/itsme_complete_ITS.fasta").exists()
    assert not (root/"03_wgs").exists()

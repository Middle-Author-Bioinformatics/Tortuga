import csv
import fcntl
import json
from pathlib import Path
import sys

import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import compact_outputs as c


def put(path,text):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(text)
    return path


@pytest.fixture
def completed(tmp_path):
    root=tmp_path/"results"
    sample=root/"SAMPLE"
    raw=put(tmp_path/"raw/SAMPLE_1.fq.gz","original reads")
    c.atomic_json(sample/"run_config.json",{"reads":[{"path":str(raw)}],"sample":"SAMPLE"})
    c.atomic_json(sample/"final/summary.json",{"sample":"SAMPLE","pipeline":"both","status":"recovered_annotated",
                                            "itsme_ok":True,"wgs_ok":False})
    c.atomic_json(root/"last_invocation.json",{"results":[{"sample":"SAMPLE","status":"recovered_annotated","wgs_ok":False}]})
    put(sample/"011_novoplasty/attempt_001/Circularized_assembly_1.fasta",">mito\nACGT\n")
    native=sample/"02_itsme/attempt_001/output"
    put(native/"final/complete_ITS.fasta",">ITS\nGGGG\n")
    put(native/"final/reported_loci.fasta",">locus\nCCCC\n")
    put(native/"final/master_summary.tsv","locus\tstatus\n1\tcomplete\n")
    put(native/"master_summary.csv","locus,status\n1,complete\n")
    put(native/"run_summary.txt","ITS successful\n")
    put(native/"final/accepted_R1.fastq.gz","discard recruited reads")
    put(native/"mapping/alignments.bam","discard mappings")
    put(native/"recruitment/round_00/frontier_baits.fasta",">bait\nAAAA\n")
    put(native/"final/assembly_graph.fastg",">graph\nATGC\n")
    put(native/"assembly/spades/assembly_graph.fastg",">graph\nATGC\n")
    put(native/"assembly/spades/K55/assembly_graph.fastg",">intermediate\nAAAA\n")
    put(native/"qc/fastp.json",'{"reads":123}')
    put(sample/"012_getorganelle/attempt_001/output/extended_1.fq.gz","discard reads")
    put(sample/"012_getorganelle/attempt_001/output/a.path_sequence.fasta",">organelle\nTTTT\n")
    put(sample/"012_getorganelle/attempt_001/output/extended_spades/assembly_graph.fastg",">organelle_graph\nGTAC\n")
    put(sample/"012_getorganelle/attempt_001/output/extended_spades/K55/assembly_graph.fastg","discard intermediate graph")
    put(sample/"03_wgs/attempt_001/output/K55/huge.bin","discard failed assembler scratch"*100)
    put(sample/"03_wgs/attempt_001/output/K55/contigs.fasta",">intermediate\nAAAA\n")
    put(sample/"03_wgs/attempt_001/output/contigs.fasta",">wgs\nATGC\n")
    put(sample/"03_wgs/attempt_001/output/assembly_graph_with_scaffolds.gfa","H\tVN:Z:1.0\n")
    put(sample/"03_wgs/attempt_001/output/spades.log","assembly timed out during K55\n")
    gff=put(sample/"05_mitos2/attempt_001/output/result.gff","##gff-version 3\n")
    c.atomic_json(sample/"final/annotation_summary.json",{"gff":str(gff)})
    seq=put(sample/"candidates/011_novoplasty_001/sequence.fasta",">candidate\nACGT\n")
    put(seq.parent/"screen.fasta",">doubled\nACGTACGT\n")
    put(seq.parent/"junction.sam","discard mapped reads")
    c.atomic_json(seq.parent/"assessment.json",{"fasta":str(seq),"source_path":str(sample/"011_novoplasty/attempt_001/Circularized_assembly_1.fasta")})
    put(seq.parent/"bold_evidence.json",'{"status":"confirmed_target"}')
    put(sample/"qc/R1.fq.gz","discard trimmed reads")
    put(sample/"qc/fastp.json",'{"passed":123}')
    put(sample/"qc/fastp.html","QC report")
    put(sample/"logs/mtgrasp.log","native failure details")
    put(root/"bold_reference/coi.00.nsq","discard BOLD index")
    put(root/"bold_reference/source.json",'{"provenance":"retained"}')
    return root,sample,raw


def test_compacts_completed_run_and_flattens_itsme_without_losing_evidence(completed):
    root,sample,raw=completed
    c.compact_run(root)
    final=sample/"final"
    assert (final/"itsme_complete_ITS.fasta").read_text()==">ITS\nGGGG\n"
    assert (final/"itsme_master_summary.csv").exists()
    assert (final/"itsme_master_summary.tsv").exists()
    assert (final/"itsme_run_summary.txt").exists()
    assert (final/"itsme_assembly_graph.fastg").exists()
    assert (final/"itsme_fastp.json").exists()
    assert not list(final.glob("itsme_*frontier*"))
    assert not list(final.glob("*intermediate*"))
    assert (final/"wgs_contigs.fasta").exists()
    assert (final/"wgs_assembly_graph_with_scaffolds.gfa").exists()
    assert (final/"mitos_result.gff").exists()
    assert (final/"getorganelle_assembly_graph.fastg").exists()
    assert not any(c.STAGE_RE.fullmatch(p.name) for p in sample.iterdir())
    assert not list(sample.rglob("*.fq.gz")) and not list(sample.rglob("*.bam"))
    assert not list(sample.rglob("*.sam"))
    assert raw.read_text()=="original reads"
    assert (sample/"qc/fastp.html").exists()
    assert (sample/"logs/wgs_spades.log").read_text().startswith("assembly timed out")
    assert (sample/"logs/mtgrasp.log").exists()
    assert not (root/"bold_reference/coi.00.nsq").exists()
    assert (root/"bold_reference/source.json").exists()
    summary=c.read_json(final/"summary.json")
    assert summary["itsme_output"]==str(final) and summary["wgs_ok"] is False
    assert c.read_json(final/"annotation_summary.json")["gff"]==str(final/"mitos_result.gff")
    candidate=c.read_json(sample/"candidates/011_novoplasty_001/assessment.json")
    assert Path(candidate["fasta"]).is_file() and Path(candidate["source_path"]).is_file()
    assert c.read_json(sample/"cleanup.json")["bytes_removed"]>1000
    c.validate_retained(sample)
    c.compact_run(root)  # Idempotent after intermediates disappear.
    assert c.read_json(root/"last_invocation.json")["results"][0]["wgs_ok"] is False


def test_cleanup_dry_run_makes_no_changes(completed):
    root,sample,raw=completed
    before={str(p):p.read_bytes() for p in root.rglob("*") if p.is_file()}
    reports=c.compact_run(root,dry_run=True)
    assert reports[0]["bytes_removed"]>0
    after={str(p):p.read_bytes() for p in root.rglob("*") if p.is_file() and p.name!=".pipeline.lock"}
    assert before==after


def test_keep_intermediates_exports_results_without_deleting(completed):
    root,sample,raw=completed
    c.compact_sample(sample,remove=False)
    assert (sample/"final/itsme_complete_ITS.fasta").exists()
    assert (sample/"qc/R1.fq.gz").exists()
    assert (sample/"03_wgs/attempt_001/output/K55/huge.bin").exists()
    assert not c.read_json(sample/"cleanup.json")["compacted"]
    c.compact_sample(sample)  # A later requested cleanup can finish compaction.
    assert not (sample/"qc/R1.fq.gz").exists()


def test_active_pipeline_lock_prevents_cleanup(completed):
    root,sample,raw=completed
    with (root/".pipeline.lock").open("w") as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        with pytest.raises(ValueError,match="still running"):
            c.compact_run(root)
    assert (sample/"qc/R1.fq.gz").exists()


def test_export_failure_never_deletes_work_files(completed,monkeypatch):
    root,sample,raw=completed
    def fail(*args,**kwargs):
        raise OSError("disk full")
    monkeypatch.setattr(c,"export_file",fail)
    with pytest.raises(OSError,match="disk full"):
        c.compact_sample(sample)
    assert (sample/"qc/R1.fq.gz").exists()
    assert (sample/"03_wgs/attempt_001/output/K55/huge.bin").exists()


def test_symlinks_never_traverse_or_delete_external_data(completed,tmp_path):
    root,sample,raw=completed
    external=put(tmp_path/"external/valuable.fasta",">valuable\nACGT\n")
    (sample/"012_getorganelle/attempt_001/external").symlink_to(external.parent,target_is_directory=True)
    (sample/"candidates/external").symlink_to(external.parent,target_is_directory=True)
    c.compact_run(root)
    assert external.read_text()==">valuable\nACGT\n"
    assert not (sample/"final/getorganelle_valuable.fasta").exists()
    assert raw.exists()


def test_protected_inputs_inside_stage_and_unknown_directories_are_untouched(completed):
    root,sample,raw=completed
    protected=put(sample/"03_wgs/attempt_001/user_original.fq.gz","original in stage")
    config=c.read_json(sample/"run_config.json")
    config["reads"].append({"path":str(protected)})
    c.atomic_json(sample/"run_config.json",config)
    unknown=put(sample/"user_data/keep.fastq.gz","unmanaged")
    c.compact_run(root)
    assert protected.read_text()=="original in stage" and unknown.read_text()=="unmanaged"


def test_colliding_native_filenames_are_preserved_separately(completed):
    root,sample,raw=completed
    put(sample/"012_getorganelle/attempt_002/output/a.path_sequence.fasta",">different\nAAAA\n")
    c.compact_sample(sample)
    files=list((sample/"final").glob("getorganelle_*a.path_sequence.fasta"))
    assert len(files)==2 and len({p.read_text() for p in files})==2


def test_missing_retained_file_is_detected_before_reusing_compacted_results(completed):
    root,sample,raw=completed
    c.compact_sample(sample)
    (sample/"final/itsme_complete_ITS.fasta").unlink()
    with pytest.raises(ValueError,match="missing or changed"):
        c.compact_sample(sample)


def test_unrecognized_directory_is_rejected_without_changes(tmp_path):
    valuable=put(tmp_path/"not_a_pipeline/reads.fq.gz","original")
    with pytest.raises(ValueError,match="Not a completed"):
        c.compact_run(valuable.parent)
    assert valuable.read_text()=="original"


def test_gff_feature_table_preserves_calls_and_decodes_gene_names(completed):
    root,sample,raw=completed
    put(sample/"05_mitos2/attempt_001/output/result.gff",
        "##gff-version 3\nseq\tMITOS2\tCDS\t2\t600\t.\t-\t0\tName=cox1;ID=gene1\n"
        "seq\tMITOS2\ttRNA\t700\t760\t.\t+\t.\tName=tRNA%20Lys\n")
    c.compact_sample(sample)
    with (sample/"final/mitos_result_features.tsv").open() as handle:
        rows=list(csv.DictReader(handle,delimiter="\t"))
    assert rows[0]["gene"]=="cox1" and rows[0]["strand"]=="-" and rows[0]["start"]=="2"
    assert rows[1]["gene"]=="tRNA Lys" and rows[1]["feature_type"]=="tRNA"
    c.validate_retained(sample)


def test_linked_metadata_is_rejected_before_deletion(completed,tmp_path):
    root,sample,raw=completed
    external=put(tmp_path/"external_metadata.tsv","precious data")
    (sample/"final/outputs_manifest.tsv").symlink_to(external)
    with pytest.raises(ValueError,match="linked metadata"):
        c.compact_sample(sample)
    assert external.read_text()=="precious data"
    assert (sample/"qc/R1.fq.gz").exists()


def test_old_reference_unavailable_itsme_moves_into_parent_final_and_preserves_status(completed):
    root,sample,raw=completed
    fallback=root/"unresolved_references_itsme"
    child=fallback/"NOREF"
    c.atomic_json(child/"run_config.json",{"sample":"NOREF","reads":[]})
    c.atomic_json(child/"final/summary.json",{"sample":"NOREF","pipeline":"itsme","status":"completed","itsme_ok":True})
    c.atomic_json(fallback/"last_invocation.json",{"results":[{"sample":"NOREF","status":"completed"}]})
    put(child/"itsme/attempt_001/output/final/complete_ITS.fasta",">NOREF_ITS\nACGT\n")
    parent=c.read_json(root/"last_invocation.json")
    parent["results"].append({"sample":"NOREF","pipeline":"both","status":"reference_unavailable","itsme_ok":True})
    c.atomic_json(root/"last_invocation.json",parent)
    c.compact_run(root)
    final=root/"NOREF/final"
    assert (final/"itsme_complete_ITS.fasta").exists()
    assert c.read_json(final/"summary.json")["status"]=="reference_unavailable"
    c.validate_retained(final.parent)
    c.compact_run(root)


def test_internal_final_product_symlink_is_materialized_before_stage_deletion(completed):
    root,sample,raw=completed
    native=sample/"02_itsme/attempt_001/output"
    target=put(native/"assembly/locus.fasta",">linked_locus\nTTTT\n")
    link=native/"final/linked_locus.fasta"
    link.symlink_to(target)
    c.compact_sample(sample)
    exported=sample/"final/itsme_linked_locus.fasta"
    assert exported.read_text()==">linked_locus\nTTTT\n" and not exported.is_symlink()
    assert not target.exists()
    c.validate_retained(sample)

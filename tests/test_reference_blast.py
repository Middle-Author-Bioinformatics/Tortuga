from pathlib import Path
import sys

import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from mito_pipeline import parser
from reference_blast import rank_hits


def args():
    return parser().parse_args([])


def reference(identifier,organism,code=5,phylum="Annelida"):
    return {"id":identifier,"path":"/refs/"+identifier+".gb","accession":identifier,
            "organism":organism,"taxonomy":["Metazoa",phylum],"length_bp":15000,"mt_code":code,"error":""}


def hit(query,subject,score=100,identity=95,start=1,end=100):
    return dict(zip("qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore".split(),
                    [query,subject,str(identity),str(end-start+1),"150","15000",str(start),str(end),"1","100","1e-30",str(score)]))


def test_sequence_support_overrides_name_based_preference():
    refs=[reference("exact","Target species"),reference("relative","Other species")]
    rows=[hit(f"read{i}/1",sid,score) for i in range(3) for sid,score in [("exact",100),("relative",150)]]
    choice,ranked=rank_hits(rows,refs,args())
    assert choice["accession"]=="relative" and choice["match_level"]=="blast_reads"
    assert choice["bitscore"]==450 and choice["supporting_templates"]==3
    assert [r["accession"] for r in ranked]==["relative","exact"]


def test_mates_and_overlapping_hsps_cannot_inflate_support():
    refs=[reference("ref","Target species")]
    rows=[hit("one/1","ref"),hit("one/2","ref"),hit("one/1","ref",score=90)]
    assert rank_hits(rows,refs,args())[0] is None
    rows += [hit("two/1","ref"),hit("three/1","ref")]
    choice,_=rank_hits(rows,refs,args())
    assert choice["supporting_templates"]==3 and choice["bitscore"]==400


def test_phylum_code_and_quality_filters_apply_before_reference_ranking():
    refs=[reference("right","Target"),reference("wrong_phylum","Other",phylum="Mollusca"),reference("wrong_code","Other",code=2)]
    rows=[hit(f"r{i}/1",r["id"],score=200 if r["id"]!="right" else 100) for i in range(3) for r in refs]
    choice,_=rank_hits(rows,refs,args(),lineage={"phylum":"Annelida"},code=5)
    assert choice["accession"]=="right"
    for row in rows:
        row["pident"]="60"
    assert rank_hits(rows,refs,args())[0] is None


def test_conflicting_near_best_translation_tables_require_override():
    refs=[reference("code5","Target",5),reference("code2","Other",2)]
    rows=[hit(f"r{i}/1",r["id"]) for i in range(3) for r in refs]
    assert rank_hits(rows,refs,args())[0] is None
    assert rank_hits(rows,refs,args(),code=5)[0]["mt_code"]==5


def test_single_contig_requires_substantial_aligned_sequence():
    refs=[reference("ref","Target")]
    assert rank_hits([hit("candidate","ref")],refs,args(),contig=True)[0] is None
    choice,_=rank_hits([hit("candidate","ref",end=1000)],refs,args(),contig=True)
    assert choice["match_level"]=="blast_contig" and choice["aligned_bp"]==1000


def test_rejected_and_oversized_reference_records_cannot_win():
    refs=[reference("good","Target"),{**reference("bad","Other"),"error":"ambiguous CDS"},
          {**reference("long","Other"),"length_bp":100000}]
    rows=[hit(f"r{i}/1",r["id"],score=100 if r["id"]=="good" else 200) for i in range(3) for r in refs]
    assert rank_hits(rows,refs,args())[0]["accession"]=="good"

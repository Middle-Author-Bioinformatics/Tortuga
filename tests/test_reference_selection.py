import csv
from pathlib import Path
import sys

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mito_pipeline import reference_data
from reference_selection import load_catalog, select_reference, taxonomy_from_row


def test_hierarchy_normalizes_rank_columns_and_species_epithets():
    lineage=taxonomy_from_row({"Genus":"Target","Species":"species","Phylum":"Annelida",
                               "Class":"Polychaeta","Order":"NA","Family":None})
    assert lineage=={"species":"Target species","genus":"Target","class":"Polychaeta","phylum":"Annelida"}
    assert taxonomy_from_row({"species":"Target species"})["species"]=="Target species"


def test_family_fallback_uses_provided_hierarchy_only(tmp_path):
    make_reference(tmp_path,"RELATIVE","Relative species")
    catalog=load_catalog(tmp_path,reference_data)
    lineage={"genus":"Missinggenus","family":"Targetidae"}
    result=select_reference(catalog,"Missinggenus species",match="hierarchy",lineage=lineage)
    assert result["match_rank"]=="family" and result["matched_taxon"]=="Targetidae"
    assert result["match_level"]=="same_family"
    with pytest.raises(ValueError,match="No suitable"):
        select_reference(catalog,"Missinggenus species",match="genus",lineage=lineage)
    with pytest.raises(ValueError,match="No suitable"):
        select_reference(catalog,"Missinggenus species",match="hierarchy",lineage={"family":"Unrelatedidae"})


def test_closer_genus_precedes_better_annotated_family_record(tmp_path):
    make_reference(tmp_path,"GENUS","Target other")
    make_reference(tmp_path,"FAMILY","Relative species",gene_count=4)
    catalog=load_catalog(tmp_path,reference_data)
    for row in catalog:
        if row["accession"]=="FAMILY":
            row["taxonomy"].remove("Target")
            row["taxonomy"].append("Relative")
    chosen=select_reference(catalog,"Target species",match="hierarchy",
                            lineage={"genus":"Target","family":"Targetidae"})
    assert chosen["accession"]=="GENUS" and chosen["match_rank"]=="genus"


def test_order_fallback_requires_explicit_limit_and_phylum_agreement(tmp_path):
    make_reference(tmp_path,"ORDER")
    catalog=load_catalog(tmp_path,reference_data)
    catalog[0]["taxonomy"].extend(["Annelida","Targetida"])
    lineage={"genus":"Missing","family":"Missingidae","order":"Targetida","phylum":"Annelida"}
    with pytest.raises(ValueError,match="No suitable"):
        select_reference(catalog,"Missing species",match="hierarchy",lineage=lineage)
    chosen=select_reference(catalog,"Missing species",match="hierarchy",lineage=lineage,min_rank="order")
    assert chosen["match_rank"]=="order"
    with pytest.raises(ValueError,match="No suitable"):
        select_reference(catalog,"Missing species",match="hierarchy",
                         lineage={**lineage,"phylum":"Mollusca"},min_rank="order")


def test_family_code_conflict_does_not_silently_choose_a_translation_table(tmp_path):
    make_reference(tmp_path,"CODE5",code=5)
    make_reference(tmp_path,"CODE2",code=2)
    catalog=load_catalog(tmp_path,reference_data)
    lineage={"genus":"Missing","family":"Targetidae"}
    with pytest.raises(ValueError,match="unique mitochondrial code"):
        select_reference(catalog,"Missing species",match="hierarchy",lineage=lineage)
    chosen=select_reference(catalog,"Missing species",mt_code=5,match="hierarchy",lineage=lineage)
    assert chosen["accession"]=="CODE5" and chosen["match_rank"]=="family"


def make_reference(folder, accession, organism="Target species", code=5,
                   organelle="mitochondrion", complete=True, gene_count=3):
    folder.mkdir(exist_ok=True)
    rec = SeqRecord(Seq("ACGT"*3000), id=accession,
                    description=f"{organism} {organelle}, {'complete genome' if complete else 'partial sequence'}")
    rec.annotations.update(molecule_type="DNA", organism=organism,
                           taxonomy=["Eukaryota","Metazoa","Targetidae","Target"],topology="circular")
    rec.features.append(SeqFeature(FeatureLocation(0,12000),type="source",qualifiers={"organelle":[organelle]}))
    for i, gene in enumerate(["COX1","ND1","CYTB","ND2"][:gene_count]):
        qualifiers={"gene":[gene],"translation":["M"+"A"*199]}
        if code is not None:
            qualifiers["transl_table"]=[str(code)]
        rec.features.append(SeqFeature(FeatureLocation(i*1000,i*1000+600),type="CDS",qualifiers=qualifiers))
    rec.features.append(SeqFeature(FeatureLocation(8000,9000),type="rRNA",qualifiers={"gene":["rrnS"]}))
    path=folder/f"{accession}.gb"
    SeqIO.write(rec,path,"genbank")
    return path


def test_exact_species_precedes_better_annotated_relative(tmp_path):
    exact=make_reference(tmp_path,"ZZ_EXACT")
    make_reference(tmp_path,"AA_RELATIVE","Target other",gene_count=4)
    catalog=load_catalog(tmp_path,reference_data)
    chosen=select_reference(catalog,"  TARGET   species ")
    assert chosen["path"]==str(exact)
    assert chosen["match_level"]=="exact_taxon" and chosen["mt_code"]==5


def test_genus_fallback_is_reported_and_can_be_disabled(tmp_path):
    make_reference(tmp_path,"RELATIVE","Target other")
    catalog=load_catalog(tmp_path,reference_data)
    assert select_reference(catalog,"Target species")["match_level"]=="same_genus"
    with pytest.raises(ValueError,match="No suitable"):
        select_reference(catalog,"Target species",match="exact")
    with pytest.raises(ValueError,match="No suitable"):
        select_reference(catalog,"Unrelated species")
    assert select_reference(catalog,"Targetidae")["match_level"]=="exact_taxon"


def test_excludes_wrong_organelle_partial_and_invalid_gene_annotations(tmp_path):
    make_reference(tmp_path,"GOOD")
    make_reference(tmp_path,"CHLORO",organelle="chloroplast")
    make_reference(tmp_path,"PARTIAL",complete=False)
    make_reference(tmp_path,"SHORT",gene_count=2)
    (tmp_path/"broken.gb").write_text("not GenBank\n")
    catalog=load_catalog(tmp_path,reference_data)
    assert len([row for row in catalog if not row["error"]])==1
    assert select_reference(catalog,"Target species")["accession"]=="GOOD"


def test_code_conflicts_require_explicit_choice_and_length_filter(tmp_path):
    make_reference(tmp_path,"CODE5",code=5)
    make_reference(tmp_path,"CODE2",code=2)
    catalog=load_catalog(tmp_path,reference_data)
    with pytest.raises(ValueError,match="unique mitochondrial code"):
        select_reference(catalog,"Target species")
    assert select_reference(catalog,"Target species",mt_code=5)["accession"]=="CODE5"
    with pytest.raises(ValueError,match="No suitable"):
        select_reference(catalog,"Target species",mt_code=5,max_bp=11000)


def test_missing_code_requires_sheet_code(tmp_path):
    make_reference(tmp_path,"NO_CODE",code=None)
    catalog=load_catalog(tmp_path,reference_data)
    with pytest.raises(ValueError,match="fill mt_code"):
        select_reference(catalog,"Target species")
    assert select_reference(catalog,"Target species",mt_code=5)["accession"]=="NO_CODE"


def test_deterministic_ties_and_incremental_cache(tmp_path):
    b=make_reference(tmp_path,"B")
    make_reference(tmp_path,"A")
    calls=[]
    def validate(record,code):
        calls.append(record.id)
        return reference_data(record,code)
    first=load_catalog(tmp_path,validate)
    assert len(calls)==2
    assert select_reference(first,"Target species")["accession"]=="A"
    load_catalog(tmp_path,validate)
    assert len(calls)==2
    with b.open("a") as handle:
        handle.write("\n")
    load_catalog(tmp_path,validate)
    assert len(calls)==3
    b.unlink()
    assert len(load_catalog(tmp_path,validate))==1


def test_real_orbiniid_record_accepts_boundary_annotations_and_rrna_products(tmp_path):
    from mito_pipeline import reference_boundary_warnings, reference_files, Sample
    source=Path(__file__).parent/"data/NC_007933.1.gb"
    record=SeqIO.read(source,"genbank")
    proteins,rnas,labels,seed=reference_data(record,5)
    assert len(proteins)==13 and {name for name,seq in rnas}=={"rrns","rrnl"}
    assert len(seed)==1563 and set(seed.upper())<=set("ACGT")
    warnings=reference_boundary_warnings(record)
    assert {row["gene"] for row in warnings}=={"cox1","cox3"}
    catalog=load_catalog(source.parent,reference_data,tmp_path/"index.sqlite")
    chosen=select_reference(catalog,"Theodisca dendritica",match="hierarchy",
                            lineage={"family":"Orbiniidae","phylum":"Annelida"})
    assert chosen["accession"]=="NC_007933.1" and chosen["match_rank"]=="family"
    assert chosen["annotation_warnings"]==warnings
    refs=reference_files(Sample("test",source,source,"Theodisca dendritica",5,source),tmp_path/"refs")
    assert refs["annotation_warnings"]==warnings
    import json
    assert json.loads((tmp_path/"refs/reference_quality.json").read_text())["annotation_warnings"]==warnings


def test_uncertain_boundaries_do_not_allow_ambiguous_sequence_or_missing_translation():
    from Bio.SeqFeature import WithinPosition
    record=SeqIO.read(Path(__file__).parent/"data/NC_007933.1.gb","genbank")
    first=next(f for f in record.features if f.type=="CDS")
    original=record.seq
    record.seq=Seq("N"+str(original)[1:])
    with pytest.raises(ValueError,match="ambiguous DNA"):
        reference_data(record,5)
    record.seq=original
    translation=first.qualifiers.pop("translation")
    with pytest.raises(ValueError,match="without an annotated translation"):
        reference_data(record,5)
    first.qualifiers["translation"]=translation
    first.location=FeatureLocation(WithinPosition(0,0,3),1563,strand=1)
    with pytest.raises(ValueError,match="unsupported uncertain coordinates"):
        reference_data(record,5)


def test_exact_boundaries_preferred_within_same_taxonomic_tier(tmp_path):
    from Bio.SeqFeature import BeforePosition
    fuzzy=make_reference(tmp_path,"A_FUZZY")
    make_reference(tmp_path,"Z_EXACT")
    record=SeqIO.read(fuzzy,"genbank")
    record.features[1].location=FeatureLocation(BeforePosition(0),600,strand=1)
    SeqIO.write(record,fuzzy,"genbank")
    catalog=load_catalog(tmp_path,reference_data)
    assert select_reference(catalog,"Target species")["accession"]=="Z_EXACT"


def test_old_rejection_cache_is_rebuilt_and_rejected_metadata_retained(tmp_path):
    import json
    import sqlite3
    path=make_reference(tmp_path,"OLD")
    catalog=load_catalog(tmp_path,reference_data)
    index=tmp_path/"tortuga_reference_index.sqlite"
    with sqlite3.connect(index) as db:
        db.execute("UPDATE entries SET version=1,data=?",(json.dumps({"path":str(path),"error":"old rejection"}),))
    catalog=load_catalog(tmp_path,reference_data)
    assert not catalog[0]["error"]
    record=SeqIO.read(path,"genbank")
    record.seq=Seq("N"+str(record.seq)[1:])
    SeqIO.write(record,path,"genbank")
    rejected=load_catalog(tmp_path,reference_data)[0]
    assert "ambiguous DNA" in rejected["error"]
    assert rejected["organism"]=="Target species" and "Targetidae" in rejected["taxonomy"]

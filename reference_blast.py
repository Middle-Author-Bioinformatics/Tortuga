"""Bounded sample-read/contig BLAST evidence for local assembly-reference choice.

This ranks eligible local genomes by sequence support. It is not a phylogenetic
closest-relative claim or independent identification of the assembled organism.
"""
from pathlib import Path
import json

FIELDS="qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore"


def eligible_references(catalog,args,lineage=None,code=None):
    from reference_selection import normalized
    phylum=(lineage or {}).get("phylum")
    return [row for row in catalog if not row.get("error")
            and args.min_genome_bp<=row["length_bp"]<=args.max_genome_bp
            and row.get("mt_code") is not None
            and (code is None or row["mt_code"]==code)
            and (not phylum or normalized(phylum) in {normalized(x) for x in row["taxonomy"]})]


def rank_hits(rows,references,args,lineage=None,code=None,contig=False):
    """Deduplicate overlapping HSPs and mate support; report all qualifying ranks."""
    eligible={r["id"]:r for r in eligible_references(references,args,lineage,code)}
    grouped={}
    for hit in rows:
        sid=hit["sseqid"].removeprefix("lcl|")
        if sid not in eligible or float(hit["pident"])<args.reference_min_identity or int(hit["length"])<args.reference_min_aligned_bp:
            continue
        grouped.setdefault((sid,hit["qseqid"]),[]).append(hit)
    scores={}
    for (sid,qid),hits in grouped.items():
        accepted=[]
        for hit in sorted(hits,key=lambda r:-float(r["bitscore"])):
            lo,hi=sorted((int(hit["qstart"]),int(hit["qend"])))
            if any(not (hi<a or lo>b) for a,b in accepted):
                continue
            accepted.append((lo,hi))
            score=scores.setdefault(sid,{"bitscore":0.,"aligned_bp":0,"identity_sum":0.,"templates":set()})
            aligned=hi-lo+1
            score["bitscore"]+=float(hit["bitscore"])
            score["aligned_bp"]+=aligned
            score["identity_sum"]+=float(hit["pident"])*aligned
            score["templates"].add(qid.rsplit("/",1)[0])
    ranked=[]
    for sid,score in scores.items():
        if len(score["templates"])<(1 if contig else args.reference_min_templates):
            continue
        if contig and score["aligned_bp"]<max(400,args.reference_min_aligned_bp):
            continue
        ranked.append({**eligible[sid],"bitscore":score["bitscore"],"aligned_bp":score["aligned_bp"],
                       "supporting_templates":len(score["templates"]),
                       "mean_identity":score["identity_sum"]/score["aligned_bp"]})
    ranked.sort(key=lambda r:(-r["bitscore"],-r["mean_identity"],-r["supporting_templates"],r["accession"],r["path"]))
    if not ranked:
        return None,[]
    top=ranked[0]
    tied=[row for row in ranked if row["bitscore"]>=top["bitscore"]*.99]
    if code is None and len({row["mt_code"] for row in tied})>1:
        return None,ranked  # No implicit translation-table guess across near ties.
    return {**top,"match_level":"blast_contig" if contig else "blast_reads",
            "match_rank":"sequence","matched_taxon":top["organism"],"eligible_candidates":len(ranked)},ranked


class ReferenceSearch:
    def __init__(self,catalog,args,runner):
        self.args,self.runner=args,runner
        self.root=args.outdir/"reference_search"
        self.root.mkdir(parents=True,exist_ok=True)
        self.references=[]
        self.database=self.root/"mitochondria"
        self.ready=False
        from Bio import SeqIO
        from mito_pipeline import atomic_json,stamp
        eligible=sorted(eligible_references(catalog,args),key=lambda r:r["path"])
        for i,row in enumerate(eligible):
            self.references.append({**row,"id":f"ref{i:06d}"})
        if not self.references:
            atomic_json(self.root/"source.json",{"references":[],"reason":"No eligible annotated animal references"})
            return
        source={"references":[{"id":r["id"],"file":stamp(r["path"])} for r in self.references]}
        self.source=source
        previous=self.root/"source.json"
        if previous.exists() and json.loads(previous.read_text())!=source:
            raise ValueError("Local BLAST reference collection changed; use a new output directory")
        atomic_json(previous,source)
        self.ready=None

    def prepare(self):
        from Bio import SeqIO
        from mito_pipeline import atomic_json
        if self.ready is not None:
            return self.ready
        fasta=self.root/"mitochondria.fasta"
        if not fasta.exists():
            with fasta.open("w") as handle:
                for row in self.references:
                    record=SeqIO.read(row["path"],"genbank")
                    handle.write(f">{row['id']}\n{record.seq}\n")
        self.ready=self.runner.run("reference_index","makeblastdb",["-in",fasta,"-dbtype","nucl","-parse_seqids","-out",self.database],
                                  self.root,self.args.validation_timeout*60,[self.database.with_suffix(".nsq")])
        return self.ready

    def query(self,path,sample,folder,lineage,code,contig=False):
        from mito_pipeline import blast_rows,atomic_json,write_table,digest_file
        import hashlib
        folder.mkdir(parents=True,exist_ok=True)
        hits=folder/"reference_hits.tsv"
        signature=hashlib.sha256(json.dumps({"query":digest_file(path),"source":getattr(self,"source",{}),
                    "lineage":lineage,"code":code,"contig":contig,"identity":self.args.reference_min_identity,
                    "alignment":self.args.reference_min_aligned_bp,"templates":self.args.reference_min_templates},sort_keys=True).encode()).hexdigest()
        selection=folder/"selection.json"
        if self.args.resume and selection.exists():
            previous=json.loads(selection.read_text())
            if previous.get("signature")==signature and previous.get("search_succeeded"):
                return previous.get("selected")
        if not self.prepare():
            atomic_json(selection,{"status":"unavailable","reason":"Reference database unavailable","signature":signature,"selected":None,"search_succeeded":not bool(self.references)})
            return None
        stage="reference_"+sample+"_"+hashlib.sha256(str(folder).encode()).hexdigest()[:12]
        ok=self.runner.run(stage,"blastn",["-query",path,"-db",self.database,"-task","blastn",
                    "-evalue","1e-10","-max_target_seqs","100","-num_threads",self.args.threads,
                    "-outfmt","6 "+FIELDS,"-out",hits],folder,self.args.validation_timeout*60,[hits])
        if self.args.dry_run:
            return None
        choice,ranked=rank_hits(blast_rows(hits) if ok else [],self.references,self.args,lineage,code,contig)
        write_table(folder/"reference_ranking.tsv",ranked,
                    ["accession","organism","path","mt_code","bitscore","aligned_bp","mean_identity","supporting_templates"])
        atomic_json(folder/"selection.json",{"status":"selected" if choice else "unavailable","selected":choice,"signature":signature,"search_succeeded":bool(ok),
                    "reason":"sequence support" if choice else "Conflicting genetic codes among near-best hits" if ranked else "No qualifying sequence support",
                    "query_kind":"contig" if contig else "bounded_read_prefix",
                    "read_pair_limit":None if contig else self.args.reference_read_pairs,
                    "min_identity":self.args.reference_min_identity,"min_aligned_bp":self.args.reference_min_aligned_bp,
                    "min_templates":1 if contig else self.args.reference_min_templates,"max_targets_per_query":100,
                    "note":"Best sequence-supported eligible local reference; not a taxonomic identification or phylogenetic closest-relative proof"})
        return choice

    def reads(self,name,pair,lineage,code):
        from mito_pipeline import fastq_records,template_id
        from itertools import islice,zip_longest
        folder=self.root/name
        folder.mkdir(exist_ok=True)
        query=folder/"read_queries.fasta"
        with query.open("w") as handle:
            for i,(left,right) in enumerate(islice(zip_longest(fastq_records(pair[0]),fastq_records(pair[1])),self.args.reference_read_pairs)):
                if left is None or right is None or template_id(left)!=template_id(right):
                    raise ValueError(f"{name}: reference-search read mates differ")
                for mate,record in enumerate((left,right),1):
                    handle.write(f">template{i}/{mate}\n{record[1].strip()}\n")
        return self.query(query,name,folder,lineage,code)

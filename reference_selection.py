"""Select validated local animal mitochondrial references by annotated taxonomy.

This uses GenBank organism/lineage names, not a sequence similarity search.
Fallback follows supplied ranks, with an explicit specificity limit.
"""
import json
from pathlib import Path
import re
import sqlite3

INDEX_VERSION = 1
SUFFIXES = {".gb", ".gbk", ".gbff", ".genbank"}
RANKS = ("species", "genus", "family", "order", "class", "phylum", "kingdom")


def normalized(name):
    return " ".join(str(name).split()).casefold()


def taxonomy_from_row(row):
    """Accept lower/title-case rank columns and species epithets or binomials."""
    fields = {str(k).strip().lower(): str(v or "").strip() for k,v in row.items() if k is not None}
    lineage = {rank: fields[rank] for rank in RANKS
               if fields.get(rank) and fields[rank].casefold() not in {"na","n/a","none"}}
    species = lineage.get("species", "")
    if species and len(species.split()) == 1:
        if lineage.get("genus"):
            lineage["species"] = lineage["genus"] + " " + species
        else:
            lineage.pop("species")
    return lineage


def describe(path, validator):
    from Bio import SeqIO
    records = list(SeqIO.parse(path, "genbank"))
    if len(records) != 1:
        raise ValueError("Reference must contain one record")
    record = records[0]
    source = [f for f in record.features if f.type == "source"]
    organelles = [v for f in source for v in f.qualifiers.get("organelle", [])]
    if not any("mitochond" in normalized(x) for x in organelles + [record.description]):
        raise ValueError("Record is not annotated as mitochondrial")
    taxonomy = record.annotations.get("taxonomy", [])
    if not {normalized(x) for x in taxonomy} & {"metazoa", "animalia"}:
        raise ValueError("Automatic selection requires an annotated animal lineage (animal_mt workflow)")
    if re.search(r"\b(partial|incomplete)\b", record.description, re.I):
        raise ValueError("Record is described as partial or incomplete")
    codes = {int(x) for f in record.features if f.type == "CDS"
             for x in f.qualifiers.get("transl_table", [])}
    if len(codes) > 1:
        raise ValueError("Reference CDS annotations specify conflicting genetic codes")
    code = next(iter(codes)) if codes else None
    proteins, rnas, _, _ = validator(record, code or 1)
    return {"path": str(path), "accession": record.id,
            "organism": record.annotations.get("organism", ""),
            "taxonomy": taxonomy, "mt_code": code,
            "length_bp": len(record.seq), "cds_count": len(proteins),
            "rrna_count": len(rnas),
            "circular": record.annotations.get("topology", "").lower() == "circular",
            "complete": "complete genome" in record.description.lower(),
            "error": ""}


def load_catalog(directory, validator, index_path=None):
    """Incremental cache: inspect new/changed files and prune removed entries."""
    directory = Path(directory).expanduser().resolve()
    if not directory.is_dir():
        raise ValueError(f"GenBank directory not found: {directory}; use --genbank-dir")
    paths = sorted(p.resolve() for p in directory.rglob("*")
                   if p.is_file() and p.suffix.lower() in SUFFIXES)
    if not paths:
        raise ValueError(f"No .gb/.gbk/.gbff/.genbank files found in {directory}")
    index_path = Path(index_path) if index_path else directory/"tortuga_reference_index.sqlite"
    index_path.parent.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(index_path, timeout=60) as db:
        db.execute("CREATE TABLE IF NOT EXISTS entries (path TEXT PRIMARY KEY, size INTEGER, mtime INTEGER, version INTEGER, data TEXT)")
        old = {row[0]: row[1:] for row in db.execute("SELECT path,size,mtime,version,data FROM entries")}
        rows = []
        print(f"Indexing {len(paths)} GenBank files (cached unchanged files are reused)", flush=True)
        for n, path in enumerate(paths, 1):
            stat = path.stat()
            cached = old.pop(str(path), None)
            if cached and cached[:3] == (stat.st_size, stat.st_mtime_ns, INDEX_VERSION):
                row = json.loads(cached[3])
            else:
                try:
                    row = describe(path, validator)
                except (ValueError, TypeError, KeyError, AttributeError, IndexError, UnicodeError) as error:
                    row = {"path": str(path), "error": str(error)}
                db.execute("INSERT OR REPLACE INTO entries VALUES (?,?,?,?,?)",
                           (str(path), stat.st_size, stat.st_mtime_ns, INDEX_VERSION, json.dumps(row)))
            rows.append(row)
            if n % 1000 == 0:
                print(f"  Indexed {n}/{len(paths)} records", flush=True)
        for path in old:
            db.execute("DELETE FROM entries WHERE path=?", (path,))
    good = sum(not r.get("error") for r in rows)
    print(f"GenBank catalog: {good} usable animal mitochondrial references; {len(rows)-good} excluded", flush=True)
    return rows


def select_reference(catalog, target, mt_code=None, match="genus", min_bp=10000, max_bp=30000,
                     lineage=None, min_rank="family"):
    """Prefer the nearest available named rank, then annotation quality."""
    query = normalized(target)
    if not query:
        raise ValueError("A target_taxon is required for automatic reference selection")
    lineage = taxonomy_from_row(lineage or {})
    genus = lineage.get("genus") or (query.split()[0] if len(query.split()) >= 2 else None)
    exact_rank = next((rank for rank,value in lineage.items() if normalized(value)==query),"target_taxon")
    tiers = [("exact_taxon", exact_rank, str(target).strip())]
    if match in {"genus","hierarchy"} and genus and normalized(genus) != query:
        tiers.append(("same_genus", "genus", genus))
    if match == "hierarchy":
        if min_rank not in {"family","order","class","phylum"}:
            raise ValueError("Reference minimum rank must be family, order, class, or phylum")
        for rank in RANKS[RANKS.index("family"):RANKS.index(min_rank)+1]:
            if lineage.get(rank) and normalized(lineage[rank]) != query:
                tiers.append((f"same_{rank}", rank, lineage[rank]))
    matches = []
    for row in catalog:
        if row.get("error") or not min_bp <= row["length_bp"] <= max_bp:
            continue
        if mt_code is not None and row["mt_code"] not in {None, mt_code}:
            continue
        names = {normalized(x) for x in row["taxonomy"] + [row["organism"]]}
        # Prevent cross-phylum choices when the input supplies a phylum.
        if lineage.get("phylum") and normalized(lineage["phylum"]) not in names:
            continue
        for tier, (level,rank,taxon) in enumerate(tiers):
            if normalized(taxon) in names:
                matches.append((row,level,rank,taxon,tier))
                break
    if not matches:
        raise ValueError(f"No suitable local mitochondrial reference for {target!r} with the requested length/code settings. "
                         "Add an appropriate annotated GenBank record, supply the family/order hierarchy with an appropriate "
                         "--reference-min-rank, or set reference_gb explicitly. Keep target_taxon as the intended identity.")
    closest = min(item[4] for item in matches)
    matches = [item for item in matches if item[4] == closest]
    if mt_code is None:
        codes = {item[0]["mt_code"] for item in matches}
        if None in codes or len(codes) != 1:
            raise ValueError(f"Cannot infer a unique mitochondrial code for {target!r}; fill mt_code in the sample sheet")
    matches.sort(key=lambda item: (-int(item[0]["complete"]),
                                  -item[0]["cds_count"], -item[0]["rrna_count"],
                                  -int(item[0]["circular"]), item[0]["accession"], item[0]["path"]))
    row, level, rank, taxon, _ = matches[0]
    return {**row, "match_level": level,"match_rank":rank,"matched_taxon":taxon,
            "eligible_candidates": len(matches)}

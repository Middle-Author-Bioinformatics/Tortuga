# Hierarchical reference selection and standardized outputs (2.1.0)

The defaults are `--reference-match hierarchy --reference-min-rank kingdom`.
Reference selection reads the organism and lineage names in local annotated
GenBank records. It does not sample reads or run BLAST to choose a reference.
BLAST remains required for gene screening and independent COI/BOLD validation.

## Matching order

Try the target taxon, genus, family, order, class, phylum, then kingdom. Stop at
the first rank with eligible records. Ranks come from the sample sheet; missing
ranks are skipped. A genus can also be obtained from a binomial target taxon.
The sample-sheet kingdom `Animalia` matches GenBank `Metazoa`, and vice versa.
An explicit `reference_gb` overrides automatic matching. `--reference-match exact`
and `--reference-match genus` remain available for stricter matching.

`--reference-min-rank family|order|class|phylum|kingdom` limits broad fallback;
kingdom is the default. A provided phylum must agree for all tiers below kingdom.
Kingdom fallback can cross phyla: its match rank is reported clearly, and such a
reference can be distant and have an unrepresentative gene complement. It is an
assembly aid, not a prediction of the sample's identity or closest relative.
The expected `target_taxon` for BOLD validation is unchanged.

Only validated animal mitochondrial records with eligible annotations and
lengths (default 10,000–30,000 bp) qualify. Explicit `mt_code` excludes conflicting
reference codes. Without an explicit code, eligible records at the closest
matching rank must specify one unique code; ambiguity triggers reference-free
recovery instead of choosing a translation table. At one rank, prefer complete
annotations, fewer uncertain boundaries, more CDS/rRNA genes, circular topology,
and then accession/path for reproducible ties. Broader matching never bypasses
annotation, length or code checks. A match is therefore not guaranteed.

The selected accession/path, reference organism, matched rank/taxon, genetic code,
annotation warnings and eligible count are retained in `reference_selection.tsv`.
Unresolved choices are explained in `reference_failures.tsv`. The full supplied
hierarchy is preserved in `resolved_samples.tsv`.

## Without a suitable custom reference

Keep the sample in the normal cohort. Skip reference-dependent NOVOPlasty,
MitoFinder and mtGrasp. GetOrganelle is the first configured method that can run
using its installed animal_mt seed/label databases without a custom `-s` seed.
Unavailable methods are skipped; WGS remains available unless disabled. If
metaSPAdes is unavailable, installed MEGAHIT can be used for reference-free WGS.
Graph-based GetOrganelle can also use its installed labels without custom genes.
ITSME still runs in combined mode.

If the GetOrganelle animal mitochondrial database is missing, prepare it once:

```bash
/home/ark/miniconda3/bin/conda run -p /home/ark/miniconda3/envs/getorganelle \
  get_organelle_config.py -a animal_mt
```

Candidates without an annotated reference are searched directly against BOLD COI
for target validation. Gene completeness stays unknown and cannot satisfy the
default all-genes policy or a positive gene-count requirement. An explicit length
policy can accept a candidate under its remaining topology/target/quality gates.
MITOS2 now annotates the best retained candidate by default, including unclosed
or unresolved candidates, when `mt_code` is known. Annotation does not change
acceptance. If the code is unavailable, the FASTA is retained and annotation is
explicitly skipped; no code is guessed. Accepted candidates without a usable
annotation code retain a nonzero batch exit.

## Outputs and deployment

Each processed sample has `ITSME/` and `MitoLooper/` folders with branch statuses,
including failures or unrequested branches. Deliverables are FASTAs, annotations,
tables and assembly graphs. Common `final/summary.json` and
`final/outputs_manifest.tsv`, QC, logs and reference provenance remain alongside.
Default cleanup exports essentials before removing heavy stage files. Active
runs cannot be compacted; unfinished samples in inactive runs are skipped.
`--keep-intermediates` preserves work files for stage retries. Compact resume
reuses finalized outcomes; use a new output directory to change analysis settings.

Merge the PR, then run:

```bash
git -C /home/ark/MAB/bin/Tortuga pull --ff-only
```

See [CONTROLLER.md](CONTROLLER.md) for the full Belly command with both branches.
Existing completed results can be standardized without reassembly:

```bash
bash /home/ark/MAB/bin/Tortuga/compact_outputs.sh \
  /home/ark/MAB/NOAA/CA_intertidal/both_results_v205
```

# One controller for ITSME, MitoLooper, or both

## Compact outputs by default (2.0.7)

Completed samples now export client deliverables into `OUT/SAMPLE/final/` and
remove generated stage work files automatically. This includes samples with
failed mitochondrial recovery/WGS: their actual status and logs are retained.
Cleanup does not rerun tools, change acceptance criteria or turn failures into
successful biological recovery.

| Retained location | Contents |
| --- | --- |
| `SAMPLE/final/itsme_*.fasta`, `.csv`, `.tsv` | Native ITSME final sequences, master/partial summaries and annotation/support tables, directly alongside mito outputs. |
| `SAMPLE/final/mitos_*` | Native MITOS2 sequences, GFF/GenBank/BED products when present, plus `mitos_result_features.tsv` generated from the GFF. |
| `SAMPLE/final/*_assembly_graph*.fastg`, `.gfa` | Final available assembly graphs from ITSME, WGS and recovery tools. |
| `SAMPLE/final/wgs_contigs.fasta`, `wgs_scaffolds.fasta` | Final WGS sequence products when produced; summary flags still identify a failed/timed-out WGS run. |
| `SAMPLE/final/outputs_manifest.tsv` | Original source path, retained path, product kind and size. |
| `SAMPLE/qc/`, `SAMPLE/logs/` | QC summaries/HTML and wrapper/native logs, including failures and timeouts. |
| `SAMPLE/candidates/`, `candidates.tsv`, `stages.*`, `references/`, `run_config.json` | Small assessment, gene/taxon evidence, provenance and checkpoint reports. Candidate FASTA paths are updated to retained final files. |
| `SAMPLE/cleanup.json`, `cleanup_removed.tsv` | Cleanup counts, discarded intermediate bytes and removed-file audit. |

Generated GetOrganelle/ITSME recruited FASTQs, trimmed/subsampled FASTQs,
SAM/BAM mappings, SPAdes K* intermediates, assembly scratch, doubled candidate
subjects, temporary junction mappings and superseded stage copies are removed.
The copied BOLD BLAST index is removed only after the cohort finishes; its
provenance remains. External input reads and shared databases are not deleted.
Final graphs and contig FASTAs may themselves be large and are intentionally
retained. Nonempty final sequences/graphs are exported; empty FASTA placeholders
and intermediate recruitment bait/round FASTAs are omitted.

Exports finish before stage deletion begins. Same-name identical products share
one retained destination; different products receive distinct names. Native
logs keep their commands and historical working paths. The cleanup does not
rewrite native annotation calls. The GFF-derived TSV retains original coordinates,
feature types, strands and attributes (GFF coordinates are 1-based inclusive).

Use `--keep-intermediates` to retain the full working tree for debugging or
stage-level retry. Final products are still exported into `final/` in that mode.
The default compact mode is intended for finished deliverables: matching
`--resume` calls reuse finalized sample results, including reported failures,
without rebuilding the BOLD index or rerunning assembly. Removed stage work cannot
be resumed; retry those analyses in a new output directory. Missing/size-changed
retained deliverables cause reuse to fail instead of silently reporting success.

### Shrink an existing completed run without rerunning analyses

After updating the Tortuga files, run this on TheBelly for the output directory
shown in your earlier command:

```bash
bash /home/ark/MAB/bin/Tortuga/compact_outputs.sh \
  /home/ark/MAB/NOAA/CA_intertidal/both_results_v205
```

An optional preview uses the same command with `--dry-run`. Specify the actual
completed output directory if your run used a different name. This also supports
older Tortuga output versions; no new pipeline run is needed. The cleanup refuses
a directory still locked by an active controller or an unrecognized directory.
For new runs no cleanup flag is needed: compaction is the default.

## Configurable stopping criteria (2.0.6)

The default remains circularity plus all reference CDS/rRNA genes and the
configured genome length range. To stop sooner, choose the topology and a rule:

| Option | Meaning |
| --- | --- |
| `--accept-topology circular` | Require established circularity (default). |
| `--accept-topology any` | Allow an unclosed/linear contig as well as a circular candidate. This does not prove a biologically linear genome. |
| `--accept-rule all-genes` | All reference CDS/rRNA genes plus the length threshold (default). |
| `--accept-rule genes` | At least `--min-detected-genes N` on one contig. |
| `--accept-rule length` | At least `--min-accepted-bp N` on one contig. |
| `--accept-rule either` | Meet the gene-count OR length threshold. Supply both thresholds. |
| `--accept-rule both` | Meet the gene-count AND length threshold. Supply both thresholds. |
| `--skip-wgs-if-accepted` | Skip WGS if read recovery already yielded a passing candidate. ITSME still runs in combined mode. |

Thresholds are batch-wide. Gene count means unique reference protein-coding
CDS plus rRNA genes passing the existing identity and query-coverage cutoffs;
it excludes tRNAs. Counts and lengths are never pooled across contigs or methods.
Every rule still requires target-confirmed COI/BOLD evidence, successful gene
screening, zero ambiguous bases and the `--min-candidate-bp` / `--max-genome-bp`
safety bounds (500 / 30,000 bp by default). Without `--accept-topology any`,
every rule also requires circularity.

`--min-accepted-bp` defaults to `--min-genome-bp` (10,000 bp). It changes the
acceptance length threshold independently of reference-database length filtering.
The `genes` rule uses the gene threshold and common safety bounds; to require
both a gene count and a minimum length, explicitly select `both`.

Examples of user-chosen thresholds, to append to an existing command:

```bash
# Accept a circular or unclosed candidate with at least 12 CDS/rRNA genes.
--accept-topology any --accept-rule genes --min-detected-genes 12

# Accept a circular or unclosed candidate of at least 14,000 bp.
--accept-topology any --accept-rule length --min-accepted-bp 14000

# Require both thresholds; avoid WGS once a candidate passes.
--accept-topology any --accept-rule both --min-detected-genes 12 --min-accepted-bp 14000 --skip-wgs-if-accepted
```

The numbers above illustrate configuration; they are not taxon-specific biological
recommendations. Larger genomes may need a larger `--max-genome-bp` as well.

A passing candidate stops the remaining mitochondrial recovery methods and
receives MITOS2 annotation. Separate WGS still runs unless `--no-assembly` or
`--skip-wgs-if-accepted` applies. If no candidate passes before WGS, that WGS
stage and assembly-based recovery remain available. Target identity and the
actual missing-gene/circularity evidence remain in the reports.

Linear/unclosed sequences are screened without duplicating the sequence across
an artificial end-to-start join. Terminal repeats are preserved unless
circularity is established. MITOS2 receives `--linear` for unclosed candidates,
as documented by the MITOS maintainers:
https://gitlab.com/Bernt/MITOS/-/blob/mitos2/README.md
This is single-contig linear screening/annotation, not a telomere or terminal
structure validation workflow, and does not assemble multipartite genomes.

`accepted_candidate.fasta` contains the candidate that met the selected policy.
Relaxed acceptances use `accepted_by_policy_annotated` or
`accepted_by_policy_annotation_failed`; missing genes remain missing and
circularity is never invented. `strict_complete_candidate` and
`complete_by_reference` report the stricter screens separately. The legacy
`circular_target_candidate.fasta` is written only when the strict circular,
all-reference-gene and genome-length checks also pass. Annotation failures
return a nonzero batch status under every policy.

`acceptance_policy.json`, `candidates.tsv`, per-candidate `assessment.json` and
sample/cohort summaries record the policy, detected genes, thresholds, failed
checks and decision reason. Console messages now print the candidate's length,
gene count, topology, target status and stopping decision.

New analysis settings require a fresh output root: use `both_results_v207`, for example.
They do not alter an already-running v2.0.5 process. Changed settings/version
cannot be resumed into old checkpoints. Existing results are retained.

## Automatic GenBank references (2.0.5)

Only `sample` and `target_taxon` are required in mitochondrial sample sheets.
Retain the full taxonomic hierarchy whenever available: `kingdom`, `phylum`,
`class`, `order`, `family`, `genus`, `species`. These columns now guide selection;
they are not merely display metadata. Original IDs, authorities, and other input
columns are preserved in `resolved_samples.tsv`. AphiaIDs are retained as AphiaIDs,
not treated as NCBI taxonomic identifiers. The included CA intertidal example now
contains all twelve real targets and the hierarchy provided in taxaIDs.xlsx.
The existing template's `reference_gb`, `mt_code`, and `seed_fasta` cells may be
left blank. An explicit reference path or genetic code still overrides selection.

The controller scans `/home/ark/MAB/NOAA/mitos_genbank` by default, configurable with
`--genbank-dir`. It reads GenBank organism/lineage annotations, not filenames or
sequence similarity. It first selects records matching the exact organism or
lineage name (case-insensitive, with whitespace normalized). If no suitable exact
match exists, selection tries the supplied genus, then the supplied family.
This is the default `--reference-match hierarchy` behavior. Missing ranks are
skipped; they are never invented. `--reference-match genus` restores the narrower
2.0.3 behavior, and `--reference-match exact` disables fallback. A genus/family
target selects within that explicitly named lineage. By default, no search
extends beyond the supplied family. `--reference-min-rank order`, `class`, or
`phylum` explicitly permits broader searches through the corresponding supplied
ranks. A supplied phylum must agree with the reference's annotated lineage.
If no reference matches within the permitted scope, that sample is reported as
`reference_unavailable`; other samples continue and the batch exits nonzero.
It does not broaden `target_taxon` to get a seed. The nearest matching rank takes
precedence over annotation-quality ties at more distant ranks. Exact-name matching
still does not resolve nomenclatural synonyms automatically, but family metadata
can provide a same-family seed despite a renamed species or genus.

Automatic candidates must identify a mitochondrial organelle and animal lineage,
pass the same CDS/rRNA validation used by recovery, and meet the configured genome
length and any supplied genetic-code filters. Partial/incomplete descriptions are
excluded. Before/after annotation boundaries are accepted as homology queries
with a warning when DNA is unambiguous and uncertain CDS features have annotated
protein translations. Common rRNA product names are recognized without `/gene`.
Uncertain boundaries cannot establish full-length gene recovery. Among matches,
complete-genome descriptions, fewer boundary warnings, larger annotated CDS/rRNA
gene sets and circular topology are preferred; accession/path settle equal ties.
These rules select a usable taxonomic reference, not a proven closest relative
or independently curated perfect genome. This workflow still uses `animal_mt`.

A blank `mt_code` is inferred only if all eligible references at the selected
match tier agree on one explicitly annotated CDS translation table. Missing or
conflicting codes require a sample-sheet override. A manually specified GenBank
reference can supply a unique code from its CDS annotations as well.

The first scan of a large collection takes time; unchanged records are reused
from `tortuga_reference_index.sqlite` in that directory. Version 2.0.5 automatically
refreshes older cached validations once; do not delete the index manually. `--reference-index` can
place this cache elsewhere. New/changed files are re-inspected and removed files
are pruned. `reference_selection.tsv` records the selected accession, organism,
match level, actual matched rank/name, code source and path. `resolved_samples.tsv`
records usable per-sample inputs with all original taxonomy and metadata columns.
Per-sample summaries and run checkpoints retain the input hierarchy as well.
Neither report changes the expected target used for the BOLD gate. For example,
a Pinnotheridae reference may seed recovery for Rathbunixa occidentalis, while
the BOLD gate still requires Rathbunixa occidentalis, not just Pinnotheridae.

For the user's current Tortuga installation and prepared databases:

```bash
bash /home/ark/MAB/bin/Tortuga/pipeline_controller.sh \
    --pipeline both \
    -i /home/ark/MAB/NOAA/CA_intertidal/raw_reads \
    --samples /home/ark/MAB/NOAA/CA_intertidal/CA_intertidal_samples_full_taxonomy.tsv \
    --genbank-dir /home/ark/MAB/NOAA/mitos_genbank \
    --bold-fasta /home/ark/MAB/NOAA/BOLD/prepared_fixed/coi.fasta \
    --bold-metadata /home/ark/MAB/NOAA/BOLD/prepared_fixed/metadata.tsv \
    --mitos-refdir /home/ark/MAB/NOAA/mitos_data \
    --mitos-refseqver refseq89m \
    --itsme-db /home/ark/databases/itsme_db \
    -o /home/ark/MAB/NOAA/CA_intertidal/both_results_v207 \
    -t 32 -m 200
```

Use `--pipeline mitolooper` to omit ITSME. New controller versions/options require
a fresh output root; existing older-version checkpoint directories are retained.

If a sample has no usable reference, `both` still runs its ITSME stage under
`OUT/SAMPLE/final/` with `itsme_` prefixes; it skips that
sample's mitochondrial recovery and WGS. `reference_selection.tsv` records the
failure and any boundary warnings, and `reference_failures.tsv` lists unresolved
samples. The cohort and invocation summaries retain `reference_unavailable` and
ITSME success/failure separately. In compact mode `--resume` reuses finalized
results; use `--keep-intermediates` for stage retries and a new output root after
a version change.

For a run allowing Rhizophysa/unclosed candidates, use `--accept-topology any`
and an appropriate rule above. Acceptance does not establish biological linearity
or full terminal structures.

| Consumer | Use of selected GenBank reference |
| --- | --- |
| MitoFinder | Passes the annotated GenBank file directly with `-r`. |
| NOVOPlasty and read GetOrganelle | Receives a COX1 nucleotide seed extracted by this wrapper. |
| Assembly GetOrganelle | Receives extracted nucleotide gene labels. |
| mtGrasp | Receives the extracted mitochondrial nucleotide genome (or a confirmed rescue partial). |
| Candidate evaluation | Uses extracted CDS proteins and rRNAs for the completeness screen and COX1 localization. |
| MITOS2 | Uses its separate `refseq89m` release; it does not receive the selected GenBank file. |
| ITSME and WGS assembly | Do not receive the selected mitochondrial GenBank file. |

Run `pipeline_controller.sh` from `base`; environment switching is automatic.
The required selector is `--pipeline itsme`, `--pipeline mitolooper`, or
`--pipeline both`. An input directory can be supplied with `-i` / `--input-dir`
or as the positional argument. Samples run sequentially.

| Selector | What runs | Reference requirements |
| --- | --- | --- |
| `itsme` | ITSME nuclear 18S/ITS/28S recovery only | Existing ITSME rRNA database root; no BOLD, mitochondrial references, MITOS models, or sample sheet required. |
| `mitolooper` | Initial mito recovery, MITOS2 for passing candidates, optional WGS assembly and mito fallback | Same mitochondrial sample sheet, BOLD and MITOS2 references as MitoLooper. ITSME is not called. |
| `both` | Initial mito recovery/annotation, then ITSME once, then WGS and post-WGS mito fallback if enabled/needed | All mitochondrial references plus the ITSME rRNA root. |

`both` uses the integrated workflow rather than launching two independent full
pipelines. ITSME therefore runs once per sample in the intended position, before
WGS. Its failure is reported while mitochondrial fallback can continue. All mode
selectors preserve per-stage environment choices and notes about other relevant
Conda installations.

## ITSME only: your current reads

```bash
bash ~/bin/MitoGenomeScanner/pipeline_controller.sh \
    --pipeline itsme \
    -i /home/ark/MAB/NOAA/CA_intertidal/raw_reads \
    -o /home/ark/MAB/NOAA/CA_intertidal/itsme_results \
    -t 60 -m 200
```

With this selector, the default rRNA database is
`/home/ark/databases/itsme_db`, as listed in your server output. Override it with
`--itsme-db /actual/path` or the `ITSME_DB` environment variable. The default
script is `~/bin/ITSME/itsme.v2.sh`; override with `--itsme-script` if needed.
No mitochondrial metadata columns are read. If `--samples` is provided in this
mode, only its `sample` column is used to select names from the FASTQ directory.

`--sensitivity 2` is the default ITSME preset; 1 and 3 are also supported. This
setting affects only ITSME, including when `--pipeline both` is selected.

## MitoLooper only or both

Use your actual reference paths in place of the capitalized placeholders:

```bash
bash ~/bin/MitoGenomeScanner/pipeline_controller.sh \
    --pipeline both \
    -i /home/ark/MAB/NOAA/CA_intertidal/raw_reads \
    --samples ~/bin/MitoGenomeScanner/samples.tsv \
    --bold-fasta /PATH/TO/bold/coi.fasta \
    --bold-metadata /PATH/TO/bold/metadata.tsv \
    --mitos-refdir /PATH/TO/mitos_data \
    --mitos-refseqver refseq89m \
    --itsme-db /home/ark/databases/itsme_db \
    -o /home/ark/MAB/NOAA/CA_intertidal/both_results \
    -t 60 -m 200
```

For mitochondrial recovery only, change `--pipeline both` to
`--pipeline mitolooper` and select an appropriate output directory. All other
mitochondrial options from `mito_looper.sh` remain available. The controller does
not require a generated configuration file or manual environment activation.

## Common controls and outputs

- `--only-sample W30877` runs one sample for a trial.
- `--resume` reuses finalized compact samples with matching inputs/settings. With
  `--keep-intermediates`, failed stages can be retried in a new attempt directory;
  directory existence alone never proves success.
- `--no-assembly` disables the separate WGS and post-WGS mito recovery stages.
  In `both` mode ITSME still runs and does its own targeted assembly.
- `--itsme-timeout 180` limits each ITSME call to 180 minutes by default.
- `--dry-run` prints planned commands in a new/empty output directory.

When `-o` is omitted, the default directory is `itsme_results`, `mito_results`, or
`both_results` according to the selector. Use a separate output directory when
switching modes; checkpoints from another mode/version are not silently reused.

ITSME-only and combined results are exported directly into `OUT/SAMPLE/final/`
with `itsme_` prefixes, including native master/partial summaries and FASTAs.
`itsme_output` in the sample/cohort summaries points there. Logs, manifests and
cleanup records retain the original working paths for provenance. Standalone
batches still continue after individual failures and return nonzero on failure.
With `--keep-intermediates`, the native deep output tree is also retained.

The old `mito_looper.sh` and legacy `--branch mito/full` flags remain available.
For the new controller, use the single `--pipeline` selector; incompatible legacy
branch combinations are rejected. No third-party assembly or ITSME installation
is bundled. The tests use synthetic tool doubles rather than real project analyses.

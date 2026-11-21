# One controller for ITSME, MitoLooper, or both

Version 2.1.0 standardizes sample ITSME/ and MitoLooper/ outputs and defaults
to hierarchical local reference selection. See
[REFERENCE_SELECTION.md](REFERENCE_SELECTION.md) for the current
algorithm, fallback safeguards and migration behavior.

## Annotation and BOLD validation (2.3.0)

Default `--annotation-policy best` annotates the best retained mitochondrial
candidate regardless of circularity, reference-gene recovery or target status.
MITOS2 uses `--linear` for unclosed sequences. Annotation is labeled exploratory
unless the candidate passed the original acceptance policy; the policy itself
is unchanged. A known genetic code remains required. `accepted` and `none` are
explicit annotation-policy alternatives.

For an existing inactive run, after updating Tortuga:

```bash
bash /home/ark/MAB/bin/Tortuga/annotate_existing.sh \
  /home/ark/MAB/NOAA/CA_intertidal/both_results_v210_hierarchy \
  --recheck-bold
```

Saved per-sample options supply the database paths, genetic codes and environment
commands. No assemblies or ITSME analyses are repeated. Omit `--recheck-bold` to
annotate without repeating COI validation. Existing unchanged annotations are
cached. The helper exports annotations, preserves old summaries/stages for audit,
updates cohort tables/dashboard, and removes only its own successful attempt
scratch and generated BOLD index. See README.md for overrides and retry controls.

BOLD starts at `--bold-max-hits 500`. A saturated near-best score boundary expands
the search up to `--bold-search-max-hits 50000`; a full result page with a clear
score gap no longer causes unconditional rejection. Incomplete near-best evidence
at the hard cap and true taxonomic conflicts stay unresolved.

ITSME badges now use its own execution outcome, not the overall mitochondrial
status. `itsme_run_status`, `itsme_partial_locus_count` and
`itsme_diagnostic_reason` distinguish execution failures from successful runs
that recovered no complete ITS-containing loci.

## Portable offline dashboard (2.2.0)

`dashboard.html` now accompanies the cohort reports. It combines a quiet sample
list with mitogenome gene maps/metrics, ITS locus maps and consensus taxonomy,
retained GFA1/FASTG graph previews, annotation tables, stage outcomes and relative
file downloads. It works by opening the file directly; move the whole result
folder to preserve downloads. New completed runs, report refresh and compaction
build it automatically.

For existing output, after merging the PR and pulling Tortuga:

```bash
bash /home/ark/MAB/bin/Tortuga/build_dashboard.sh \
  /home/ark/MAB/NOAA/CA_intertidal/both_results_v210_hierarchy
```

This only writes dashboard.html and refuses active runs. See DASHBOARD.md for
preview bounds, supported graph formats, evidence interpretation and validation.

## ITS lengths and taxonomy in cohort reports (2.1.1)

`cohort_summary.csv` includes `itsme_complete_locus_count`, matching JSON arrays
of locus IDs, full locus lengths, ITS1–5.8S–ITS2 lengths, consensus taxa and ranks,
and `itsme_reporting_status`. `cohort_its_loci.csv` gives one row per reported
complete locus, with individual region lengths and the full consensus lineage.
Taxonomy comes from ITSME evidence; expected target names are not predictions.

Existing completed results need only a reporting refresh, after updating Tortuga:

```bash
git -C /home/ark/MAB/bin/Tortuga pull --ff-only

bash /home/ark/MAB/bin/Tortuga/refresh_reports.sh \
  /home/ark/MAB/NOAA/CA_intertidal/both_results_v210_hierarchy
```

This does not run analysis or cleanup. It refuses active runs and skips unfinished
or stale sample summaries. New combined and ITSME-only runs include these reports
automatically, with or without compact cleanup.

## Compact outputs by default (2.1.0)

Completed samples export client deliverables into `OUT/SAMPLE/ITSME/` and
`OUT/SAMPLE/MitoLooper/` and
remove generated stage work files automatically. This includes samples with
failed mitochondrial recovery/WGS: their actual status and logs are retained.
Cleanup does not rerun tools, change acceptance criteria or turn failures into
successful biological recovery.

| Retained location | Contents |
| --- | --- |
| `SAMPLE/ITSME/*.fasta`, `.csv`, `.tsv` | Native ITSME final sequences, master/partial summaries and annotation/support tables, in the same sample tree as mito outputs. |
| `SAMPLE/MitoLooper/mitos_*` | Native MITOS2 sequences, GFF/GenBank/BED products when present, plus `mitos_result_features.tsv` generated from the GFF. |
| `SAMPLE/MitoLooper/*_assembly_graph*.fastg`, `.gfa`; `SAMPLE/ITSME/*graph*` | Final available assembly graphs from ITSME, WGS and recovery tools. |
| `SAMPLE/MitoLooper/wgs_contigs.fasta`, `wgs_scaffolds.fasta` | Final WGS sequence products when produced; summary flags still identify a failed/timed-out WGS run. |
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
stage-level retry. Final products are still exported into the branch directories in that mode.
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
output directory if your run used a different name. This also supports older
Tortuga output versions; no new pipeline run is needed. A missing batch-level
`last_invocation.json` no longer blocks cleanup of individual finished samples.
Only sample directories with a fresh terminal `final/summary.json` are compacted;
unfinished, unstarted, planned and stale-summary samples are skipped. The shared
BOLD index is retained unless full batch completion is established.
The cleanup still refuses the entire operation while a controller holds that
output directory's lock. It never cleans a different run directory. A missing
path and an unrecognized folder now produce distinct messages.
For new runs no cleanup flag is needed: compaction is the default.

## Circle-first recovery (default, 2.5.0)

Read methods run sequentially until the first circular candidate is recovered.
That circle stops all later mitochondrial assemblers, including WGS. In combined
mode ITSME still runs independently. Taxonomy and reference-gene checks are
recorded but do not prevent circular recovery from stopping the method loop.

If the read methods yield no circle and `--no-assembly` is absent, metaSPAdes runs
without a deadline by default, followed by GetOrganelle graph recovery and then
MitoFinder contig recovery only as needed. Stop at the first circle in this phase
as well. `--no-assemble` is an alias of `--no-assembly`.

If no circle is found, retain the single candidate with the most expected
CDS/rRNA genes, then greatest length. Without gene-screen evidence, this selects
the longest candidate. `--linear-ranking length` reverses those priorities.
Ambiguous-base count breaks a remaining tie. Taxonomic confirmation does not
override this fallback ordering. Annotation remains enabled for the selected
linear candidate by default.

`--assembly-timeout 0` means unlimited; positive values impose an explicit
deadline in minutes. Other method/validation limits retain their own settings.
Use a new output directory for this workflow; rebuilding a historical dashboard
does not rerun methods or change which sequence that run selected.

## Optional legacy acceptance stopping policy

Select `--recovery-policy acceptance` to restore the earlier gated stopping
workflow. In that mode, the default requires circularity plus all reference
CDS/rRNA genes and the configured genome length range. To stop sooner, choose
the topology and a rule:

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
--recovery-policy acceptance --accept-topology any --accept-rule genes --min-detected-genes 12

# Accept a circular or unclosed candidate of at least 14,000 bp.
--recovery-policy acceptance --accept-topology any --accept-rule length --min-accepted-bp 14000

# Require both thresholds; avoid WGS once a candidate passes.
--recovery-policy acceptance --accept-topology any --accept-rule both --min-detected-genes 12 --min-accepted-bp 14000 --skip-wgs-if-accepted
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

New analysis settings require a fresh output root: use `both_results_v210`, for example.
They do not alter an already-running v2.0.5 process. Changed settings/version
cannot be resumed into old checkpoints. Existing results are retained.

## Hierarchical references (2.1.0)

The defaults are `--reference-match hierarchy --reference-min-rank kingdom`.
Try the target, genus, family, order, class, phylum, then kingdom; skip missing
ranks and stop at the first eligible match. An explicit `reference_gb` overrides
this search. BLAST reference selection and its read-sampling flags are removed.
Annotation, genome length and genetic-code checks still apply. Conflicting codes
at the closest rank require an explicit `mt_code`; the pipeline then uses
reference-independent recovery instead of guessing. See REFERENCE_SELECTION.md.

Merge the PR, then update the existing Belly checkout:

```bash
git -C /home/ark/MAB/bin/Tortuga pull --ff-only
```

Run both branches from scratch into a new output directory:

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
  --sensitivity 2 \
  --validation-timeout 120 \
  -t 32 -m 200 \
  -o /home/ark/MAB/NOAA/CA_intertidal/both_results_v210_hierarchy
```

Run `pipeline_controller.sh` from `base`; environment switching is automatic.
The required selector is `--pipeline itsme`, `--pipeline mitolooper`, or
`--pipeline both`. An input directory can be supplied with `-i` / `--input-dir`
or as the positional argument. Samples run sequentially.

| Selector | What runs | Reference requirements |
| --- | --- | --- |
| `itsme` | ITSME nuclear 18S/ITS/28S recovery only | Existing ITSME rRNA database root; no BOLD, mitochondrial references, MITOS models, or sample sheet required. |
| `mitolooper` | Initial mito recovery, MITOS2 for the best retained candidate, optional WGS assembly and mito fallback | Same mitochondrial sample sheet, BOLD and MITOS2 references as MitoLooper. ITSME is not called. |
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
  --sensitivity 2 \
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

ITSME-only and combined results are exported directly into `OUT/SAMPLE/ITSME/`
with native filenames, including native master/partial summaries and FASTAs.
`itsme_output` in the sample/cohort summaries points there. Logs, manifests and
cleanup records retain the original working paths for provenance. Standalone
batches still continue after individual failures and return nonzero on failure.
With `--keep-intermediates`, the native deep output tree is also retained.

The old `mito_looper.sh` and legacy `--branch mito/full` flags remain available.
For the new controller, use the single `--pipeline` selector; incompatible legacy
branch combinations are rejected. No third-party assembly or ITSME installation
is bundled. The tests use synthetic tool doubles rather than real project analyses.


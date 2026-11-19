# Tortuga 2.3.1

Every sample has **ITSME/** and **MitoLooper/** results with explicit branch
statuses, including disabled or failed branches. Native ITSME FASTAs, tables and
graphs are in ITSME/; mitochondrial sequences, annotations and graphs are in
MitoLooper/. Shared summaries/QC/logs remain alongside them. Compact cleanup is
still the default; --keep-intermediates preserves work for debugging.

**Automatic reference choice follows annotated taxonomy**: target, genus, family,
order, class, phylum, then kingdom, stopping at the first eligible match.
Only supplied ranks are used (genus can also be taken from a binomial target).
BLAST reference selection has been removed. Explicit reference overrides remain
available. If no eligible reference exists, GetOrganelle uses its installed
animal_mt database; ITSME and WGS still run normally.

See [REFERENCE_SELECTION.md](REFERENCE_SELECTION.md) for rank matching,
genetic-code safeguards, reports, reference-independent recovery and migration.
Reference choice does not change the intended target for COI/BOLD validation.
The default stopping policy still requires circularity and all reference genes;
relaxed gene/length rules remain explicit options documented in CONTROLLER.md.

**Unified controller:** `pipeline_controller.sh --pipeline itsme|mitolooper|both`
selects standalone ITSME, mitochondrial recovery, or the integrated workflow.
Use one of the three literal values (not the `|` notation). See `CONTROLLER.md`
for commands. ITSME-only mode needs no mitochondrial databases or metadata.

**TheBelly users:** run `bash mito_looper.sh ...` directly, including from `base`.
The launcher enters the existing `itsme` environment for Python and runs each
software stage in its selected Conda environment automatically. No configuration
generation or manual activation is required. Read `BELLY_SETUP.md` for the command
and reference arguments. The generic installation instructions below are for other hosts.

The first existing executable in the environment priority lists wins. Those lists
and comments about relevant alternatives are in `configure_belly.py` beside
`SPECS`; the launcher also notes the main alternatives. GetOrganelle tries
`getorganelle`, then `mito-getorganelle`; mtGrasp tries `mito-mtgrasp`, then `mtgrasp`.
Selections and alternatives are printed and written to `environment_selection.tsv`
and `.json` in the output directory. `--config` remains an optional explicit override.

Sequential mitochondrial recovery updated for the SCCWRP proposal. The uploaded
shell pipeline has been retained under `legacy/`; the new controller is
`mito_pipeline.py`, launched through either `mito_looper.sh` or `mito_pipeline.sh`.
There is no GNU parallel dispatch: both samples and methods run sequentially.

## What runs

| Stage | Behavior |
| --- | --- |
| Read preparation | Paired fastp QC by default; `--no-qc` uses the input reads directly. Mates are retained as pairs; no overlap merging. |
| 1. Mitochondria from reads | NOVOPlasty, GetOrganelle, MitoFinder, then mtGrasp. Validate each output before trying the next method. Stop remaining read methods when a candidate meets the configured policy and target gate; the default requires circularity and all reference genes. Annotate that candidate with MITOS2 immediately. |
| 2. ITS/rDNA | `--branch full` calls the separately installed ITSME v1.9 shell pipeline on original reads. Default `--branch mito` omits nuclear rDNA. |
| 3. WGS assembly | metaSPAdes, with `--only-assembler`, a memory setting, wall-clock timeout, and optional percentage of reads. MEGAHIT is an alternative. This stage still runs after successful read recovery unless `--no-assembly` or `--skip-wgs-if-accepted` applies. |
| 4. Mitochondria from assembly | If stage 1 failed to produce a passing candidate, GetOrganelle processes the metaSPAdes FASTG graph; MitoFinder then processes contigs if needed. New passing candidates go to MITOS2; otherwise annotate the best retained candidate after recovery finishes. |
| Reporting | Per-sample candidate/evidence tables, all method logs and statuses, best contiguous sequence, and a cohort CSV. |

The proposal's optional ABySS/SOAPdenovo2 benchmarking is not implemented; the
supported whole-genome assemblers are metaSPAdes and MEGAHIT. The latter does not
produce a supported assembly graph in this controller, so it goes straight to
MitoFinder for stage 4.

## Install

The controller uses Python 3.10+ and Biopython. A separate orchestration environment
avoids changing your working ITSME environment:

```bash
cd ~/bin/MitoGenomeScanner
mamba env create -f environment.yml
conda activate mito_workflow
```

The included environment installs the controller dependencies, BLAST+, fastp,
minimap2, SPAdes and MEGAHIT. Install the remaining assemblers separately, or
point the JSON command configuration at existing installations:

```bash
# Skip these if the environments already exist.
mamba create -n getorganelle -c conda-forge -c bioconda getorganelle
conda run -n getorganelle get_organelle_config.py --add-organelle animal_mt
mamba create -n mtgrasp -c conda-forge -c bioconda mtgrasp
mamba create -n mitos2 -c conda-forge -c bioconda 'mitos>=2'
```

NOVOPlasty is a Perl script: obtain the official distribution and set its actual
path/version in `config.json`. For MitoFinder, its official Singularity image
avoids its older Python dependencies. The example configuration uses
`singularity run`; change it to your native executable or Apptainer prefix if
appropriate. Bind every filesystem containing reads, references and outputs into
the container. `/home:/home` is only a starting example.

```bash
cp config.example.json config.json
# Edit NOVOPlasty/MitoFinder paths and any environment names in config.json.
python check_tools.py --config config.json
```

All command prefixes are JSON arrays, with `~` and environment variables expanded.
No shell evaluation or shell command interpolation is used. Run the controller
from the `mito_workflow` environment: unlisted tools are resolved from its PATH.
The help probes also save tool help/version text in `tool_checks.json`. An
installed tool returning nonzero for `-h` needs inspection, even if usable.
Probe success checks launchability, not database validity or biological performance.

## MITOS output-directory fix (2.3.1)

The MITOS launcher now creates the `-o` directory before starting runmitos.py.
This fixes `runmitos.py: error: no such directory .../output` in both new runs and
saved-candidate annotation. For an existing run where BOLD completed but MITOS
failed, merge the fix, update Tortuga and retry **annotation only**:

```bash
git -C /home/ark/MAB/bin/Tortuga pull --ff-only

bash /home/ark/MAB/bin/Tortuga/annotate_existing.sh \
  /home/ark/MAB/NOAA/CA_intertidal/both_results_v210_hierarchy
```

No `--recheck-bold` or `--force` is needed: failed annotations are retried in a
fresh attempt directory. Assemblies, ITSME and previous BOLD searches are not
repeated. Failed attempts are retained; successful unchanged annotations remain
cached. Cohort reports and dashboard are refreshed after annotation.

## Annotate retained candidates by default (2.3.0)

MITOS2 now annotates the best retained mitochondrial candidate even when it is
unclosed, contains ambiguous bases, has missing reference genes or has unresolved
COI identity. Unclosed candidates receive `--linear`; known circular candidates
keep circular mode. Annotation requires a known mitochondrial genetic code.
`--annotation-policy best` is the default; `accepted` restores accepted-only
annotation, and `none` explicitly disables it. This does not relax biological
acceptance or convert an annotated partial/unresolved sequence into a complete
confirmed mitogenome. `annotation_status`, `annotation_scope` and the annotated
candidate ID are recorded separately. `annotation_input.fasta` and MITOS products
are retained in MitoLooper/.

For your existing Belly results, after merging the PR:

```bash
git -C /home/ark/MAB/bin/Tortuga pull --ff-only

bash /home/ark/MAB/bin/Tortuga/annotate_existing.sh \
  /home/ark/MAB/NOAA/CA_intertidal/both_results_v210_hierarchy \
  --recheck-bold
```

This uses the saved candidate FASTAs, genetic codes, Conda commands and database
paths. It does not run an assembler, QC, ITSME, or reference selection. It repeats
BOLD validation using the saved COI database inputs and then performs MITOS
annotation. It reuses saved COX1 protein-hit coordinates when available, including
circular-origin spans; otherwise it queries the whole retained candidate against
the COI-only database. Omit `--recheck-bold` for annotation only. Successful unchanged
annotations are cached; `--force` repeats them. Unknown codes require an explicit
`--mt-code` or per-sample saved code. Changed database locations can be supplied
with `--mitos-refdir`, `--mitos-refseqver`, `--bold-fasta`, `--bold-metadata`;
`--only-sample` restricts processing. The active-run lock is respected.

With `--recheck-bold`, the existing selected candidate is reassessed under the
original run's acceptance policy; no thresholds are loosened. Previously published
acceptance files are archived if the updated evidence no longer passes. Original
sample summaries and stage records are retained under postprocess_history/.
MITOS outputs are exported before deleting only that attempt's scratch, and the
helper's generated BOLD index is removed after processing. `--keep-intermediates`
retains helper work for debugging. Reports and dashboard are refreshed. A BOLD
recheck failure is reported with a nonzero exit but does not block annotation.

ITSME execution and recovery are also separate: `itsme_run_status` reports whether
the program completed, failed or timed out, while complete/partial locus counts
report biological recovery. A successful run with zero complete loci is not
labeled an ITSME execution failure. `itsme_diagnostic_reason` records available
stage failure reasons or points to partials and the native run summary.

## Portable HTML dashboard (2.2.0)

Open `dashboard.html` in the results directory for a sample-centered view of
mitogenome sequences and annotations, ITS loci and consensus taxonomy, retained
assembly graphs, QC and stage outcomes. It works offline and travels with the
whole results folder. Gene/region maps, graph previews and display data are
embedded; original downloads use relative links. New completed runs and report
refreshes generate it automatically. No additional production libraries are needed.

Build it for an existing run without rerunning analysis or cleanup:

```bash
bash /home/ark/MAB/bin/Tortuga/build_dashboard.sh \
  /home/ark/MAB/NOAA/CA_intertidal/both_results_v210_hierarchy
```

See [DASHBOARD.md](DASHBOARD.md) for deployment, graph limits and interpretation.
Missing products stay explicitly unavailable; the display does not change the
underlying taxonomy or completeness calls.

## ITS results in the cohort summary (2.1.1)

`cohort_summary.csv` now includes `itsme_complete_locus_count`,
`itsme_locus_ids`, `itsme_its_lengths_bp`, `itsme_locus_lengths_bp`,
`itsme_taxa`, `itsme_taxonomy_ranks` and `itsme_reporting_status`.
The list-valued cells are JSON arrays in matching locus order. ITS length means
ITS1–5.8S–ITS2, distinct from the full rDNA locus length. It is computed from
inclusive ITSx coordinates, with matching extracted ITS FASTA as a fallback.
Unknown lengths or assignments remain unknown; a missing report is not zero loci.
Taxa are the deepest ITSME consensus assignments, not supplied targets or seed
reference names. Rank labels come from explicit consensus rank fields; older
reports with only lineage text retain the taxon with an unknown rank.

`cohort_its_loci.csv` is the easier table for reviewing multiple loci: one row per
reported complete locus, with full locus length, ITS cassette length, separate
ITS1/5.8S/ITS2 lengths, consensus taxon/rank/lineage and taxonomy status. Partials
and chimeras are not counted as complete loci; they remain in native ITSME outputs.
New runs create both tables automatically in combined and standalone ITSME mode.

To update existing completed results after merging the PR:

```bash
git -C /home/ark/MAB/bin/Tortuga pull --ff-only

bash /home/ark/MAB/bin/Tortuga/refresh_reports.sh \
  /home/ark/MAB/NOAA/CA_intertidal/both_results_v210_hierarchy
```

This command only refreshes reporting files. It does not launch analyses, move
assemblies or perform cleanup. It refuses an active controller and skips stale
or unfinished sample summaries. It can read current compact, legacy prefixed
final, and retained native ITSME reports. Missing native evidence requires no
invented taxonomic assignment and is reported explicitly.

## References and sample identities

Reference-dependent tools need a suitable mitochondrial reference; GetOrganelle
can bootstrap from its installed database without a custom reference. Target
validation still needs the intended sample identity, which cannot be inferred
from a sample ID such as W30877.
Create your sample sheet once:

```bash
python mito_pipeline.py ~/NOAA/CA_intertidal/raw_reads \
    --make-sample-sheet samples.tsv
```

Fill the tab-separated columns:

| Column | Meaning |
| --- | --- |
| `sample` | FASTQ stem, e.g. `W30877`. |
| `target_taxon` | Expected species or other chosen rank, spelled as in BOLD metadata. Confirmation is only at this requested taxon, not automatically at species level. |
| `mt_code` | Optional explicit NCBI mitochondrial translation table. Inferred only when eligible references at the closest matching rank agree on a code. Without a reference it remains unknown and annotation requires an explicit override. |
| `reference_gb` | Optional path override for ONE annotated mitochondrial genome. Blank/missing uses local hierarchical taxonomy matching by default. Relative override paths are resolved against the sheet directory. |
| `seed_fasta` | Optional mitochondrial/COI nucleotide seed. If blank, the controller derives COX1 from `reference_gb`. |
| `kingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species` | Recommended full input hierarchy. Rank column names may be lower or title case. `species` may contain the epithet or full binomial. Missing ranks remain blank. |

The reference must have uniquely named CDS/rRNA genes, including COX1,
at least three CDS genes and at least one rRNA. The expected CDS/rRNA gene set is
derived from that reference; it does not assume every animal has exactly 13 CDS
genes. An unrepresentative or incomplete reference weakens this screen. Repeated
gene names, ambiguous reference gene DNA or multiple GenBank records are rejected
for explicit curation. Before/after (`<`/`>`) annotation boundaries are allowed as
homology queries when the extracted DNA is unambiguous; uncertain CDS boundaries
also require an explicit GenBank protein translation. They are reported in
`references/reference_quality.json` and the final summary, and cannot establish
full-length gene boundaries. Other uncertain coordinate types remain rejected.
rRNA product labels `s-rRNA`/`l-rRNA`, 12S/16S and small/large subunit labels are
recognized when `/gene` is absent. This version does not model duplicated genes
or multipartite mitochondrial genomes. Single-contig linear/unclosed candidates
are annotated by default; accepting them requires `--accept-topology any`. This is not validation
of telomeres or biological terminal structures.

An unavailable initial local reference now triggers custom-reference-free
mitochondrial recovery in the same sample tree. In combined mode ITSME still
runs; WGS remains available. Initial selection diagnostics are retained in
reference_failures.tsv. See REFERENCE_SELECTION.md for bootstrap assessment
when no gene annotation reference is available.

From GenBank the controller derives: a mitochondrial FASTA for mtGrasp, COX1 seed,
CDS translations for completeness screening, and CDS/rRNA nucleotide labels for
GetOrganelle assembly extraction. No rDNA database from ITSME substitutes for
these mitochondrial references.

## BOLD validation

Use a local BOLD COI reference export, with taxonomic metadata and provenance. The
controller does not upload sample data to BOLD and does not call an undocumented
online identification API. This is a BLAST comparison against BOLD-derived records,
not the hosted BOLD identification engine.

Download a suitable public BOLD TSV export or reference snapshot, then normalize
its COI records:

```bash
python prepare_bold.py /path/to/bold_export.tsv.gz \
    --source 'BOLD dataset URL/version; downloaded YYYY-MM-DD' \
    -o /path/to/bold_coi
```

The helper recognizes `nuc`/`nucleotides`/`nucleotide_sequence`/`sequence`,
`processid`/`process_id`, and `marker_code`/`marker`/`markercode`. Explicit column
overrides are available in `--help`. Taxonomy comes from standard kingdom through
subspecies columns. It retains eligible COI/COI-5P records of at least 400 bp and
logs skipped records. It does not download the reference dataset.

Alternatively supply an existing FASTA and a normalized metadata TSV with exactly
matching FASTA identifiers and these fields:

`sequence_id`, `processid`, `taxon`, `lineage`, `source`.

Lineage is semicolon-separated, with full taxon names. `source` should identify
the BOLD dataset/release and retrieval date. `examples/bold_metadata.tsv` is a
format example containing fictional placeholder data, not a reference database.

Include plausible alternative taxa and contaminants. A target-only reference set
cannot exclude other organisms that are absent from that set. References can be
incomplete, mislabeled or unable to separate species; matching COI does not verify
every other base in the assembled contig or exclude NUMTs/mitochondrial chimeras.

COX1 is located by translated comparison to the curated reference. Its nucleotide
sequence is searched against the BOLD COI database. Defaults require >=97%
identity, >=400 aligned bp and >=80% of the shorter marker/barcode sequence.
Every qualifying hit within 2% of the best bit score must agree with the requested
taxon. Missing metadata, conflicting near-best hits, a higher-scoring hit failing
quality gates leaves the candidate unresolved. Reaching the subject limit alone
no longer rejects the candidate: a score boundary below the near-best window
allows the returned near-best evidence to be evaluated. If the window is still
truncated, the search expands tenfold from `--bold-max-hits` (default 500) up to
`--bold-search-max-hits` (default 50,000). A tied/truncated window at the hard cap
remains unresolved. Search limits, boundary scores and attempts are saved in
`bold_evidence.json`; taxonomic conflicts and missing metadata still block confirmation.

These are configurable starting heuristics, not calibrated species boundaries.
Per-candidate raw hits and BOLD process IDs, taxonomy, source and decision reason
are retained. Species-level expectations without an adequate BOLD match remain
unresolved; broaden `target_taxon` only when that is scientifically appropriate.

## Run your mitochondrial branch

After filling `samples.tsv`, setting the command paths, and preparing references:

```bash
conda activate mito_workflow

bash ~/bin/MitoGenomeScanner/mito_looper.sh \
    ~/NOAA/CA_intertidal/raw_reads \
    --samples samples.tsv \
    --config config.json \
    --bold-fasta /path/to/bold_coi/coi.fasta \
    --bold-metadata /path/to/bold_coi/metadata.tsv \
    --mitos-refdir /path/to/MITOS2-reference-root \
    --mitos-refseqver refseq-release-directory \
    --branch mito \
    -o ~/NOAA/CA_intertidal/mito_results \
    -t 60 -m 200
```

Use your actual MITOS2 reference release folder name, not the literal placeholder.
Download/extract its reference data separately using the MITOS2 project's reference
data instructions. The specified release must exist under `--mitos-refdir`.

This command discovers all sample pairs, so an additional `for` loop is unnecessary.
It accepts your `W30877_1.fq.gz` / `W30877_2.fq.gz` names as well as common Illumina
`_R1_001.fastq.gz` / `_R2_001.fastq.gz`, `_R1.fastq.gz` / `_R2.fastq.gz`,
`_R1.fq.gz` / `_R2.fq.gz`, and `_1.fastq.gz` / `_2.fastq.gz` names. Duplicate
sample stems or missing mates fail explicitly. Merge sequencing lanes beforehand.

For only initial mitochondrial recovery, add `--no-assembly`. This disables the
separate WGS and post-WGS recovery stages. Initial recovery methods still assemble
reads internally; MitoFinder may itself invoke MEGAHIT on the supplied reads.

To run the full proposal including ITSME between initial mitochondrial recovery
and WGS, replace `--branch mito` with:

```bash
--branch full --itsme-db "$HOME/../databases/itsme_db"
```

The default ITSME command prefix uses `~/bin/ITSME/itsme.v2.sh` in its `itsme`
environment. ITSME failures are recorded but do not suppress the mitochondrial
fallback stages. `--no-assembly` can be combined with either branch.

## Runtime and restart controls

| Option | Default / effect |
| --- | --- |
| `--method-timeout` | 60 minutes per NOVOPlasty/GetOrganelle/MitoFinder/mtGrasp call. |
| `--assembly-timeout` | 120 minutes per WGS assembler process group. |
| `--annotation-timeout` | 120 minutes for MITOS2. |
| `--validation-timeout` | 20 minutes per BLAST or junction-read verification call. |
| `--itsme-timeout` | 180 minutes for the ITSME process group in full mode. |
| `--assembly-percent` | 100. Change to e.g. `10` for WGS read subsetting only. |
| `--sampling-seed` | 17. Reproducible Bernoulli template sampling; actual percentage may differ slightly. |
| `--assembler` | `metaspades`; `megahit` is supported as an alternative. |
| `--skip-method` | Repeat to disable individual recovery tools explicitly. |
| `--resume` / `--skip-existing` | Reuse finalized compact samples; with `--keep-intermediates`, reuse successful matching stages. |
| `--keep-intermediates` | Preserve full stage work for debugging/retries; default exports essentials and deletes scratch files. |
| `--only-sample` | Run one sample from the sheet, preserving other saved cohort entries. |
| `--dry-run` | Show planned command arrays without launching external tools. |
| `--env-root` | `/home/ark/miniconda3/envs`; select existing environments automatically under this root. |
| `--conda` | Optional Conda executable override; the matching Miniconda executable is preferred. |

Timeouts are per invocation, not a guaranteed total run ceiling. Child processes in
the owned process group receive TERM and then KILL if needed. Native process
groups are supported; remotely managed/container-daemon jobs are not guaranteed
to follow the same cleanup behavior. Input scans, Python validation/subsetting,
and reference preparation are not whole-run time limited. NOVOPlasty has no
general thread-count flag; `-t` is passed only to tools that support it.

`-m` is provided to NOVOPlasty, MitoFinder and WGS assemblers, not enforced as a
global memory limit for every child process. Fastp uses at most 16 threads.
These controls are starting settings for benchmarking, not optimized parameters.

Successful checkpoints store commands and input fingerprints. Failed tool
directories can be retained with `--keep-intermediates` and retried in new `attempt_NNN` directories. Changing
reference content, biological options, input FASTQ stamps, or command configuration
requires a new output directory. Ordinary `--resume` can be used for an interrupted
run with unchanged settings. Completed compact samples are reused as a whole;
removed intermediates cannot be resumed. Retry finalized analyses in a new root. Interrupted stages are rerun, not resumed internally.
Do not run two controllers into the same output directory: an advisory lock blocks
concurrent writers. Use separate output directories for independent jobs.

Dry runs require a new/empty output directory and create a plan/checkpoint tree with hypothetical outputs and no biological
assessment. Use a separate dry-run output directory, or `--resume` when subsequently
executing that same plan. Decisions dependent on recovery cannot be predicted by
the dry run, so the plan shows the methods that could be needed.

## Circularity, completeness, rescue and fallback output

Circularity is accepted from a GetOrganelle closed graph output, a NOVOPlasty
`Circularized_assembly` output, or mtGrasp's explicit circular FASTA label.
An `incomplete` path filename never establishes circularity. MitoFinder's
single-contig circularization log is usable only with additional raw-read junction
support. Exact terminal repeats alone also require that support: at least three
distinct templates with MAPQ >=20, a continuous alignment >=20 bp on both sides
of the join, and at least two alignment start/strand combinations. Alignments
compete against the full candidate to reduce internal-repeat artifacts.

The reference gene screen requires every expected CDS and rRNA to have >=80%
single-HSP query coverage, at >=30% amino-acid identity / >=70% rRNA nucleotide
identity. Established circular candidates are doubled to handle origin-spanning
genes. Unclosed candidates are screened once with their sequence ends preserved.
The default requires all reference genes and 10-30 kb length; optional gene-count,
length, OR and AND rules are described in `CONTROLLER.md`. Every policy still
requires target confirmation and zero ambiguous bases. All these
thresholds are adjustable. This is evidence for reference gene recovery, not a
proof of full tRNA content, intact reading frames, or perfect assembly. MITOS2
annotation is retained for gene-boundary/tRNA review; its success is separate from
the recovery gate and is not claimed to independently confirm every gene.

Only noncircular candidates with `confirmed_target` COI status are eligible as
sample-specific rescue seeds. The strongest eligible partial can seed subsequent
GetOrganelle-from-reads/mtGrasp attempts. GetOrganelle's assembly mode has no
`-s` read-seed argument; verified partial sequence is appended to its nucleotide
label set instead. MitoFinder uses the annotated GenBank reference, not an
unsupported seed flag. Every rescue source is recorded in `rescue_seed.json`.

When no candidate passes, the best single contiguous candidate is still exported.
Confirmed-target sequences rank ahead of unresolved/non-target candidates; gene
recovery, ambiguity and length break ties. Contigs are not concatenated across
gaps to imitate completeness. If no candidate sequence exists, no FASTA is
fabricated. Unresolved, non-target and partial candidates retain their true status.

## Outputs

| Path | Contents |
| --- | --- |
| `cohort_summary.csv` | One row per sample: mitochondrial recovery and execution flags, plus reported ITSME locus counts, ITS/locus lengths and consensus taxa/ranks. |
| `cohort_its_loci.csv` | One row per reported complete ITSME locus: ITS cassette/region lengths, full locus length and consensus taxonomy. |
| `environment_selection.tsv`, `.json` | Selected stage environments and other matching installations; generated automatically when no explicit config is supplied. |
| `cohort_gene_recovery.csv` | Per-gene reference coverage and identity for each sample's selected candidate. |
| `SAMPLE/MitoLooper/best_sequence.fasta` | Best single contiguous candidate, potentially partial or unresolved; always read its summary. |
| `acceptance_policy.json` | Active stopping rule, topology, length thresholds and gene-count definition. |
| `SAMPLE/MitoLooper/accepted_candidate.fasta` | Candidate that met the selected stopping policy; may be incomplete or unclosed under relaxed settings. |
| `SAMPLE/MitoLooper/circular_target_candidate.fasta` | Written only for candidates passing circularity, reference gene recovery and COI target checks. |
| `SAMPLE/final/summary.json` | Selection, evidence, statuses and caveats. |
| `SAMPLE/MitoLooper/annotation_summary.json` | MITOS2 GFF feature counts and reported gene labels. |
| `SAMPLE/final/gene_recovery.csv` | CDS/rRNA homology recovery against the curated reference; CDS lengths are amino acids, rRNA lengths are nucleotides. |
| `SAMPLE/candidates.tsv` | All assessed candidate paths, gene recovery and taxonomic/circularity decisions. |
| `SAMPLE/candidates/*/bold_hits.tsv` | Raw COI BLAST alignments. |
| `SAMPLE/candidates/*/bold_evidence.json` | BOLD process IDs, taxa, provenance and decision reason. |
| `SAMPLE/candidates/*/assessment.json` | Per-candidate reference-gene, ambiguity and circularity evidence. |
| `SAMPLE/stages.tsv`, `stages.json` | Every attempted, cached, skipped, failed or timed-out stage. |
| `SAMPLE/logs/` | Exact commands and native program output. |
| `SAMPLE/MitoLooper/wgs_sampling.json` | Requested and actual assembly subset sizes, if sampling was enabled. |
| `SAMPLE/MitoLooper/mitos_*` | MITOS2 final annotation/sequence files and a GFF-derived feature TSV. |
| `SAMPLE/ITSME/*` | ITSME final FASTAs, tables and graphs, in standalone or combined mode. |
| `SAMPLE/final/outputs_manifest.tsv` | Retained output inventory and original source paths. |
| `SAMPLE/cleanup.json`, `cleanup_removed.tsv` | Compaction statistics and removed intermediate paths. |

`recovered_annotated` means the defined recovery gates passed and a nonempty MITOS2
GFF was produced. `recovered_annotation_failed` means recovery passed but annotation
failed. Relaxed-policy winners use `accepted_by_policy_annotated` or
`accepted_by_policy_annotation_failed`; those names do not claim complete or
circular mitogenomes. Other statuses include `partial_target_confirmed`, `unresolved_candidate`,
`no_candidate` and `failed`. A failed requested WGS/ITSME execution stage or attempted annotation produces
a nonzero batch exit status while subsequent samples continue. Recovery-method
failures followed by successful fallback remain visible without invalidating that
successful recovery. Compaction preserves these statuses and logs even when
failed-stage scratch is removed. Without a recovery winner, method/validation failures produce
a nonzero exit status. Biological non-recovery with successfully executed tools
is distinguished from tool execution failure.

## Testing and practical limits

Run the bundled tests with `python -m pytest -q tests` (install pytest separately).
They exercise the controller using synthetic data and CLI doubles: method ordering,
early stopping, annotation timing, WGS fallbacks, target conflicts, safe seed reuse,
resume rules, paired subsampling, junction checks, no-assembly mode and child-process
timeouts. They do not benchmark real assembler recovery or validate tool/database
compatibility on TheBelly. Start with one sample (`--only-sample W30877`) and inspect
native logs and biological outputs before the complete batch.

## Primary documentation used for tool interfaces

- NOVOPlasty and its config format: https://github.com/ndierckx/NOVOPlasty
- GetOrganelle and assembly graph input: https://github.com/Kinggerm/GetOrganelle
- MitoFinder: https://github.com/RemiAllio/MitoFinder
- mtGrasp: https://github.com/BirolLab/mtGrasp
- MITOS2: https://gitlab.com/Bernt/MITOS/-/tree/mitos2
- BOLD exports/API documentation: https://boldsystems.org/data/api/
- BOLD reference packages: https://boldsystems.org/data/data-packages/
- Minimap2: https://github.com/lh3/minimap2

The shipped package contains orchestration code and your original legacy scripts;
third-party programs, reference databases and sample biological metadata are not
bundled.


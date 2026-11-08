# MitoGenomeScanner 2.0.2

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
| 1. Mitochondria from reads | NOVOPlasty, GetOrganelle, MitoFinder, then mtGrasp. Validate each output before trying the next method. Stop remaining read methods after a circular, target-confirmed candidate passes the reference gene screen. Annotate that candidate with MITOS2 immediately. |
| 2. ITS/rDNA | `--branch full` calls the separately installed ITSME v1.9 shell pipeline on original reads. Default `--branch mito` omits nuclear rDNA. |
| 3. WGS assembly | metaSPAdes, with `--only-assembler`, a memory setting, wall-clock timeout, and optional percentage of reads. MEGAHIT is an alternative. This stage still runs after successful read recovery unless `--no-assembly` is set. |
| 4. Mitochondria from assembly | If stage 1 failed to produce a passing candidate, GetOrganelle processes the metaSPAdes FASTG graph; MitoFinder then processes contigs if needed. New passing candidates go to MITOS2. |
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

## References and sample identities

The recovery tools need suitable mitochondrial references; target validation also
needs the intended sample identity. These cannot be inferred from `W30877`.
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
| `mt_code` | Correct NCBI mitochondrial translation table for that organism. Avoid assigning one code to all marine invertebrates without checking. |
| `reference_gb` | ONE annotated, curated mitochondrial reference genome in GenBank format, preferably from a close relative. Relative paths are resolved against the sheet directory. |
| `seed_fasta` | Optional mitochondrial/COI nucleotide seed. If blank, the controller derives COX1 from `reference_gb`. |

The reference must have uniquely named, complete CDS/rRNA genes, including COX1,
at least three CDS genes and at least one rRNA. The expected CDS/rRNA gene set is
derived from that reference; it does not assume every animal has exactly 13 CDS
genes. An unrepresentative or incomplete reference weakens this screen. Repeated
gene names, fuzzy feature boundaries, ambiguous reference genes or multiple
GenBank records are rejected for explicit curation. This first version does not
model duplicated genes or multipartite/linear mitochondrial genomes.

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
quality gates, or reaching the hit ceiling leaves the candidate unresolved.

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
| `--resume` / `--skip-existing` | Reuse successful matching stages; directory existence alone never proves success. |
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
directories are retained and retried in new `attempt_NNN` directories. Changing
reference content, biological options, input FASTQ stamps, or command configuration
requires a new output directory. Ordinary `--resume` can be used for an interrupted
run with unchanged settings. Interrupted stages are rerun, not resumed internally.
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
identity. Screening a doubled candidate handles genes across the circular origin.
Candidates also need 10-30 kb length by default and zero ambiguous bases. All these
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
| `cohort_summary.csv` | One row per sample: recovery status, winning method, reference taxon, genetic code, failure/annotation/WGS flags. |
| `environment_selection.tsv`, `.json` | Selected stage environments and other matching installations; generated automatically when no explicit config is supplied. |
| `cohort_gene_recovery.csv` | Per-gene reference coverage and identity for each sample's selected candidate. |
| `SAMPLE/final/best_sequence.fasta` | Best single contiguous candidate, potentially partial or unresolved; always read its summary. |
| `SAMPLE/final/circular_target_candidate.fasta` | Written only for candidates passing circularity, reference gene recovery and COI target checks. |
| `SAMPLE/final/summary.json` | Selection, evidence, statuses and caveats. |
| `SAMPLE/final/annotation_summary.json` | MITOS2 GFF feature counts and reported gene labels. |
| `SAMPLE/final/gene_recovery.csv` | CDS/rRNA homology recovery against the curated reference; CDS lengths are amino acids, rRNA lengths are nucleotides. |
| `SAMPLE/candidates.tsv` | All assessed candidate paths, gene recovery and taxonomic/circularity decisions. |
| `SAMPLE/candidates/*/bold_hits.tsv` | Raw COI BLAST alignments. |
| `SAMPLE/candidates/*/bold_evidence.json` | BOLD process IDs, taxa, provenance and decision reason. |
| `SAMPLE/candidates/*/assessment.json` | Per-candidate reference-gene, ambiguity and circularity evidence. |
| `SAMPLE/stages.tsv`, `stages.json` | Every attempted, cached, skipped, failed or timed-out stage. |
| `SAMPLE/logs/` | Exact commands and native program output. |
| `SAMPLE/03_wgs/attempt_*/sampling.json` | Requested and actual assembly subset sizes, if sampling was enabled. |
| `SAMPLE/05_mitos2/attempt_*/output/` | MITOS2 annotation files; exact path also in summary JSON. |
| `SAMPLE/02_itsme/attempt_*/output/` | Separate ITSME results in full mode. |

`recovered_annotated` means the defined recovery gates passed and a nonempty MITOS2
GFF was produced. `recovered_annotation_failed` means recovery passed but annotation
failed. Other statuses include `partial_target_confirmed`, `unresolved_candidate`,
`no_candidate` and `failed`. A failed requested WGS/ITS stage or annotation produces
a nonzero batch exit status while subsequent samples continue. Recovery-method
failures followed by successful fallback remain visible without invalidating that
successful recovery. Without a recovery winner, method/validation failures produce
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

# Mitochondrial Genome Assembly and QC Pipeline

A restart-safe pipeline for assembling and quality-assessing mitochondrial genomes from paired-end Illumina FASTQ files. For each sample, reads are quality-trimmed and merged, then assembled in parallel by two complementary approaches — de novo assembly with GetOrganelle and reference-guided assembly with mtGrasp — followed by per-sample summarization of assembly completeness and cohort-level aggregation.

---

## Overview

```
fastqs/
├── <sample>_R1_001.fastq.gz
├── <sample>_R2_001.fastq.gz
└── <sample>.fasta              ← reference mitogenome (required for mtGrasp)
        │
  mito_pipeline.sh  ◄──────────────────────────────── main entry point
        │
        ├── [1] QC_reads.sh  (per sample)
        │         ↓
        │   Trimmomatic PE → FLASH2 overlap merge → quality trim → pigz
        │   → qc/<sample>/<sample>.notCombined.final_1P.gz  (paired)
        │   → qc/<sample>/<sample>.notCombined.final_2P.gz  (paired)
        │   → qc/<sample>/<sample>_reads.U.fq.gz            (unpaired pool)
        │
        ├── [2] GetOrganelle  (per sample, de novo)
        │         ↓
        │   get_organelle_from_reads.py -F animal_mt
        │   → getorganelle/<sample>-mitos/
        │       ├── *.complete.graph1.1.path_sequence.fasta
        │       ├── *.hits.csv
        │       └── get_org.log.txt
        │
        ├── [3] mtGrasp  (per sample, reference-guided)
        │         ↓
        │   mtgrasp.py -r <sample>.fasta -m <mt_code>
        │   → mtgrasp/<sample>/
        │
        ├── [4] mito2summary.py  (per sample, per assembler)
        │         ↓
        │   parse hits CSV + assembled FASTA + log
        │   → summaries/getorganelle/<sample>_run_summary.csv
        │   → summaries/getorganelle/<sample>_gene_lengths_long.csv
        │   → summaries/mtgrasp/<sample>_run_summary.csv
        │   → summaries/mtgrasp/<sample>_gene_lengths_long.csv
        │
        └── [5] aggregate_mito.sh  (cohort level)
                  ↓
            concatenate per-sample CSVs across all samples
            → cohort_getorganelle/cohort_run_summary.csv
            → cohort_getorganelle/cohort_gene_lengths_long.csv
            → cohort_mtgrasp/cohort_run_summary.csv
            → cohort_mtgrasp/cohort_gene_lengths_long.csv
```

Each sample is **restart-safe**: `--skip-existing` detects which stages have already completed and resumes from the furthest successful point.

---

## Repository Contents

| File | Description |
|---|---|
| `mito_pipeline.sh` | Main pipeline entry point; orchestrates all stages for every sample in a FASTQ directory |
| `QC_reads.sh` | Adapter trimming (Trimmomatic PE), read overlap merging (FLASH2), quality trimming, and gzip compression |
| `getorgLooperParallel.sh` | Standalone GetOrganelle looper with optional subsampling; superseded by `mito_pipeline.sh` for standard runs |
| `mtgrasp.sh` | Standalone mtGrasp wrapper; superseded by `mito_pipeline.sh` for standard runs |
| `mito2summary.py` | Parses a single GetOrganelle or mtGrasp output directory into a run summary CSV and a long-format gene completeness CSV |
| `aggregate_mito.sh` | Loops over all per-sample summary directories, calls `mito2summary.py` for each, and concatenates outputs into cohort-level tables |

---

## Dependencies

### Conda Environments

The pipeline uses three named conda environments. Create them before running.

#### `getorganelle` — GetOrganelle de novo assembly

```bash
conda create -n getorganelle -c bioconda -c conda-forge getorganelle
conda activate getorganelle
# Download the animal mitochondrial reference database (required once):
get_organelle_config.py --add-organelle animal_mt
# Verify:
get_organelle_from_reads.py --version
```

#### `mtgrasp` — Reference-guided mitogenome assembly

```bash
conda create -n mtgrasp -c bioconda -c conda-forge python=3.9
conda activate mtgrasp
pip install mtgrasp
# Verify:
mtgrasp.py --version
```

mtGrasp requires the following tools in the active environment:

```bash
conda install -c bioconda -c conda-forge \
    bwa samtools minimap2 seqtk ntlink abyss
```

#### `mito_qc` — Read QC tools

```bash
conda create -n mito_qc -c bioconda -c conda-forge \
    trimmomatic flash2 pigz
conda activate mito_qc
# Verify:
trimmomatic -version
flash2 --version
```

### Additional Tools (in `$PATH` or active environment)

| Tool | Version | Purpose | Install |
|---|---|---|---|
| Trimmomatic | ≥ 0.39 | Illumina PE adapter and quality trimming | `conda install -c bioconda trimmomatic` |
| FLASH2 | ≥ 2.2 | Paired-read overlap merging | `conda install -c bioconda flash2` |
| pigz | ≥ 2.6 | Parallel gzip compression | `conda install -c conda-forge pigz` |
| GNU parallel | any | Parallel sample processing (`-j > 1`) | `conda install -c conda-forge parallel` |
| Python | ≥ 3.8 | `mito2summary.py` and `aggregate_mito.sh` | system or conda |

### Trimmomatic Adapter File

`QC_reads.sh` requires the Trimmomatic PE adapter FASTA. Update the adapter path in `QC_reads.sh` to match your installation:

```bash
# Typical path in a conda environment:
ADAPTERS=$(dirname $(which trimmomatic))/../share/trimmomatic/adapters/TruSeq3-PE-2.fa

# Or with a manual installation:
ADAPTERS=/path/to/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa
```

### Path Configuration

The following paths are hardcoded and must be updated before running on a new system:

| Script | Variable / location | Default value |
|---|---|---|
| `mito_pipeline.sh` | `TRIMMOMATIC_PATH` (line 36) | `/home/ark/bin/Trimmomatic-0.39` |
| `mito_pipeline.sh` | `MITO2SUMMARY` | resolved from PATH; override with `--mito2summary` |
| `mito_pipeline.sh` | `AGGREGATE` | resolved from PATH; override with `--aggregate` |
| `QC_reads.sh` | adapter FASTA path | hardcoded; update before first use |

---

## Input Data

### FASTQ Files

Place all sample FASTQ files in a single flat directory. Files must follow the Illumina naming convention:

```
fastqs/
├── SampleA_R1_001.fastq.gz
├── SampleA_R2_001.fastq.gz
├── SampleB_R1_001.fastq.gz
└── SampleB_R2_001.fastq.gz
```

The sample name is parsed as everything before `_R1_001.fastq.gz`. Do not use underscores in positions that would conflict with this suffix.

### Reference Mitogenomes (mtGrasp only)

For each sample that will be processed by mtGrasp, place a reference mitogenome FASTA in the same directory as the FASTQs (or a separate directory specified with `--ref-dir`):

```
fastqs/
├── SampleA_R1_001.fastq.gz
├── SampleA_R2_001.fastq.gz
└── SampleA.fasta              ← reference for SampleA
```

If no `<sample>.fasta` is found for a given sample, mtGrasp is skipped for that sample with a warning; GetOrganelle still runs normally.

### Expected Gene Lengths (optional)

A two-column CSV (`Gene,Average_Length_bp`) used by `mito2summary.py` to compute per-gene completeness fractions. If not provided, raw lengths are reported without a completeness percentage:

```
Gene,Average_Length_bp
COX1,1545
COX2,684
COX3,789
ND1,957
...
```

Pass this file with `--expected-csv`.

---

## Pipeline SOP

### Step 0 — Setup

```bash
# Add scripts to PATH:
export PATH=$PATH:/path/to/this/repo
chmod +x mito_pipeline.sh QC_reads.sh aggregate_mito.sh

# Update TRIMMOMATIC_PATH in mito_pipeline.sh:
# TRIMMOMATIC_PATH="/path/to/your/Trimmomatic-0.39"

# Verify tools are available:
conda activate getorganelle && get_organelle_from_reads.py --version
conda activate mtgrasp      && mtgrasp.py --version
conda activate mito_qc      && trimmomatic -version && flash2 --version
```

---

### Step 1 — Run the Pipeline

```bash
bash mito_pipeline.sh \
    --mt-code 2 \
    -o mito_results \
    -t 16 \
    fastqs/
```

#### Arguments

| Flag | Required | Description |
|---|---|---|
| `<fastq_dir>` | Yes | Directory containing `*_R1_001.fastq.gz` / `*_R2_001.fastq.gz` pairs |
| `-o, --outdir DIR` | No | Output directory (default: `mito_out`) |
| `-t, --threads N` | No | Threads per tool per sample (default: 12) |
| `-j, --jobs N` | No | Parallel samples via GNU parallel (default: 1) |
| `--skip-existing` | No | Skip a stage if its output directory already exists |
| `-v, --verbose` | No | Print extra progress messages |

#### QC options

| Flag | Description |
|---|---|
| `--no-qc` | Skip `QC_reads.sh`; raw reads are passed directly to assemblers |
| `--qc-minlen N` | Trimmomatic `SLIDINGWINDOW` quality threshold (default: 30) |
| `--trimmomatic-path P` | Path to Trimmomatic installation directory |

#### GetOrganelle options

| Flag | Description |
|---|---|
| `--no-getorg` | Skip GetOrganelle assembly |
| `--getorg-rounds N` | Maximum extension rounds (default: 30) |

#### mtGrasp options

| Flag | Description |
|---|---|
| `--no-mtgrasp` | Skip mtGrasp assembly |
| `--mt-code N` | NCBI translation table number (required when mtGrasp is active) |
| `--ref-dir DIR` | Directory containing `<sample>.fasta` references (default: `<fastq_dir>`) |

#### Post-processing options

| Flag | Description |
|---|---|
| `--mito2summary PATH` | Full path to `mito2summary.py` (default: searched in PATH) |
| `--aggregate PATH` | Full path to `aggregate_mito.sh` (default: searched in PATH) |
| `--expected-csv FILE` | CSV of expected gene lengths passed to `mito2summary.py` |

#### NCBI Translation Table Reference

| Code | Taxon |
|---|---|
| 2 | Vertebrate mitochondria |
| 5 | Invertebrate mitochondria |
| 9 | Echinoderm / flatworm mitochondria |
| 13 | Ascidian mitochondria |
| 14 | Flatworm mitochondria |

See [NCBI Genetic Codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) for the full list.

---

### Step 2 — Read QC (`QC_reads.sh`)

Called automatically per sample; can also be run standalone.

```bash
bash QC_reads.sh \
    <threads> \
    <sample_name> \
    <workdir> \
    <R1.fastq.gz> <R2.fastq.gz> \
    trimmomatic \
    <trimmomatic_dir> \
    flash2 \
    pigz \
    <sliding_window_quality>
```

`QC_reads.sh` performs three sequential operations:

1. **Trimmomatic PE** — removes adapters and low-quality bases; writes paired and unpaired read sets
2. **FLASH2** — merges overlapping read pairs into extended fragments; non-overlapping pairs are retained separately
3. **Quality trim** — applies a second Trimmomatic pass to merged reads; all outputs are gzip-compressed with pigz

#### Outputs

| File | Description |
|---|---|
| `<sample>.notCombined.final_1P.gz` | Trimmed R1 reads with surviving mate (non-overlapping) |
| `<sample>.notCombined.final_2P.gz` | Trimmed R2 reads with surviving mate (non-overlapping) |
| `<sample>_reads.U.fq.gz` | Unpaired reads (singletons from Trimmomatic + FLASH2 merged reads) |

The paired outputs are passed to both assemblers. The unpaired pool is passed to GetOrganelle only via `-u`.

---

### Step 3 — De Novo Assembly (`GetOrganelle`)

Called automatically per sample; the output directory is `getorganelle/<sample>-mitos/`.

GetOrganelle is invoked as:

```
get_organelle_from_reads.py \
    -1 <R1> -2 <R2> -u <unpaired> \
    -F animal_mt \
    -o getorganelle/<sample>-mitos/ \
    -t <threads> \
    -R 30 \
    --max-reads inf \
    --disentangle-time-limit 200000000 \
    --reduce-reads-for-coverage inf
```

#### Outputs

| File | Description |
|---|---|
| `*.complete.graph1.1.path_sequence.fasta` | Complete circular assembly (preferred by downstream steps) |
| `*.contigs.graph1.1.path_sequence.fasta` | Incomplete / fragmented assembly (used if complete not available) |
| `*.hits.csv` | Per-read alignment hits table; input to `mito2summary.py` |
| `get_org.log.txt` | Full GetOrganelle run log; input to `mito2summary.py` |

If GetOrganelle produces multiple graph paths, the first complete path sequence is used by `mito2summary.py`.

---

### Step 4 — Reference-Guided Assembly (`mtGrasp`)

Called automatically per sample when a `<sample>.fasta` reference is found; the output directory is `mtgrasp/<sample>/`.

mtGrasp is invoked as:

```
mtgrasp.py \
    -r1 <R1> -r2 <R2> \
    -r  <sample>.fasta \
    -m  <mt_code> \
    -o  mtgrasp/<sample>/ \
    -t  <threads>
```

If no reference FASTA is found for a sample, mtGrasp is skipped for that sample with a warning and GetOrganelle output is still processed normally.

---

### Step 5 — Per-Sample Summary (`mito2summary.py`)

Called automatically after each assembler, once per assembler per sample. Can also be run standalone on any GetOrganelle or compatible output directory.

```bash
python mito2summary.py \
    -c <hits.csv> \
    -f <assembly.fasta> \
    -l <get_org.log.txt> \
    -o <output_prefix> \
    [--expected-csv <gene_lengths.csv>]
```

#### Arguments

| Flag | Required | Description |
|---|---|---|
| `-c` | Yes | GetOrganelle hits CSV |
| `-f` | Yes | Assembled FASTA (complete path sequence preferred) |
| `-l` | Yes | GetOrganelle log file (`get_org.log.txt`) |
| `-o` | Yes | Output prefix (directory + base name, no extension) |
| `--expected-csv` | No | CSV of expected gene lengths (`Gene,Average_Length_bp`) |

#### Outputs

| File | Description |
|---|---|
| `<prefix>_run_summary.csv` | One row per sample: assembly stats, total mitochondrial bp, number of fragments, run metadata |
| `<prefix>_gene_lengths_long.csv` | Long-format table: one row per gene per sample, with observed and (if `--expected-csv` provided) expected length and completeness fraction |

---

### Step 6 — Cohort Aggregation (`aggregate_mito.sh`)

Called automatically after all per-sample processing is complete. Can also be run standalone on any directory of per-sample summary CSVs.

```bash
bash aggregate_mito.sh \
    -o <output_dir> \
    [--expected-csv <gene_lengths.csv>] \
    [--suffix ""] \
    <summaries_directory>
```

#### Arguments

| Flag | Required | Description |
|---|---|---|
| `<summaries_directory>` | Yes | Directory containing per-sample `*_run_summary.csv` and `*_gene_lengths_long.csv` files |
| `-o, --outdir DIR` | No | Output directory for cohort tables (default: `mito_cohort_out`) |
| `--mito2summary PATH` | No | Full path to `mito2summary.py` if not in PATH |
| `--expected-csv FILE` | No | CSV of expected gene lengths; forwarded to `mito2summary.py` |
| `--suffix STR` | No | Suffix filter for sample subdirectories (default: `-mitos`; pass `""` for all) |
| `-j, --jobs N` | No | Parallel samples via GNU parallel (default: 1) |
| `--skip-existing` | No | Skip samples whose per-sample CSVs already exist |

#### Outputs

| File | Description |
|---|---|
| `<outdir>/cohort_run_summary.csv` | All per-sample run summaries concatenated into a single table |
| `<outdir>/cohort_gene_lengths_long.csv` | All per-sample gene length tables concatenated into a single table |
| `<outdir>/aggregate_mito.log` | Full aggregation log |

---

## Complete Usage Examples

### Full run, vertebrate mitogenomes, 8 parallel samples

```bash
bash mito_pipeline.sh \
    --mt-code 2 \
    -o mito_results \
    -t 16 \
    -j 8 \
    --expected-csv vertebrate_mt_gene_lengths.csv \
    fastqs/
```

### GetOrganelle only (no reference required)

```bash
bash mito_pipeline.sh \
    --no-mtgrasp \
    -o getorg_results \
    -t 12 \
    fastqs/
```

### mtGrasp only, references in a separate directory, skip QC

```bash
bash mito_pipeline.sh \
    --no-getorg \
    --no-qc \
    --mt-code 5 \
    --ref-dir /path/to/references/ \
    -o mtgrasp_results \
    -t 16 \
    fastqs/
```

### Resume an interrupted run

```bash
bash mito_pipeline.sh \
    --mt-code 2 \
    --skip-existing \
    -o mito_results \
    -t 16 \
    fastqs/
```

### Standalone cohort aggregation after a completed run

```bash
bash aggregate_mito.sh \
    --suffix "" \
    --expected-csv vertebrate_mt_gene_lengths.csv \
    -o mito_results/cohort_getorganelle \
    mito_results/summaries/getorganelle/

bash aggregate_mito.sh \
    --suffix "" \
    --expected-csv vertebrate_mt_gene_lengths.csv \
    -o mito_results/cohort_mtgrasp \
    mito_results/summaries/mtgrasp/
```

---

## Output Structure

After a complete run with both assemblers:

```
mito_results/
├── qc/
│   └── <sample>/
│       ├── <sample>.notCombined.final_1P.gz
│       ├── <sample>.notCombined.final_2P.gz
│       └── <sample>_reads.U.fq.gz
├── getorganelle/
│   └── <sample>-mitos/
│       ├── *.complete.graph1.1.path_sequence.fasta
│       ├── *.hits.csv
│       └── get_org.log.txt
├── mtgrasp/
│   └── <sample>/
│       └── (mtGrasp output files)
├── summaries/
│   ├── getorganelle/
│   │   ├── <sample>_run_summary.csv
│   │   └── <sample>_gene_lengths_long.csv
│   └── mtgrasp/
│       ├── <sample>_run_summary.csv
│       └── <sample>_gene_lengths_long.csv
├── cohort_getorganelle/
│   ├── cohort_run_summary.csv        ← primary deliverable
│   ├── cohort_gene_lengths_long.csv  ← primary deliverable
│   └── aggregate_mito.log
├── cohort_mtgrasp/
│   ├── cohort_run_summary.csv
│   ├── cohort_gene_lengths_long.csv
│   └── aggregate_mito.log
├── logs/
│   ├── <sample>.qc.log
│   ├── <sample>.getorg.log
│   ├── <sample>.mtgrasp.log
│   └── <sample>.mito2summary.log
└── mito_pipeline.log                 ← full pipeline log
```

---

## Practical Notes

- **Conda activation in subshells:** `mito_pipeline.sh` calls `eval "$(conda shell.bash hook)"` before activating each environment. If your conda installation requires `source ~/.bashrc` or a custom init path, update the activation blocks in `run_getorganelle()` and `run_mtgrasp()` accordingly.
- **GetOrganelle multiple paths:** When GetOrganelle resolves multiple graph paths, it writes numbered FASTA files (`.graph1.1`, `.graph1.2`, etc.). The pipeline and `mito2summary.py` use the first complete path sequence found. Review the `get_org.log.txt` if an assembly looks unexpectedly small or fragmented.
- **mtGrasp reference quality:** mtGrasp performs best with a conspecific or closely related reference (ideally > 90% identity to the target). A highly divergent reference will reduce mapping efficiency and may produce a fragmented assembly. GetOrganelle output from the same or a related sample is a good reference source.
- **Unpaired reads routing:** The unpaired read pool (`_reads.U.fq.gz`) is passed to GetOrganelle via `-u` but not to mtGrasp, which requires strictly paired input. This is intentional.
- **`mito2summary.py` sample name parsing:** The script parses the sample name by stripping the `_getOrganelleReads_NOseed` suffix pattern from the hits CSV filename. For mtGrasp outputs, this pattern will not match and the full filename is used. Verify that output filenames are being parsed correctly if sample names appear malformed in the cohort tables.
- **Parallel jobs and conda:** When using `-j > 1`, each sample's conda activation runs in a GNU parallel subshell. Ensure conda is initialized in the shell environment before calling the pipeline (`conda init bash` or equivalent) and that `~/.bashrc` is sourced.
- **Disk space:** GetOrganelle can write several GB of intermediate read files per sample. By default these are retained in `getorganelle/<sample>-mitos/`. If disk space is a concern, add `--remove-duplicates` to the GetOrganelle call or clean intermediate files after the run.

---

## References

- Jin J-J et al. GetOrganelle: a fast and versatile toolkit for accurate de novo assembly of organelle genomes. *Genome Biology* (2020). doi: [10.1186/s13059-020-02154-5](https://doi.org/10.1186/s13059-020-02154-5)
- Coombe L et al. mtGrasp: reference-guided mitochondrial genome assembly. *Bioinformatics Advances* (2023). doi: [10.1093/bioadv/vbad171](https://doi.org/10.1093/bioadv/vbad171)
- Bolger AM et al. Trimmomatic: a flexible trimmer for Illumina sequence data. *Bioinformatics* (2014). doi: [10.1093/bioinformatics/btu170](https://doi.org/10.1093/bioinformatics/btu170)
- Magoč T and Salzberg SL. FLASH: fast length adjustment of short reads to improve genome assemblies. *Bioinformatics* (2011). doi: [10.1093/bioinformatics/btr507](https://doi.org/10.1093/bioinformatics/btr507)

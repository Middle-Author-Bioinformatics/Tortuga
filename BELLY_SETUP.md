# Run on TheBelly: automatic environment switching

Use the existing installations. No environment creation, manual activation, or
config-generation step is needed for normal execution.

## Run

After extracting the package to `~/bin/MitoGenomeScanner`, run from `base`:

```bash
cd ~/bin/MitoGenomeScanner

bash mito_looper.sh ~/NOAA/CA_intertidal/raw_reads \
    --samples samples.tsv \
    --bold-fasta /PATH/TO/bold/coi.fasta \
    --bold-metadata /PATH/TO/bold/metadata.tsv \
    --mitos-refdir /PATH/TO/mitos_data \
    --mitos-refseqver refseq89m \
    --branch mito --only-sample W30877 \
    -o ~/NOAA/CA_intertidal/mito_results \
    -t 60 -m 200
```

Replace capitalized reference paths with actual locations and use your installed
MITOS2 release name. After checking the first sample, remove `--only-sample W30877`
and add `--resume` to run the remaining samples sequentially. The controller already
loops over the input directory.

Add `--no-assembly` to omit the separate WGS and post-WGS recovery stages. Initial
mitochondrial tools can still assemble reads internally. `--assembly-percent 10`
gives the WGS stage an approximately 10% paired subset. Keep biological settings
unchanged between a trial and its resumed batch.

## Automatic choices and notes about alternatives

The launcher starts the Python controller in `/home/ark/miniconda3/envs/itsme`.
Each stage then runs through `conda run -p` in the first matching environment,
without changing the user's parent shell.

| Component | First choice | Other relevant installation |
| --- | --- | --- |
| Controller, fastp, minimap2 | `itsme` | None selected from the supplied inventory |
| metaSPAdes | `itsme` | `mito-spades` |
| BLAST+ | `itsme` if its executables exist | `mito-blast`; expected choice with your supplied package list |
| NOVOPlasty | `mito-novoplasty` | Actual Perl script detected from bin/package metadata |
| GetOrganelle | `getorganelle` | `mito-getorganelle` |
| MitoFinder | `mito-mitofinder` | None selected from the supplied inventory |
| mtGrasp | `mito-mtgrasp` | `mtgrasp` |
| MITOS2 | `mito-mitos2` | `mitos2` if subsequently present |
| Optional MEGAHIT | `mito-megahit` | `itsme` if subsequently present |
| Optional nuclear rDNA | `itsme` | Calls `~/bin/ITSME/itsme.v2.sh` |

The source comments beside `SPECS` in `configure_belly.py` note the alternatives;
the launcher also has comments for the main duplicates. The first environment
containing the executable wins. Found alternatives are noted in stdout and in
`environment_selection.tsv` and `.json` under the output directory. No prompts
occur. A failed chosen program is reported as a tool failure; a different package
version is not silently retried.

The absolute prefixes prevent duplicate environment names in Miniforge from
selecting a different installation. No environments are installed, merged or
removed. The old Snakemake entry point is not used; mtGrasp retains its own
internal dependencies.

For another host, `MITO_ENV_ROOT=/other/conda/envs` selects the environment root
for both the shell launcher and tool selection. `MITO_CONDA_EXE=/path/to/conda`
overrides Conda. `--env-root`/`--conda` override the controller's tool selection.
`MITO_PYTHON` can explicitly override the controller interpreter. An explicit
`--config existing.json` remains supported for deliberate manual command mappings;
it bypasses automatic tool selection.

## Required metadata and references

Environment switching does not replace biological reference data. Copy and fill
the sample sheet:

```bash
cp examples/CA_intertidal_samples.tsv samples.tsv
```

Fill `target_taxon`, `mt_code`, and `reference_gb` for all twelve samples. The
optional `seed_fasta` may remain blank: COX1 from the GenBank reference supplies the
initial seed. Multiple samples can share a suitable reference.

| Resource | Requirement |
| --- | --- |
| Mitochondrial GenBank reference | One curated, annotated mitochondrial genome appropriate for each sample, recorded in `samples.tsv`. |
| GetOrganelle `animal_mt` resources | Its mitochondrial seed/label databases; these may already be installed. |
| MITOS2 models | An extracted release such as `refseq89m` below `--mitos-refdir`. |
| BOLD COI FASTA and metadata TSV | Matching sequence IDs, taxonomy and provenance for target validation. |
| ITSME references, full branch only | `/home/ark/databases/itsme_db`, already shown in your server output. |

The rRNA resources `itsme2_db` and `itsme2_silva138_2` do not replace the
mitochondrial GenBank reference, BOLD library, or MITOS2 models. `prepare_bold.py`
normalizes supported public BOLD TSV exports; see README. MITOS2 reference downloads
are provided at https://zenodo.org/records/4284483; reuse existing copies if present.

For the optional nuclear rDNA step, use `--branch full` with
`--itsme-db /home/ark/databases/itsme_db`. ITSME's own process must be able to find
its dependencies, including BLAST+. Selecting `mito-blast` for the mitochondrial
controller does not install BLAST+ into `itsme`.

## Optional diagnostics

These are not required setup steps:

```bash
conda activate /home/ark/miniconda3/envs/itsme
python check_tools.py
```

`configure_belly.py` also remains as an optional audit/manual-config exporter,
including a bounded search for likely reference locations. Normal workflow execution
calls its selection functions internally and needs no generated config.

Use a new output directory for a run started with an older pipeline version;
changed version/environment settings are not silently treated as equivalent
checkpoints. Tests exercise orchestration with synthetic tools, not actual
assembler recovery on the project data.

# Run on TheBelly: automatic environment switching

For the new unified entry point, use `pipeline_controller.sh` with
`--pipeline itsme`, `--pipeline mitolooper`, or `--pipeline both`; see
`CONTROLLER.md`. The direct mitochondrial launcher below remains supported.

Use the existing installations. No environment creation, manual activation, or
config-generation step is needed for normal execution.

For a new installation, follow **[Download and prepare the databases](DATABASES.md)**
before running. It includes download sources, preparation commands and the exact
resource paths used below. Skip steps for databases already prepared on TheBelly.

## Run

From the cloned Tortuga repository, run from `base`:

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
    --only-sample W30877 \
    -o /home/ark/MAB/NOAA/CA_intertidal/both_results_new \
    -t 32 -m 200
```

Change `--pipeline both` to `--pipeline mitolooper` for mitochondria only, or
`--pipeline itsme` for nuclear rDNA only. The mitochondrial and ITSME paths shown
match the current Belly setup. After checking the first sample, remove `--only-sample W30877`
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
| Optional nuclear rDNA | `itsme` | Calls `/home/ark/MAB/bin/ITSME/itsme.v2.sh` |

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

Fill `target_taxon` for each sample. Leave `reference_gb` blank to select from
`/home/ark/MAB/NOAA/mitos_genbank` automatically, or set `--genbank-dir` to another collection.
Leave `mt_code` blank to infer it when matching GenBank CDS annotations agree on
one explicit code; fill it when the controller reports missing/conflicting codes.
The optional `seed_fasta` may remain blank: COX1 from the selected GenBank reference
supplies the initial seed. Multiple samples can share a suitable reference.
Retain the full kingdom-to-species hierarchy from source metadata. The default
selector uses the target, genus, family, order, class, phylum, then kingdom,
stopping at the first eligible match. `--reference-min-rank` can restrict how far
it searches. Selection reports the
actual matching rank and keeps the expected BOLD target unchanged.

| Resource | Requirement |
| --- | --- |
| Mitochondrial GenBank references | Local annotated collection selected by `target_taxon`, or a per-sample `reference_gb` override. Selection is reported in `reference_selection.tsv`. |
| GetOrganelle `animal_mt` resources | Its mitochondrial seed/label databases; these may already be installed. |
| MITOS2 models | An extracted release such as `refseq89m` below `--mitos-refdir`. |
| BOLD COI FASTA and metadata TSV | Matching sequence IDs, taxonomy and provenance for target validation. |
| ITSME references, full branch only | `/home/ark/databases/itsme_db`, already shown in your server output. |

The rRNA resources `itsme2_db` and `itsme2_silva138_2` do not replace the
mitochondrial GenBank reference, BOLD library, or MITOS2 models. Follow
[DATABASES.md](DATABASES.md) for downloading/extracting MITOS models, preparing
single-record GenBank references, downloading the correct BOLD TSV and producing
matching COI FASTA/taxonomy files. Reuse existing prepared copies if present.

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


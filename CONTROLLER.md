# One controller for ITSME, MitoLooper, or both

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
    -i ~/NOAA/CA_intertidal/raw_reads \
    -o ~/NOAA/CA_intertidal/itsme_results \
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
    -i ~/NOAA/CA_intertidal/raw_reads \
    --samples ~/bin/MitoGenomeScanner/samples.tsv \
    --bold-fasta /PATH/TO/bold/coi.fasta \
    --bold-metadata /PATH/TO/bold/metadata.tsv \
    --mitos-refdir /PATH/TO/mitos_data \
    --mitos-refseqver refseq89m \
    --itsme-db /home/ark/databases/itsme_db \
    -o ~/NOAA/CA_intertidal/both_results \
    -t 60 -m 200
```

For mitochondrial recovery only, change `--pipeline both` to
`--pipeline mitolooper` and select an appropriate output directory. All other
mitochondrial options from `mito_looper.sh` remain available. The controller does
not require a generated configuration file or manual environment activation.

## Common controls and outputs

- `--only-sample W30877` runs one sample for a trial.
- `--resume` reuses completed stages with matching inputs/settings. Failed stages
  are retried in a new attempt directory; directory existence alone never proves success.
- `--no-assembly` disables the separate WGS and post-WGS mito recovery stages.
  In `both` mode ITSME still runs and does its own targeted assembly.
- `--itsme-timeout 180` limits each ITSME call to 180 minutes by default.
- `--dry-run` prints planned commands in a new/empty output directory.

When `-o` is omitted, the default directory is `itsme_results`, `mito_results`, or
`both_results` according to the selector. Use a separate output directory when
switching modes; checkpoints from another mode/version are not silently reused.

ITSME-only results are at `OUT/SAMPLE/itsme/attempt_NNN/output/`, with the native
`master_summary.csv`, FASTAs, coverage and other files written there by ITSME.
The actual path is also recorded in `OUT/cohort_summary.csv` and
`OUT/SAMPLE/final/summary.json`. Logs and stage checkpoints are alongside them.
Its standalone batch continues after individual sample failures and returns a
nonzero status if any sample fails. The mitochondrial/both layouts are described
in README, with ITSME under `02_itsme/attempt_NNN/output/` in combined runs.

The old `mito_looper.sh` and legacy `--branch mito/full` flags remain available.
For the new controller, use the single `--pipeline` selector; incompatible legacy
branch combinations are rejected. No third-party assembly or ITSME installation
is bundled. The tests use synthetic tool doubles rather than real project analyses.

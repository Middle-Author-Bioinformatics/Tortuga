# One controller for ITSME, MitoLooper, or both

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
    -o /home/ark/MAB/NOAA/CA_intertidal/both_results_v205 \
    -t 32 -m 200
```

Use `--pipeline mitolooper` to omit ITSME. New controller versions/options require
a fresh output root; existing older-version checkpoint directories are retained.

If a sample has no usable reference, `both` still runs its ITSME stage under
`OUT/unresolved_references_itsme/SAMPLE/itsme/attempt_NNN/output/`; it skips that
sample's mitochondrial recovery and WGS. `reference_selection.tsv` records the
failure and any boundary warnings, and `reference_failures.tsv` lists unresolved
samples. The cohort and invocation summaries retain `reference_unavailable` and
ITSME success/failure separately. Resume with `--resume` to reuse completed
independent ITSME attempts; use a new output root after a version change.

The Rhizophysa sample still requires separate linear-genome evaluation; this
release changes reference selection, not the circularity acceptance criterion.

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

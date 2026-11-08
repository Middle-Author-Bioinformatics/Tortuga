# 2.0.2 - 2026-11-07

- Added `pipeline_controller.sh` with `--pipeline itsme`, `mitolooper`, or `both`.
- Standalone ITSME mode processes FASTQ pairs without mitochondrial references,
  BOLD, MITOS2 models, or a required sample sheet. Its native outputs are retained.
- Combined mode uses the existing integrated ordering so ITSME runs once per sample.
- Added `-i/--input-dir`, mode-specific default output folders, ITSME sensitivity,
  matching checkpoints, failure continuation, and controller usage documentation.
- Shared the automatic Python/Conda launcher across entry points and limited
  ITSME-only environment reports to its selected tool.
- Added seven selector/standalone batch regression tests.

# 2.0.1 - 2026-11-07

- Made environment selection automatic during execution. No configuration-generation
  step or manual activation is required for the standard TheBelly launcher.
- The launcher starts its Python controller in `itsme`; each stage runs in its own
  selected Conda environment through an absolute prefix.
- The first available installation wins. Code comments and runtime reports record
  other relevant environments, including duplicate GetOrganelle/mtGrasp installations.
- Kept explicit `--config` overrides and the optional diagnostic/setup helper.
- Added tests for per-stage switching, launch from an unactivated shell, automatic
  selection without a config file, and reporting missing executables.

# 2.0.0 - 2026-11-07

## Existing-environment setup update

- Added `configure_belly.py` to discover existing executables and write a config
  using absolute Conda environment prefixes, including actual NOVOPlasty entry points.
- Added bounded discovery of likely existing reference locations and an audit report.
- Added `BELLY_SETUP.md` with server-specific install/reuse requirements and trial/batch commands.
- Added four setup-helper tests; the core workflow is unchanged.

- Replaced the old two-method looper with a sequential proposal-driven controller.
- Added NOVOPlasty, GetOrganelle, MitoFinder and mtGrasp read-recovery stages.
- Added COI comparisons to local BOLD-derived records, provenance and explicit
  target/unresolved/non-target decisions.
- Added reference CDS/rRNA recovery screens and method-specific circularity evidence.
- Added raw-read junction verification for unsupported circularity candidates.
- Added MITOS2 annotation after passing recovery, including feature summaries.
- Added optional standalone ITSME integration and a mitochondrial-only branch.
- Added bounded metaSPAdes/MEGAHIT WGS assembly, paired read subsetting and
  GetOrganelle/MitoFinder post-assembly fallback.
- Added verified partial rescue seeds; no unconfirmed taxon is reused as a seed.
- Added support for NOAA `_1.fq.gz` / `_2.fq.gz` names, per-sample genetic codes,
  failed-attempt retention, checkpoints and cohort summaries.
- Added explicit `--no-assembly`, per-stage deadlines, dry runs and tool probes.
- Kept the uploaded scripts unchanged under `legacy/`; removed reliance on the
  missing `aggregate_mito.sh` helper by writing cohort reports in the controller.
- Added synthetic integration/regression tests. Real assembler performance on
  the project dataset remains to be benchmarked.

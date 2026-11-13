# 2.0.8 - 2026-11-12

- Fix an overly strict cleanup prerequisite: missing batch last_invocation.json
  no longer prevents compacting individually completed samples in an inactive
  or interrupted batch. Recognize ordinary Tortuga sample checkpoints instead.
- Preserve the active-controller lock guard. Skip unfinished/planned samples,
  unstarted manifest entries and stale summaries from interrupted retries.
- Retain the shared BOLD index unless the batch completion record and individual
  completion checks establish that no work remains. Do not invent completion.
- Support parent-final export of finished independent ITSME results even without
  the enclosing batch record; retain reference-unavailable status.
- Give distinct missing-path, unrecognized-folder and active-controller messages.
- Add seven regression cases for partial batches, locks without completion
  records, stale/planned/unstarted samples, independent ITSME and missing paths.

# 2.0.7 - 2026-11-12

- Compact results by default after each completed sample. Export ITSME final
  FASTAs, tables and graphs directly into the sample final directory with itsme_
  prefixes, alongside mitochondrial products. Export final WGS products and
  native MITOS2 annotations; add a GFF-derived feature TSV.
- Delete generated recruited/QC/subsampled reads, mapping files, assembler K*
  scratch, temporary candidate subjects and duplicated stage trees after exports
  complete. Remove the copied BOLD index after the cohort finishes, retaining
  provenance, QC summaries, native/wrapper logs and assessment evidence.
- Preserve input files/shared databases; avoid traversing symlinks, refuse active
  run cleanup, preserve conflicting products and retain failed/timed-out status.
- Add compact_outputs.sh for existing completed runs and --keep-intermediates for
  debugging/retries. Compact --resume reuses verified finalized outputs without
  rebuilding BOLD or rerunning tools; removed stage work requires a fresh run.
- Record retained-file mappings and cleanup audits. Keep annotation/candidate
  pointers current and expose reference-unavailable ITSME products in parent final.
- Add regression coverage for cleanup/export behavior, existing-run migration,
  input/link/lock protection, export failure, collisions, annotation TSVs and
  compact resume in combined and standalone modes.

# 2.0.6 - 2026-11-12

- Add configurable stopping policies: circularity required or optional; all
  reference genes, a minimum detected gene count, a minimum contig length,
  either threshold or both. Default acceptance remains strict.
- Count unique qualifying CDS+rRNA genes per contig, without pooling across
  methods; keep target confirmation, ambiguity and length safety gates.
- Record acceptance policy/checks/reasons in console, candidate and cohort
  reports. Label relaxed winners accepted_by_policy separately from strict
  recovery, and retain missing-gene and circularity evidence.
- Export accepted_candidate.fasta for policy winners. Preserve the legacy strict
  circular export only when its original criteria pass. Annotate unclosed
  candidates with MITOS2 --linear; report annotation failures as nonzero.
- Establish topology before homology screening. Screen unclosed contigs once
  with ends preserved, and double/trim terminal overlap only for established
  circular candidates. Keep raw-read junction evidence in the reports.
- Add --skip-wgs-if-accepted while preserving independent ITSME in combined mode.
- Add 21 regression cases covering thresholds, logic, linear treatment, target
  rejection, single-contig counts, fallback stopping, annotation and safety checks.

# 2.0.5 - 2026-11-12

- Distinguish uncertain annotation boundaries from ambiguous nucleotide bases.
  Permit before/after boundaries for homology queries with unambiguous gene DNA;
  require explicit translations for uncertain CDS annotations. Report warnings
  in reference selection, reference quality and sample summaries. Reject genuinely
  ambiguous genes and unsupported uncertain coordinate types.
- Recognize rRNA product labels s-rRNA/l-rRNA and small/large subunit aliases when
  gene qualifiers are absent. NC_007933.1 now yields 13 CDS and both rRNAs.
- Prefer fewer boundary warnings within the same taxonomic match tier. Refresh
  old catalog validations automatically; retain metadata for rejected records.
- Continue other samples after individual reference errors. In both mode, run
  independent ITSME for unresolved samples; report reference_unavailable and a
  nonzero batch status, without attempting their mitochondrial recovery or WGS.
- Added seven regression cases covering the public NC_007933.1 record, invalid
  DNA/annotations, cache refresh, reference ranking and batch continuation.

# 2.0.4 - 2026-11-10

- Restored the full supplied taxonomy and source identifiers in the twelve-sample
  CA intertidal sheet. Updated the template to include rank columns.
- Added default hierarchy selection: target, genus, then supplied family; optional
  broader fallback stops at the order/class/phylum selected by --reference-min-rank.
- Require agreement with a supplied phylum. Prefer the closest available named
  rank before annotation quality, and retain exact/genus-only modes.
- Record matched rank/name, preserve all input columns in resolved_samples.tsv,
  and retain normalized input hierarchy in checkpoints and sample summaries.
- Keep the original target for BOLD identification even when a family-level seed
  is chosen. Genetic-code ambiguity still requires an explicit override.
- Added six regression tests for hierarchy ranking, scope limits, phylum mismatch,
  code conflicts, metadata preservation and family-seeded target confirmation.

# 2.0.3 - 2026-11-09

- Added automatic local GenBank reference selection from target_taxon; explicit
  reference_gb overrides remain supported. Default collection: ~/NOAA/mitos_genbank.
- Match exact organism/lineage first, then allow a recorded same-genus fallback;
  --reference-match exact disables that fallback. No silent cross-genus fallback.
- Validate animal mitochondrial identity and CDS/rRNA annotations; apply genome
  length and supplied genetic-code filters before deterministic selection.
- Infer a blank mt_code only when eligible reference CDS annotations agree on one
  explicit genetic code. Report choices in reference_selection.tsv and resolved_samples.tsv.
- Cache unchanged GenBank entries in an incremental SQLite index, including
  exclusion reasons, to avoid repeatedly parsing the full collection.

# BOLD preparation fix - 2026-11-08

- Raised the TSV parser field limit from Python's default 128 KiB to 2 GiB minus
  one byte for large BOLD export fields. Added a regression test using a 200 KB
  field and concise reporting of CSV parsing errors.

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

# 2.4.2 - 2026-11-19

- Keep retained source graphs downloadable when the selected sequence cannot be
  mapped onto them. Distinguish an unavailable selected path from a missing
  source graph, without claiming that an unmapped graph produced the sequence.
- Rename the dashboard's filtered, branch-specific file list to `Result files`.
  Its count excludes hidden duplicate/intermediate FASTAs and is not the total
  number of files retained for the sample.

# 2.4.1 - 2026-11-19

- Recover observed COI reference names from existing BOLD lineages when the
  terminal taxon is the literal `None`, `null` or another missing-value marker.
  Apply the same normalization during database preparation and validation.
- Label dashboard taxa as reference matches or shared reference lineages,
  separately from confirmation of the sample's expected identity. Existing
  alignments and metadata are sufficient to rebuild the display.

# 2.4.0 - 2026-11-19

- Fix MITOS coding-gene display when RNA features coexist; support CDS parents,
  split features, circular-origin coordinates and stranded sequence extraction.
  Add gene/region/node sequence popovers with copy and FASTA download.
- Separate canonical 13 CDS/two rRNA counts from additional reference genes and
  distinguish annotation from homology screening. Show the expected hierarchy.
- Follow reported ITS/mitogenome source paths, reduce graphs to contributing
  nodes, highlight shared loci, label taxonomic assignment scope and reconcile
  sequence lengths with overlaps. Preserve reduced graphs/provenance before
  compact cleanup. Hide duplicate FASTA aliases and unrelated graph inventories.
- Support 150-bp COI alignments/reference preparation by default, retain explicit
  threshold overrides and record thresholds in BOLD evidence. Add conservative
  match-confidence grades without substituting expected names for observations.

# Earlier updates

- Add a database installation guide with download/preparation commands for ITSME,
  GetOrganelle animal_mt, MITOS refseq89m (including its checksum), annotated NCBI
  mitochondrial GenBank records and BOLD COI/taxonomy. Link it from installation
  instructions and update Belly paths and hierarchy defaults.

- Restore observed COI reference names when older evidence lacks scores or the
  confirmation filter retained no hits. Read saved alignments and selectively
  resolve names from the configured BOLD metadata in one streaming pass; add a
  --bold-metadata override. Embed names and match metrics without rerunning BOLD
  or changing confirmation/assembly acceptance.
- Remove technical warning paragraphs, raw stage reasons and annotation-policy
  commentary from all dashboard views and expandable panels. Deduplicate feature
  warnings in the saved display data; retain native reports and failure-log links.

- Simplify dashboard captions and metric notes. Keep technical target-check
  reasons and exploratory annotation scope in downloadable evidence files.
  Distinguish observed top COI reference names from confirmed target identity;
  do not repeat the expected name as an unresolved prediction.
- Link failed/timed-out Run record statuses directly to retained native logs.
  Preserve relative downloads after relocating results; leave missing logs
  unlinked and refuse paths outside the sample.

# 2.3.1 - 2026-11-19

- Create MITOS's output directory before invoking runmitos.py. Fix the native
  argument-validation failure "no such directory .../output" in ordinary and
  existing-candidate annotation, with circular and unclosed inputs.
- Make the MITOS CLI double require a pre-existing output directory instead of
  silently creating one, and add a regression case for each topology. Confirm
  the regressions fail before the fix. Document annotation-only retry without
  repeating successful BOLD searches or requiring --force.

# 2.3.0 - 2026-11-18

- Annotate the best retained mitochondrial candidate by default, independently of
  circularity, gene completeness and target acceptance. Use MITOS linear mode
  for unclosed candidates and label exploratory annotations. Preserve all
  original acceptance criteria; do not guess unknown genetic codes.
- Replace unconditional BOLD hit-limit rejection with near-best score-boundary
  checks and adaptive tenfold expansion from 500 up to 50,000 subjects. Keep
  genuinely truncated/tied evidence at the hard cap unresolved, and record
  search attempts, scores and decision reasons.
- Add annotate_existing.sh to annotate saved candidates without reassembly or
  ITSME. Optional --recheck-bold repeats COI validation under the original policy.
  Cache unchanged annotations, retain prior summaries/stages, refresh output
  inventories/cohort/dashboard and remove only owned postprocessing scratch.
  BOLD search failure does not prevent annotation.
- Separate ITSME execution status from complete/partial locus recovery counts.
  Display branch-specific ITSME badges and available diagnostic reasons, so
  successful zero-complete-locus runs are not presented as execution failures.

# 2.2.0 - 2026-11-17

- Add an offline dashboard.html that travels with the result directory, with a
  sample list and focused Mitogenome, ITS/rDNA and Run record views. Embed all
  display data/styles/scripts; retain relative links to native deliverables.
- Display mitochondrial sequence evidence, MITOS gene/annotation maps and tables,
  ITSME loci/lengths/consensus taxonomy, read QC, stage outcomes and cohort length
  comparisons. Keep detailed tables collapsed and separate target expectations
  from predicted taxonomy and COI target validation.
- Render bounded real GFA1/FASTG graph previews with pan/zoom, oriented links,
  segment metadata and original downloads. Label WGS/method graph scope, missing
  graphs, preview limits and unsupported GFA2 explicitly; never fabricate graphs.
- Add build_dashboard.sh for existing inactive runs. Automatically regenerate
  after report refresh, including new completed runs and compaction. Preserve
  active-run guards, scientific outcomes, input files and heavy work files.
- Test relocation/offline embedding, graph parsing, annotation coordinates,
  absent/empty/stale results, lock protection and DOM interactions. No new
  production environment dependencies.

# 2.1.1 - 2026-11-16

- Include reported ITSME locus counts, ITS cassette and full locus lengths, and
  consensus taxa/ranks in the cohort summary for combined and ITSME-only runs.
  Preserve multiple loci as aligned arrays rather than choosing one taxon.
- Add cohort_its_loci.csv with one row per complete reported locus, separate
  ITS1/5.8S/ITS2 lengths, full consensus lineage and taxonomy status. Missing
  results are distinguished from zero recovery; partials/chimeras are excluded.
- Add refresh_reports.sh for existing completed output. It reads saved native or
  compact ITSME evidence and refreshes reports without analysis or cleanup,
  respecting the active-run lock and skipping stale/unfinished summaries.
- Retain ITS reporting during output compaction, migration and compact resume.

# 2.1.0 - 2026-11-14

- Standardize every processed sample with ITSME/ and MitoLooper/ branch folders
  and explicit statuses, even when a branch is disabled, fails or has no reference.
  Keep shared provenance/QC/logs and support migration of older compact flat output.
- Default automatic reference selection to annotated hierarchy matching through
  kingdom. Stop at the closest eligible supplied rank, preserve length/code/gene
  checks and report the matched rank. Support Animalia/Metazoa kingdom aliases.
  Remove BLAST reference selection and read-prefix/contig selection searches.
- Recover missing-reference samples with GetOrganelle's installed animal_mt
  defaults; continue ITSME/WGS normally without diverting sample directories.
- If no annotation reference exists, permit direct candidate COI/BOLD validation
  while keeping gene completeness unknown. Never guess an unprovided genetic
  code; export accepted sequence and report unavailable annotation explicitly.
- Keep compact cleanup, active-run protection and input preservation. Retain
  reference provenance and legacy search cleanup for older completed outputs.

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


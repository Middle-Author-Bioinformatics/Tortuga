# Validation - 2026-11-12

Controller version: 2.1.0.

- 113 tests passed using Python 3.12, Biopython and pytest, including automatic
  environment selection and per-stage switching with simulated Conda/tool entry points.
- Covered sequential method order, immediate MITOS2 annotation after read recovery,
  optional ITSME before WGS, post-WGS graph/contig fallback, early stopping,
  BOLD taxonomic conflicts and hit limits, target-confirmed rescue seeds,
  annotation failures, unchanged-run caching and changed-option rejection,
  fresh retries for failed tools, deterministic paired subsetting,
  circularity/junction evidence, dry runs, BOLD export normalization and
  killing timed-out process groups including owned children.
- Python compilation, Bash launcher syntax, and CLI help/version checks passed.
- Setup-helper tests cover absolute Conda prefix mapping, fallback environment
  selection, NOVOPlasty discovery from package metadata and bounded reference scans.
- Automatic-routing tests cover launch without manual activation, first-choice
  selection with alternative notes, command execution under different stage prefixes,
  pipeline startup without a config file and reporting missing executables.
- Controller tests cover all three selectors, exactly one ITSME call in combined
  mode, no mitochondrial prerequisites for standalone ITSME, standalone retries and
  caching, failure continuation across samples, selector validation and the shell entry point.
- Original uploaded scripts are preserved byte-for-byte under `legacy/`.
- Automatic reference tests cover exact species precedence, reported same-genus
  fallback, strict matching, named higher lineages, unrelated-taxon rejection,
  organelle/annotation/partial/length filtering, conflicting and missing codes,
  deterministic ties, changed/removed-file cache refresh and a taxon-only sample
  sheet end-to-end through simulated recovery with selection reports.
- BOLD preparation tests include a 200 KB TSV field exceeding Python's old limit.
- Hierarchy tests cover normalized source rank columns/species epithets,
  same-family fallback, closer-genus precedence, explicit order-scope limits,
  phylum disagreement, family-level code conflicts, complete metadata retention
  and successful target confirmation when the seed reference is a different genus.

- Real-record regression: NCBI NC_007933.1 (Orbinia latreillii) supplies 13 CDS,
  two rRNAs and a 1,563 bp unambiguous COX1 seed, with COX1/COX3 boundary warnings.
- New safety cases retain ambiguous-base rejection, require annotated translation
  for uncertain CDS boundaries, reject unsupported uncertain coordinates, prefer
  exact-boundary references within a tier and refresh old validation caches.
- Cohort tests verify continued recovery after reference failure, independent
  ITSME in combined mode, honest failure reporting, no BOLD index when no references
  are usable and unchanged independent-ITSME resume caching.

- Policy tests cover default strict acceptance, linear all-gene recovery, gene
  thresholds, exact length boundaries, OR/AND truth tables, per-contig unique
  counts, non-target/conflicting-target rejection, ambiguous bases, length bounds,
  failed validation, invalid options and annotation failures under relaxed rules.
- Integration cases verify early stopping in read and post-WGS recovery,
  optional WGS skipping with exactly one independent ITSME call, strict versus
  relaxed output labels/files, MITOS2 linear mode, unclosed end preservation
  and single versus doubled homology subjects.

- Compaction tests cover export of native ITSME products, final graphs, WGS
  contigs and MITOS annotations; deletion of generated reads/mappings/K* scratch
  and copied BOLD indexes; preservation of QC, failure logs, provenance and inputs;
  dry-run nonmutation, repeat cleanup, conflicting filenames, missing retained
  outputs, active locks, linked metadata and external symlink protection.
- Combined and standalone integrations verify default compaction, native ITSME
  results at parent final, failed-stage status retention, and resume without new
  tool calls or a rebuilt BOLD database. Older reference-unavailable batches are
  migrated into the same parent final layout without changing biological status.
- GFF-to-TSV tests retain original coordinates, strand/type and decoded names.

- Partial-batch cleanup tests verify that missing batch records still allow
  finished samples to be compacted while unfinished/planned/stale samples and
  shared indexes remain untouched. Active locks block cleanup even without a
  batch record; unstarted manifest entries prevent shared-index removal.
- Independent finished ITSME results can be promoted without fabricating a batch
  completion record; missing paths receive a separate diagnostic.

The integration tests use CLI doubles and synthetic references, not the actual
external assembly programs. No project FASTQs were supplied here, and real
NOVOPlasty/GetOrganelle/MitoFinder/mtGrasp/MITOS2 jobs were not run. Reference
compatibility, recovery sensitivity and runtime require a trial on TheBelly.
The README gives the necessary tool probes and a single-sample trial command.


Version 2.1.0 additionally tests default hierarchical fallback through kingdom,
first matching rank preference, explicit fallback limits, Animalia/Metazoa
aliases, phylum boundaries and genetic-code conflicts. The CLI rejects removed
BLAST reference-selection options. Integrations cover no-reference GetOrganelle
defaults, direct COI validation with unknown completeness, genetic-code-dependent
annotation, unavailable-tool WGS fallback, both branch directories/statuses,
early QC failure and migration of preexisting compact output without reassembly.
Native assembler jobs and real project FASTQs were not run; the integration
programs remain CLI doubles.

Version 2.1.1 tests ITS cassette versus full-locus length, inclusive/reverse
coordinates, matching extracted FASTA fallback, multiple loci with distinct taxa,
reconstructed graph-source taxonomy mapping, unknown ranks and missing evidence,
partial/chimera exclusion, zero versus unavailable counts, retained mitochondrial
columns and files, active locks and stale samples. Integration coverage checks
both standalone and combined modes with and without compact cleanup.

Version 2.2.0 adds offline dashboard tests for embedded assets and moved-folder
links, preserved input files, GFA1/FASTG orientation/self loops and bounded previews,
GFA2 download-only behavior, GFF bounds/multiple contigs, missing/empty/stale
results, HTML/script-like text, external links and active-run locks. Report refresh
and existing pipeline integrations cover automatic generation. A jsdom-only test
checks filtering, sample/branch/locus changes, features, graphs, zoom controls,
comparison and missing-result states. Local Chromium could not launch, and the
cloud browser blocked file:// previews; visual browser rendering was not verified.

Version 2.3.0 tests default annotation of unclosed/unresolved candidates without
promotion, explicit accepted-only/disabled policies, missing-code handling,
clear BOLD score boundaries versus truly truncated ties, adaptive expansion and
hard caps, saved-output annotation/revalidation without assembler/ITSME calls,
annotation caching, active locks and continued annotation after BOLD failure.
ITSME tests distinguish successful zero-complete-locus runs, partial counts and
actual execution timeouts. Native MITOS2 jobs on real Belly data were not run.

Version 2.3.1 reproduces the reported native MITOS precondition by making the
CLI double reject missing -o directories. Circular and unclosed regression cases
fail before the patch and pass once the shared launcher creates the directory.
The existing saved-output tests exercise the same shared launcher and retry path.

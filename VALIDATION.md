# Validation - 2026-11-10

Controller version: 2.0.4.

- 46 tests passed using Python 3.12, Biopython and pytest, including automatic
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

The integration tests use CLI doubles and synthetic references, not the actual
external assembly programs. No project FASTQs were supplied here, and real
NOVOPlasty/GetOrganelle/MitoFinder/mtGrasp/MITOS2 jobs were not run. Reference
compatibility, recovery sensitivity and runtime require a trial on TheBelly.
The README gives the necessary tool probes and a single-sample trial command.

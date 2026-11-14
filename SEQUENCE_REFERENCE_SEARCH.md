# Sequence-based reference choice and stable branch results (2.1.0)

Every processed sample now has `ITSME/` and `MitoLooper/` directories. Each has
`status.json`, including when its branch was not requested, could not run or
failed. Native ITSME final sequences/tables/graphs are exported into `ITSME/` with
native filenames. Mitochondrial candidates, annotations and final assembly graphs
are in `MitoLooper/`. Shared QC/logs/provenance remain beside those folders;
`final/summary.json` and `final/outputs_manifest.tsv` remain the common metadata
entry points. Historical flat final products are migrated with compatibility
links when compacting old completed runs. Retained-output manifests point to real
files, not those compatibility links.

The native ITSME command/API is unchanged. These directories are the controller's
standard result layout, not extra execution of either branch. `--pipeline itsme`
still runs ITSME only; `MitoLooper/status.json` says `not_requested`. In mito-only
runs the ITSME status says `not_requested`. No-reference samples use the same
sample tree as every other sample, rather than being diverted into a separate
ITSME-only batch.

## Reference selection

`--reference-match blast` is now the default. If an old command explicitly says
`--reference-match hierarchy`, change that value to `blast` to enable the new
selector. Explicit `reference_gb` sample-sheet overrides still take precedence.
The old `exact`, `genus` and `hierarchy` modes remain available intentionally.

The sequence-based path:

1. Validate and incrementally catalog the local annotated animal mitochondrial
   GenBank records. Annotation, genome-length and explicit CDS translation-table
   requirements still apply; BLAST does not make a bad reference valid.
2. Sample a bounded deterministic prefix of paired reads: at most 10,000 pairs
   per sample by default. This is a bounded read prefix, not uniform sampling of
   the entire sequencing library. Reads remain on the server.
3. BLAST the queries against an indexed FASTA of eligible local mitochondrial
   genomes. No remote NCBI BLAST or silent database downloads are performed.
4. Require at least 80% identity, at least 80 aligned bases per HSP and support
   from at least three read templates. Mates count as one template. Overlapping
   HSPs for the same query/reference cannot inflate the summed score.
5. Apply a provided phylum and explicit mitochondrial code as eligibility filters.
   Within that eligible pool, rank by summed nonoverlapping-HSP bit score, then
   aligned-length-weighted identity, supporting template count and deterministic
   accession/path ties. Names do not determine the winning rank in BLAST mode.
   If near-best scores (within 1%) imply different genetic codes, require an
   explicit code or continue to bootstrap recovery; do not guess a code.

Options: `--reference-read-pairs`, `--reference-min-identity`,
`--reference-min-aligned-bp`, `--reference-min-templates`. The BLAST output retains
at most 100 target sequences per query. The result is the best sequence-supported
eligible local assembly reference under those settings, not proof of the true
phylogenetically closest relative. Large, mixed or divergent samples can have
insufficient or misleading read-prefix evidence; target validation remains
independent. Full input taxonomy, especially phylum, should be retained in the
sample sheet. Without a provided phylum, there is no inferred cross-phylum filter.

Reports retain `reference_selection.tsv`, the sampled-read ranking and selection
JSON in `reference_search/SAMPLE/`, the selected accession/path/code, scores,
alignment lengths, identity and template support. `source.json` fingerprints the
reference collection. The top-level table records the initial read-based choice;
if a reference is first resolved after assembly, that choice is recorded in the
sample's `references/sequence_selection.json`, candidate search evidence and final
summary. Reference selection does not relabel the expected target organism.

## No suitable custom reference

The sample is still processed. Reference-dependent NOVOPlasty, MitoFinder and
mtGrasp attempts are skipped until a usable annotated reference exists.
GetOrganelle is the first configured recovery tool that can run without a custom
reference. Its `-s` option is omitted so it uses its installed `animal_mt` seed
and label databases. If unavailable/disabled/unsuccessful, WGS remains enabled
unless explicitly disabled. An unavailable WGS assembler can fall through from
metaSPAdes to installed MEGAHIT. Graph-based GetOrganelle can likewise use its
installed labels without a custom `--genes` file.

GetOrganelle is de novo assembly using database-guided recruitment/extraction;
this mode is custom-reference-free, not database-free. Its maintainers document
the default seeds and database setup here:
https://github.com/Kinggerm/GetOrganelle

If the existing GetOrganelle environment lacks its animal mitochondrial database,
prepare it on TheBelly using the installed environment:

```bash
/home/ark/miniconda3/bin/conda run -p /home/ark/miniconda3/envs/getorganelle \
  get_organelle_config.py -a animal_mt
```

After a bootstrap contig is recovered, BLAST it against the same local reference
collection. Contig selection requires at least 400 aligned query bases and uses
nonoverlapping HSP support; it does not require three read templates. A successful
choice supplies the reference genes and translation table for normal assessment
and subsequent reference-dependent recovery attempts.

If no annotated local reference can be selected even then, compare the candidate
DNA directly against BOLD COI references. This can validate the target identity,
but it does not establish recovery of all mitochondrial genes. Gene completeness
remains unknown (`gene_screen_available=false`, `genes_expected=null`) and cannot
satisfy `all-genes` or a positive gene-count criterion. A user-selected length
policy can still accept a target-confirmed candidate under its other gates.
MITOS2 can annotate it if `mt_code` is explicitly provided; otherwise export the
accepted FASTA, skip annotation with an explanation and report
`accepted_by_policy_annotation_unavailable` with a nonzero batch status.

ITSME still runs in combined mode regardless of reference availability. Skipped,
failed and unavailable stages remain visible in stage/status reports. Biological
non-recovery is distinct from a program execution error.

## Deployment and existing runs

Merge the PR and update the existing checkout:

```bash
git -C /home/ark/MAB/bin/Tortuga pull --ff-only
```

Use a fresh output directory for new 2.1.0 analyses; older checkpoints cannot be
silently reused under changed selection rules. Existing completed results can be
reorganized without reassembly with the usual cleanup command:

```bash
bash /home/ark/MAB/bin/Tortuga/compact_outputs.sh \
  /home/ark/MAB/NOAA/CA_intertidal/both_results_v205
```

The active-controller lock guard still applies. Cleanup skips unfinished samples
and retains shared indexes when whole-batch completion is uncertain. Default
compaction removes generated BLAST indices, query prefixes and raw search hits
only after the cohort finishes, preserving rankings and selection provenance.
`--keep-intermediates` retains work files for debugging. Successful unchanged
reference searches are cached; execution failures remain retryable.

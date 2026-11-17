# A report that travels with the results

Tortuga writes **dashboard.html** in the cohort output directory. Open it directly
in a modern browser. All display data, CSS, JavaScript, gene maps and bounded graph
previews are embedded: no server, account, CDN, internet access or browser plugins
are needed. Original sequence, annotation, graph and table downloads use relative
links. Move or share the **whole output directory** to keep those links working.
No genomic sequences or raw reads are embedded in the HTML.

The first version uses a quiet sample list and one workspace, rather than a grid
of metric cards. It has three views:

- **Mitogenome:** accepted assembly versus best retained candidate, length,
  topology, reference-gene recovery, recovery method, MITOS feature map and table,
  COI/BOLD target evidence, and final sequence/annotation downloads.
- **ITS / rDNA:** a locus selector, full locus and ITS cassette lengths,
  SSU–ITS1–5.8S–ITS2–LSU region map, predicted consensus taxon/rank and lineage,
  complete-locus table and sequence downloads. Multiple loci stay separate.
- **Run record:** available fastp read/base-quality metrics, individual stage
  outcomes and durations, retained logs, and the saved sample summary. Timeouts
  remain visible even if another method recovered a sequence.

A display switch in each assembly view opens the retained assembly graph with
pan/zoom, file selection, segment lengths/coverage and original graph downloads.
**Compare sample lengths** opens an on-demand cohort comparison, using separate
mitochondrial and ITS scales. Feature tables, file inventories and COI records
are collapsed initially to keep the main view simple. Keyboard navigation,
search by sample/taxon, narrow-screen CSS and print styles are included.

## Build for your existing Belly run

Merge the dashboard PR and update Tortuga:

```bash
git -C /home/ark/MAB/bin/Tortuga pull --ff-only

bash /home/ark/MAB/bin/Tortuga/build_dashboard.sh \
  /home/ark/MAB/NOAA/CA_intertidal/both_results_v210_hierarchy
```

This creates:

```text
/home/ark/MAB/NOAA/CA_intertidal/both_results_v210_hierarchy/dashboard.html
```

Copy the result directory to your computer and double-click **dashboard.html**.
An optional `--title "California intertidal"` changes the displayed collection
name. Building the report reads saved outputs and replaces only dashboard.html;
it does not run assemblers, annotations, classification, or cleanup. It refuses
an active controller. If the report is a stale snapshot after an interrupted
sample retry, old recovery claims are hidden and a warning is shown.

New completed runs generate the dashboard automatically when cohort reports are
refreshed. The existing `refresh_reports.sh` and completed-output compaction also
rebuild it. Report-generation errors warn without turning a completed assembly
into a failed analysis; the standalone builder reports errors directly. No new
production dependencies are required beyond the existing Python environment.

## How to interpret the views

- An **accepted assembly** met the selected run policy. A **best candidate** may
  be partial, unclosed or unresolved. Read acceptance/circularity/gene evidence
  separately; reference-gene recovery does not establish every tRNA or exclude
  all NUMTs/chimeras. Policy acceptance is not a universal completeness claim.
- ITS predicted taxonomy comes from native ITSME consensus reports. Unknown
  taxa/ranks stay unknown. Expected target names are labeled as expectations.
- Mitochondrial taxonomy remains the existing COI/BOLD **target check**, with
  supporting reference records when retained. It does not independently infer
  the deepest possible taxonomy, and the assembly reference is not a taxonomic
  prediction. This dashboard does not change scientific inference.
- A gene map is a coordinate display of retained MITOS GFF annotations, shown
  unrolled even for circular assemblies. It is separate from the assembly graph.
  Feature records are not assumed to be distinct genes. Multiple annotation
  sequence IDs disable the overlaid map; the table and GFF remain available.
- Graphs are the actual retained GFA1/FASTG method graphs; no graph is synthesized
  from a FASTA or inferred from circularity. A method graph can contain alternate
  and off-target paths. WGS graphs are explicitly labeled as whole-genome graphs,
  not as the selected mitochondrial path. Node positions/widths are schematic;
  segment metadata gives actual recorded lengths. Link orientations are retained.
- For portability and responsiveness, each graph preview is limited to the first
  240 segments and up to 720 stored oriented links. Files over 100 MiB are offered
  as downloads instead. GFA2 is downloadable but not rendered in this first
  version. Truncation and unsupported formats are labeled explicitly. Full graph
  files are never altered. Annotation previews cap at 3,000 feature records.
- Missing graph, QC, annotation, sequence or taxonomy evidence is shown as
  unavailable, not as zero recovery or as an invented result. Failed annotation
  outputs are labeled for review. Downloads stay available where retained.

## Validation

Backend tests cover relative-link relocation, input preservation, current and
legacy compact products, graph orientations/self loops/caps, GFA2 refusal,
annotation bounds/multiple contigs, absent outputs, empty FASTAs, stale summaries,
active locks, report refresh and safe JSON embedding. Existing pipeline tests
exercise automatic generation in standalone/combined and compact/retained modes.

A DOM-only test exercises sample search/selection, branch/locus switches, gene
selection, GFA/FASTG displays, zoom buttons, comparison, timeouts and missing
results, without network access. `jsdom` is an optional development test dependency,
not required on TheBelly. For visual review, generate illustrative data with
`tests/dashboard_fixture.py` and open the resulting report in a real browser.
This environment blocked local browser rendering; the initial release was checked
through data and DOM tests, not a completed visual browser review.

# A report that travels with the results

Tortuga writes **dashboard.html** in the cohort output directory. Open it directly
in a modern browser. All display data, CSS, JavaScript, gene maps and bounded graph
previews are embedded: no server, account, CDN, internet access or browser plugins
are needed. Original sequence, annotation, graph and table downloads use relative
links. Move or share the **whole output directory** to keep those links working.
Selected gene/ITS-region sequences and reduced graph-node sequences are embedded
for copying and FASTA download. Raw reads and whole-genome assemblies are not embedded.

The `Result files` count belongs to the currently selected branch and includes
only its displayed deliverables. Duplicate FASTA aliases, intermediates and logs
are not counted there. Four mitogenome files can be the complete displayed set:
the selected FASTA, gene-recovery CSV, MITOS GFF and feature TSV. Different counts
do not indicate that fewer assembly methods ran.

In the assembly-graph view, `Source graphs` provides downloads of available native
graphs independently of selected-path mapping. `Selected path unavailable` means
source graphs exist but the report could not establish the displayed sequence's
path. `Source graph unavailable` means no applicable source graph was found by
the report inventory. Either can occur without an assembly failure. The report
does not infer a replacement graph from FASTA or substitute a different method's
graph. Deleted graphs cannot be restored by rebuilding HTML.

COI reference names are displayed independently of the identity cutoff used to
confirm the expected target. Version 2.4.1 recovers names from existing lineages
when a prepared BOLD record ends with a missing-value marker such as `None`.
Rebuild the dashboard with the existing metadata to repair the display; no database
download, preparation, BLAST, assembly or annotation rerun is needed for this fix.
Reference matches and shared reference lineages are observed database evidence,
not automatic confirmation of the expected species.

To change the validation cutoff on saved candidates, use `annotate_existing.sh`
with `--recheck-bold --coi-min-bp 150 --coi-identity 95`. This repeats COI validation
and refreshes reports while reusing unchanged MITOS annotations. It does not rerun
assembly or ITSME. The default identity cutoff remains 97%; a lower cutoff can
admit more distant reference matches and does not force target confirmation.

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
  remain visible even if another method recovered a sequence. Click a failed or
  timed-out status to download its native log; a status stays plain text when
  no log was retained. Log links survive moving the entire result directory.

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
name. Building the report reads saved outputs, writes selected graph/provenance
assets and replaces dashboard.html. It does not run assemblers, annotations,
classification or cleanup. It refuses
an active controller. If the report is a stale snapshot after an interrupted
sample retry, old recovery claims are hidden and a warning is shown.

New completed runs generate the dashboard automatically when cohort reports are
refreshed. The existing `refresh_reports.sh` and completed-output compaction also
rebuild it. Report-generation errors warn without turning a completed assembly
into a failed analysis; the standalone builder reports errors directly. No new
production dependencies are required beyond the existing Python environment.

## Updated assessment labels (2.3.0)

The MITOS map can now represent an exploratory annotation of an unclosed or
unresolved candidate. Its scope is saved in the evidence files; annotation does
not establish target identity or acceptance. BOLD reasons are retained in the saved evidence and
cohort summary. In the ITS view, the badge reports **ITSME execution**, independent
of the overall mitochondrial outcome; complete/partial locus counts remain
separate. A successfully completed ITSME run with no complete loci is not a
program failure.

The client-facing report shows measurements, taxonomy, maps and downloads without
developer commentary or warning paragraphs, including inside expandable panels.
Technical reasons and annotation scope remain in the downloadable native reports.
**Taxonomic evidence** contains a concise target agreement result and reference
records. Failed/timed-out stages retain their log download links.

**COI reference match** uses the highest-scoring saved alignment names, including
ties within the bounded 12-record preview, independently of whether the target
confirmation filter accepted them. Percent identity and aligned length accompany
the name. Known names in older evidence without scores remain visible as reference
matches, without claiming a ranking. Confirmed checks show **COI-supported identity**.
Expected sample names are never substituted for absent observed reference names.

If confirmation saved no qualifying hits, the builder reads `bold_hits.tsv` and
resolves the selected IDs from the BOLD metadata path saved in `run_config.json`.
It scans each metadata TSV once per report and stores only the needed records;
it does not load the full reference database into memory or rerun BLAST. An explicit
path can be supplied when the saved path is unavailable:

```bash
bash /home/ark/MAB/bin/Tortuga/build_dashboard.sh \
  /home/ark/MAB/NOAA/CA_intertidal/both_results_v210_hierarchy \
  --bold-metadata /home/ark/MAB/NOAA/BOLD/prepared_fixed/metadata.tsv
```

Names and match metrics are embedded in `dashboard.html` for offline sharing.
Building the dashboard preserves the saved confirmation and acceptance decisions.

## How to interpret the views

### Graph-to-sequence view (2.4.0)

The graph view follows the older ITSME graph-to-sequence dashboard: choose a locus,
highlight its path, distinguish shared nodes, inspect taxonomy and copy sequences.
Only graphs associated with reported loci or the chosen mitogenome are offered.
Native `locus_source_map.tsv`, `contigs.paths`, original manifest paths and exact
sequence reconstruction distinguish nodes from different source assemblies.
Collapsed alternative paths are not treated as the chosen representative.

The path table gives segment orientation, length, incoming overlap and added
bases. Its length equation appears only when the reconstructed path contains the
final FASTA sequence (or its reverse complement), accounting for overlaps and
trimmed ends. Discontinuous paths retain their separate pieces; unknown overlaps
are not assumed to be zero. A FASTG overlap can be inferred only when recorded
links and the inferred common overlap reconstruct the final sequence exactly.
Short segments required by a selected path are kept; unrelated short graphs are
excluded. Distinct reported ITS taxa are not collapsed into one arbitrary locus.

New compact runs preserve `ITSME/selected_graphs/` and
`MitoLooper/selected_graphs/` before stage cleanup. These contain reduced GFA
graphs and an index of path/sequence evidence. Full source graph exports may be
removed once every reported sequence in that branch has a preserved subgraph;
`--keep-intermediates` preserves native work. Old missing source/path files cannot
always be recreated from FASTA. No synthetic assembly graph is drawn as a substitute.

Node labels/popovers show observed locus-consensus or assembly-COI associations,
with the assignment scope identified. They are not independent BLAST classifications
of every node. Nodes with no evidence remain unassigned; expected taxonomy is never
projected onto graph nodes. Gene, ITS-region and node popovers provide Copy,
FASTA download and an NCBI BLAST link; clicking the link opens the BLAST form and
does not automatically submit sample sequences.

### Gene counts and COI support

MITOS `gene` entries for known coding genes are displayed as CDS even when RNA
features coexist. Explicit CDS children take precedence over parent gene records;
multipart CDS parts are joined and minus strands reverse-complemented. Circular
features extending over the origin are split on the map and extracted across
both ends. Invalid/out-of-range features remain excluded without on-screen
developer warnings; native annotations remain downloadable.

Core counts show the canonical 13 CDS and two rRNAs separately, with additional
genes listed explicitly. The label distinguishes MITOS annotation counts from
reference-screen counts. Reference-specific recovery remains in its own table;
a sponge reference with `atp9` can legitimately have 16 CDS+rRNA genes total.

The expected identity includes the supplied ranked hierarchy. Observed COI names
come only from saved alignments/metadata; even a stored confirmed-target flag
does not substitute the expected name when retained reference evidence is absent.
The evidence grade is conservative and uncalibrated: High requires at least 99%
identity, 400 aligned bp, 90% shorter-sequence coverage and complete near-best
metadata/agreement. Moderate requires at least 97%, 150 bp and 80% coverage.
Mixed indicates different near-best reference taxa; Limited indicates incomplete
metadata, score-window truncation or unavailable scores. No matches is Unavailable;
other matches are Low. The deepest common retained lineage is shown only when
the evaluated near-best window is complete. These grades are not probabilities
or replacements for the pipeline's explicit target-confirmation policy.

Download inventories collapse FASTA files with identical sequence content even
when headers differ, and show the selected graphs in the graph workspace instead
of listing every intermediate graph. Report rebuilding does not delete existing
native FASTAs.

- An **accepted assembly** met the selected run policy. A **best candidate** may
  be partial, unclosed or unresolved. Read acceptance/circularity/gene evidence
  separately; reference-gene recovery does not establish every tRNA or exclude
  all NUMTs/chimeras. Policy acceptance is not a universal completeness claim.
- ITS predicted taxonomy comes from native ITSME consensus reports. Unknown
  taxa/ranks stay unknown. Expected target names are labeled as expectations.
- Mitochondrial taxonomy combines displayed observed COI reference records with
  the conservative common-lineage summary above. The saved COI/BOLD target check
  remains separate. Assembly reference names are not taxonomic predictions.
- A gene map is a coordinate display of retained MITOS GFF annotations, shown
  unrolled even for circular assemblies. It is separate from the assembly graph.
  Feature records are not assumed to be distinct genes. Multiple annotation
  sequence IDs disable the overlaid map; the table and GFF remain available.
- Graphs are the actual retained GFA1/FASTG method graphs; no graph is synthesized
  from a FASTA or inferred from circularity. A method graph can contain alternate
  and off-target paths. WGS graphs are explicitly labeled as whole-genome graphs,
  not as the selected mitochondrial path. Node positions/widths are schematic;
  segment metadata gives actual recorded lengths. Link orientations are retained.
- Source-graph reading is limited to 100 MiB per file. Only selected-path nodes
  are embedded; paths are not truncated to the first arbitrary nodes of a source
  graph. GFA1 and FASTG are supported; GFA2 is not rendered. Annotation parsing
  caps at 3,000 source feature records. Existing source files are not altered by
  report rebuilding.
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

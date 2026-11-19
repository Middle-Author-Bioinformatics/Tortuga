# Download and prepare the databases

Run these setup steps once, after installing the software. Reuse databases that
are already prepared; the existing Belly paths below can be used directly.
Commands use full paths for TheBelly. On another host, substitute your own paths
and Conda environment prefixes.

The download commands use `curl`, `wget`, `tar`, `bzip2` support and `md5sum`.
The GenBank preparation step uses Biopython from the controller environment.

## What each pipeline needs

| Resource | Used by | Controller argument |
| --- | --- | --- |
| ITSME SILVA seeds, NCBI rRNA/ITS BLAST databases and taxonomy | `itsme` and `both` | `--itsme-db /home/ark/databases/itsme_db` |
| Annotated mitochondrial GenBank records | Reference selection, NOVOPlasty seeds, MitoFinder, mtGrasp and reference-gene screening | `--genbank-dir /home/ark/MAB/NOAA/mitos_genbank` |
| GetOrganelle `animal_mt` seed and label databases | GetOrganelle, including recovery without a suitable custom reference | Installed with GetOrganelle's configuration command |
| MITOS reference models | Mitochondrial annotation | `--mitos-refdir /home/ark/MAB/NOAA/mitos_data --mitos-refseqver refseq89m` |
| BOLD COI FASTA and matching taxonomy TSV | COI target checking and reference-match reporting | `--bold-fasta .../coi.fasta --bold-metadata .../metadata.tsv` |

`--pipeline itsme` needs only the ITSME databases. `mitolooper` needs the
mitochondrial resources; `both` needs both sets. Without an eligible custom
GenBank reference, Tortuga can start with GetOrganelle's installed `animal_mt`
resources. MITOS still requires its models and a known mitochondrial genetic code.

The mitochondrial GenBank collection, MITOS models and BOLD library are separate
resources. Nucleotide/protein FASTAs alone do not replace annotated GenBank files.
The `itsme2_db` and `itsme2_silva138_2` layouts are not replacements for the
`--itsme-db` layout used by the ITSME shell pipeline.

## 1. ITSME databases

Use the setup helper in the separately installed
[ITSME repository](https://github.com/Middle-Author-Bioinformatics/ITSME).
If that repository is not installed yet:

```bash
git clone https://github.com/Middle-Author-Bioinformatics/ITSME.git \
  /home/ark/MAB/bin/ITSME
```

Download and prepare the databases from the environment containing BLAST+:

```bash
/home/ark/miniconda3/bin/conda run --no-capture-output \
  -p /home/ark/miniconda3/envs/itsme \
  bash /home/ark/MAB/bin/ITSME/setup_db.sh \
  --db-dir /home/ark/databases/itsme_db \
  --threads 8
```

The helper downloads the SILVA eukaryotic 18S/28S seed FASTAs, NCBI
`SSU_eukaryote_rRNA`, `ITS_eukaryote_sequences`, `LSU_eukaryote_rRNA`, `taxdb`,
and the NCBI taxonomy dump. It extracts the BLAST databases and taxonomy files.
BLAST+'s `update_blastdb` must be available in this environment. Existing seed
files are reused unless `--force` is requested; the NCBI downloader checks its
database files for updates.

Expected layout:

```text
/home/ark/databases/itsme_db/
  silva-euk-18s-id95.fasta
  silva-euk-28s-id98.fasta
  NCBI_rRNA_BLAST/
    SSU_eukaryote_rRNA.*
    ITS_eukaryote_sequences.*
    LSU_eukaryote_rRNA.*
    taxdb.*
  taxdump/
    nodes.dmp
    names.dmp
    merged.dmp
```

No separate manual `makeblastdb` step is needed for these downloaded databases.

## 2. GetOrganelle animal mitochondrial databases

Install the seed and label databases in the GetOrganelle environment Tortuga will
use. On TheBelly, `getorganelle` is checked before `mito-getorganelle`:

```bash
/home/ark/miniconda3/bin/conda run --no-capture-output \
  -p /home/ark/miniconda3/envs/getorganelle \
  get_organelle_config.py --add animal_mt
```

Run setup as the same user that runs Tortuga. GetOrganelle uses its per-user
configuration directory by default; a configured `GETORG_PATH` changes that
location and must be the same during setup and execution. GetOrganelle's helper
downloads and initializes its resources; Tortuga does not build these itself.
For download problems, use the upstream
[local-file initialization instructions](https://github.com/Kinggerm/GetOrganelleDB).

## 3. MITOS annotation models

For this animal mitochondrial workflow, download the pinned **refseq89m** release
from the [MITOS reference-data record](https://zenodo.org/records/4284483).
Skip this step if `/home/ark/MAB/NOAA/mitos_data/refseq89m` is already extracted.

```bash
mkdir -p /home/ark/MAB/NOAA/mitos_data

curl --fail --location --retry 3 \
  https://zenodo.org/records/4284483/files/refseq89m.tar.bz2 \
  --output /home/ark/MAB/NOAA/mitos_data/refseq89m.tar.bz2

printf '%s\n' 'bb6325e27612e61a6995e88e8ffcecf8  /home/ark/MAB/NOAA/mitos_data/refseq89m.tar.bz2' \
  | md5sum --check && \
tar -xjf /home/ark/MAB/NOAA/mitos_data/refseq89m.tar.bz2 \
  -C /home/ark/MAB/NOAA/mitos_data
```

The checksum is published with that release. The resulting directory must be:

```text
/home/ark/MAB/NOAA/mitos_data/refseq89m/
```

Pass its **parent directory** as `--mitos-refdir`, and `refseq89m` as
`--mitos-refseqver`. Do not pass the archive or repeat `refseq89m` in the parent
path. This archive contains annotation models/databases, not the per-sample
GenBank references used for assembly and reference-gene screening.

## 4. Annotated mitochondrial reference genomes

If `/home/ark/MAB/NOAA/mitos_genbank` already contains your downloaded single-record
`.gb` files, reuse it. To prepare a collection from scratch, download the genomic
GenBank flat files from the official
[NCBI RefSeq mitochondrial release directory](https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/).
Use **`.genomic.gbff.gz`**, which includes sequence, features and taxonomy.

```bash
mkdir -p /home/ark/MAB/NOAA/refseq_mito_download

wget --continue --recursive --level=1 --no-parent --no-directories \
  --accept 'mitochondrion.*.genomic.gbff.gz' \
  --directory-prefix /home/ark/MAB/NOAA/refseq_mito_download \
  https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/
```

The download pattern includes all numbered genomic GenBank volumes, without
downloading the separate protein/FASTA releases. Save the download date and release
information with your database records; this NCBI directory is updated over time.

Tortuga expects **one annotated genome per file**. Split the downloaded volumes
and retain animal records using Biopython from the existing controller environment:

```bash
/home/ark/miniconda3/envs/itsme/bin/python - <<'PY'
import gzip
import re
from pathlib import Path
from Bio import SeqIO

source = Path('/home/ark/MAB/NOAA/refseq_mito_download')
destination = Path('/home/ark/MAB/NOAA/mitos_genbank')
archives = sorted(source.glob('mitochondrion.*.genomic.gbff.gz'))
if not archives:
    raise SystemExit('No genomic GenBank archives found; check the download.')
destination.mkdir(parents=True, exist_ok=True)
written = existing = 0
for archive in archives:
    with gzip.open(archive, 'rt') as handle:
        for record in SeqIO.parse(handle, 'genbank'):
            if not {'Metazoa', 'Animalia'}.intersection(record.annotations.get('taxonomy', [])):
                continue
            if not re.fullmatch(r'[A-Za-z0-9_.-]+', record.id):
                raise ValueError('Unexpected accession: ' + record.id)
            output = destination / (record.id + '.gb')
            if output.exists():
                existing += 1
                continue
            temporary = output.with_suffix('.gb.tmp')
            SeqIO.write(record, temporary, 'genbank')
            temporary.replace(output)
            written += 1
print(f'Prepared {written} animal GenBank records; kept {existing} existing files.')
PY
```

This preserves the annotations; renaming a FASTA to `.gb` would not. Existing
accession files are left in place. For a clean later release, prepare a new
database directory instead of mixing reference snapshots.

At first execution, Tortuga creates `tortuga_reference_index.sqlite` automatically
and caches unchanged records. It filters records for animal mitochondrial
annotation, usable CDS/rRNA features, genetic code and configured length bounds;
not every downloaded record will be eligible. Reference names are matched up the
supplied taxonomy hierarchy. No BLAST reference-selection index is needed.

NOVOPlasty receives a derived COX1 seed, MitoFinder the annotated GenBank genome,
mtGrasp its nucleotide sequence, and the gene screen its CDS/rRNA queries.
Tortuga derives these products automatically; you do not need to build separate
mitochondrial protein or nucleotide BLAST databases manually.

## 5. BOLD COI references and taxonomy

Open the official [BOLD data packages page](https://boldsystems.org/data/data-packages/).
Select a public reference snapshot and download its **TSV containing barcode
sequences and taxonomy**, not just its FASTA or JSON package descriptor. Follow
any account/sign-in steps required by BOLD's download page. Extract the TSV if it
is inside an archive; `prepare_bold.py` also reads a gzip-compressed TSV directly.

To reproduce the current Belly library, select **02-Oct-2026**
([snapshot page](https://bench.boldsystems.org/index.php/datapackage?id=BOLD_Public.02-Oct-2026))
and place the extracted TSV at:

```text
/home/ark/MAB/NOAA/BOLD/BOLD_Public.02-Oct-2026.tsv
```

Prepare it with:

```bash
/home/ark/miniconda3/bin/conda run --no-capture-output \
  -p /home/ark/miniconda3/envs/itsme \
  python /home/ark/MAB/bin/Tortuga/prepare_bold.py \
  /home/ark/MAB/NOAA/BOLD/BOLD_Public.02-Oct-2026.tsv \
  --source 'BOLD public snapshot 02-Oct-2026; https://bench.boldsystems.org/index.php/datapackage?id=BOLD_Public.02-Oct-2026' \
  -o /home/ark/MAB/NOAA/BOLD/prepared_fixed
```

Use an empty/new output directory. **Skip preparation if your existing
`prepared_fixed` directory already contains the following three files:**

```text
/home/ark/MAB/NOAA/BOLD/prepared_fixed/coi.fasta
/home/ark/MAB/NOAA/BOLD/prepared_fixed/metadata.tsv
/home/ark/MAB/NOAA/BOLD/prepared_fixed/provenance.json
```

The helper selects eligible COI records, creates matching FASTA/metadata IDs,
retains taxonomy/source information and records input provenance. Its TSV field
limit already supports large BOLD exports. The package's JSON descriptor is not
the `metadata.tsv` expected by Tortuga; that TSV is generated by this helper.

Use both files from the **same preparation**. Do not combine the older 2025 BOLD
FASTA/BLAST files with 2026 metadata. Tortuga builds its own BOLD BLAST index from
`coi.fasta` on the first run, so no manual `makeblastdb` is required.
Allow space for the downloaded/extracted export, prepared FASTA/metadata and the
temporary BLAST index. The existing full prepared Belly library alone is about
16 GB; preparation is not a small download or an in-place conversion.

## Ready to run

Use these resource arguments with the controller commands in
[CONTROLLER.md](CONTROLLER.md):

```bash
--genbank-dir /home/ark/MAB/NOAA/mitos_genbank \
--bold-fasta /home/ark/MAB/NOAA/BOLD/prepared_fixed/coi.fasta \
--bold-metadata /home/ark/MAB/NOAA/BOLD/prepared_fixed/metadata.tsv \
--mitos-refdir /home/ark/MAB/NOAA/mitos_data \
--mitos-refseqver refseq89m \
--itsme-db /home/ark/databases/itsme_db
```

These are arguments to append to a controller invocation, not a standalone shell
command. Fill the sample sheet's intended identities and available taxonomy
before running. The reference index and BOLD BLAST index are generated by Tortuga;
database setup does not infer sample identities.

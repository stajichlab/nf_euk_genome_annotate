# Example sample sheets

`samples.csv` is read by `funannotate.nf`. Columns:

| Column | Required | Notes |
|---|---|---|
| `SPECIES` | yes | species name; combined with `STRAIN` into the output tag |
| `STRAIN` | recommended | isolate/strain; first `;`-token used; blank allowed |
| `ASMID` | yes | assembly id. With NCBI input it is the accession used to resolve the genome; with local FASTA input it is just the unique key for naming/caching (`input_clean_genomes/<ASMID>.fa`) |
| `LOCUSTAG` | yes | locus-tag prefix for funannotate |
| `BUSCO_LINEAGE` | yes | BUSCO odb lineage (e.g. `saccharomycetes_odb10`, `mucoromycota_odb10`) |
| `TRANSL_TABLE` | yes | NCBI translation table (`1` for most fungi) |
| `NCBI_TAXONID` | yes | taxid — used by `GENOME_CLEAN` (taxonkit → phylum) for FCS-GX purge |
| `GENOME` | input-mode dependent | empty → NCBI_ASM resolution by `ASMID`; set → use this local FASTA |

The two input modes can be **mixed in one sheet** (per row).

## 1. NCBI-derived input — [`samples_ncbi.csv`](samples_ncbi.csv)

`GENOME` is left empty, so each genome is resolved from the NCBI_ASM source dir:

```
<--source>/<ASMID>/<ASMID>_genomic.fna.gz
```

Set `--source` (default in `conf/profile_annotate.config`) to your NCBI_ASM mirror:

```bash
sbatch nextflow/run_annotate.sh \
    --samples nextflow/examples/samples_ncbi.csv \
    --source /bigdata/stajichlab/shared/projects/1KFG/2026/NCBI_fungi/source/NCBI_ASM
```

## 2. Local FASTA-only input — [`samples_fasta.csv`](samples_fasta.csv)

`GENOME` points directly at an assembly FASTA (`.fa`/`.fna`, gzipped or plain).
Relative paths resolve against the **launch dir** (where you run nextflow /
submit the sbatch); absolute paths are used as-is. `--source` is irrelevant here.

```bash
sbatch nextflow/run_annotate.sh \
    --samples nextflow/examples/samples_fasta.csv
```

`NCBI_TAXONID` is still required (the FCS-GX cleaning step needs it); use the
closest valid taxid for the organism even for unpublished assemblies.

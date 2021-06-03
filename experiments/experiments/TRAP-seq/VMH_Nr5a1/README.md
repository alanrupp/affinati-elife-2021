# VMH *Nr5a1-Cre* TRAP-seq

Analysis of TRAP-seq data from dissection of ventromedial hypothalamus (VMH) of *Nr5a1-Cre* mice.

## Files
* **counts/*ReadsPerGene.out.tab**: directory of count files per sample from `STAR`
* **analysis.R**: `R` script to find gene enrichment
* **genes.csv**: gene_id to gene_name mapper using Ensembl annotation file (Mus_musculus.GRCm38.92.gtf.gz)
* **metadata.csv**: sample annotation file

## Pipeline
1. `R -e source('analysis.R')`

## Output
* **results.csv**: Enrichment and significance results from `DESeq2`

# VMH *Lepr* TRAP-seq

Looking for genes enriched in *Lepr*-expressing cells of the ventromedial hypothalamus (VMH). To do so, we performed TRAP-seq on LepRb-Cre mice from microdissected VMH. This folder contains the raw FASTQ files and all scripts necessary to find genes enriched in the bead over supernatant.

## Data
The **data** folder contains all the necessary data files for generating the analysis output.
* **Run_1607**: the file from the UM sequencing core containing all FASTQ files and their related metadata
* **Cre.fa** & **Cre.gtf**: *Cre* sequence as its own chromosome in FASTA format and *Cre* gene information in GTF format
* **Gfp_L10a.fa** & **Gfp_L10a.gtf**: *EGFP:L10a* sequence as its own chromosome in FASTA format and *EGFP:L10a* gene information in GTF format
* **Mus_musculus.GRCm38.100.gtf.gz**: Ensembl GTF file for mouse genome
* **Mus_musculus.GRCm38.dna.primary_assembly.fa.gz**: Ensembl FASTA file for mouse genome

## Scripts
The **scripts** folder contains the files necessary for analyzing the files in the **data** folder.
* **mapping.py**: : This will run FastQC for quality control (output: **FASTQC**), filter low quality reads, generate a custom indexed genome with `STAR` containing *Cre* and *EGFP:L10a*, then map the reads to that genome and save the gene counts in a folder called **counts**. It will also create a **tx2gene.csv** file to map Ensembl gene_id to gene_name.
* **analysis.Rmd**: This will analyze the resulting count files from **counts** and determine genes that are enriched in bead over sup, creating a summary PDF (**analysis.pdf**) and a CSV file containing the summarized data (**enrichment.csv**).

## Pipeline
From the command line with the top directory (same directory as this `README.md` file) as the working directory, run:
1. `python scripts/mapping.py --directory data/Run_1607/mmyers --custom_reference`
2. `R -e "rmarkdown::render('scripts/analysis.Rmd', output_dir = 'results')"`

## Output
The **results** folder contains the relevant outputs:
* **analysis.pdf**: R markdown file of analysis of count files
* **enrichment.csv**: Ouput from `DESeq2` containing expression level in bead and sup and enrichment (and significance) of all genes.

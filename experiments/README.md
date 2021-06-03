# Experiments
This directory contains all the data and analyses from the manuscript.

## Organization
* **snRNA-seq**: all snRNA-seq experiments
  * **GSE93374**: data from [Campbell et al. 2017](https://www.nature.com/articles/nn.4495) for the `campbell` analysis
  * **NovaA.+** and **Run.+**: raw matrix files from `cellranger`
  * **ypx3sw2f7c.2**: data from [Kim et al. 2019](https://doi.org/10.17632/ypx3sw2f7c.2) for the `integration` analysis
  * `mouse_to_macaque.tsv.gz`: ortholog mapper from Ensembl
* **analysis**
  * **campbell**: integrating our dataset with data from [Campbell et al. 2017](https://www.nature.com/articles/nn.4495)
  * **celltypes**: generating cell type markers from ESmu scores from published hypothalamus datasets
  * **conservation**: generating a combined mouse & macaque VMH atlas
  * **integration**: integrating our dataset with data from [Kim et al. 2019](https://doi.org/10.1016/j.cell.2019.09.020)
  * **kim**: analyzing data from [Kim et al. 2019](https://doi.org/10.1016/j.cell.2019.09.020)
  * **macaque**: atlas of macaque (*Macaca mulatta*) VMH populations
  * **mouse**: atlas of mouse (*Mus musculus*) VMH populations
* **TRAP-seq**
  * **comparison**: script for comparing results from TRAP-seq experiments
  * **Hypo_Lepr**: *Lepr-Cre* TRAP-seq analysis from hypothalamus
  * **VMH_Lepr**: *Lepr-Cre* TRAP-seq analysis from VMH
  * **VMH_Nr5a1**: *Nr5a1-Cre* TRAP-seq analysis from VMH
  * **VMH_Slc17a6+Lepr**: *Slc17a6-Flpo*;*Lepr-Cre* TRAP-seq from VMH

## Requirements
The required R and python libraries to run this portion of the analysis
* `R`
  * `ggdendro`
  * `ggrepel`
  * `kableExtra`
  * `knitr`
  * `rcartocolor`
  * `reticulate`
  * `scran`
  * `SCTransform`
  * `Seurat`
  * `tidyverse`
  * `wesanderson`
* `python`
  * `CELLEX`
  * `scrublet`

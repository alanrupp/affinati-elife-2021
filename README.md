# Affinati et al. 2021
This folder contains all the data and analysis required to generate the results from the Affinati et al. manuscript (full citation below).

## Organization
See the subdirectories for more information
* **experiments**: the raw data and analysis code
  * **snRNA-seq**: 10x data and analyses
  * **TRAP-seq**: TRAP-seq data and analyses
* **figures**: scripts and images to generate the figures; figures are in OpenDocument and Inkscape SVG formats
* **tools**: common functions that are used throughout the analysis
  * **dropseq3**: functions and data for analyzing snRNA-seq data

## Requirements
All analysis was run with `R 3.6.3` and `python 3.8`. Additional libraries are required for different experiments and are listed separately in each experiment's folder.

## Order
The analysis needs to be run in a particular order to generate all the files that are inputs to subsequent steps. This should follow the order of the paper:
1. Analyze TRAP-seq data: each directory has instructions on how to go from FASTQ to count files to enrichment values
  * generate the *Nr5a1* VMH TRAP-seq dataset: `experiments/TRAP-seq/VMH_Nr5a1`
  * generate the *Lepr* VMH TRAP-seq dataset: `experiments/TRAP-seq/VMH_Lepr`
  * generate the *Lepr* hypothalamus TRAP-seq dataset: `experiments/TRAP-seq/Hypo_Lepr`
  * generate the *Slc17a6*/*Lepr* VMH TRAP-seq dataset: `experiments/TRAP-seq/VMH_Slc17a6+Lepr`
  * compare the enrichment of VMH *Lepr*-specific genes in each dataset: `experiments/TRAP-seq/comparison`
2. Analyze the mouse snRNA-seq data: navigate to each directory to run/knit the R/Rmd code and save output
  * generate cell-type modules from published data: `experiments/snRNA-seq/modules.R`
  * generate the mouse dataset: `experiments/snRNA-seq/mouse/mouse.Rmd`
  * analyze the Kim et al. dataset: `experiments/snRNA-seq/kim/kim.Rmd`
  * compare to Kim et al. dataset: `experiments/snRNA-seq/integration/integration.Rmd`
  * compare to Campbell et al. dataset: `experiments/snRNA-seq/campbell/campbell.Rmd`
3. Analyze macaque snRNA-seq data: navigate to each directory to knit the Rmd document and save output
  * generate the macaque dataset: `experiments/snRNA-seq/macaque/macaque.Rmd`
  * generate the conserved mouse+macaque dataset: `experiments/snRNA-seq/conservation/conservation.Rmd`
4. Generate figure panels: `figures/figures.Rmd`

## Citation
Affinati et al., Cross-Species analysis defines the conservation of anatomically-segregated VMH neuron populations. *eLife*. 2021 May 21;10:e69065. doi: [10.7554/eLife.69065](https://doi.org/10.7554/elife.69065).

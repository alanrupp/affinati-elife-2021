# Kim et al. data

This is the dataset from [Kim et al.](https://data.mendeley.com/datasets/ypx3sw2f7c/2) of single-cell RNA-seq in the VMH. I've organized the data into this folder and copied their description of the data (**Description of this data**). Data came as a series of `*.zip` files (see below **Experiment data files**) with file names corresponding to treatment. I wrote a script (`preprocess.py`) that generates a `samples.csv` file to link the original file name that contains treatment info with the sample date (the output from unzipping), unzipped everything, then gzipped the individual files in each subfolder.

## Description of this data
Single-cell RNA-seq (scRNA-seq) data collected from the ventral part of ventromedial hypothalamus (VMH) are listed below (after alignment process). You can find the detailed information for these data in the following paper, "Multimodal Analysis of Cell Types in a Hypothalamic Node Controlling Social Behavior" (doi: 10.1016/j.cell.2019.09.020). Briefly, we performed two independent scRNA-seq methods, SMART-seq & 10x (Chromium v2 chemistry), and all 10x data (each dataset was collected from 2-4 brains) were collected for Act-seq (behavior + scRNA-seq), whereas SMART-seq data included Retro-seq data (retrograde-labeling + scRNA-seq; see the detail in the metadata).

## Experiment data files
* 10x_VMH_Female_Control_1.zip
* 10x_VMH_Female_Mating_not_receptive_1.zip
* 10x_VMH_Female_Mating_not_receptive_2.zip
* 10x_VMH_Female_Mating_not_receptive_3.zip
* 10x_VMH_Female_Plain_1.zip
* 10x_VMH_Male_Aggression_1.zip
* 10x_VMH_Male_Aggression_2.zip
* 10x_VMH_Male_Aggression_3.zip
* 10x_VMH_Male_Control_1.zip
* 10x_VMH_Male_Control_2.zip
* 10x_VMH_Male_Control_3.zip
* 10x_VMH_Male_M-F_CI_dangled_1.zip
* 10x_VMH_Male_M-F_CI_pencil_cup_1.zip
* 10x_VMH_Male_M-M_CI_dangled_1.zip
* 10x_VMH_Male_M-M_CI_dangled_2.zip
* 10x_VMH_Male_M-M_CI_pencil_cup_1.zip
* 10x_VMH_Male_Mating_1.zip
* 10x_VMH_Male_Others_1.zip
* 10x_VMH_Male_Others_2.zip
* 10x_VMH_Male_Others_3.zip
* 10x_VMH_Male_Others_4.zip
* 10x_VMH_Male_Plain_1.zip
* 10x_VMH_Male_Plain_2.zip
* 10x_VMH_Male_Social_Fear_Group_housed_1.zip
* 10x_VMH_Male_Social_Fear_Singly-housed_1.zip
* 10x_VMH_Male_Social_Fear_Singly-housed_2.zip
* 10x_VMH_metadata.csv
* SMART-seq_VMH_cpm.rda.zip
* SMART-seq_VMH_metadata.csv

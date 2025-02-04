# TDA1_proj

Transcriptomic analysis of *TDA1* deletion in *Saccharomyces cerevisiae* during the diauxic shift. This repository contains the code and intermediate data necessary to reproduce the analysis, figure, and table. The .fastq files can be found at the European Nucleotide Archive (ENA) under the project ID "PRJEB59812".

Last updated: 2025-02-04

Connected paper can be found at:

## Project information

![Figure1](./plots/Fig3.svg){width="100%"}

The biological role of *TDA1* in *Saccharomyces cerevisiae* was investigated by conducting three independent transcriptomic experiments. A reference strain (BY4741) and *tda1∆* strain were cultivated in minimal media and were sampled for RNA-seq at 9 hours (pre-diauxic shift/log-phase) and 26 hours (post-diauxic shift/PDS-phase). Growth data for the samples in these three experiments can be found in the connected paper.

## Directory structure

```         
TDA1_proj
|   README.md
|   LICENSE
|   TDA_flask_exp1.Rmd (Script to analyse data from experiment 1)
|   TDA_flask_exp2.Rmd (Script to analyse data from experiment 2)
|   TDA_flask_exp3.Rmd (Script to analyse data from experiment 3)
|   Union_TDA.Rmd (Script to analyse combined data from all experiments)
|   TDA_growth.Rmd (Generates growth curve from separate growth experiment)
|
└───data
│   └───db (Contains data from databases, e.g. gene names, GO-sets, TF target)
│   └───Exp1 (Contains raw count data from Experiment 1)
│   └───Exp2 (Contains raw count data from Experiment 2)
│   └───Exp3 (Contains raw count data from Experiment 3)
│   └───Growth (Contains measurements from separate growth experiment)
|
└───plots
|   (Contains all plots found in connected paper)
|   └───Unedited (Contains unedited plots directly generated from code)
|
└───r_objects (Intermediate r objects to transfer data between Rmd scripts)
|
└───scripts
|   supp.R (Contains functions made/used in this project)
|
└───tables (Supplementary tables)
```

## Analysis

The differential gene expression analysis is performed for each expriment separately. TDA_flask_exp1.Rmd, TDA_flask_exp2.Rmd, and TDA_flask_exp3.Rmd is run to generate the results tables from DESeq2 analysis, which is then saved in "/r_objects". Union_TDA.Rmd is then run to analyse and compare the output from each experiment. The output of Union_TDA.Rmd are found in "./plots" and "./tables".

TDA_growth.Rmd is run to generate a smoothed growth curve using the data found in "/data/Growth".

## Citation

If you make use of the data or functions used in this project, we would appreciate a citation of the original paper.

# COVIDs vs. ARDS files

## comet.genes.csv
A matrix of gene counts for the subjects included in this study. This can also be found on GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163426

## covidardsmetadata-final.csv
This file corresponds to Table 1 in the manuscript. There are three key differences. First, to maintain data privacy, patient ages >89 are censored. One subject was older than 89, and their age was replaced with 89 in this table. As this is the greatest value in the distribution, it has no effect on the reported median or Wilcoxon test comparing the ages in each cohort. Second, two subjects did not have recorded PF ratios and were included on the basis of SF ratios < 315. These subjects are identified with an asterisk in the version of this table that is included in the manuscript, but are reported as missing for this table to prevent any confusion in downstream analyses. Third, the table here also includes the SARS-CoV-2 viral load quantified in each sample from the metagenomic sequencing (reads-per-million, rpm). These values were generated with the pipeline here: https://github.com/czbiohub/sc2-msspe-bioinfo.

## COVID vs ARDS.Rmd
This R notebook can be used to reproduce the differential expression analysis and heatmaps in the manuscript. 

## rnaseqfunctions.r
Functions to facilitate analysis and visualization of RNASeq data.

## NP_TA_IFN_viral_load_regression.R
This code performs robust regression of interferon stimulated gene expression against SARS-CoV-2 viral load (rpm) in the tracheal aspirate (TA) samples reported in this study as compared to nasopharyngeal (NP) swab samples we [previously reported](https://doi.org/10.1038/s41467-020-19587-y).
Also included in the repo are the NP swab gene counts & metadata and the interferon geneset we used.

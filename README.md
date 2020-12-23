# COVIDs vs. ARDS files

## Metadata file
This CSV file corresponds to Table 1 in the manuscript. There are two key differences. First, to maintain data privacy, patient ages >89 are censored. One subject was older than 89, and their age was replaced with 89 in this table. As this is the greatest value in the distribution, it has no effect on the reported median or Wilcoxon test comparing the ages in each cohort. Second, two subjects did not have recorded PF ratios and were included on the basis of SF ratios < 315. These subjects are identified with an asterisk in the version of this table that is included in with the manuscript, but are reported as missing for this table to prevent any confusion in downstream analyses. 

## COVID vs ARDS.Rmd
This R notebook can be used to reproduce the differential expression analysis and heatmaps in the manuscript. 

## comet.genes.csv
A matrix of gene counts for the subjects included in this study. This can also be found on GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163426

## rnaseqfunctions.r
Functions to facilitate analysis and visualization of RNASeq data

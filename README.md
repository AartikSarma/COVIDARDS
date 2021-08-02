# COVIDs vs. ARDS files

## comet.genes.csv
A matrix of gene counts for the subjects included in this study. This can also be found on GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163426

## covidardsmetadata.csv
This metadata file is the source data for Table 1 in the manuscript. There are three key differences. First, to maintain data privacy, patient ages >89 are censored. One subject was older than 89, and their age was replaced with 89 in this table. As this is the greatest value in the distribution, it has no effect on the reported median or Wilcoxon test comparing the ages in each cohort. Second, two subjects did not have recorded PF ratios and were included on the basis of SF ratios < 315. These subjects are identified with an asterisk in the version of this table that is included in the manuscript, but are reported as missing here to prevent any confusion in downstream analyses. Third, the table here also includes the SARS-CoV-2 viral load quantified in each sample from the metagenomic sequencing (reads-per-million, rpm). These values were generated with the pipeline here: https://github.com/czbiohub/sc2-msspe-bioinfo.

## COVID vs ARDS.Rmd
This R notebook can be used to reproduce the differential expression analysis, heatmaps, and PCA plots in the manuscript. Results of differential expression analyses can found in the .deseq.csv files.  

## rnaseqfunctions.r
Functions to facilitate analysis and visualization of RNASeq data.

## NP_TA_IFN_viral_load_regression.R
This code performs robust regression of interferon (IFN) stimulated gene expression against SARS-CoV-2 viral load (rpm) in the tracheal aspirate (TA) samples reported in this study as compared to nasopharyngeal (NP) swab samples we [previously reported](https://doi.org/10.1038/s41467-020-19587-y).
Also included in the repo are the NP swab gene counts & metadata and the interferon gene set used.

## SingleCell_COVIDARDS_Data_Ingest_QC_Annotated.ipynb
This code ingests all single cell data used in the ARDS study. This data was generated using the 5' 10x Single Cell V(D)Jv1.1 method. Data was ingested, QC'd, Batch Corrected using Harmony https://github.com/immunogenomics/harmony and Annotated using the Scanpy pipeline https://scanpy.readthedocs.io/en/stable/. Cell type composition was also analyzed within COVID19 patients within the ARDS cohort. 

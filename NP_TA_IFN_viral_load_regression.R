library(dplyr)
library(tidyr)
library(tibble)
library(limma)
library(robustbase)

TA_counts_path <- "comet.genes.csv"
TA_metadata_path <- "covidardsmetadata-final.csv"
NP_counts_path <- "NP_swab_gene_counts.csv"
NP_metadata_path <- "NP_swab_metadata.csv"
ifn_geneset_path <- "IFN_geneset.txt"

## Read in TA gene counts

TA_counts <- read.csv(file = TA_counts_path,
                      stringsAsFactors = F,
                      header = T,
                      check.names = F)

# Parse the gene ID to gene symbol mapping
TA_counts %>% 
  select(ensembl_gene_id, gene_symbol) %>% 
  column_to_rownames("ensembl_gene_id") ->
  gene_names
# And remove the gene symbol column from the counts table
TA_counts %>% 
  select(-gene_symbol) %>% 
  column_to_rownames("ensembl_gene_id") ->
  TA_counts

## Read in TA sample metadata

# Calculate log10 reads-per-million (rpm) mapping to SARS-CoV-2 (SC2),
# adding a value of 0.1 to samples with rpm less than 0.1
TA_samples <- read.csv(TA_metadata_path,
                       stringsAsFactors = F,
                       header = T,
                       check.names = F) %>% 
  mutate("log10_rpm" = ifelse(SC2_rpm < 0.1, log10(SC2_rpm+0.1), log10(SC2_rpm)))

# Subset to genes with more than 10 counts in at least 20% of samples
keep <- rowSums(TA_counts > 10) > 0.2*ncol(TA_counts)
TA_counts_sub <- TA_counts[keep,]

# Quantile normalize with limma
design_matrix <- model.matrix(~1, data = TA_samples)
TA_vwts <- voom(TA_counts_sub,
                design=design_matrix,
                normalize.method="quantile",
                plot=T)

## Read in IFN geneset

ifn_geneset_sym <- read.table(ifn_geneset_path,
                              stringsAsFactors = F,
                              header = F)[,1]
# Convert gene symbols to gene IDs
gene_names %>% 
  rownames_to_column("gene_id") %>% 
  filter(gene_symbol %in% ifn_geneset_sym) %>% 
  pull(gene_id) ->
  ifn_geneset_id

## Robust regression of IFN genes on viral load in TA samples

TA_rpm_gene_robs <- as.data.frame(t(sapply(ifn_geneset_id, function(g) {
  TA_samples %>%
    mutate("norm_counts" = TA_vwts$E[g, ]) %>%
    filter(`ARDS group` == "COVID-ARDS") ->
    reg_df
  rpm_gene_rob <- summary(lmrob(norm_counts ~ log10_rpm,
                                data = reg_df,
                                setting = "KS2014"))
  c("TA_intercept" = rpm_gene_rob$coefficients[1,1],
    "TA_viral_load_slope" = rpm_gene_rob$coefficients[2,1],
    "TA_adjR2" = rpm_gene_rob$adj.r.squared,
    "TA_viral_load_pval" = rpm_gene_rob$coefficients[2,4]
    )
})))
# Adjust for multiple testing
TA_rpm_gene_robs[,"TA_viral_load_padj"] <- p.adjust(TA_rpm_gene_robs$TA_viral_load_pval,
                                                        method="BH")
## Read in NP gene counts

NP_counts <- read.csv(file = NP_counts_path,
                      stringsAsFactors = F,
                      header = T,
                      check.names = F,
                      row.names = 1)

## Read in NP sample metadata

# Calculate log10 reads-per-million (rpm) mapping to SARS-CoV-2 (SC2),
# adding a value of 0.1 to samples with rpm less than 0.1
NP_swab_samples <- read.csv(NP_metadata_path,
                            stringsAsFactors = F,
                            header = T,
                            check.names = F) %>% 
  mutate("log10_rpm" = ifelse(SC2_rpm < 0.1, log10(SC2_rpm+0.1), log10(SC2_rpm)))

# Subset to genes with more than 10 counts in at least 20% of samples
keep <- rowSums(NP_counts > 10) > 0.2*ncol(NP_counts)
NP_counts_sub <- NP_counts[keep,]

# Quantile normalize with limma
design_matrix <- model.matrix(~1, data = NP_swab_samples)
NP_vwts <- voom(NP_counts_sub,
                design=design_matrix,
                normalize.method="quantile",
                plot=T)

## Robust regression of IFN genes on viral load in TA samples

NP_rpm_gene_robs <- as.data.frame(t(sapply(ifn_geneset_id, function(g) {
  NP_swab_samples %>%
    mutate("norm_counts" = NP_vwts$E[g, ]) %>%
    filter(viral_status == "SC2") ->
    reg_df
  rpm_gene_rob <- summary(lmrob(norm_counts ~ log10_rpm,
                                data = reg_df,
                                setting = "KS2014"))
  c("NP_intercept" = rpm_gene_rob$coefficients[1,1],
    "NP_viral_load_slope" = rpm_gene_rob$coefficients[2,1],
    "NP_adjR2" = rpm_gene_rob$adj.r.squared,
    "NP_viral_load_pval" = rpm_gene_rob$coefficients[2,4]
    )
})))
# Adjust for multiple testing
NP_rpm_gene_robs[,"NP_viral_load_padj"] <- p.adjust(NP_rpm_gene_robs$NP_viral_load_pval,
                                                        method="BH")

## Combine NP and TA stats

comb_rpm_gene_robs <- cbind(TA_rpm_gene_robs, NP_rpm_gene_robs)
comb_rpm_gene_robs[,"gene"] <- gene_names[rownames(comb_rpm_gene_robs), "gene_symbol"]
comb_rpm_gene_robs %>%
  select("gene", everything()) ->
  comb_rpm_gene_robs

## Output combined table

comb_rpm_gene_robs %>%
  write.csv("TA_NP_SC2_rpm_regression.csv")

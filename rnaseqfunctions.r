library(dendsort)
library(pheatmap)
library(viridis)
library(dendextend)

#Create a new experiment folder
new.exp <- function(experimentname){
  dir.create(experimentname)
  setwd(paste("./", experimentname, sep= ""))
  dir.create("Rdafiles")
  dir.create("deseq")
  dir.create("gsea")
  dir.create("ipaoutput")
  dir.create("plots")
  setwd("..")
}


deseqresults <- function(dds, #A DESeq object
                         methods = c("ashr"), #Select method for shrinking contrasts
                         writefiles = TRUE #If TRUE, write the DESeq result table to a CSV file
){
  
  stopifnot(class(dds) == "DESeqDataSet")
  
  #Extract the experiment name from dds object
  expname <- deparse(substitute(dds)) %>% str_remove(".dds")
  dir.create(expname) 
  setwd(expname)
  paste("Current experiment", expname, "...") %>% print
  
  #Write the transcript count and metadata used for this experiment in the folder
  dds %>% colData %>% as.data.frame %>%
    write_csv(file = paste(expname,".subjects.csv", sep = ""))
  dds %>% counts %>% as.data.frame %>%  rownames_to_column("ensembl_gene_id") %>%
    write_csv(file = paste(expname,".genes.csv", sep = ""))
  resnames <- resultsNames(dds)[-1]
  
  
  #Extract results for each coef using selected methods
  for(coefname in resnames){
    resultfolder <- coefname %>% str_remove(".+?_") %>% str_remove_all("_") 
    new.exp(resultfolder)
    setwd(resultfolder)
    for (method in methods){
      resultname <- paste(expname, ".", resultfolder, ".", method, sep ="")
      paste("Applying", method, "to", expname, "for coef", resultfolder, "...") %>% print
      print(dds)
      print(method)
      print(coefname)
      foo <- lfcShrink(dds, type = method, coef = coefname)
      foo %<>%
        as.data.frame %>%
        rownames_to_column("ensembl_gene_id") %>%
        left_join(ensembltohgnc) %>%
        mutate(ihwP = ihw(pvalue ~ baseMean,  data = ., alpha = 0.1) %>% adj_pvalues())
      assign(x=resultname, value= foo, envir = .GlobalEnv)
      if(writefiles){write_csv(foo, file = paste("./deseq/", resultname, ".deseq.csv", sep = ""))}
    }
    setwd("..") #End contrast specific results
  }
  setwd("..") #Return to main directory
}


###########################
#### Heatmap functions ####
###########################

sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

#Create a vector of deciles of gene expression
#Can use this for make color scales for heatmaps
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

#Function to facilitate rotating clustered heatmaps
pheatmap.rotate <- function(input.mat, #Matrix of gene expression
                            uselastdend = F, #Once the function has been run once, can set this option to TRUE to use the prior dendrogram
                            cluster_cols = TRUE,
                            cluster_rows = F,
                            rotate = T, #If TRUE, prompts to rotate dendrogram
                            clustermethod = "euclidean", #Cluster method to use for columns
                            cutree_cols = F,
                            ...){
  
  mat_cluster_rows <- sort_hclust(hclust(dist((input.mat))))
  mat_cluster_cols <- sort_hclust(hclust(dist(t(input.mat), method = clustermethod)))
  
  if(cluster_cols == FALSE){
    mat_cluster_cols <- FALSE
  }
  
  if(uselastdend == TRUE){
    mat_cluster_cols <- rotated
  }
  
  pheatmap(mat = input.mat, cluster_cols = mat_cluster_cols, cluster_rows = cluster_rows, ...)
  
  if(rotate == T){ 
    click <- (toupper(readline(prompt="Rotate dendrogram? (Y/N)")) == "Y")
  }
  else{
    click = FALSE
  }
  
  if(click == TRUE) {
    rotated <<- mat_cluster_cols %>%
      as.dendrogram %>%
      click_rotate(continue = TRUE) %>% 
      as.hclust
    
    pheatmap.rotate(input.mat, uselastdend = T, cluster_cols = cluster_cols, cluster_rows = cluster_rows, rotate = T, cutree_cols = cutree_cols, ...)
    
  }
  
}
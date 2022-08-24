################################################################################
## AUTHOR: NAROA LEGARRA MARCOS
## DATE: 10/03/22

## Functions used to process and perform the analysis of sc RNAseq and FACS-bulk 
## RNAseq data
################################################################################

wd <- "~/Desktop/TFM" #set the working directory
pseudobulk <- "/DATA/pseudo/" #path to the folder to store the pseudobulks
single_cell <- "/DATA/single_cell/"

############################## LOAD LIBRARIES ##################################
library(readr)
library(tidyverse)
library(stringr)
library(biomaRt)
library(gplots)
library(org.Hs.eg.db)
library(RColorBrewer)
library(pheatmap)
library(reshape)
library(compiler)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(ggfortify)
library(ggvenn)
library(hues)
library(circlize)
library(EnsDb.Hsapiens.v79)
library(AnnotationDbi)
library(clusterProfiler)
library(enrichplot)
library(ggnewscale)
library(DOSE)
library(gridExtra)
library(scuttle)
library(FactoMineR)

# FUNCTIONS ################################################################################

sparsity <- function(matrix){
  #####################################################################################
  #Function that calculates the sparsity of a matrix
  
  # INPUT VARIABLES
  #matrix <- a matrix 
  
  # OUTPUT VARIABLES
  # sparsity <- the sparsity of the matrix
  #####################################################################################
  return(sum(matrix==0)/(ncol(matrix)*nrow(matrix)))
}

cell_num <- function(matrix,name){
  #####################################################################################
  #Function that calculates the number of each cell type in the matrix
  
  # INPUT VARIABLES
  #matrix <- sc sequencing matrix
  
  # OUTPUT VARIABLES
  # cell_num <- a tibble with the number of HSC, MEP, GMP and LMPP in the matrix
  #####################################################################################
  #measure the number of columns of each cell
  HSC_sc <- ncol(matrix[,str_detect(colnames(matrix),"^HSC")]) 
  MEP_sc <- ncol(matrix[,str_detect(colnames(matrix),"^MEP")])
  GMP_sc <- ncol(matrix[,str_detect(colnames(matrix),"^GMP")])
  cell_num <- tibble(HSC=HSC_sc,MEP=MEP_sc,GMP=GMP_sc) #construct a tibble with the number of each cell type
  saveRDS(cell_num,paste0(wd,pseudobulk,name,"_cellnum.rds")) #save 
}

assign_labels <- function(matrix,labels){
  #####################################################################################
  #Function that assings the cell type label to a matrix count with cell tags
  
  # INPUT VARIABLES
  #matrix <- matrix with counts and sequecing tags
  #labels <- cell type for each sequencing tag
  
  # OUTPUT VARIABLES
  # counts <- matrix with counts and cell type labels
  #####################################################################################
  
  counts <- matrix@assays$RNA@counts #extract the count matrix
  tags <- as.data.frame(str_split(counts@Dimnames[[2]],"-",simplify=T)) #extract the sequencing tags in the count matrix
  tags2 <- str_split(labels$X1,"_",simplify=T) #extract the sequencing tags in the labels matrix
  labels$X1 <- tags2[,1] #add the sequenicng tags as a new column in the labels matrix
  tags_full <- left_join(tags,labels,by=c("V2"="X1")) #join the labels matrix to the count matrix tags
  counts@Dimnames[[2]] <- tags_full$X2 #add the corresponding cell labels to the count matrix
  counts2 <- as.data.frame(counts) #convert count matrix to data.frame
  counts2$labels <- rownames(counts2) #add rownames (gene names) as a new column names "labels"
  colnames(counts2) <- make.unique(colnames(counts2),sep="_") #as cell types are repeated, add a nuber suffix
  colnames(counts2)[is.na(colnames(counts2))] <- "Unknown" #replace NAs by "Unknown" 
  return(counts2)
}

pseudosc_mean <- function(matrix,suffix){
  ##############################################################################
  #Function that generates pseudosc matrix from single cell data
  
  # INPUT VARIABLES
  #matrix <- single cell matrix
  #suffix <- condition (senior or mds)
  
  # OUTPUT VARIABLES
  # pseudo_sc <- pseudo single cell matrix
  # HSC_sc_lc <- low genes expressed genes (<5 counts in total)
  # MEP_sc_lc <- low genes expressed genes (<5 counts in total)
  # GMP_sc_lc <- low genes expressed genes (<5 counts in total)
  #####################################################################################
  genes <- matrix$labels #extract the gene names
  
  #extract data for HSC, MEP, GMP and LMPP
  HSC_sc <- matrix[,str_detect(colnames(matrix),"^HSC")] 
  MEP_sc <- matrix[,str_detect(colnames(matrix),"^MEP")]
  GMP_sc <- matrix[,str_detect(colnames(matrix),"^GMP")]
  
  #analyse low expressed genes
  HSC_sc_lc <- genes[apply(HSC_sc,1,sum)<=5]
  MEP_sc_lc <- genes[apply(MEP_sc,1,sum)<=5]
  GMP_sc_lc <- genes[apply(GMP_sc,1,sum)<=5]

  #convert to numeric (just in case)
  HSC_sc <- apply(HSC_sc,1,as.numeric)
  MEP_sc <- apply(MEP_sc,1,as.numeric)
  GMP_sc <- apply(GMP_sc,1,as.numeric)
  
  #calculate mean and generate pseudosingle cell matrix (multiply by 1e5 and convert to integer because
  # the analysis needs input as integers and we only care about proportions)
  HSC_sc_mean <- as.integer(apply(HSC_sc,2,mean)*10000)
  MEP_sc_mean <- as.integer(apply(MEP_sc,2,mean)*10000)
  GMP_sc_mean <- as.integer(apply(GMP_sc,2,mean)*10000)
  pseudo_sc <- data.frame(genes,HSC_sc_mean,MEP_sc_mean,GMP_sc_mean)
  
  colnames(pseudo_sc) <- str_c(colnames(pseudo_sc),suffix)
  return(list(pseudo_sc=pseudo_sc,HSC_sc_lc=HSC_sc_lc,
              MEP_sc_lc=MEP_sc_lc,GMP_sc_lc=GMP_sc_lc))
}

pseudosc_sum <- function(matrix,suffix){
  ##############################################################################
  #Function that generates pseudosc matrix from single cell data
  
  # INPUT VARIABLES
  #matrix <- single cell matrix
  #suffix <- condition (senior or mds)
  
  # OUTPUT VARIABLES
  # pseudo_sc <- pseudo single cell matrix
  # HSC_sc_lc <- low genes expressed genes (<5 counts in total)
  # MEP_sc_lc <- low genes expressed genes (<5 counts in total)
  # GMP_sc_lc <- low genes expressed genes (<5 counts in total)
  ##############################################################################
  genes <- matrix$labels
  
  #extract data for HSC, MEP, GMP and LMPP
  HSC_sc <- matrix[,str_detect(colnames(matrix),"^HSC")]
  MEP_sc <- matrix[,str_detect(colnames(matrix),"^MEP")]
  GMP_sc <- matrix[,str_detect(colnames(matrix),"^GMP")]
  
  #analyze low expressed genes
  HSC_sc_lc <- genes[apply(HSC_sc,1,sum)<=5]
  MEP_sc_lc <- genes[apply(MEP_sc,1,sum)<=5]
  GMP_sc_lc <- genes[apply(GMP_sc,1,sum)<=5]
  
  #convert to numeric (just in case)
  HSC_sc <- apply(HSC_sc,1,as.numeric)
  MEP_sc <- apply(MEP_sc,1,as.numeric)
  GMP_sc <- apply(GMP_sc,1,as.numeric)
  
  #calculate mean and generate pseudosingle cell matrix (multiply by 1e5 and convert to integer because
  # the analysis needs input as integers and we only care about proportions)
  HSC_sc_sum <- as.integer(apply(HSC_sc,2,sum))
  MEP_sc_sum <- as.integer(apply(MEP_sc,2,sum))
  GMP_sc_sum <- as.integer(apply(GMP_sc,2,sum))
  pseudo_sc <- data.frame(genes,HSC_sc_sum,MEP_sc_sum,GMP_sc_sum)
  
  colnames(pseudo_sc) <- str_c(colnames(pseudo_sc),suffix)
  return(list(pseudo_sc=pseudo_sc,HSC_sc_lc=HSC_sc_lc,
              MEP_sc_lc=MEP_sc_lc,GMP_sc_lc=GMP_sc_lc))
}

pseudosc_rep_sum <- function(matrix,suffix){
  genes <- matrix$labels
  ##############################################################################
  #Function that generates pseudosc matrix from single cell data
  
  # INPUT VARIABLES
  #matrix <- single cell matrix
  #suffix <- condition (senior or mds)
  
  # OUTPUT VARIABLES
  # pseudo_sc <- pseudo single cell matrix
  # HSC_sc_lc <- low genes expressed genes (<5 counts in total)
  # MEP_sc_lc <- low genes expressed genes (<5 counts in total)
  # GMP_sc_lc <- low genes expressed genes (<5 counts in total)
  ##############################################################################
  
  # extract data for HSC, MEP, GMP and LMPP
  HSC_sc <- matrix[,str_detect(colnames(matrix),"^HSC")]
  MEP_sc <- matrix[,str_detect(colnames(matrix),"^MEP")]
  GMP_sc <- matrix[,str_detect(colnames(matrix),"^GMP")]
  
  # convert to numeric (just in case)
  HSC_sc <- apply(HSC_sc,1,as.numeric)
  MEP_sc <- apply(MEP_sc,1,as.numeric)
  GMP_sc <- apply(GMP_sc,1,as.numeric)
  
  # Replicate 1 ################################################################
  
  # select a sample
  HSC_sc_sample <- HSC_sc[sample(1:nrow(HSC_sc),0.2*nrow(HSC_sc),replace=F),]
  MEP_sc_sample <- MEP_sc[sample(1:nrow(MEP_sc),0.2*nrow(MEP_sc),replace=F),]
  GMP_sc_sample <- GMP_sc[sample(1:nrow(GMP_sc),0.2*nrow(GMP_sc),replace=F),]
  
  # analyse low expressed gene sin the sample
  HSC_sc_lc1 <- apply(HSC_sc_sample,2,sum)<=5
  MEP_sc_lc1 <- apply(MEP_sc_sample,2,sum)<=5
  GMP_sc_lc1 <- apply(GMP_sc_sample,2,sum)<=5
  
  #calculate mean and generate pseudo single cell matrix by the sum of all the counts
  HSC_sc_sum <- as.integer(apply(HSC_sc_sample,2,sum))
  MEP_sc_sum <- as.integer(apply(MEP_sc_sample,2,sum))
  GMP_sc_sum <- as.integer(apply(GMP_sc_sample,2,sum))
  pseudo_sc_1 <- data.frame(genes,HSC_sc_sum,MEP_sc_sum,GMP_sc_sum)
  
  # Replicate 2 ################################################################
  
  HSC_sc_sample <- HSC_sc[sample(1:nrow(HSC_sc),0.2*nrow(HSC_sc),replace=F),]
  MEP_sc_sample <- MEP_sc[sample(1:nrow(MEP_sc),0.2*nrow(MEP_sc),replace=F),]
  GMP_sc_sample <- GMP_sc[sample(1:nrow(GMP_sc),0.2*nrow(GMP_sc),replace=F),]
  
  HSC_sc_lc2 <- apply(HSC_sc_sample,2,sum)<=5
  MEP_sc_lc2 <- apply(MEP_sc_sample,2,sum)<=5
  GMP_sc_lc2 <- apply(GMP_sc_sample,2,sum)<=5
  
  HSC_sc_sum <- as.integer(apply(HSC_sc_sample,2,sum))
  MEP_sc_sum <- as.integer(apply(MEP_sc_sample,2,sum))
  GMP_sc_sum <- as.integer(apply(GMP_sc_sample,2,sum))
  pseudo_sc_2 <- data.frame(genes,HSC_sc_sum,MEP_sc_sum,GMP_sc_sum)
  
  # Replicate 3 ################################################################
  
  HSC_sc_sample <- HSC_sc[sample(1:nrow(HSC_sc),0.2*nrow(HSC_sc),replace=F),]
  MEP_sc_sample <- MEP_sc[sample(1:nrow(MEP_sc),0.2*nrow(MEP_sc),replace=F),]
  GMP_sc_sample <- GMP_sc[sample(1:nrow(GMP_sc),0.2*nrow(GMP_sc),replace=F),]
  
  HSC_sc_lc3 <- apply(HSC_sc_sample,2,sum)<=5
  MEP_sc_lc3 <- apply(MEP_sc_sample,2,sum)<=5
  GMP_sc_lc3 <- apply(GMP_sc_sample,2,sum)<=5
  
  HSC_sc_sum <- as.integer(apply(HSC_sc_sample,2,sum))
  MEP_sc_sum <- as.integer(apply(MEP_sc_sample,2,sum))
  GMP_sc_sum <- as.integer(apply(GMP_sc_sample,2,sum))
  pseudo_sc_3 <- data.frame(genes,HSC_sc_sum,MEP_sc_sum,GMP_sc_sum)
  ##############################################################################
  
  #Join all replicates
  pseudo_sc <- full_join(pseudo_sc_1,pseudo_sc_2,by="genes",suffix=c("_A","_B"))
  pseudo_sc <- full_join(pseudo_sc,pseudo_sc_3,by="genes",suffix=c(" ","_C"))
  
  #Join low expressed genes in each replicate
  HSC_sc_lc <- genes[(HSC_sc_lc1+HSC_sc_lc2+HSC_sc_lc3)==3]
  MEP_sc_lc <- genes[(MEP_sc_lc1+MEP_sc_lc2+MEP_sc_lc3)==3]
  GMP_sc_lc <- genes[(GMP_sc_lc1+GMP_sc_lc2+GMP_sc_lc3)==3]
  
  #add suffix
  colnames(pseudo_sc) <- str_c(colnames(pseudo_sc),suffix)
  
  return(list(pseudo_sc=pseudo_sc,HSC_sc_lc=HSC_sc_lc,
              MEP_sc_lc=MEP_sc_lc,GMP_sc_lc=GMP_sc_lc))
}

pseudosc_rep_mean <- function(matrix,suffix){
  genes <- matrix$labels
  ##############################################################################
  #Function that generates pseudosc matrix from single cell data
  
  # INPUT VARIABLES
  #matrix <- single cell matrix
  #suffix <- condition (senior or mds)
  
  # OUTPUT VARIABLES
  # pseudo_sc <- pseudo single cell matrix
  # HSC_sc_lc <- low genes expressed genes (<5 counts in total)
  # MEP_sc_lc <- low genes expressed genes (<5 counts in total)
  # GMP_sc_lc <- low genes expressed genes (<5 counts in total)
  ##############################################################################
  
  #extract data for HSC, MEP, GMP and LMPP
  HSC_sc <- matrix[,str_detect(colnames(matrix),"^HSC")]
  MEP_sc <- matrix[,str_detect(colnames(matrix),"^MEP")]
  GMP_sc <- matrix[,str_detect(colnames(matrix),"^GMP")]
  
  #convert to numeric (just in case)
  HSC_sc <- apply(HSC_sc,1,as.numeric)
  MEP_sc <- apply(MEP_sc,1,as.numeric)
  GMP_sc <- apply(GMP_sc,1,as.numeric)
  
  # Replicate 1 ################################################################
  
  # select the saples
  HSC_sc_sample <- HSC_sc[sample(1:nrow(HSC_sc),0.2*nrow(HSC_sc),replace=F),]
  MEP_sc_sample <- MEP_sc[sample(1:nrow(MEP_sc),0.2*nrow(MEP_sc),replace=F),]
  GMP_sc_sample <- GMP_sc[sample(1:nrow(GMP_sc),0.2*nrow(GMP_sc),replace=F),]
  
  # analyse low expressed genes
  HSC_sc_lc1 <- apply(HSC_sc_sample,2,sum)<=5
  MEP_sc_lc1 <- apply(MEP_sc_sample,2,sum)<=5
  GMP_sc_lc1 <- apply(GMP_sc_sample,2,sum)<=5
  
  #calculate mean and generate pseudosingle cell matrix (multiply by 1e5 and convert to integer because
  # the analysis needs input as integers and we only care about proportions)
  HSC_sc_mean <- as.integer(apply(HSC_sc_sample,2,mean)*10000)
  MEP_sc_mean <- as.integer(apply(MEP_sc_sample,2,mean)*10000)
  GMP_sc_mean <- as.integer(apply(GMP_sc_sample,2,mean)*10000)
  pseudo_sc_1 <- data.frame(genes,HSC_sc_mean,MEP_sc_mean,GMP_sc_mean)
  
  # Replicate 2 ################################################################
  
  HSC_sc_sample <- HSC_sc[sample(1:nrow(HSC_sc),0.2*nrow(HSC_sc),replace=F),]
  MEP_sc_sample <- MEP_sc[sample(1:nrow(MEP_sc),0.2*nrow(MEP_sc),replace=F),]
  GMP_sc_sample <- GMP_sc[sample(1:nrow(GMP_sc),0.2*nrow(GMP_sc),replace=F),]
  
  HSC_sc_lc2 <- apply(HSC_sc_sample,2,sum)<=5
  MEP_sc_lc2 <- apply(MEP_sc_sample,2,sum)<=5
  GMP_sc_lc2 <- apply(GMP_sc_sample,2,sum)<=5
  
  HSC_sc_mean <- as.integer(apply(HSC_sc_sample,2,mean)*10000)
  MEP_sc_mean <- as.integer(apply(MEP_sc_sample,2,mean)*10000)
  GMP_sc_mean <- as.integer(apply(GMP_sc_sample,2,mean)*10000)
  pseudo_sc_2 <- data.frame(genes,HSC_sc_mean,MEP_sc_mean,GMP_sc_mean)
  
  # Replicate 3 ################################################################
  
  HSC_sc_sample <- HSC_sc[sample(1:nrow(HSC_sc),0.2*nrow(HSC_sc),replace=F),]
  MEP_sc_sample <- MEP_sc[sample(1:nrow(MEP_sc),0.2*nrow(MEP_sc),replace=F),]
  GMP_sc_sample <- GMP_sc[sample(1:nrow(GMP_sc),0.2*nrow(GMP_sc),replace=F),]
  
  HSC_sc_lc3 <- apply(HSC_sc_sample,2,sum)<=5
  MEP_sc_lc3 <- apply(MEP_sc_sample,2,sum)<=5
  GMP_sc_lc3 <- apply(GMP_sc_sample,2,sum)<=5
  
  HSC_sc_mean <- as.integer(apply(HSC_sc_sample,2,mean)*10000)
  MEP_sc_mean <- as.integer(apply(MEP_sc_sample,2,mean)*10000)
  GMP_sc_mean <- as.integer(apply(GMP_sc_sample,2,mean)*10000)
  pseudo_sc_3 <- data.frame(genes,HSC_sc_mean,MEP_sc_mean,GMP_sc_mean)
  ##############################################################################
  
  #Join all replicates
  pseudo_sc <- full_join(pseudo_sc_1,pseudo_sc_2,by="genes",suffix=c("_A","_B"))
  pseudo_sc <- full_join(pseudo_sc,pseudo_sc_3,by="genes",suffix=c(" ","_C"))
  
  #Join low expressed genes in each replicate
  HSC_sc_lc <- genes[(HSC_sc_lc1+HSC_sc_lc2+HSC_sc_lc3)==3]
  MEP_sc_lc <- genes[(MEP_sc_lc1+MEP_sc_lc2+MEP_sc_lc3)==3]
  GMP_sc_lc <- genes[(GMP_sc_lc1+GMP_sc_lc2+GMP_sc_lc3)==3]
  
  #add suffix
  colnames(pseudo_sc) <- str_c(colnames(pseudo_sc),suffix)
  
  return(list(pseudo_sc=pseudo_sc,HSC_sc_lc=HSC_sc_lc,
              MEP_sc_lc=MEP_sc_lc,GMP_sc_lc=GMP_sc_lc))
}

keep_protein_coding <- function(matrix){
  ##############################################################################
  #Function that keeps only the protein coding genes in the sc count matrix
  
  # INPUT VARIABLES
  #matrix <- sc count matrix
  
  # OUTPUT VARIABLES
  # matrix <- count matrix with only protein coding genes
  ##############################################################################
  
  gene_types = mapIds(org.Hs.eg.db, 
                      keys=matrix$labels, 
                      column="GENETYPE",
                      keytype="SYMBOL",
                      multiVals="first") #obtain the gene type
  gene_types[is.na(gene_types)] <- "NAN" #replace NA by "NAN"
  keep <- gene_types=="protein-coding" #extract indexfor protein coding genes
  matrix <- matrix[keep,] #only keep protein coding genes
  return(matrix)
}

keep_protein_coding_bulk <- function(matrix){
  ##############################################################################
  #Function that keeps only the protein coding genes in the FACS-bulk count matrix
  
  # INPUT VARIABLES
  #matrix <- FACS-bulk count matrix
  
  # OUTPUT VARIABLES
  # matrix <- FACS-bulk matrix with only protein coding genes
  ##############################################################################
  gene_types = mapIds(org.Hs.eg.db, 
                      keys=matrix$genes, 
                      column="GENETYPE",
                      keytype="ENSEMBL",
                      multiVals="first")
  gene_types[is.na(gene_types)] <- "NAN"
  keep <- gene_types=="protein-coding"
  matrix <- matrix[keep,]
  return(matrix)
}

keep_cells <- function(matrix,name){
  ##############################################################################
  # Function that saves the single cell matrix only containing cells of interest
  # with the corresponding labels
  
  # INPUT VARIABLES
  # matrix <- sc counts matrix
  # name <- identifier of the matrix
  
  #####################################################################################
  genes <- matrix$labels #extract the gene names
  
  #extract data for HSC, MEP, GMP and LMPP
  HSC_sc <- matrix[,str_detect(colnames(matrix),"^HSC")] 
  MEP_sc <- matrix[,str_detect(colnames(matrix),"^MEP")]
  GMP_sc <- matrix[,str_detect(colnames(matrix),"^GMP")]
  sc <- data.frame(HSC_sc,MEP_sc,GMP_sc)
  colnames(sc) <- c(rep("HSC",ncol(HSC_sc)),rep("MEP",ncol(MEP_sc)),rep("GMP",ncol(GMP_sc)))
  rownames(sc) <- genes
  write.table(sc,file = (paste0(wd,single_cell,name,".csv")))
}


make_pseudo <- function(matrix,name,suffix){
  ##############################################################################
  # Function that makes the pseudobulk matrices using the mean, sum, mean with 
  # replicas and sum with replicas and saves them
  
  # INPUT VARIABLES
  # matrix <- sc counts matrix
  # name <- identifier of the matrix
  # suffix <- suffix indicating the condition ("_mds" or "_senior")
  ##############################################################################
  
  #generate the pseudobulk matrices with the corresponding function
  ps_mean <- pseudosc_mean(matrix,suffix)
  ps_sum <- pseudosc_sum(matrix,suffix)
  ps_rep_mean <- pseudosc_rep_mean(matrix,suffix)
  ps_rep_sum <- pseudosc_rep_sum(matrix,suffix)
  
  #save
  saveRDS(ps_mean,paste0(wd,pseudobulk,"pseudo_",name,"_mean.rds")) 
  saveRDS(ps_sum,paste0(wd,pseudobulk,"pseudo_",name,"_sum.rds")) 
  saveRDS(ps_rep_mean,paste0(wd,pseudobulk,"pseudo_",name,"_rep_mean.rds")) 
  saveRDS(ps_rep_sum,paste0(wd,pseudobulk,"pseudo_",name,"_rep_sum.rds")) 
}
preprocess_mds <- function(path_to_raw,name,suffix){
  ##############################################################################
  # Read single cell data and calculate the pseudosc and number of cells in each cell type (HSC, MEP, GMP)
  # in mds data
  
  # INPUT VARIABLES
  # path_to_raw <- path to the rwcounts matrix
  # name <- identifier of the matrix
  # suffix <- suffix indicating the condition ("_mds" or "_senior")
  ##############################################################################
  rawcounts <- read_delim(path_to_raw,delim = "\t", escape_double = FALSE,trim_ws = TRUE, skip = 1) #read raw counts
  rawcounts <-  keep_protein_coding(rawcounts) #keep only protein coding genes
  keep_cells(rawcounts,name) 
  pseudo_mds1 <- make_pseudo(rawcounts,name,suffix) #make the pseudobulk matrices
  cell_num(rawcounts,name) #calculate the number of each cell type and save 
}

preprocess_senior <- function(path_to_raw,path_to_labels,name,suffix){
  ##############################################################################
  # Read single cell data and calculate the pseudosc and number of cells in each cell type (HSC, MEP, GMP)
  # in senior data
  
  # INPUT VARIABLES
  # path_to_raw <- path to the rawcounts matrix
  # path_to_labels <- path to the cell labels
  # name <- identifier of the matrix
  # suffix <- suffix indicating the condition ("_mds" or "_senior")
  ##############################################################################
  senior <- read_rds(path_to_raw)
  senior_labels <- read_delim(path_to_labels,
                              delim = "\t", escape_double = FALSE,
                              col_names = FALSE, trim_ws = TRUE, skip = 1)
  rawcounts_senior <- assign_labels(senior,senior_labels)
  rawcounts_senior <- keep_protein_coding(rawcounts_senior)
  keep_cells(rawcounts_senior,name) 
  pseudo_senior <- make_pseudo(rawcounts_senior,name,suffix)
  cell_num(rawcounts_senior,name)
}

join_pseudo <- function(mds1,mds3,mds5,mds10,senior1,senior2,senior3,name){
  ##############################################################################
  # Function that makes the pseudobulk matrices using the mean, sum, mean with 
  # replicas and sum with replicas and saves them
  
  # INPUT VARIABLES
  # mds1,mds3,mds5,mds10,senior1,senior2,senior3 <- pseudobulk matrices
  # name <- mathod used to calculate the pseudobulk matrices ("sum","mean","rep_sum","rep_mean")
  ##############################################################################
  
  pseudo_all <- inner_join(mds1$pseudo_sc,mds3$pseudo_sc,by=c("genes_mds"),suffix=c("_1","_3")) 
  pseudo_all <- inner_join(pseudo_all,mds5$pseudo_sc,by=c("genes_mds"),suffix=c("","_5")) 
  pseudo_all <- inner_join(pseudo_all,mds10$pseudo_sc,by=c("genes_mds"),suffix=c("","_10")) 
  pseudo_all <- inner_join(pseudo_all,senior1$pseudo_sc,by=c("genes_mds"="genes_senior"),suffix=c("","_1")) 
  pseudo_all <- inner_join(pseudo_all,senior2$pseudo_sc,by=c("genes_mds"="genes_senior"),suffix=c("","_2")) 
  pseudo_all <- inner_join(pseudo_all,senior3$pseudo_sc,by=c("genes_mds"="genes_senior"),suffix=c("","_3")) 
  
  pseudo_all <- as.data.frame(pseudo_all) # convert sc matrix to a data frame
  pseudo_all_genes <- pseudo_all$genes_mds #save the original gene names
  
  pseudo_all$genes= mapIds(org.Hs.eg.db, #obtain the ensemble id
                           keys=pseudo_all$genes_mds, 
                           column="ENSEMBL",
                           keytype="SYMBOL",
                           multiVals="first")
  pseudo_all$genes[is.na(pseudo_all$genes_mds)] <- pseudo_all_genes[is.na(pseudo_all$genes)] #fill the genes without ensemble id with the gene symbol
  pseudo_all$genes_mds <- NULL #delete column
  saveRDS(pseudo_all,paste0(wd,pseudobulk,"pseudo_all_",name,".rds")) #save for further analysis 
  
  zero_genes <- list(
    HSC_mds_lc = unique(c(mds1$HSC_sc_lc,mds3$HSC_sc_lc,mds5$HSC_sc_lc,mds10$HSC_sc_lc)),
    MEP_mds_lc = unique(c(mds1$MEP_sc_lc,mds3$MEP_sc_lc,mds5$MEP_sc_lc,mds10$MEP_sc_lc)),
    GMP_mds_lc = unique(c(mds1$GMP_sc_lc,mds3$GMP_sc_lc,mds5$GMP_sc_lc,mds10$GMP_sc_lc)),
    
    HSC_senior_lc = unique(c(senior1$HSC_sc_lc,senior2$HSC_sc_lc,senior3$HSC_sc_lc)),
    MEP_senior_lc = unique(c(senior1$MEP_sc_lc,senior2$MEP_sc_lc,senior3$MEP_sc_lc)),
    GMP_senior_lc = unique(c(senior1$GMP_sc_lc,senior2$GMP_sc_lc,senior3$GMP_sc_lc))
  )
  saveRDS(zero_genes,paste0(wd,pseudobulk,"zero_genes_",name,".rds")) #save for further analysis
}

# Select the samples corresponding to each cell type and condition
select_samples <- function(dataframe,regex){ 
  #####################################################################################
  #Function that selects all samples from each cell type/condition
  
  # INPUT VARIABLES
  #dataframe <- dataframe with counts
  #regex <- regex expression to identify the labels
  
  # OUTPUT VARIABLES
  # subset <- matrix with counts corresponding to all the samples that correspond to the regex label
  #####################################################################################
  subset <- as.matrix(dataframe[,str_detect(colnames(dataframe),regex)])
  rownames(subset) <- dataframe[,"genes"] 
  return(subset)
}


DEA <- function(matrix_1,matrix_2,condition_1,condition_2,name_exp){
  #####################################################################################
  #Function that calculates the differential expression analysis between two conditions
  
  # INPUT VARIABLES
  #matrix_1 <- reference condition matrix with counts
  #matrix_2 <- alternative condition matrix with counts
  #condition 1 <- reference condition
  #condition_2 <- alternative condition
  #name_exp <- name of the experiment
  
  # OUTPUT VARIABLES
  #_results_filtered <- results filtered by padj<0.05
  #_results_upregulated <- upregulated genes with their logFC
  #_results_downregulated <- downregulated genes with their logFC
  #norm <- normalized counts of the original matrix
  #sparsity <- sparsity of each matrix
  #PCAplots <- PCA plots before and after normalization
  #####################################################################################
  
  #Bind both matrices and generate the design matrix
  matrix <- cbind(matrix_1,matrix_2)  
  design_matrix <- data.frame(Condition = factor(c(rep(condition_1,dim(matrix_1)[2]), rep(condition_2,dim(matrix_2)[2]))))
  
  # Assign sample names to the desing matrix
  rownames(design_matrix) <- colnames(matrix)
  
  # Set 'Condition_1' as reference 
  design_matrix$Condition <- relevel(design_matrix$Condition, ref = condition_1)
  
  #check that rownames in design matrix and colnames in count matrix are equal
  stopifnot(all(rownames(design_matrix) %in% colnames(matrix)),
            all(rownames(design_matrix) == colnames(matrix)))
  
  #eliminate low count genes
  matrix <- matrix[apply(matrix,1,max)>50,]
  
  #generate the dds object
  dds <- DESeqDataSetFromMatrix(countData = matrix, 
                                colData = design_matrix,
                                design = ~ Condition) 
  
  # run differential expressin analysis with DESEQ
  dds <- DESeq(dds)
  
  normalized_counts <- counts(dds,normalized=TRUE) #save normalize counts
  
  results_dds <- results(dds) #save results
  
  #write.table(results_dds,
  #            file = paste0(wd,"/results/",name_exp,"_results.txt"), quote = F,
  #           sep = "\t", row.names = F, col.names = T)
  
  #we will apply FDR cutoff=0.05 for padj value
  results_filtered <- results(dds, alpha = 0.05) #FDR cutoff for padj value
  
  #filter genes by padj < 0.05 and generate dataframe with the GENE (gene name) and LFC
  results_filtered <- results_filtered[which(results_dds$padj < 0.05),]
  filtered <- data.frame("GENE"=rownames(results_filtered), 
                         "LFC"=results_filtered[,2])
  
  # obtain upregulated gene subset (LFC>0) and generate dataframe with the GENE (gene name) and LFC
  upreg <- subset(results_filtered, log2FoldChange > 0)
  results_upregulated <- data.frame("GENE"=rownames(upreg), 
                                    "LFC"=upreg[,2])
  
  # obtain upregulated gene subset (LFC<0) and generate dataframe with the GENE (gene name) and LFC
  downreg <- subset(results_filtered, log2FoldChange < 0)
  results_downregulated <- data.frame("GENE"=rownames(downreg), 
                                      "LFC"=downreg[,2])
  
  #calculate sparsity of each matrix 
  sparsity <- paste0(condition_1,"-> ",as.character(round(sparsity(matrix_1),3)),"  ",condition_2,"-> ",as.character(round(sparsity(matrix_2),3)))
  
  # Generate the PCA plots of raw and nromalized counts
  pr.out_raw <- prcomp(t(matrix)) #principal component analysis of the raw counts
  pr.res_raw <- as.data.frame(pr.out_raw$x) #convert the component matrix into a dataframe
  pr.res_raw$cell <- factor(c(rep(condition_1,dim(matrix_1)[2]), rep(condition_2,dim(matrix_2)[2]))) #assing the cell names
  plot1 <- pr.res_raw %>% ggplot(aes(x=PC1,y=PC2,col=cell))+ 
    geom_point()+
    labs(title=paste0("Raw counts ",name_exp," ",as.character(nrow(pr.res_raw))),col="CELL TYPE")+
    scale_color_manual(values=c("#ff9dc5","#82dde0"))+
    theme_light()+
    theme(plot.title = element_text(size=12)) # plot the PCA
  
  #perfom the same procedure with the normalized data
  pr.out_norm <- prcomp(t(normalized_counts)) 
  pr.res_norm <- as.data.frame(pr.out_norm$x)
  pr.res_norm$cell <- factor(c(rep(condition_1,dim(matrix_1)[2]), rep(condition_2,dim(matrix_2)[2])))
  plot2 <- pr.res_raw %>% ggplot(aes(x=PC1,y=PC2,col=cell))+
    geom_point()+
    labs(title=paste0("Normalized counts ",name_exp," ",as.character(nrow(pr.res_norm))),col="CELL TYPE")+
    scale_color_manual(values=c("#ff9dc5","#82dde0"))+
    theme_light()+
    theme(plot.title = element_text(size=12))
  
  PCAplots <- cowplot::plot_grid(plot1,plot2) #join both plots
  
  return(list(results=results_dds,filtered=filtered,
              upregulated=results_upregulated,
              downregulated=results_downregulated,
              sparsity=sparsity,
              norm=normalized_counts,
              PCAplots=PCAplots
  ))
}


common_results <- function(bulk,sc,only_bulk_genes,zero_genes1,zero_genes2,markersA,markersB){
###################################################################################  
  #Function that compares the results from two DEA results
  
  # INPUT VARIABLES
  # bulk <- DEA results from comparing two bulk samples
  # sc <- DEA results from comparing two sc samples
  # only bulk genes <- gene list with genes sequenced by bulk bulk but not by sc
  # zero genes1 <-  gene list with genes that are low expressed in the cell type 1
  # zero genes2 <-  gene list with genes that are low expressed in the cell type 2
  
  # OUTPUT VARIABLES
  # fused_gene_list <- fused list of significant gene sin both DEA and its corresponding LFC
  # DEG_num_comparison <- tibble with the number of DEG in each experiment
  # plot_DEG_num_comparison <- plot of DEG in each experiment
  # title
  # LFC_graph=LFC_graph
  # LFC_table=LFC_table
  # bulk_sparsity=bulk_sparsity
  # sc_sparsity=sc_sparsity
  # bulk_normcounts=bulk_normcounts
  # sc_normcounts=sc_normcounts
  # sc_PCAplots=sc_PCAplots
  # bulk_PCAplots=bulk_PCAplots
  # CC_enrichment_common=p1
  # BP_enrichment_common=p2
  # MF_enrichment_common=p3
  # CC_enrichment_sc=p4
  # BP_enrichment_sc=p5
  # MF_enrichment_sc=p6
  # CC_enrichment_bulk=p7
  # BP_enrichment_bulk=p8
  # MF_enrichment_bulk=p9))
################################################################################### 
  markers <- mapIds(org.Hs.eg.db, #gene types only in bulk RNAseq
                    keys=c(markersA,markersB), 
                    column="ENSEMBL",
                    keytype="SYMBOL",
                    multiVals="first")
  
  bulk_sparsity <- bulk$sparsity
  sc_sparsity <- sc$sparsity
  bulk_normcounts=bulk$norm
  sc_normcounts=sc$norm
  sc_PCAplots=sc$PCAplots
  bulk_PCAplots=bulk$PCAplots
  bulk <- bulk$filtered
  sc <- sc$filtered
  
  # Fuse the gene lists of both results
  fused_gene_list <- full_join(bulk,sc,by=c("GENE"),suffix=c("Bulk","Sc"))
  common_genes <- intersect(bulk$GENE,sc$GENE)
  common_n <- length(common_genes)
  common_markers <- length(intersect(common_genes,markers))
  common_LFC <- fused_gene_list[complete.cases(fused_gene_list),]
  
  bulk_genes_unique=fused_gene_list$GENE[is.na(fused_gene_list$LFCSc)] #genes only detected in the bulk RNAseq
  bulk_genes_LFC=fused_gene_list$LFCBulk[is.na(fused_gene_list$LFCSc)]
  only_bulk_n=length(intersect(bulk_genes_unique,only_bulk_genes))
  
  sc_genes_unique=fused_gene_list$GENE[is.na(fused_gene_list$LFCBulk)] #genes only detected in the sc RNAseq
  sc_genes_LFC=fused_gene_list$LFCSc[is.na(fused_gene_list$LFCBulk)]
  
 LFC_table <- tibble(Condition = factor(c(rep("common genes sc",common_n),rep("common genes bulk",common_n),rep("scRNAseq unique",length(sc_genes_LFC)),rep("FACS-bulk RNAseq unique",length(bulk_genes_LFC)))),
        LFC=abs(c(common_LFC$LFCSc,common_LFC$LFCBulk,sc_genes_LFC,bulk_genes_LFC)))
  
 LFC_graph <- LFC_table %>% ggplot(aes(x=LFC,col=Condition,fill=Condition)) +
    geom_density(kernel = "gaussian",alpha=0.3)+
    labs(y="Density",x="Absolute value of LFC",fill="")+
    scale_color_manual(values=c("#fae980","#c5b5fb","#82dde0","#ff9dc5"))+
    scale_fill_manual(values=c("#fae980","#c5b5fb","#82dde0","#ff9dc5"))+
    guides(col="none")+
    theme_light()
  
  zero_genes <- mapIds(org.Hs.eg.db, #gene types only in bulk RNAseq
                       keys=unique(c(zero_genes1,zero_genes2)), 
                       column="ENSEMBL",
                       keytype="SYMBOL",
                       multiVals="first")
  zero_genes_only_sc <- length(intersect(sc_genes_unique,zero_genes))
  
  # Generate a tibble and a plot with the gene proportions detected by each experiment
  DEG_num_comparison <- tibble(sample=c("All","Bulk","Sc"),
                               common_markers=c(common_markers,common_markers,common_markers),
                               common=c(common_n,common_n,common_n),
                               unique_genes=c(nrow(fused_gene_list),(nrow(bulk)-only_bulk_n),(nrow(sc)-zero_genes_only_sc)),
                               only_bulk=c(0,nrow(bulk),0),
                               low_expressed=c(0,0,nrow(sc)))
  title <- paste0("BULK_",bulk_sparsity,"
                    PSEUDO_SC_",sc_sparsity)                           
  
  bulk_genes_all <- bulk$GENE #select all significant genes detected by bulk
  sc_genes_all <- sc$GENE #select all significant genes detected by sc
  
  #GO analysis of common genes
  ENTREZID = mapIds(org.Hs.eg.db, #obtain the ENTREZ ID 
                    keys=common_genes, 
                    column="ENTREZID",
                    keytype="ENSEMBL",
                    multiVals="first") 
  genes = as.character(ENTREZID) #set as character
  # Obtain GO terms with cluster profiler
  ggo_cc <- clusterProfiler::groupGO(gene     = genes,
                                     OrgDb    = org.Hs.eg.db,
                                     ont      = "CC",   
                                     level    = 3,
                                     readable = TRUE)
  # GO enrichment analysis
  ego <- clusterProfiler::enrichGO(gene          = genes,
                                   OrgDb         = org.Hs.eg.db,
                                   ont           = "CC",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.01,
                                   qvalueCutoff  = 0.01, 
                                   readable      = TRUE)
  p1 <- ego@result$ID
  # Obtain GO terms with cluster profiler
  ggo_cc <- clusterProfiler::groupGO(gene     = genes,
                                     OrgDb    = org.Hs.eg.db,
                                     ont      = "BP",   
                                     level    = 3,
                                     readable = TRUE)
  # GO enrichment analysis
  ego <- clusterProfiler::enrichGO(gene          = genes,
                                   OrgDb         = org.Hs.eg.db,
                                   ont           = "BP",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.01,
                                   qvalueCutoff  = 0.01, 
                                   readable      = TRUE)
  p2 <- ego@result$ID
  # Obtain GO terms with cluster profiler
  ggo_cc <- clusterProfiler::groupGO(gene     = genes,
                                     OrgDb    = org.Hs.eg.db,
                                     ont      = "MF",     #MF , BP CC
                                     level    = 3,
                                     readable = TRUE)
  # GO enrichment analysis
  ego <- clusterProfiler::enrichGO(gene          = genes,
                                   OrgDb         = org.Hs.eg.db,
                                   ont           = "MF",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.01,
                                   qvalueCutoff  = 0.01, 
                                   readable      = TRUE)
  p3 <- ego@result$ID
  
  #GO analysis of sc genes
  ENTREZID = mapIds(org.Hs.eg.db,
                    keys=sc_genes_all, 
                    column="ENTREZID",
                    keytype="ENSEMBL",
                    multiVals="first") 
  genes = as.character(ENTREZID) 
  # Obtain GO terms with cluster profiler
  ggo_cc <- clusterProfiler::groupGO(gene     = genes,
                                     OrgDb    = org.Hs.eg.db,
                                     ont      = "CC",     #MF , BP CC
                                     level    = 3,
                                     readable = TRUE)
  # GO enrichment analysis
  ego <- clusterProfiler::enrichGO(gene          = genes,
                                   OrgDb         = org.Hs.eg.db,
                                   ont           = "CC",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.01,
                                   qvalueCutoff  = 0.01, 
                                   readable      = TRUE)
  p4 <- ego@result$ID
  ggo_cc <- clusterProfiler::groupGO(gene     = genes,
                                     OrgDb    = org.Hs.eg.db,
                                     ont      = "BP",     #MF , BP CC
                                     level    = 3,
                                     readable = TRUE)
  # GO enrichment analysis
  ego <- clusterProfiler::enrichGO(gene          = genes,
                                   OrgDb         = org.Hs.eg.db,
                                   ont           = "BP",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.01,
                                   qvalueCutoff  = 0.01, 
                                   readable      = TRUE)
  p5 <- ego@result$ID
  
  # Obtain GO terms with cluster profiler
  ggo_cc <- clusterProfiler::groupGO(gene     = genes,
                                     OrgDb    = org.Hs.eg.db,
                                     ont      = "MF",     #MF , BP CC
                                     level    = 3,
                                     readable = TRUE)
  # GO enrichment analysis
  ego <- clusterProfiler::enrichGO(gene          = genes,
                                   OrgDb         = org.Hs.eg.db,
                                   ont           = "MF",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.01,
                                   qvalueCutoff  = 0.01, 
                                   readable      = TRUE)
  p6 <- ego@result$ID
  
  #GO analysis of bulk genes
  ENTREZID = mapIds(org.Hs.eg.db, #obtain the ENTREZ ID 
                    keys=bulk_genes_all, 
                    column="ENTREZID",
                    keytype="ENSEMBL",
                    multiVals="first") 
  genes = as.character(ENTREZID) 
  # Obtain GO terms with cluster profiler
  ggo_cc <- clusterProfiler::groupGO(gene     = genes,
                                     OrgDb    = org.Hs.eg.db,
                                     ont      = "CC",     #MF , BP CC
                                     level    = 3,
                                     readable = TRUE)
  # GO enrichment analysis
  ego <- clusterProfiler::enrichGO(gene          = genes,
                                   OrgDb         = org.Hs.eg.db,
                                   ont           = "CC",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.01,
                                   qvalueCutoff  = 0.01, 
                                   readable      = TRUE)
  p7 <- ego@result$ID
  # Obtain GO terms with cluster profiler
  ggo_cc <- clusterProfiler::groupGO(gene     = genes,
                                     OrgDb    = org.Hs.eg.db,
                                     ont      = "BP",     #MF , BP CC
                                     level    = 3,
                                     readable = TRUE)
  # GO enrichment analysis
  ego <- clusterProfiler::enrichGO(gene          = genes,
                                   OrgDb         = org.Hs.eg.db,
                                   ont           = "BP",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.01,
                                   qvalueCutoff  = 0.01, 
                                   readable      = TRUE)
  p8 <- ego@result$ID
  # Obtain GO terms with cluster profiler
  ggo_cc <- clusterProfiler::groupGO(gene     = genes,
                                     OrgDb    = org.Hs.eg.db,
                                     ont      = "MF",     #MF , BP CC
                                     level    = 3,
                                     readable = TRUE)
  # GO enrichment analysis
  ego <- clusterProfiler::enrichGO(gene          = genes,
                                   OrgDb         = org.Hs.eg.db,
                                   ont           = "MF",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.01,
                                   qvalueCutoff  = 0.01, 
                                   readable      = TRUE)
  p9 <- ego@result$ID
  
  return(list(fused_gene_list=fused_gene_list,
              DEG_num_comparison=DEG_num_comparison,
              title=title,
              LFC_graph=LFC_graph,
              LFC_table=LFC_table,
              bulk_sparsity=bulk_sparsity,
              sc_sparsity=sc_sparsity,
              bulk_normcounts=bulk_normcounts,
              sc_normcounts=sc_normcounts,
              sc_PCAplots=sc_PCAplots,
              bulk_PCAplots=bulk_PCAplots,
              CC_enrichment_common=p1,
              BP_enrichment_common=p2,
              MF_enrichment_common=p3,
              CC_enrichment_sc=p4,
              BP_enrichment_sc=p5,
              MF_enrichment_sc=p6,
              CC_enrichment_bulk=p7,
              BP_enrichment_bulk=p8,
              MF_enrichment_bulk=p9))
}

GO_comparison <- function(pseudotype,cells,sample){
  ###################################################################################  
  #Function the GO term enrichment result
  
  # INPUT VARIABLES
  # pseudotype <- "pseudobulk (mean)" or "pseudobulk (sum)
  # cells <-  "cell1_cell2_condition"
  
  # OUTPUT VARIABLES
  # GO enrichment comparison table
  ################################################################################### 

  common_CC <- intersect(sample$CC_enrichment_sc,sample$CC_enrichment_bulk)
  n_common_CC <- length(common_CC)
  n_bulk_CC <- length(sample$CC_enrichment_bulk)
  n_sc_CC <- length(sample$CC_enrichment_sc)
  common_BP <- intersect(sample$BP_enrichment_sc,sample$BP_enrichment_bulk)
  n_common_BP <- length(common_BP)
  n_bulk_BP <- length(sample$BP_enrichment_bulk)
  n_sc_BP <- length(sample$BP_enrichment_sc)
  common_MF <- intersect(sample$MF_enrichment_sc,sample$MF_enrichment_bulk)
  n_common_MF <- length(common_MF)
  n_bulk_MF <- length(sample$MF_enrichment_bulk)
  n_sc_MF <- length(sample$MF_enrichment_sc)
  
  GO_num_comparison <- tibble(pseudo=rep(pseudotype,6),
                              celltypes=rep(cells,6),
                              experiment=rep(c("FACS-bulk DEGS","sc DEGS"),3),
                              GO=c("CC","CC","BP","BP","MF","MF"),
                              common=c(n_common_CC,n_common_CC,n_common_BP,n_common_BP,n_common_MF,n_common_MF),
                              total=c(n_bulk_CC,n_sc_CC,n_bulk_BP,n_sc_BP,n_bulk_MF,n_sc_MF))
  return(GO_num_comparison)}

marker_genes <- function(bulk,sc,markerA,markerB){
  #####################################################################################
  #Function that compares the results from two DEA resultswith respect to marker genes
  
  # INPUT VARIABLES
  # bulk <- DEA results from comparing two bulk samples
  # sc <- DEA results from comparing two sc samples
  # markerA <- marker genes of cell type 1
  # markerB <- marker genes of cell type 2
  
  # OUTPUT VARIABLES
  # bulk_A <- A marker genes detected by bulk
  # bulk_B <- B marker genes detected by bulk
  # sc_A <- A marker genes detected by sc
  # sc_B <- B marker genes detected by sc
  # percentages <-  percentage of the marker genes detected by each bulk and sc
  #####################################################################################
  
  bulk_down <- mapIds(org.Hs.eg.db, #gene types only in bulk RNAseq
                      keys=bulk$downregulated$GENE, 
                      column="SYMBOL",
                      keytype="ENSEMBL",
                      multiVals="first")
  bulk_A_perc <- mean(markerA %in% bulk_down)
  bulk_A <- markerA[markerA %in% bulk_down]
  
  bulk_up <- mapIds(org.Hs.eg.db, #gene types only in bulk RNAseq
                    keys=bulk$upregulated$GENE, 
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first")
  bulk_B_perc <- mean(markerB %in% bulk_up)
  bulk_B <- markerB[markerB %in% bulk_up]
  
  sc_down <- mapIds(org.Hs.eg.db, #gene types only in bulk RNAseq
                    keys=sc$downregulated$GENE, 
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first")
  sc_A_perc <- mean(markerA %in% sc_down)
  sc_A <- markerA[markerA %in% sc_down]
  
  sc_up <- mapIds(org.Hs.eg.db, #gene types only in bulk RNAseq
                  keys=sc$upregulated$GENE, 
                  column="SYMBOL",
                  keytype="ENSEMBL",
                  multiVals="first")
  sc_B_perc <- mean(markerB %in% sc_up)
  sc_B <- markerB[markerB %in% sc_up]
  
  percentages <- tibble(bulk=c(bulk_A_perc,bulk_B_perc),
                        sc=c(sc_A_perc,sc_B_perc))
  
  return(list(percentages=percentages,bulk_A=bulk_A,
              bulk_B=bulk_B,sc_A=sc_A,sc_B=sc_B))}

analyse_sc_bulk <- function(sc,bulk,zero_genes){
  ##############################################################################
  # Function that selects the samples, performs the DEA analysis between each pair of cells
  # and compares the results
  
  # INPUT VARIABLES
  # bulk <- bulk matrix
  # sc <- pseudobulk sc matrix
  # zero_genes <- list with the genes that are zero for all cells in each of the pseudobulk matrices
  
  # OUTPUT VARIABLES
  # cell1_cell2_bulkvssc_condition <- result of comparison between the DEA analysis of bulk vs pseudobulk for those two cell types
  # HSC_MEP_bulkvssc_mds_markers <- result of the comparison at the level of marker genes
  # cell_experiment_condition <- matrix with the samples corresponding to this configuration
  ##############################################################################
  
  # Join all the data and select only common genes in both datasets ###############
  bulk_sc <- full_join(bulk,sc,by=c("genes")) #join bulk and sc data by genes
  bulk_sc[is.na(bulk_sc)] <- 0
  only_bulk_genes <- bulk$genes[as.logical(abs((bulk$genes %in% sc$genes)-1))]
  
  # Select samples
  HSC_sc_mds <- select_samples(bulk_sc,"^HSC.*sc.*mds")
  MEP_sc_mds <- select_samples(bulk_sc,"^MEP.*sc.*mds")
  GMP_sc_mds <- select_samples(bulk_sc,"^GMP.*sc.*mds")
  
  HSC_sc_senior <- select_samples(bulk_sc,"^HSC.*sc.*senior")
  MEP_sc_senior <- select_samples(bulk_sc,"^MEP.*sc.*senior")
  GMP_sc_senior <- select_samples(bulk_sc,"^GMP.*sc.*senior")
  
  HSC_bulk_mds <- select_samples(bulk_sc,"^SMD.*HSC")
  MEP_bulk_mds <- select_samples(bulk_sc,"^SMD.*MEP")
  GMP_bulk_mds <- select_samples(bulk_sc,"^SMD.*GMP")
  
  HSC_bulk_old <- select_samples(bulk_sc,"MO.*HSC")
  MEP_bulk_old <- select_samples(bulk_sc,"MO.*MEP")
  GMP_bulk_old <- select_samples(bulk_sc,"MO.*GMP")
  
  # SPECIFY MARKERS
  
  HSC_markers <- c("CRHBP","HOPX","KYT","CD34","AVP","PRSS2","MLLT3","IDS","BTS2")
  GMP_markers <- c("CSF3R","CTSG","PRTN3","MPO","CFD","CSTA","CST7")
  MEP_markers <- c("NFE2","HFS1","TAL1","FCER1A","PBX1","PLEK","DAD1","IGSF10")
  
  ####################  DIFFERENTIAL EXPRESSION ANALYSIS  ########################
  ################################################################################
  
  #DIFFERENTIAL EXPRESSION ANALYSIS bulk_vs_bulk vs pseudosc_vs_pseudosc #########
  
  # MDS
  bulk_HSCvsMEP_mds <- DEA(HSC_bulk_mds,MEP_bulk_mds,"HSC","MEP","bulk_HSCvsMEP")
  sc_HSCvsMEP_mds <- DEA(HSC_sc_mds,MEP_sc_mds,"HSC","MEP","sc_HSCvsMEP")
  HSC_MEP_bulkvssc_mds <- common_results(bulk_HSCvsMEP_mds,sc_HSCvsMEP_mds,only_bulk_genes,zero_genes$HSC_mds_lc,zero_genes$MEP_mds_lc,HSC_markers,MEP_markers)
  HSC_MEP_bulkvssc_mds_markers <- marker_genes(bulk_HSCvsMEP_mds,sc_HSCvsMEP_mds,HSC_markers,MEP_markers)
  
  bulk_HSCvsGMP_mds <- DEA(HSC_bulk_mds,GMP_bulk_mds,"HSC","GMP","bulk_HSCvsGMP")
  sc_HSCvsGMP_mds <- DEA(HSC_sc_mds,GMP_sc_mds,"HSC","GMP","sc_HSCvsGMP")
  HSC_GMP_bulkvssc_mds <- common_results(bulk_HSCvsGMP_mds,sc_HSCvsGMP_mds,only_bulk_genes,zero_genes$HSC_mds_lc,zero_genes$GMP_mds_lc,HSC_markers,GMP_markers)
  HSC_GMP_bulkvssc_mds_markers <- marker_genes(bulk_HSCvsGMP_mds,sc_HSCvsGMP_mds,HSC_markers,GMP_markers)
  
  bulk_MEPvsGMP_mds <- DEA(MEP_bulk_mds,GMP_bulk_mds,"MEP","GMP","bulk_MEPvsGMP")
  sc_MEPvsGMP_mds <- DEA(MEP_sc_mds,GMP_sc_mds,"MEP","GMP","sc_MEPvsGMP")
  MEP_GMP_bulkvssc_mds <- common_results(bulk_MEPvsGMP_mds,sc_MEPvsGMP_mds,only_bulk_genes,zero_genes$MEP_mds_lc,zero_genes$GMP_mds_lc,MEP_markers,GMP_markers)
  MEP_GMP_bulkvssc_mds_markers <- marker_genes(bulk_MEPvsGMP_mds,sc_MEPvsGMP_mds,MEP_markers,GMP_markers)
  
  
  # Senior
  bulk_HSCvsMEP_senior <- DEA(HSC_bulk_old,MEP_bulk_old,"HSC","MEP","bulk_HSCvsMEP")
  sc_HSCvsMEP_senior <- DEA(HSC_sc_senior,MEP_sc_senior,"HSC","MEP","sc_HSCvsMEP")
  HSC_MEP_bulkvssc_senior <- common_results(bulk_HSCvsMEP_senior,sc_HSCvsMEP_senior,only_bulk_genes,zero_genes$HSC_senior_lc,zero_genes$MEP_senior_lc,HSC_markers,MEP_markers)
  HSC_MEP_bulkvssc_senior_markers <- marker_genes(bulk_HSCvsMEP_senior,sc_HSCvsMEP_senior,HSC_markers,MEP_markers)
  
  bulk_HSCvsGMP_senior <- DEA(HSC_bulk_old,GMP_bulk_old,"HSC","GMP","bulk_HSCvsGMP")
  sc_HSCvsGMP_senior <- DEA(HSC_sc_senior,GMP_sc_senior,"HSC","GMP","sc_HSCvsGMP")
  HSC_GMP_bulkvssc_senior <- common_results(bulk_HSCvsGMP_senior,sc_HSCvsGMP_senior,only_bulk_genes,zero_genes$HSC_senior_lc,zero_genes$GMP_senior_lc,HSC_markers,GMP_markers)
  HSC_GMP_bulkvssc_senior_markers <- marker_genes(bulk_HSCvsGMP_senior,sc_HSCvsGMP_senior,HSC_markers,GMP_markers)
  
  bulk_MEPvsGMP_senior <- DEA(MEP_bulk_old,GMP_bulk_old,"MEP","GMP","bulk_MEPvsGMP")
  sc_MEPvsGMP_senior <- DEA(MEP_sc_senior,GMP_sc_senior,"MEP","GMP","sc_MEPvsGMP")
  MEP_GMP_bulkvssc_senior <- common_results(bulk_MEPvsGMP_senior,sc_MEPvsGMP_senior,only_bulk_genes,zero_genes$MEP_senior_lc,zero_genes$GMP_senior_lc,MEP_markers,GMP_markers)
  MEP_GMP_bulkvssc_senior_markers <- marker_genes(bulk_MEPvsGMP_senior,sc_MEPvsGMP_senior,MEP_markers,GMP_markers)
  
  return(list(HSC_MEP_bulkvssc_mds=HSC_MEP_bulkvssc_mds,
              HSC_MEP_bulkvssc_mds_markers=HSC_MEP_bulkvssc_mds_markers,
              HSC_GMP_bulkvssc_mds=HSC_GMP_bulkvssc_mds,
              HSC_GMP_bulkvssc_mds_markers=HSC_GMP_bulkvssc_mds_markers,
              MEP_GMP_bulkvssc_mds=MEP_GMP_bulkvssc_mds,
              MEP_GMP_bulkvssc_mds_markers=MEP_GMP_bulkvssc_mds_markers,
              HSC_MEP_bulkvssc_senior=HSC_MEP_bulkvssc_senior,
              HSC_MEP_bulkvssc_senior_markers=HSC_MEP_bulkvssc_senior_markers,
              HSC_GMP_bulkvssc_senior=HSC_GMP_bulkvssc_senior,
              HSC_GMP_bulkvssc_senior_markers=HSC_GMP_bulkvssc_senior_markers,
              MEP_GMP_bulkvssc_senior=MEP_GMP_bulkvssc_senior,
              MEP_GMP_bulkvssc_senior_markers=MEP_GMP_bulkvssc_senior_markers,
              only_bulk_genes=only_bulk_genes,
              HSC_sc_mds=HSC_sc_mds,MEP_sc_mds=MEP_sc_mds,GMP_sc_mds=GMP_sc_mds,
              HSC_sc_senior=HSC_sc_senior,MEP_sc_senior=MEP_sc_senior,GMP_sc_senior=GMP_sc_senior,
              HSC_bulk_mds=HSC_bulk_mds,MEP_bulk_mds=MEP_bulk_mds,GMP_bulk_mds=GMP_bulk_mds,
              HSC_bulk_old=HSC_bulk_old,MEP_bulk_old=MEP_bulk_old,GMP_bulk_old=GMP_bulk_old))
}

compare_DEG_graph <- function(mean,sum,rep_mean,rep_sum){
  ##############################################################################
  # Function that compared the DEG for the different techniques used to make the pseudobulk matrix
  
  # INPUT VARIABLES
  # mean <- results of the analyse_sc_bulk function between bulk and psuedobulk (mean) for two cel types
  # sum <- results of the analyse_sc_bulk function between bulk and psuedobulk (sum) for two cel types
  # rep_sum <- results of the analyse_sc_bulk function between bulk and psuedobulk (rep_mean) for two cel types
  # rep_mean <- results of the analyse_sc_bulk function between bulk and psuedobulk (rep_mean) for two cel types
  
  # OUTPUT VARIABLES
  # data <- data comparing the DEG detected by each method
  # plot <- plot that comapres the DEG detected by each method
  ##############################################################################
  data <- rbind(mean$DEG_num_comparison,
                sum$DEG_num_comparison,
                rep_mean$DEG_num_comparison,
                rep_sum$DEG_num_comparison)
  data$sample <- c("All DEG (mean)","Bulk DEG (mean)","PSEUDO-SC DEG (mean)",
                   "All DEG (sum)","Bulk DEG (sum)","PSEUDO-SC DEG (sum)",
                   "All DEG (mean+replicas)","Bulk DEG (mean+replicas)","SC DEG (mean+replicas)",
                   "All DEG (sum+replicas)","Bulk DEG (sum+replicas)","SC DEG (sum+replicas)")
  title <- paste0("BULK: ",mean$bulk_sparsity,"
SC_MEAN: ",mean$sc_sparsity,"
SC_SUM: ",sum$sc_sparsity,"
SC_REP_MEAN: ",rep_mean$sc_sparsity,"
SC_REP_SUM: ",rep_sum$sc_sparsity)
  plot <- data %>% 
    ggplot() + 
    geom_col(aes(x=sample,y=low_expressed,fill="Low expressed genes"))+
    geom_col(aes(x=sample,y=only_bulk,fill="Only bulk genes"))+
    geom_col(aes(x=sample,y=unique_genes,fill="Unique genes"))+
    geom_col(aes(x=sample,y=common,fill="Common genes"))+
    geom_col(aes(x=sample,y=common_markers,fill="Marker genes"))+
    theme(axis.text.x=element_text(angle = -70, hjust = 0))+
    labs(y="DEG number",x=" ",fill="",title=title)+
    scale_fill_manual(values=c("#ff9dc5","#82dde0","#d24835","#fae980","#c5b5fb"))+
    theme(plot.title = element_text(size=8))
  return(list(data=data,plot=plot))
}

compare_DEG_graph2 <- function(mean,sum){
  ##############################################################################
  # Function that compared the DEG for the different techniques used to make the pseudobulk matrix
  # (without using the replicates)
  
  # INPUT VARIABLES
  # mean <- results of the analyse_sc_bulk function between bulk and psuedobulk (mean) for two cel types
  # sum <- results of the analyse_sc_bulk function between bulk and psuedobulk (sum) for two cel types
  
  # OUTPUT VARIABLES
  # data <- data comparing the DEG detected by each method
  # plot <- plot that comapres the DEG detected by each method
  ##############################################################################
  data <- rbind(mean$DEG_num_comparison,
                sum$DEG_num_comparison)
  data$sample <- c("All DEG (mean)","Bulk DEG (mean)","PSEUDO-SC DEG (mean)",
                   "All DEG (sum)","Bulk DEG (sum)","PSEUDO-SC DEG (sum)")
  title <- paste0("BULK: ",mean$bulk_sparsity,"
SC_MEAN: ",mean$sc_sparsity,"
SC_SUM: ",sum$sc_sparsity)
  plot <- data %>% 
    ggplot() + 
    geom_col(aes(x=sample,y=low_expressed,fill="Low expressed genes"))+
    geom_col(aes(x=sample,y=only_bulk,fill="Only bulk genes"))+
    geom_col(aes(x=sample,y=unique_genes,fill="Unique genes"))+
    geom_col(aes(x=sample,y=common,fill="Common genes"))+
    theme_light()+
    theme(axis.text.x=element_text(angle = -70, hjust = 0))+
    labs(y="DEG number",x=" ",fill="",title=title)+
    scale_fill_manual(values=c("#fae980","#c5b5fb","#82dde0","#ff9dc5"))+
    theme(plot.title = element_text(size=8))
  return(list(data=data,plot=plot))
}

compare_DEG_graph3 <- function(results,method){
  ##############################################################################
  # Function that compared the DEG for the different techniques used to make the pseudobulk matrix
  # (without using the replicates)
  
  # INPUT VARIABLES
  # results <- comparison to be analized
  # method <- type pf method to construct pseudobulk matrix
  
  # OUTPUT VARIABLES
  # data <- data comparing the DEG detected by each method
  # plot <- plot that comapres the DEG detected by each method
  ##############################################################################
  data <- results$DEG_num_comparison
  data$sample <- c("All DEG","Bulk DEG",paste0("PSEUDO-SC DEG ",method))
  plot <- data %>% 
    ggplot() + 
    geom_col(aes(x=sample,y=low_expressed,fill="Low expressed genes"))+
    geom_col(aes(x=sample,y=only_bulk,fill="Only bulk genes"))+
    geom_col(aes(x=sample,y=unique_genes,fill="Unique genes"))+
    geom_col(aes(x=sample,y=common,fill="Common genes"))+
    theme_light()+
    theme(axis.text.x=element_text(angle = -70, hjust = 0))+
    labs(y="DEG number",x=" ",fill="")+
    scale_fill_manual(values=c("#fae980","#c5b5fb","#82dde0","#ff9dc5"))+
    theme(plot.title = element_text(size=8))
  return(list(data=data,plot=plot))
}

process_results_markers <- function(results,name){
  ##############################################################################
  # Function that processes the results of the analyse_sc_bulk function comparing the marker genes
  
  # INPUT VARIABLES
  # reults <- reults of the analyse_sc_bulk function for an specific method
  
  # OUTPUT VARIABLES
  # total percentages <- tibble with all the percentages detected in each comparison
  ##############################################################################
  total_percentages <- rbind(results$HSC_MEP_bulkvssc_mds_markers$percentages,results$HSC_MEP_bulkvssc_senior_markers$percentages,
                             results$HSC_GMP_bulkvssc_mds_markers$percentages,results$HSC_GMP_bulkvssc_senior_markers$percentages,
                             results$MEP_GMP_bulkvssc_mds_markers$percentages,results$MEP_GMP_bulkvssc_senior_markers$percentages)
  total_percentages$COMPARISON <- c("HSC VS. MEP mds","HSC VS. MEP mds",
                                    "HSC VS.MEP senior","HSC VS.MEP senior",
                                    "HSC VS. GMP mds","HSC VS. GMP mds",
                                    "HSC VS. GMP senior","HSC VS. GMP senior",
                                    "MEP VS. GMP mds","MEP VS. GMP mds",
                                    "MEP VS. GMP senior","MEP_ VS. GMP senior")
  total_percentages$EXPERIMENT <- c(rep(name,12))
  total_percentages <- total_percentages %>% pivot_longer(cols=c(1:2),names_to="category") 
  return(total_percentages)
}

top_genes <- function(bulk,start,end,n,path){
  ##############################################################################
  # Function that processes the results of the analyse_sc_bulk function comparing the marker genes
  
  # INPUT VARIABLES
  # bulk <- bulk RNAseq matrix
  # start <- position to start extracting genes in each cell type
  # end <- position to end extracting genes in each cell type
  # n <- number of genes to extract
  # path <- path to save the top genes in csv format
  
  # OUTPUT VARIABLES
  # heatmap <- heatmap with the log(counts+1) the selected genes
  ##############################################################################
  
  bulk <- cbind("genes"=rownames(bulk),bulk) # generate a column with rownames
  bulk <- keep_protein_coding_bulk(bulk) # only keep protein coding genes
  
  HSC_bulk <- select_samples(bulk,"HSC") #select HSC
  MEP_bulk <- select_samples(bulk,"MEP") #select MEP
  GMP_bulk <- select_samples(bulk,"GMP") #select GMP
  
  bulk_MEPvsHSC <- DEA(MEP_bulk,HSC_bulk,"MEP","HSC","bulk_MEPvsHSC") #DEA between two cell types
  bulk_MEPvsHSC <- bulk_MEPvsHSC$upregulated %>% arrange(desc(LFC)) #arrange by descending LFC
  bulk_GMPvsHSC <- DEA(GMP_bulk,HSC_bulk,"GMP","HSC","bulk_GMPvsHSC")
  bulk_GMPvsHSC <- bulk_GMPvsHSC$upregulated %>% arrange(desc(LFC))
  
  bulk_HSCvsMEP <- DEA(HSC_bulk,MEP_bulk,"HSC","MEP","bulk_HSCvsMEP")
  bulk_HSCvsMEP <- bulk_HSCvsMEP$upregulated %>% arrange(desc(LFC))
  bulk_GMPvsMEP <- DEA(GMP_bulk,MEP_bulk,"GMP","MEP","bulk_GMPvsMEP")
  bulk_GMPvsMEP <- bulk_GMPvsMEP$upregulated %>% arrange(desc(LFC))
  
  bulk_MEPvsGMP <- DEA(MEP_bulk,GMP_bulk,"MEP","GMP","bulk_MEPvsGMP")
  bulk_MEPvsGMP <- bulk_MEPvsGMP$upregulated %>% arrange(desc(LFC))
  bulk_HSCvsGMP <- DEA(HSC_bulk,GMP_bulk,"HSC","GMP","bulk_HSCvsGMP")
  bulk_HSCvsGMP <- bulk_HSCvsGMP$upregulated %>% arrange(desc(LFC))
  
  HSC <- inner_join(bulk_MEPvsHSC,bulk_GMPvsHSC,by="GENE") #join the common genes
  MEP <- inner_join(bulk_HSCvsMEP,bulk_GMPvsMEP,by="GENE")
  GMP <- inner_join(bulk_MEPvsGMP,bulk_HSCvsGMP,by="GENE")
  
  top_genes <- rbind(HSC[start:end,],MEP[start:end,],GMP[start:end,]) #select the top genes
  bulk_top_genes <- inner_join(top_genes,bulk,by=c("GENE"="genes")) #inner join with the bulk
  bulk_top_genes$GENE <- NULL #eiminate innecessary columns
  bulk_top_genes$LFC.x <- NULL
  bulk_top_genes$LFC.y <- NULL
  bulk_top_genes <- as.matrix(bulk_top_genes) #convert to matrix
  row.names(bulk_top_genes) <- c(rep("HSC_gene",n),rep("MEP_gene",n),rep("GMP_gene",n)) #add rownames with the cell type
  bulk_top_genes <- log2(bulk_top_genes+1) #calculate the log (+1 to avoid log 0)
  plot <- heatmap.2(bulk_top_genes,dendrogram = "none",trace="none",Rowv = "none",Colv="none",
                    cexRow = 0.5,cexCol = 0.5) #plot the heatmap to check
  
  top_genes_symbol = unname(mapIds(org.Hs.eg.db, 
                                   keys=top_genes$GENE, 
                                   column="SYMBOL",
                                   keytype="ENSEMBL",
                                   multiVals="first")) #obtain the gene symbol
  top_genes_symbol <- cbind(top_genes_symbol,category =c(rep("HSC_gene",n),rep("MEP_gene",n),rep("GMP_gene",n))) #addthe cell type
  write.csv(top_genes_symbol,file=path) #save csv file
  return(list(plot <- plot,top_genes_symbol <- top_genes_symbol))
}


save_tsv_pseudo_senior <- function(rds1,rds2,rds3,tsv){
  ##############################################################################
  # Function that saves the pseudobulk matrix in tsv file with HSC, MEP and GMP cells
  # INPUT VARIABLES
  # rds <- name of the rds file without extension
  
  ##############################################################################
  matrix1 <- readRDS(paste0(wd,pseudobulk,rds1,".rds"))$pseudo_sc
  matrix2 <- readRDS(paste0(wd,pseudobulk,rds2,".rds"))$pseudo_sc
  matrix3 <- readRDS(paste0(wd,pseudobulk,rds3,".rds"))$pseudo_sc
  pseudo <- inner_join(matrix1,matrix2,by=c("genes_senior"),suffix=c("_1","_2")) 
  pseudo <- inner_join(pseudo,matrix3,by=c("genes_senior"),suffix=c("","_3")) 
  write.table(pseudo,
              file=paste0(wd,pseudobulk,tsv,".tsv"),
              sep="\t")
}

save_tsv_pseudo_mds <- function(rds1,rds2,rds3,rds4,tsv){
  ##############################################################################
  # Function that saves the pseudobulk matrix in tsv file with HSC, MEP and GMP cells
  # INPUT VARIABLES
  # rds <- name of the rds file without extension
  
  ##############################################################################
  matrix1 <- readRDS(paste0(wd,pseudobulk,rds1,".rds"))$pseudo_sc
  matrix2 <- readRDS(paste0(wd,pseudobulk,rds2,".rds"))$pseudo_sc
  matrix3 <- readRDS(paste0(wd,pseudobulk,rds3,".rds"))$pseudo_sc
  matrix4 <- readRDS(paste0(wd,pseudobulk,rds4,".rds"))$pseudo_sc
  pseudo <- inner_join(matrix1,matrix2,by=c("genes_mds"),suffix=c("_1","_2")) 
  pseudo <- inner_join(pseudo,matrix3,by=c("genes_mds"),suffix=c("","_3")) 
  pseudo <- inner_join(pseudo,matrix4,by=c("genes_mds"),suffix=c("","_4")) 
  write.table(pseudo,
              file=paste0(wd,pseudobulk,tsv,".tsv"),
              sep="\t")
}


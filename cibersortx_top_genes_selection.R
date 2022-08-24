################################################################################
## AUTHOR: NAROA LEGARRA MARCOS
## DATE: 01/06/22

## Selects the top differentially expressed genes for each type of cell vs the other two
## using FACS-bulk RNAseq data
################################################################################


#LOAD AND PROCESS DATA ######################################################################
#Bulk RNA data (senior and mds)-------------------------------------------------
load(paste0(wd,"DATA/bulkrna/MDSCounts_HSC_progenitors_02-10-20.Rdata")) #read counts data
counts <- counts[apply(counts, 1, function(x) !all(x==0)),] #remove genes with 0 counts in all samples

################################################################################
#select the healthy HSC, GMP and MEP cells from the bulk data
bulk_senior <-as.data.frame(counts[,c(1,147,158,162,166,184,188,193, #HSC elderly
                               11,142,144,148,151,155,159,163,167,185,189,194,197, #MEP elderly
                               141,143,150,153,157,161,165,183,187,191,192,196)]) #GMP elderly
#obtain the 150 top DEG for each celltype
top_genes_senior_150 <- top_genes(bulk_senior,1,150,150,paste0(wd,"DATA/top_genes/top_genes_senior_150.csv"))
#obtain the 280 top DEG for each celltype after the top 150 
top_genes_senior_after150 <- top_genes(bulk_senior,151,430,280,paste0(wd,"DATA/top_genes/top_genes_senior_after150.csv"))

#################################################################################
#select the mds HSC, GMP and MEP cells from the bulk data
bulk_mds <-as.data.frame(counts[,c(204,210,212,215,220,223,227,233,277,288,302,322,349,353,360,433, #HSC MLD 
                               259,262,216,219,224,228,236,278,289,307,323,338,350,354,361,434,438,378,#MEP MLD
                               258,261,214,218,222,226,235,276,287,301,321,337,348,352,359,436)]) #GMP MLD
#obtain the 150 top DEG for each celltype
top_genes_mds_150 <- top_genes(bulk_mds,1,150,150,paste0(wd,"DATA/top_genes/top_genes_mds_150.csv"))
#obtain the 300 top DEG for each celltype after the top 150 
top_genes_mds_after150 <- top_genes(bulk_mds,151,450,300,paste0(wd,"DATA/top_genes/top_genes_mds_after150.csv"))

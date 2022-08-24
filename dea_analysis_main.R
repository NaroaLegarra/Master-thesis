################################################################################
## AUTHOR: NAROA LEGARRA MARCOS
## DATE: 10/03/22

## Performs the differential gene expression analysis between cell types in 
## FACS-bulk RNAseq and pseudobulk RNAseq in the CUN datasets and compares results
################################################################################

wd <- "~/Desktop/TFM" #set the working directory
pseudobulk <- "/DATA/pseudo/" #path to the folder to store the pseudobulks
figures <- "/FIGURES"

###################################### LOAD DATA ###############################

# SC RNA data ------------------------------------------------------------------
pseudo_all_mean <- readRDS(paste0(wd,pseudobulk,"pseudo_all_mean.rds")) 
zero_genes_mean <- readRDS(paste0(wd,pseudobulk,"zero_genes_mean.rds")) 

pseudo_all_sum <- readRDS(paste0(wd,pseudobulk,"pseudo_all_sum.rds"))
zero_genes_sum <- readRDS(paste0(wd,pseudobulk,"zero_genes_sum.rds")) 

pseudo_all_rep_mean <- readRDS(paste0(wd,pseudobulk,"pseudo_all_rep_mean.rds")) 
zero_genes_rep_mean <- readRDS(paste0(wd,pseudobulk,"zero_genes_rep_mean.rds")) 

pseudo_all_rep_sum <- readRDS(paste0(wd,pseudobulk,"pseudo_all_rep_sum.rds")) 
zero_genes_rep_sum <- readRDS(paste0(wd,pseudobulk,"zero_genes_rep_sum.rds"))

################################## PERFORM ANALYSIS ############################

#All bulk samples
bulk <- readRDS(paste0(wd,"/DATA/bulkrna/bulk_all.rds")) 
results_mean <- analyse_sc_bulk(pseudo_all_mean,bulk,zero_genes_mean)
saveRDS(results_mean,file=paste0(wd,"/DATA/DEA_results/results_mean_all_bulk.rds"))
results_sum <- analyse_sc_bulk(pseudo_all_sum,bulk,zero_genes_sum)
saveRDS(results_sum,file=paste0(wd,"/DATA/DEA_results/results_sum_all_bulk.rds"))

#Reduced bulk samples (repeat for each 5 replicas)
bulk <- readRDS(paste0(wd,"/DATA/bulkrna/bulk_1.rds")) 
results_mean <- analyse_sc_bulk(pseudo_all_mean,bulk,zero_genes_mean)
saveRDS(results_mean,file=paste0(wd,"/DATA/DEA_results/results_mean_4_rep1_bulk.rds"))
results_sum <- analyse_sc_bulk(pseudo_all_sum,bulk,zero_genes_sum)
saveRDS(results_sum,file=paste0(wd,"/DATA/DEA_results/results_sum_4_rep1_bulk.rds"))

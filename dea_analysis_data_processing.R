################################################################################
## AUTHOR: NAROA LEGARRA MARCOS
## DATE: 10/03/22

## Process the sc RNAseq and FACS-bulk RNAseq data from 
## the CUN dataset data for further analysis
################################################################################

wd <- "~/Desktop/TFM" #set the working directory
pseudobulk <- "/DATA/pseudo/" #path to the folder to store the pseudobulks
single_cell <- "/DATA/single_cell/"
sc_mds <- "/DATA/mds_scrna" #path to the folder with mds sc data
sc_senior <- "/DATA/healthy_scrna" #path to the folder with senior sc data

################################ SC RNA data ###################################
# SC RNA data (mds)-------------------------------------------------------------
#Read single cell data and calculate the pseudosc
preprocess_mds(paste0(wd,sc_mds,"/rawcounts_mds1.txt"),
               "mds1","_mds")
preprocess_mds(paste0(wd,sc_mds,"/rawcounts_mds3.txt"),
               "mds3","_mds")
preprocess_mds(paste0(wd,sc_mds,"/rawcounts_mds5.txt"),
               "mds5","_mds")
preprocess_mds(paste0(wd,sc_mds,"/rawcounts_mds10.txt"),
               "mds10","_mds")

# SC RNA data (senior)----------------------------------------------------------
#Read single cell data and calculate the pseudosc
preprocess_senior(paste0(wd,sc_senior,"/GSM5460411_seurat_obj_norm.rds"),
                  paste0(wd,sc_senior,"/senior1_manual_labels.txt"),
                  "senior1","_senior")
preprocess_senior(paste0(wd,sc_senior,"/GSM5460412_seurat_obj_norm.rds"),
                  paste0(wd,sc_senior,"/senior2_manual_labels.txt"),
                  "senior2","_senior")
preprocess_senior(paste0(wd,sc_senior,"/GSM5460413_seurat_obj_norm.rds"),
                  paste0(wd,sc_senior,"/senior3_manual_labels.txt"),
                  "senior3","_senior")

############################# Join alll cell num data ##########################
# Read data
cells_mds1 <- readRDS(paste0(wd,pseudobulk,"pseudo_mds1_cellnum.rds"))
cells_mds3 <- readRDS(paste0(wd,pseudobulk,"pseudo_mds3_cellnum.rds"))
cells_mds5 <- readRDS(paste0(wd,pseudobulk,"pseudo_mds5_cellnum.rds"))
cells_mds10 <- readRDS(paste0(wd,pseudobulk,"pseudo_mds10_cellnum.rds"))
cells_senior1 <- readRDS(paste0(wd,pseudobulk,"pseudo_senior1_cellnum.rds"))
cells_senior2 <- readRDS(paste0(wd,pseudobulk,"pseudo_senior2_cellnum.rds"))
cells_senior3 <- readRDS(paste0(wd,pseudobulk,"pseudo_senior3_cellnum.rds"))

cell_num <- rbind(cells_mds1,cells_mds3,cells_mds5,cells_mds10,
                  cells_senior1,cells_senior2,cells_senior3) #bind data
cell_num$sample <- c("mds1","mds3","mds5","mds10",
                     "senior1","senior2","senior3") #add sample names
write.csv(cell_num,paste0(wd,pseudobulk,"/cell_num.csv")) #save

# Read all pseudosc data---------------------------------------------------------
pseudo_mds1_mean <- readRDS(paste0(wd,pseudobulk,"pseudo_mds1_mean.rds"))
pseudo_mds1_sum <- readRDS(paste0(wd,pseudobulk,"pseudo_mds1_sum.rds"))
pseudo_mds1_rep_mean <- readRDS(paste0(wd,pseudobulk,"pseudo_mds1_rep_mean.rds"))
pseudo_mds1_rep_sum <- readRDS(paste0(wd,pseudobulk,"pseudo_mds1_rep_sum.rds"))

pseudo_mds3_mean <- readRDS(paste0(wd,pseudobulk,"pseudo_mds3_mean.rds"))
pseudo_mds3_sum <- readRDS(paste0(wd,pseudobulk,"pseudo_mds3_sum.rds"))
pseudo_mds3_rep_mean <- readRDS(paste0(wd,pseudobulk,"pseudo_mds3_rep_mean.rds"))
pseudo_mds3_rep_sum <- readRDS(paste0(wd,pseudobulk,"pseudo_mds3_rep_sum.rds"))

pseudo_mds5_mean <- readRDS(paste0(wd,pseudobulk,"pseudo_mds5_mean.rds"))
pseudo_mds5_sum <- readRDS(paste0(wd,pseudobulk,"pseudo_mds5_sum.rds"))
pseudo_mds5_rep_mean <- readRDS(paste0(wd,pseudobulk,"pseudo_mds5_rep_mean.rds"))
pseudo_mds5_rep_sum <- readRDS(paste0(wd,pseudobulk,"pseudo_mds5_rep_sum.rds"))

pseudo_mds10_mean <- readRDS(paste0(wd,pseudobulk,"pseudo_mds10_mean.rds"))
pseudo_mds10_sum <- readRDS(paste0(wd,pseudobulk,"pseudo_mds10_sum.rds"))
pseudo_mds10_rep_mean <- readRDS(paste0(wd,pseudobulk,"pseudo_mds10_rep_mean.rds"))
pseudo_mds10_rep_sum <- readRDS(paste0(wd,pseudobulk,"pseudo_mds10_rep_sum.rds"))

pseudo_senior1_mean <- readRDS(paste0(wd,pseudobulk,"pseudo_senior1_mean.rds"))
pseudo_senior1_sum <- readRDS(paste0(wd,pseudobulk,"pseudo_senior1_sum.rds"))
pseudo_senior1_rep_mean <- readRDS(paste0(wd,pseudobulk,"pseudo_senior1_rep_mean.rds"))
pseudo_senior1_rep_sum <- readRDS(paste0(wd,pseudobulk,"pseudo_senior1_rep_sum.rds"))

pseudo_senior2_mean <- readRDS(paste0(wd,pseudobulk,"pseudo_senior2_mean.rds"))
pseudo_senior2_sum <- readRDS(paste0(wd,pseudobulk,"pseudo_senior2_sum.rds"))
pseudo_senior2_rep_mean <- readRDS(paste0(wd,pseudobulk,"pseudo_senior2_rep_mean.rds"))
pseudo_senior2_rep_sum <- readRDS(paste0(wd,pseudobulk,"pseudo_senior2_rep_sum.rds"))

pseudo_senior3_mean <- readRDS(paste0(wd,pseudobulk,"pseudo_senior3_mean.rds"))
pseudo_senior3_sum <- readRDS(paste0(wd,pseudobulk,"pseudo_senior3_sum.rds"))
pseudo_senior3_rep_mean <- readRDS(paste0(wd,pseudobulk,"pseudo_senior3_rep_mean.rds"))
pseudo_senior3_rep_sum <- readRDS(paste0(wd,pseudobulk,"pseudo_senior3_rep_sum.rds"))


#Select genes that are common and join all matrices-----------------------------
## MEAN
join_pseudo(pseudo_mds1_mean,pseudo_mds3_mean,pseudo_mds5_mean,pseudo_mds10_mean,
            pseudo_senior1_mean,pseudo_senior2_mean,pseudo_senior3_mean,
            "mean")
## SUM
join_pseudo(pseudo_mds1_sum,pseudo_mds3_sum,pseudo_mds5_sum,pseudo_mds10_sum,
            pseudo_senior1_sum,pseudo_senior2_sum,pseudo_senior3_sum,
            "sum")

## REP MEAN
join_pseudo(pseudo_mds1_rep_mean,pseudo_mds3_rep_mean,pseudo_mds5_rep_mean,pseudo_mds10_rep_mean,
            pseudo_senior1_rep_mean,pseudo_senior2_rep_mean,pseudo_senior3_rep_mean,
            "rep_mean")

## REP SUM
join_pseudo(pseudo_mds1_rep_sum,pseudo_mds3_rep_sum,pseudo_mds5_rep_sum,pseudo_mds10_rep_sum,
            pseudo_senior1_rep_sum,pseudo_senior2_rep_sum,pseudo_senior3_rep_sum,
            "rep_sum")

###################### Bulk RNA data (senior and mds) ##########################

load(paste0(wd,"/DATA/bulkrna/MDSCounts_HSC_progenitors_02-10-20.Rdata")) # load data
counts <- counts[apply(counts, 1, function(x) !all(x==0)),] #remove genes with 0 counts in all samples

## Obtain the samples corresponding to the cells HSC, MEP, GMP and CMP (Old and MLD) and convert it to a data frame
bulk <-as.data.frame(counts[,c(1,147,158,162,166,184,188,193, #HSC elderly
                               204,210,212,215,220,223,227,233,277,288,302,322,349,353,360,433,437, #HSC MLD + SLD (last one)
                               11,142,144,148,151,155,159,163,167,185,189,194,197, #MEP elderly
                               259,262,216,219,224,228,236,278,289,307,323,338,350,354,361,434,438,378,#MEP MLD
                               141,143,150,153,157,161,165,183,187,191,192,196, #GMP elderly
                               258,261,214,218,222,226,235,276,287,301,321,337,348,352,359,436)]) #GMP MLD

bulk <- cbind("genes"=rownames(bulk),bulk) # generate a column with rownames
bulk <- keep_protein_coding_bulk(bulk) # only keep protein coding genes
saveRDS(bulk,paste0(wd,"/DATA/bulkrna/bulk_all.rds")) #save data
write_csv(bulk,paste0(wd,"/DATA/bulkrna/bulk_all.csv"))

bulk1 <-as.data.frame(counts[,c(1,147,158, #HSC elderly
                               204,210,212,215, #HSC MLD
                               11,142,144, #MEP elderly
                               259,262,216,219,#MEP MLD
                               141,143,150, #GMP elderly
                               258,261,214,218)]) #GMP MLD
#bulk2 <-as.data.frame(counts[,c(162,166,184, #HSC elderly
                               220,223,227,233, #HSC MLD
                               148,151,155, #MEP elderly
                               224,228,236,278,#MEP MLD
                               153,157,161, #GMP elderly
                               222,226,235,276)]) #GMP MLD
#bulk3 <-as.data.frame(counts[,c(184,188,193, #HSC elderly
                                353,360,433,437, #HSC MLD
                                189,194,197, #MEP elderly
                                361,434,438,378,#MEP MLD
                                191,192,196, #GMP elderly
                                348,352,359,436)]) #GMP MLD
#bulk4 <-as.data.frame(counts[,c(1,162,188, #HSC elderly
                               277,288,302,322, #HSC MLD
                               159,163,167, #MEP elderly
                               289,307,323,338,#MEP MLD
                               165,183,187, #GMP elderly
                               287,301,321,337)]) #GMP MLD
#bulk5 <-as.data.frame(counts[,c(147,166,193, #HSC elderly
                               204,220,277,349, #HSC MLD
                               11,148,159, #MEP elderly
                               259,224,289,350,#MEP MLD
                               141,153,165, #GMP elderly
                               258,222,287,348)]) #GMP MLD
bulk <- cbind("genes"=rownames(bulk),bulk) # generate a column with rownames
bulk <- keep_protein_coding_bulk(bulk) # only keep protein coding genes
saveRDS(bulk,paste0(wd,"/DATA/bulkrna/bulk_1.rds")) #save data
write_csv(bulk,paste0(wd,"/DATA/bulkrna/bulk_1.csv"))

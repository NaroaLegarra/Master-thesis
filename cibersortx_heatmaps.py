######################################################################################################################
## AUTHOR: NAROA LEGARRA MARCOS
## DATE: 10/03/22

## Script to make heatmaps comparing reference and CIBERSORTx predicted results
######################################################################################################################

wd = "/home/nlegarra/data/jfuente/digital_cytometry/UpsamplingSingleNN" #set the working directory
figures = "/home/nlegarra/data/jfuente/digital_cytometry/UpsamplingSingleNN/reports/figures" #set folder to store heatmaps

################################################### IMPORT LIBRARIES ################################################
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

################################################### DEFINE FUNCTIONS ################################################

def preprocess_matrix(matrix):
    #####################################################################################
    #Function that processes the cibersort predicted matrix
  
    # INPUT VARIABLES
    # matrix <- sc sequencing matrix
  
    # OUTPUT VARIABLES
    # proccessed matrix
    #####################################################################################
    matrix=matrix.set_index("GeneSymbol")
    matrix=np.log(matrix+1)
    return(matrix)

def heatmap(pred_matrix,ref_matrix,not_matrix,path,title):
    #####################################################################################
    #Function that generates heatmaps for cibersor predicted matrix, refrence matrix and incorrectly predicted cells

    # INPUT VARIABLES
    # pred_matrix <- cibersort predicted matrix
    # ref_matrix <- refrence matrix
    # not matrix <– incorrectly predicted cells (this matrix should not have expression
    # path <- path for saving
    # title <- title

    # OUTPUT VARIABLES
    # heatmap
    #####################################################################################
    plt.figure(figsize=(15, 9), dpi=100)
    plt.subplot(1,3,1)
    sns.heatmap(ref_matrix,cmap="viridis",vmin=-1,vmax=10,yticklabels=False,cbar=False,xticklabels=True).get_figure()
    plt.title(title,fontsize=15)
    plt.xticks(fontsize=8)
    plt.ylabel("Gene expression log10",fontsize=18)
    plt.subplot(1,3,2)
    sns.heatmap(pred_matrix,cmap="viridis",vmin=-1,vmax=10,yticklabels=False,cbar=False,xticklabels=True).get_figure()
    plt.title("CIBERSORTx predicted GEP",fontsize=15)
    plt.xticks(fontsize=8)
    plt.ylabel("")
    plt.subplot(1,3,3)
    sns.heatmap(not_matrix,cmap="viridis",vmin=-1,vmax=10,yticklabels=False,xticklabels=True).get_figure()
    plt.xticks(fontsize=4)
    plt.ylabel("")
    plt.title("Incorrectly predicted cells",fontsize=15)
    plt.tight_layout()
    plt.savefig(path)
    plt.clf()

def plot_heatmap(path_HSC,path_GMP,path_MEP,ref,genes,path_for_heat,title):
    #####################################################################################
    #Function that processes the data and generates the heatmap
  
    # INPUT VARIABLES
    # path_HSC <- path to HSC cibersortx predicted matrix
    # path_GMP <- path to GMP cibersortx predicted matrix
    # path_MEP <- path to MEP cibersortx predicted matrix
    # ref <- reference matrix
    # genes <- genes used for CIBERSORTx high resolution prediction
    # path_for_heat <- path for saving the hetamap
    # title <- title
  
    # OUTPUT VARIABLES
    # heatmap
    #####################################################################################
    HSC=preprocess_matrix(pd.read_table(path_HSC))
    filter_col1 = [col for col in HSC if 'HSC' in col]
    filter_col2 = [col for col in HSC if not 'HSC' in col]
    not_HSC=HSC[filter_col2]
    HSC=HSC[filter_col1]
    GMP=preprocess_matrix(pd.read_table(path_GMP))
    filter_col1 = [col for col in GMP if 'GMP' in col]
    filter_col2 = [col for col in GMP if not 'GMP' in col]
    not_GMP=GMP[filter_col2]
    GMP=GMP[filter_col1]
    MEP=preprocess_matrix(pd.read_table(path_MEP))
    filter_col1 = [col for col in MEP if 'MEP' in col]
    filter_col2 = [col for col in MEP if not 'MEP' in col]
    not_MEP=MEP[filter_col2]
    MEP=MEP[filter_col1]
    [i for i in genes if i in HSC.index.values]
    genes=genes['top_genes_symbol'].to_numpy()
    genes=listcomp(genes.astype(str),HSC.index.values.astype(str))
    genes=listcomp(genes,ref.index.values.astype(str))
    ref=ref.loc[genes]
    ref=np.log(ref+1)
    predicted=pd.concat([HSC,MEP,GMP],axis=1).fillna(0)
    predicted=predicted.loc[genes]
    not_matrix=pd.concat([not_HSC.add_suffix("_HSC"),not_MEP.add_suffix("_MEP"),not_GMP.add_suffix("_GMP")],axis=1).fillna(0)
    not_matrix=not_matrix.loc[genes]
    heatmap(predicted,ref,not_matrix,path_for_heat,title)

def order_matrix(matrix):
    #####################################################################################
    #Function that orders the pseudobulk matrix
  
    # INPUT VARIABLES
    # matrix <- pseudobulk matrix
  
    # OUTPUT VARIABLES
    # ordered_matrix <- ordered pseudobulk matrix
    #####################################################################################
    filter_col1 = [col for col in matrix if col.startswith("HSC")]
    filter_col2 = [col for col in matrix if col.startswith("MEP")]
    filter_col3 = [col for col in matrix if col.startswith("GMP")]
    ordered_matrix=pd.concat([matrix[filter_col1],matrix[filter_col2],matrix[filter_col3]],axis=1)
    return(ordered_matrix)

def listcomp(k, d):
    #####################################################################################
    #Function that intersects k with d maintaining the order in k
    # INPUT VARIABLES
    # k <- vector
    # d <- vector
  
    # OUTPUT VARIABLES
    # k elements present in d
    #####################################################################################
    return [i for i in k if i in d]

##################################################################### READ DATA ##################################################################
bulk=pd.read_csv("/home/nlegarra/data/jfuente/digital_cytometry/UpsamplingSingleNN/data/raw/bulkrna/healthy_mds/healthy_mds_bulk.csv")
bulk=bulk.set_index("genes")
bulk_mds=bulk[['SMD19957_HSC', 'SMD20793_HSC','SMD21949_HSC','SMD22929_HSC','SMD27762_HSC','SMD28381_HSC', 
            'SMD19957_MEP','SMD20793_MEP','SMD21949_MEP','SMD22929_MEP','SMD27762_MEP','SMD28381_MEP', 
            'SMD19957_GMP','SMD20793_GMP','SMD21949_GMP','SMD22929_GMP','SMD27762_GMP','SMD28381_GMP']]
bulk_healthy=bulk[["MOS6_HSC","MOS11_HSC","MOS17_19268_HSC","MOS19_HSC","MOS2_HSC","MOS3_HSC","MOS4_HSC","MOS7_HSC",
            "MOS1_MEP","MOS10_MEP","MOS11_MEP","MOS12_MEP","MOS16_19267_MEP","MOS17_19268_MEP","MOS19_MEP","MOS2_MEP","MOS3_MEP","MOS4_MEP","MOS7_MEP","MOS9_MEP_PET16674",
            "MOS1_GMP","MOS10_GMP","MOS12_GMP","MOS16_19267_GMP","MOS17_19268_GMP","MOS19_GMP","MOS2_GMP","MOS3_GMP","MOS4_GMP","MOS6_GMP","MOS7_GMP","MOS9_GMP_PET16674"]]

top_150_healthy=pd.read_csv(wd+"/data/processed/healthy_mds/top_genes/top_150_genes_healthy.csv",sep=';')
top_150_mds=pd.read_csv(wd+"/data/processed/healthy_mds/top_genes/top_150_genes_mds.csv")

top_300_healthy=pd.read_csv(wd+"/data/processed/healthy_mds/top_genes/top_300_after_top150_genes_healthy.csv",sep=";")
top_300_mds=pd.read_csv(wd+"/data/processed/healthy_mds/top_genes/top_300_after_top150_genes_mds.csv",sep=";")

pseudo_healthy_mean=order_matrix(pd.read_csv(wd+"/data/processed/healthy_mds/pseudo_bulk/pseudo_senior_mean.tsv",sep='\t').set_index("genes_senior"))
pseudo_healthy_sum=order_matrix(pd.read_csv(wd+"/data/processed/healthy_mds/pseudo_bulk/pseudo_senior_sum.tsv",sep='\t').set_index("genes_senior"))
pseudo_mds_mean=order_matrix(pd.read_csv(wd+"/data/processed/healthy_mds/pseudo_bulk/pseudo_mds_mean.tsv",sep='\t').set_index("genes_mds"))
pseudo_mds_sum=order_matrix(pd.read_csv(wd+"/data/processed/healthy_mds/pseudo_bulk/pseudo_mds_sum.tsv",sep='\t').set_index("genes_mds"))

########################################################### HEATMAP TOP 150 DEG ###################################################################

plot_heatmap((wd+"/outputs/validation/healthy1/GEPInference/GEPInference_Healthy_SMode_150MarkerGenes_3CT/CIBERSORTxHiRes_Job75_HSC_Window12.txt"),
            (wd+"/outputs/validation/healthy1/GEPInference/GEPInference_Healthy_SMode_150MarkerGenes_3CT/CIBERSORTxHiRes_Job75_GMP_Window12.txt"),
            (wd+"/outputs/validation/healthy1/GEPInference/GEPInference_Healthy_SMode_150MarkerGenes_3CT/CIBERSORTxHiRes_Job75_MEP_Window12.txt"),
            bulk_healthy,top_150_healthy,"/home/nlegarra/data/nlegarra/DigitalCytometry_Validation/FIGURES/heatmap_healthy_top150.pdf","FACS-bulk RNAseq healthy top 150")

plot_heatmap((wd+"/outputs/validation/mds1/GEPInference/NoBatchCorrection_GEPInference_top150MarkerGenes_3CT/CIBERSORTxHiRes_Job79_HSC_Window9.txt"),
                (wd+"/outputs/validation/mds1/GEPInference/NoBatchCorrection_GEPInference_top150MarkerGenes_3CT/CIBERSORTxHiRes_Job79_GMP_Window9.txt"),
                (wd+"/outputs/validation/mds1/GEPInference/NoBatchCorrection_GEPInference_top150MarkerGenes_3CT/CIBERSORTxHiRes_Job79_MEP_Window9.txt"),
                bulk_mds,top_150_healthy,"/home/nlegarra/data/nlegarra/DigitalCytometry_Validation/FIGURES/heatmap_mds_top150.pdf","FACS-bulk RNAseq mds top 150")

########################################################### HEATMAP TOP 300 AFTER 150 DEG ###################################################################

plot_heatmap("/home/nlegarra/data/nlegarra/DigitalCytometry_Validation/DATA/top_300_genes/CIBERSORTxHiRes_Job82_HSC_Window12.txt",
                "/home/nlegarra/data/nlegarra/DigitalCytometry_Validation/DATA/top_300_genes/CIBERSORTxHiRes_Job82_GMP_Window12.txt",
                "/home/nlegarra/data/nlegarra/DigitalCytometry_Validation/DATA/top_300_genes/CIBERSORTxHiRes_Job82_MEP_Window12.txt",
                bulk_healthy,top_300_healthy,"/home/nlegarra/data/nlegarra/DigitalCytometry_Validation/FIGURES/heatmap_healthy_top300.pdf","FACS-bulk RNAseq healthy top 151-450")

plot_heatmap("/home/nlegarra/data/nlegarra/DigitalCytometry_Validation/DATA/top_300_genes/CIBERSORTxHiRes_Job80_HSC_Window9.txt",
                "/home/nlegarra/data/nlegarra/DigitalCytometry_Validation/DATA/top_300_genes/CIBERSORTxHiRes_Job80_GMP_Window9.txt",
                "/home/nlegarra/data/nlegarra/DigitalCytometry_Validation/DATA/top_300_genes/CIBERSORTxHiRes_Job80_MEP_Window9.txt",
                bulk_mds,top_300_mds,"/home/nlegarra/data/nlegarra/DigitalCytometry_Validation/FIGURES/heatmap_mds_top300.pdf","FACS-bulk RNAseq mds top 151-450")

########################################################### HEATMAP TOP 150 MEAN ###################################################################

plot_heatmap("/home/nlegarra/data/nlegarra/DigitalCytometry_Validation/DATA/CIBERSORTx/GEPinference/HEALTHY_MEAN_TOP_150_NO_NORM/CIBERSORTxHiRes_Job21_HSC_Window5.txt",
                "/home/nlegarra/data/nlegarra/DigitalCytometry_Validation/DATA/CIBERSORTx/GEPinference/HEALTHY_MEAN_TOP_150_NO_NORM/CIBERSORTxHiRes_Job21_GMP_Window5.txt",
                "/home/nlegarra/data/nlegarra/DigitalCytometry_Validation/DATA/CIBERSORTx/GEPinference/HEALTHY_MEAN_TOP_150_NO_NORM/CIBERSORTxHiRes_Job21_MEP_Window5.txt",
                pseudo_healthy_mean,top_150_healthy,"/home/nlegarra/data/nlegarra/DigitalCytometry_Validation/FIGURES/heatmap_healthy_top150_pseudomean.pdf","pseudo-bulk (mean) healthy top 150")

plot_heatmap("/home/nlegarra/data/nlegarra/DigitalCytometry_Validation/DATA/CIBERSORTx/GEPinference/MDS_MEAN_TOP_150_NO_NORM/CIBERSORTxHiRes_Job18_HSC_Window6.txt",
                "/home/nlegarra/data/nlegarra/DigitalCytometry_Validation/DATA/CIBERSORTx/GEPinference/MDS_MEAN_TOP_150_NO_NORM/CIBERSORTxHiRes_Job18_GMP_Window6.txt",
                "/home/nlegarra/data/nlegarra/DigitalCytometry_Validation/DATA/CIBERSORTx/GEPinference/MDS_MEAN_TOP_150_NO_NORM/CIBERSORTxHiRes_Job18_MEP_Window6.txt",
                pseudo_mds_mean,top_150_mds,"/home/nlegarra/data/nlegarra/DigitalCytometry_Validation/FIGURES/heatmap_mds_top150_pseudomean.pdf","pseudo-bulk (mean) mds top 150")

########################################################### HEATMAP TOP 150 SUM ###################################################################

plot_heatmap("/home/nlegarra/data/nlegarra/DigitalCytometry_Validation/DATA/CIBERSORTx/GEPinference/HEALTHY_SUM_TOP_150_NO_NORM/CIBERSORTxHiRes_Job22_HSC_Window5.txt",
                "/home/nlegarra/data/nlegarra/DigitalCytometry_Validation/DATA/CIBERSORTx/GEPinference/HEALTHY_SUM_TOP_150_NO_NORM/CIBERSORTxHiRes_Job22_GMP_Window5.txt",
                "/home/nlegarra/data/nlegarra/DigitalCytometry_Validation/DATA/CIBERSORTx/GEPinference/HEALTHY_SUM_TOP_150_NO_NORM/CIBERSORTxHiRes_Job22_MEP_Window5.txt",
                pseudo_healthy_sum,top_150_healthy,"/home/nlegarra/data/nlegarra/DigitalCytometry_Validation/FIGURES/heatmap_healthy_top150_pseudosum.pdf","pseudo-bulk (sum) healthy top 150")

plot_heatmap("/home/nlegarra/data/nlegarra/DigitalCytometry_Validation/DATA/CIBERSORTx/GEPinference/MDS_SUM_TOP_150_S_MODE/CIBERSORTxHiRes_Job20_HSC_Window6.txt",
                "/home/nlegarra/data/nlegarra/DigitalCytometry_Validation/DATA/CIBERSORTx/GEPinference/MDS_SUM_TOP_150_S_MODE/CIBERSORTxHiRes_Job20_GMP_Window6.txt",
                "/home/nlegarra/data/nlegarra/DigitalCytometry_Validation/DATA/CIBERSORTx/GEPinference/MDS_SUM_TOP_150_S_MODE/CIBERSORTxHiRes_Job20_MEP_Window6.txt",
                pseudo_mds_sum,top_150_mds,"/home/nlegarra/data/nlegarra/DigitalCytometry_Validation/FIGURES/heatmap_mds_top150_pseudosum.pdf","pseudo-bulk (sum) mds top 150")

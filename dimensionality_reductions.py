######################################################################################################################
## AUTHOR: NAROA LEGARRA MARCOS
## DATE: 10/03/22

## Script to make PCA+TSNE in sc, pseudobulk and pseudobulk replicas with the CUN dataset, Granja dataset and Processed Granja dataset. 
######################################################################################################################

wd = "/home/nlegarra/data/nlegarra/DigitalCytometry_Validation" #set the working directory
figures = "/home/nlegarra/data/nlegarra/DigitalCytometry_Validation/FIGURES" #set folder to store heatmaps

############################################# IMPORT LIBRARIES ##########################################
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import pickle5
import numpy as np
import pandas as pd
import random
import re

labels=["HSC","MEP","GMP"]

def plot_reduction(sc,sc_labels,pseudomean,pseudosum,pseudorepmean,pseudorepsum,file,dim1,dim2):
    #####################################################################################
    # Function generates the plot with the first two components of the dimensionality reduction
  
    # INPUT VARIABLES
    # sc <- matrix with the reduction of the sc data
    # sc_labels <- vector with cell type labels in sc matrix
    # pseudomean <- matrix with the reduction of the pseudobulk (mean) data
    # pseudosum <- matrix with the reduction of the pseudobulk (sum) data
  
    # OUTPUT VARIABLES
    # plot
    #####################################################################################
    colormap=["#8d70c9","#6ca74d","#c8588c"] 
    colormap9=["#8d70c9","#6ca74d","#c8588c","#8d70c9","#6ca74d","#c8588c","#8d70c9","#6ca74d","#c8588c"] 
    figura=plt.figure(figsize=(10,5))
    ((axs1,axs2))=figura.subplots(1,2)
    i=0
    for cell_type in labels:
        axs1.scatter(sc[sc_labels==cell_type,0],sc[sc_labels==cell_type,1],label=cell_type,c=colormap[i],s=8,marker=".",alpha=0.1)
        i+=1
    axs1.scatter(pseudorepmean[:,0],pseudorepmean[:,1],marker='+',s=8,c=colormap9,label='Pseudobulk replicas')
    axs1.scatter(pseudomean[:,0],pseudomean[:,1],marker='*',s=8,c=colormap,label='Pseudobulk')
    axs1.set_title('Pseudobulk (mean) + sc',fontsize=12)
    axs1.tick_params(axis='both', which='major', labelsize=4)
    axs1.tick_params(axis='both', which='minor', labelsize=4)
    axs1.set_xlabel(dim1,fontsize=8)
    axs1.set_ylabel(dim2,fontsize=8)
    i=0
    for cell_type in labels:
        axs2.scatter(sc[sc_labels==cell_type,0],sc[sc_labels==cell_type,1],label=cell_type,c=colormap[i],s=8,marker=".",alpha=0.1)
        i+=1
    axs2.scatter(pseudorepsum[:,0],pseudorepsum[:,1],marker='+',s=5,c=colormap9,label='Pseudobulk replicas')
    axs2.scatter(pseudosum[:,0],pseudosum[:,1],marker='*',s=8,c=colormap,label='Pseudobulk')
    axs2.set_title('Pseudobulk (sum) + sc',fontsize=12)
    axs2.tick_params(axis='both', which='major', labelsize=4)
    axs2.tick_params(axis='both', which='minor', labelsize=4)
    axs2.set_xlabel(dim1,fontsize=8)
    axs2.set_ylabel(dim2,fontsize=8)
    axs2.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=8,markerscale=2)
    plt.tight_layout()
    plt.savefig(file)
    plt.clf()

def plot_reduction_all(all_matrix,sc_matrix,name):
    #####################################################################################
    # Function generates the dimensionality reductions and saves the plots with the first 
    # two components of the dimensionality reduction
  
    # INPUT VARIABLES
    # all_matrix <- fused matrix
    # sc_matrix <- matrix with only sc data
    # name <- name of the experiment
  
    #####################################################################################
    sc_labels=sc_matrix.columns.values
    sc_labels=np.array([label[0:3] for label in sc_labels])
    # TSNE REDUCTION -----------------------------------------------------------------------------------------
    # Perform the TSNE REDUCTION
    tsne_reduction = TSNE(n_components=2,init='random').fit_transform(all_matrix.T)
    sc_tsne=tsne_reduction[0:sc_matrix.shape[1],:]
    pseudobulk_mean_tsne=tsne_reduction[sc_matrix.shape[1]:sc_matrix.shape[1]+3,:]
    pseudobulk_sum_tsne=tsne_reduction[sc_matrix.shape[1]+3:sc_matrix.shape[1]+6,:]
    pseudobulk_rep_mean_tsne=tsne_reduction[sc_matrix.shape[1]+6:sc_matrix.shape[1]+15,:]
    pseudobulk_rep_sum_tsne=tsne_reduction[sc_matrix.shape[1]+15:sc_matrix.shape[1]+24,:]
    # PCA REDUCTION ------------------------------------------------------------------------------------------
    pca_reduction = PCA().fit_transform(all_matrix.T) # perform PCA 
    sc_pca=pca_reduction[0:sc_matrix.shape[1],:]
    pseudobulk_mean_pca=pca_reduction[sc_matrix.shape[1]:sc_matrix.shape[1]+3,:]
    pseudobulk_sum_pca=pca_reduction[sc_matrix.shape[1]+3:sc_matrix.shape[1]+6,:]
    pseudobulk_rep_mean_pca=pca_reduction[sc_matrix.shape[1]+6:sc_matrix.shape[1]+15,:]
    pseudobulk_rep_sum_pca=pca_reduction[sc_matrix.shape[1]+15:sc_matrix.shape[1]+24,:]
    # PCA+TSNE REDUCTION -------------------------------------------------------------------------------------
    tsne_pca_reduction = TSNE(n_components=2,init='random').fit_transform(pca_reduction[:,0:100]) # perform TSNE reduction in the first 5 dimensions
    sc_tsne_pca=tsne_pca_reduction[0:sc_matrix.shape[1],:]
    pseudobulk_mean_tsne_pca=tsne_pca_reduction[sc_matrix.shape[1]:sc_matrix.shape[1]+3,:]
    pseudobulk_sum_tsne_pca=tsne_pca_reduction[sc_matrix.shape[1]+3:sc_matrix.shape[1]+6,:]
    pseudobulk_rep_mean_tsne_pca=tsne_pca_reduction[sc_matrix.shape[1]+6:sc_matrix.shape[1]+15,:]
    pseudobulk_rep_sum_tsne_pca=tsne_pca_reduction[sc_matrix.shape[1]+15:sc_matrix.shape[1]+24,:]
    # PLOT ---------------------------------------------------------------------------------------------------
    plot_reduction(sc_tsne,sc_labels,pseudobulk_mean_tsne,pseudobulk_sum_tsne,pseudobulk_rep_mean_tsne,pseudobulk_rep_sum_tsne,figures+'/TSNE_'+name+'.pdf','tSNE1','tSNE2')
    plot_reduction(sc_pca,sc_labels,pseudobulk_mean_pca,pseudobulk_sum_pca,pseudobulk_rep_mean_pca,pseudobulk_rep_sum_pca,figures+'/PCA_'+name+'.pdf','PC1','PC2')
    plot_reduction(sc_tsne_pca,sc_labels,pseudobulk_mean_tsne_pca,pseudobulk_sum_tsne_pca,pseudobulk_rep_mean_tsne_pca,pseudobulk_rep_sum_tsne_pca,figures+'/TSNE+PCA_'+name+'.pdf','tSNE1','tSNE2')

def replicate3(matrix,percentage):
   #####################################################################################
    # Function generates the plot with the first two components of the dimensionality reduction
  
    # INPUT VARIABLES
    # matrix <– matrix from which to make replicas
    # percentage <- percentage of cells for subsampling
  
    # OUTPUT VARIABLES
    # reps_mean <- 3 pseudobulk replicates made by the mean
    # reps_sum <- 3 pseudobulk replicates made by the sum
    #####################################################################################
    rep1=matrix.iloc[:,random.sample(range(0, matrix.shape[1]),round(matrix.shape[1]*percentage))]
    rep1_mean=np.mean(np.array(rep1),axis=1)
    rep1_sum=np.sum(np.array(rep1),axis=1)
    rep2=matrix.iloc[:,random.sample(range(0, matrix.shape[1]),round(matrix.shape[1]*percentage))]
    rep2_mean=np.mean(np.array(rep2),axis=1)
    rep2_sum=np.sum(np.array(rep2),axis=1)
    rep3=matrix.iloc[:,random.sample(range(0, matrix.shape[1]),round(matrix.shape[1]*percentage))]
    rep3_mean=np.mean(np.array(rep3),axis=1)
    rep3_sum=np.sum(np.array(rep3),axis=1)
    reps_mean=np.vstack((rep1_mean,rep2_mean,rep3_mean))
    reps_sum=np.vstack((rep1_sum,rep2_sum,rep3_sum))
    return reps_mean,reps_sum

############################################ CUN DATASET #############################################

# Read and process pseudobulks

pseudo_senior_mean=pd.read_table(wd+"/DATA/pseudo/pseudo_senior_mean.tsv",sep="\t").fillna(0) #read pseudobulk data
pseudo_senior_mean=pseudo_senior_mean.set_index("genes_senior") #set genes as index
pseudo_senior_mean=(pseudo_senior_mean-pseudo_senior_mean.min())/(pseudo_senior_mean.max()-pseudo_senior_mean.min()) #min-max normalization
pseudo_senior1_mean=pseudo_senior_mean.iloc[:,0:3] #select senior1 sample
pseudo_senior2_mean=pseudo_senior_mean.iloc[:,3:6] #select senior2 sample
pseudo_senior3_mean=pseudo_senior_mean.iloc[:,6:9] #select senior3 sample

pseudo_senior_sum=pd.read_table(wd+"/DATA/pseudo/pseudo_senior_sum.tsv",sep="\t").fillna(0)
pseudo_senior_sum=pseudo_senior_sum.set_index("genes_senior")
pseudo_senior_sum=(pseudo_senior_sum-pseudo_senior_sum.min())/(pseudo_senior_sum.max()-pseudo_senior_sum.min())
pseudo_senior1_sum=pseudo_senior_sum.iloc[:,0:3]
pseudo_senior2_sum=pseudo_senior_sum.iloc[:,3:6]
pseudo_senior3_sum=pseudo_senior_sum.iloc[:,6:9]

pseudo_mds_mean=pd.read_table(wd+"/DATA/pseudo/pseudo_mds_mean.tsv",sep="\t").fillna(0)
pseudo_mds_mean=pseudo_mds_mean.set_index("genes_mds")
pseudo_mds_mean=(pseudo_mds_mean-pseudo_mds_mean.min())/(pseudo_mds_mean.max()-pseudo_mds_mean.min())
pseudo_mds1_mean=pseudo_mds_mean.iloc[:,0:3]
pseudo_mds2_mean=pseudo_mds_mean.iloc[:,3:6]
pseudo_mds3_mean=pseudo_mds_mean.iloc[:,6:9]
pseudo_mds4_mean=pseudo_mds_mean.iloc[:,9:12]

pseudo_mds_sum=pd.read_table(wd+"/DATA/pseudo/pseudo_mds_sum.tsv",sep="\t").fillna(0)
pseudo_mds_sum=pseudo_mds_sum.set_index("genes_mds")
pseudo_mds_sum=(pseudo_mds_sum-pseudo_mds_sum.min())/(pseudo_mds_sum.max()-pseudo_mds_sum.min())
pseudo_mds1_sum=pseudo_mds_sum.iloc[:,0:3]
pseudo_mds2_sum=pseudo_mds_sum.iloc[:,3:6]
pseudo_mds3_sum=pseudo_mds_sum.iloc[:,6:9]
pseudo_mds4_sum=pseudo_mds_sum.iloc[:,9:12]

pseudo_senior_rep_mean=pd.read_table(wd+"/DATA/pseudo/pseudo_senior_rep_mean_0.2.tsv",sep="\t").fillna(0) #read pseudobulk data
pseudo_senior_rep_mean=pseudo_senior_rep_mean.set_index("genes_senior") #set genes as index
pseudo_senior_rep_mean=(pseudo_senior_rep_mean-pseudo_senior_rep_mean.min())/(pseudo_senior_rep_mean.max()-pseudo_senior_rep_mean.min()) #min-max normalization
pseudo_senior1_rep_mean=pseudo_senior_rep_mean.iloc[:,0:9] #select senior1 sample
pseudo_senior2_rep_mean=pseudo_senior_rep_mean.iloc[:,9:18] #select senior2 sample
pseudo_senior3_rep_mean=pseudo_senior_rep_mean.iloc[:,18:27] #select senior3 sample

pseudo_senior_rep_sum=pd.read_table(wd+"/DATA/pseudo/pseudo_senior_rep_sum_0.2.tsv",sep="\t").fillna(0)
pseudo_senior_rep_sum=pseudo_senior_rep_sum.set_index("genes_senior")
pseudo_senior_rep_sum=(pseudo_senior_rep_sum-pseudo_senior_rep_sum.min())/(pseudo_senior_rep_sum.max()-pseudo_senior_rep_sum.min())
pseudo_senior1_rep_sum=pseudo_senior_rep_sum.iloc[:,0:9]
pseudo_senior2_rep_sum=pseudo_senior_rep_sum.iloc[:,9:18]
pseudo_senior3_rep_sum=pseudo_senior_rep_sum.iloc[:,18:27]

pseudo_mds_rep_mean=pd.read_table(wd+"/DATA/pseudo/pseudo_mds_rep_mean_0.2.tsv",sep="\t").fillna(0)
pseudo_mds_rep_mean=pseudo_mds_rep_mean.set_index("genes_mds")
pseudo_mds_rep_mean=(pseudo_mds_rep_mean-pseudo_mds_rep_mean.min())/(pseudo_mds_rep_mean.max()-pseudo_mds_rep_mean.min())
pseudo_mds1_rep_mean=pseudo_mds_rep_mean.iloc[:,0:9]
pseudo_mds2_rep_mean=pseudo_mds_rep_mean.iloc[:,9:18]
pseudo_mds3_rep_mean=pseudo_mds_rep_mean.iloc[:,18:27]
pseudo_mds4_rep_mean=pseudo_mds_rep_mean.iloc[:,27:36]

pseudo_mds_rep_sum=pd.read_table(wd+"/DATA/pseudo/pseudo_mds_rep_sum_0.2.tsv",sep="\t").fillna(0)
pseudo_mds_rep_sum=pseudo_mds_rep_sum.set_index("genes_mds")
pseudo_mds_rep_sum=(pseudo_mds_rep_sum-pseudo_mds_rep_sum.min())/(pseudo_mds_rep_sum.max()-pseudo_mds_rep_sum.min())
pseudo_mds1_rep_sum=pseudo_mds_rep_sum.iloc[:,0:9]
pseudo_mds2_rep_sum=pseudo_mds_rep_sum.iloc[:,9:18]
pseudo_mds3_rep_sum=pseudo_mds_rep_sum.iloc[:,18:27]
pseudo_mds4_rep_sum=pseudo_mds_rep_sum.iloc[:,27:36]

# Read sc

senior1=pd.read_table(wd+"/DATA/single_cells/senior1.csv",sep=" ").fillna(0) #read sc data
senior1=(senior1-senior1.min())/(senior1.max()-senior1.min()) #min-max normalization
senior1_all=pd.concat([senior1, pseudo_senior1_mean,pseudo_senior1_sum,pseudo_senior1_rep_mean,pseudo_senior1_rep_sum], axis=1,join="inner") #concatenate all data
plot_reduction_all(senior1_all,senior1,'senior1_0.2') #generate the reductions and the plots

senior2=pd.read_table(wd+"/DATA/single_cells/senior2.csv",sep=" ").fillna(0)
senior2=(senior2-senior2.min())/(senior2.max()-senior2.min())
senior2_all=pd.concat([senior2, pseudo_senior2_mean, pseudo_senior2_sum,pseudo_senior2_rep_mean,pseudo_senior2_rep_sum], axis=1,join="inner")
plot_reduction_all(senior2_all,senior2,'senior2_0.2')

senior3=pd.read_table(wd+"/DATA/single_cells/senior3.csv",sep=" ").fillna(0)
senior3=(senior3-senior3.min())/(senior3.max()-senior3.min())
senior3_all=pd.concat([senior3, pseudo_senior3_mean, pseudo_senior3_sum,pseudo_senior3_rep_mean,pseudo_senior3_rep_sum], axis=1,join="inner")
plot_reduction_all(senior3_all,senior3,'senior3_0.2')

mds1=pd.read_table(wd+"/DATA/single_cells/mds1.csv",sep=" ").fillna(0)
mds1=(mds1-mds1.min())/(mds1.max()-mds1.min())
mds1_all=pd.concat([mds1, pseudo_mds1_mean, pseudo_mds1_sum,pseudo_mds1_rep_mean,pseudo_mds1_rep_sum], axis=1, join="inner")
plot_reduction_all(mds1_all,mds1,'mds1_0.2')

mds2=pd.read_table(wd+"/DATA/single_cells/mds3.csv",sep=" ").fillna(0)
mds2=(mds2-mds2.min())/(mds2.max()-mds2.min())
mds2_all=pd.concat([mds2, pseudo_mds2_mean, pseudo_mds2_sum,pseudo_mds2_rep_mean,pseudo_mds2_rep_sum], axis=1, join="inner")
plot_reduction_all(mds2_all,mds2,'mds3_0.2')

mds3=pd.read_table(wd+"/DATA/single_cells/mds5.csv",sep=" ").fillna(0)
mds3=(mds3-mds3.min())/(mds3.max()-mds3.min())
mds3_all=pd.concat([mds3, pseudo_mds3_mean, pseudo_mds3_sum,pseudo_mds3_rep_mean,pseudo_mds3_rep_sum], axis=1, join="inner")
plot_reduction_all(mds3_all,mds3,'mds5_0.2')

mds4=pd.read_table(wd+"/DATA/single_cells/mds10.csv",sep=" ").fillna(0)
mds4=(mds4-mds4.min())/(mds4.max()-mds4.min())
mds4_all=pd.concat([mds4, pseudo_mds4_mean, pseudo_mds4_sum,pseudo_mds4_rep_mean,pseudo_mds4_rep_sum], axis=1, join="inner")
plot_reduction_all(mds4_all,mds4,'mds10_0.2')

################################################ GRANJA ####################################################################

granja=pickle5.load(open(wd+'/DATA/single_cells/granja_scrna.pickle','rb')) #read sc data
granja.columns=granja.iloc[0,:] #set cell types as column index
cells=['01_HSC','25_NK','12_CD14.Mono.2','05_CMP.LMPP','08_GMP.Neut'] #selected cell typess

granja_HSC=granja.loc[:,np.array((granja.iloc[0,:]=='01_HSC').reset_index(drop=True))] #select cells of each cell type in the sc
granja_NK=granja.loc[:,np.array((granja.iloc[0,:]=='25_NK').reset_index(drop=True))]
granja_Mono=granja.loc[:,np.array((granja.iloc[0,:]=='12_CD14.Mono.2').reset_index(drop=True))]
granja_CMP=granja.loc[:,np.array((granja.iloc[0,:]=='05_CMP.LMPP').reset_index(drop=True))]
granja_GMP=granja.loc[:,np.array((granja.iloc[0,:]=='08_GMP.Neut').reset_index(drop=True))]

#paste to obtain the sc with the selected cells
sc=pd.concat([granja_HSC.iloc[1:,:],granja_NK.iloc[1:,:],granja_Mono.iloc[1:,:],granja_CMP.iloc[1:,:],granja_GMP.iloc[1:,:]],axis=1).astype(int)
sc=(sc-sc.min())/(sc.max()-sc.min()) #min max normalization
#paste labels to obtain final labels
sc_labels=pd.concat([granja_HSC.iloc[0,:],granja_NK.iloc[0,:],granja_Mono.iloc[0,:],granja_CMP.iloc[0,:],granja_GMP.iloc[0,:]],axis=0)

#generate pseudobulk for each cell type (sum)
granja_HSC_sum=np.sum(np.array(granja_HSC.iloc[1:,:]).astype(int),axis=1) 
granja_NK_sum=np.sum(np.array(granja_NK.iloc[1:,:]).astype(int),axis=1)
granja_Mono_sum=np.sum(np.array(granja_Mono.iloc[1:,:]).astype(int),axis=1)
granja_CMP_sum=np.sum(np.array(granja_CMP.iloc[1:,:]).astype(int),axis=1)
granja_GMP_sum=np.sum(np.array(granja_GMP.iloc[1:,:]).astype(int),axis=1)

pseudosum=np.vstack((granja_HSC_sum,granja_NK_sum,granja_Mono_sum,granja_CMP_sum,granja_GMP_sum)) #stack pseudobulk
pseudosum=pd.DataFrame(pseudosum.T) #transpose because the sc has different orientation
pseudosum=(pseudosum-pseudosum.min())/(pseudosum.max()-pseudosum.min())  #min max normalization
pseudosum.columns=cells #add column labels
pseudosum.index=sc.index #add gene names

#generate pseudobulk for each cell type (mean)
granja_HSC_mean=np.mean(np.array(granja_HSC.iloc[1:,:]).astype(int),axis=1)
granja_NK_mean=np.mean(np.array(granja_NK.iloc[1:,:]).astype(int),axis=1)
granja_Mono_mean=np.mean(np.array(granja_Mono.iloc[1:,:]).astype(int),axis=1)
granja_CMP_mean=np.mean(np.array(granja_CMP.iloc[1:,:]).astype(int),axis=1)
granja_GMP_mean=np.mean(np.array(granja_GMP.iloc[1:,:]).astype(int),axis=1)

pseudomean=np.vstack((granja_HSC_mean,granja_NK_mean,granja_Mono_mean,granja_CMP_mean,granja_GMP_mean))
pseudomean=pd.DataFrame(pseudomean.T)
pseudomean=(pseudomean-pseudomean.min())/(pseudomean.max()-pseudomean.min())
pseudomean.columns=cells
pseudomean.index=sc.index

sc_all=pd.concat([sc,pseudosum,pseudomean],axis=1) #concatenate all the matrices
pca_reduction = PCA(n_components=100).fit_transform(sc_all.T) # perform PCA 
tsne_pca_reduction = TSNE(n_components=2,init='random').fit_transform(pca_reduction) # perform TSNE reduction in the first 100 dimensions
#separate each matrix
sc_tsne_pca=tsne_pca_reduction[0:sc.shape[1],:]
pseudobulk_sum_tsne_pca_sum=tsne_pca_reduction[sc.shape[1]:sc.shape[1]+5,:]
pseudobulk_sum_tsne_pca_mean=tsne_pca_reduction[sc.shape[1]+5:sc.shape[1]+10,:]

colormap=["#8d70c9","#6ca74d","#c8588c","#49adad","#cc5a43"] #set colormap

figura=plt.figure(figsize=(10,5)) #open figure
((axs1,axs2))=figura.subplots(1,2) #set axis

i=0 #set index as 0
#plot sc
for cell_type in cells:
   axs1.scatter(sc_tsne_pca[sc_labels==cell_type,0],sc_tsne_pca[sc_labels==cell_type,1],label=cell_type,c=colormap[i],s=8,marker=".",alpha=0.1)
   i+=1

#plot pseudobulk
axs1.scatter(pseudobulk_sum_tsne_pca_mean[:,0],pseudobulk_sum_tsne_pca_mean[:,1],marker='*',s=25,c=colormap,label='Pseudobulk')
axs1.set_title('Pseudobulk (mean) + sc',fontsize=12)
axs1.tick_params(axis='both', which='major', labelsize=4)
axs1.tick_params(axis='both', which='minor', labelsize=4)
axs1.set_xlabel('tSNE1',fontsize=8)
axs1.set_ylabel('tSNE2',fontsize=8)

i=0
for cell_type in cells:
   axs2.scatter(sc_tsne_pca[sc_labels==cell_type,0],sc_tsne_pca[sc_labels==cell_type,1],label=cell_type,c=colormap[i],s=8,marker=".",alpha=0.1)
   i+=1

axs2.scatter(pseudobulk_sum_tsne_pca_sum[:,0],pseudobulk_sum_tsne_pca_sum[:,1],marker='*',s=25,c=colormap,label='Pseudobulk')
axs2.set_title('Pseudobulk (sum) + sc',fontsize=12)
axs2.tick_params(axis='both', which='major', labelsize=4)
axs2.tick_params(axis='both', which='minor', labelsize=4)
axs2.set_xlabel('tSNE1',fontsize=8)
axs2.set_ylabel('tSNE2',fontsize=8)
axs2.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=8,markerscale=2)
plt.tight_layout()
plt.savefig(figures+"/tSNE+PCA_granja.pdf")
plt.clf()

##################################################################################################################

percentage=0.001 #set subsamplig percentage

#make replicas
granja_HSC_rep_mean,granja_HSC_rep_sum=replicate3(granja_HSC.iloc[1:,:].astype(int),percentage)
granja_NK_rep_mean,granja_NK_rep_sum=replicate3(granja_NK.iloc[1:,:].astype(int),percentage)
granja_Mono_rep_mean,granja_Mono_rep_sum=replicate3(granja_Mono.iloc[1:,:].astype(int),percentage)
granja_CMP_rep_mean,granja_CMP_rep_sum=replicate3(granja_CMP.iloc[1:,:].astype(int),percentage)
granja_GMP_rep_mean,granja_GMP_rep_sum=replicate3(granja_GMP.iloc[1:,:].astype(int),percentage)

#concatenate replicas
pseudomean_rep=pd.DataFrame(np.vstack((granja_HSC_rep_mean,granja_NK_rep_mean,granja_Mono_rep_mean,granja_CMP_rep_mean,granja_GMP_rep_mean)).T)
pseudomean_rep=(pseudomean_rep-pseudomean_rep.min())/(pseudomean_rep.max()-pseudomean_rep.min()) #min max noralization
pseudomean_rep.index=sc.index #set gene names

pseudosum_rep=pd.DataFrame(np.vstack((granja_HSC_rep_sum,granja_NK_rep_sum,granja_Mono_rep_sum,granja_CMP_rep_sum,granja_GMP_rep_sum)).T)
pseudosum_rep=(pseudosum_rep-pseudosum_rep.min())/(pseudosum_rep.max()-pseudosum_rep.min())
pseudosum_rep.index=sc.index

sc_all=pd.concat([sc,pseudosum,pseudomean,pseudomean_rep,pseudosum_rep],axis=1).fillna(0)
pca_reduction = PCA(n_components=100).fit_transform(sc_all.T) # perform PCA 
tsne_pca_reduction = TSNE(n_components=2,init='random').fit_transform(pca_reduction) # perform TSNE reduction in the first 100 dimensions
sc_tsne_pca=tsne_pca_reduction[0:sc.shape[1],:]
pseudobulk_tsne_pca_sum=tsne_pca_reduction[sc.shape[1]:sc.shape[1]+5,:]
pseudobulk_tsne_pca_mean=tsne_pca_reduction[sc.shape[1]+5:sc.shape[1]+10,:]
pseudobulk_tsne_pca_repsum=tsne_pca_reduction[sc.shape[1]+10:sc.shape[1]+25,:]
pseudobulk_tsne_pca_repmean=tsne_pca_reduction[sc.shape[1]+25:sc.shape[1]+40,:]

colormap=["#8d70c9","#6ca74d","#c8588c","#49adad","#cc5a43"]

#set colormap in order for the replicas
colormap_rep=["#8d70c9","#8d70c9","#8d70c9","#6ca74d","#6ca74d","#6ca74d","#c8588c","#c8588c","#c8588c","#49adad","#49adad","#49adad","#cc5a43","#cc5a43","#cc5a43"]

figura=plt.figure(figsize=(10,5))
((axs1,axs2))=figura.subplots(1,2)
i=0
for cell_type in cells:
   axs1.scatter(sc_tsne_pca[sc_labels==cell_type,0],sc_tsne_pca[sc_labels==cell_type,1],label=cell_type,c=colormap[i],s=8,marker=".",alpha=0.05)
   i+=1

axs1.scatter(pseudobulk_tsne_pca_mean[:,0],pseudobulk_tsne_pca_mean[:,1],marker='*',s=25,c=colormap,label='Pseudobulk')
axs1.scatter(pseudobulk_tsne_pca_repmean[:,0],pseudobulk_tsne_pca_repmean[:,1],marker='^',s=25,c=colormap_rep,label='Pseudobulk replicas 0.001')
axs1.set_title('Pseudobulk (mean) + sc',fontsize=12)
axs1.tick_params(axis='both', which='major', labelsize=4)
axs1.tick_params(axis='both', which='minor', labelsize=4)
axs1.set_xlabel('tSNE1',fontsize=8)
axs1.set_ylabel('tSNE2',fontsize=8)

i=0
for cell_type in cells:
   axs2.scatter(sc_tsne_pca[sc_labels==cell_type,0],sc_tsne_pca[sc_labels==cell_type,1],label=cell_type,c=colormap[i],s=8,marker=".",alpha=0.05)
   i+=1

axs2.scatter(pseudobulk_tsne_pca_sum[:,0],pseudobulk_tsne_pca_sum[:,1],marker='*',s=25,c=colormap,label='Pseudobulk')
axs2.scatter(pseudobulk_tsne_pca_repsum[:,0],pseudobulk_tsne_pca_repsum[:,1],marker='^',s=25,c=colormap_rep,label='Pseudobulk replicas 0.001')
axs2.set_title('Pseudobulk (sum) + sc',fontsize=12)
axs2.tick_params(axis='both', which='major', labelsize=4)
axs2.tick_params(axis='both', which='minor', labelsize=4)
axs2.set_xlabel('tSNE1',fontsize=8)
axs2.set_ylabel('tSNE2',fontsize=8)
axs2.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=8,markerscale=2)
plt.tight_layout()
plt.savefig(figures+"/tSNE+PCA_granja_0.001.pdf")
plt.clf()

################################### PROCCESSED GRANJA ###############################################################################

granja_processed=pickle5.load(open(wd+'/DATA/single_cells/scrna_granja_Scanpy_LS_sNN.pickle','rb')) #read sc data
granja_processed_labels=pd.read_table(wd+'/DATA/single_cells/scrna_granja_Scanpy_LS_labels.txt') #read sc labels
granja_processed.columns=granja_processed_labels.iloc[:,1] #set labels as columns

cells=['Natural killer cell','Plasmacytoid dendritic cell','Dendritic cell','B cell','CD14+ monocyte','CD4+ T cell']

granjap_NK=granja_processed.loc[:,np.array((granja_processed_labels.iloc[:,1]=='Natural killer cell').reset_index(drop=True))]
granjap_pl_dendri=granja_processed.loc[:,np.array((granja_processed_labels.iloc[:,1]=='Plasmacytoid dendritic cell').reset_index(drop=True))]
granjap_dendri=granja_processed.loc[:,np.array((granja_processed_labels.iloc[:,1]=='Dendritic cell').reset_index(drop=True))]
granjap_Bcell=granja_processed.loc[:,np.array((granja_processed_labels.iloc[:,1]=='B cell').reset_index(drop=True))]
granjap_Mono=granja_processed.loc[:,np.array((granja_processed_labels.iloc[:,1]=='CD14+ monocyte').reset_index(drop=True))]
granjap_Tcell=granja_processed.loc[:,np.array((granja_processed_labels.iloc[:,1]=='CD4+ T cell').reset_index(drop=True))]

sc=pd.concat([granjap_NK,granjap_pl_dendri,granjap_dendri,granjap_Bcell,granjap_Mono,granjap_Tcell],axis=1)
sc=(sc-sc.min())/(sc.max()-sc.min())
sc_labels=sc.columns

granjap_NK_sum=np.sum(np.array(granjap_NK),axis=1)
granjap_pl_dendri_sum=np.sum(np.array(granjap_pl_dendri),axis=1)
granjap_dendri_sum=np.sum(np.array(granjap_dendri),axis=1)
granjap_Bcell_sum=np.sum(np.array(granjap_Bcell),axis=1)
granjap_Mono_sum=np.sum(np.array(granjap_Mono),axis=1)
granjap_Tcell_sum=np.sum(np.array(granjap_Tcell),axis=1)

pseudosum=np.vstack((granjap_NK_sum,granjap_pl_dendri_sum,granjap_dendri_sum,granjap_Bcell_sum,granjap_Mono_sum,granjap_Tcell_sum))
pseudosum=pd.DataFrame(pseudosum.T)
pseudosum=(pseudosum-pseudosum.min())/(pseudosum.max()-pseudosum.min())
pseudosum.columns=cells
pseudosum.index=sc.index

granjap_NK_mean=np.mean(np.array(granjap_NK),axis=1)
granjap_pl_dendri_mean=np.mean(np.array(granjap_pl_dendri),axis=1)
granjap_dendri_mean=np.mean(np.array(granjap_dendri),axis=1)
granjap_Bcell_mean=np.mean(np.array(granjap_Bcell),axis=1)
granjap_Mono_mean=np.mean(np.array(granjap_Mono),axis=1)
granjap_Tcell_mean=np.mean(np.array(granjap_Tcell),axis=1)

pseudomean=np.vstack((granjap_NK_mean,granjap_pl_dendri_mean,granjap_dendri_mean,granjap_Bcell_mean,granjap_Mono_mean,granjap_Tcell_mean))
pseudomean=pd.DataFrame(pseudomean.T)
pseudomean=(pseudomean-pseudomean.min())/(pseudomean.max()-pseudomean.min())
pseudomean.columns=cells
pseudomean.index=sc.index

sc_all=pd.concat([sc,pseudosum,pseudomean],axis=1)
pca_reduction = PCA(n_components=50).fit_transform(sc_all.T) # perform PCA 
tsne_pca_reduction = TSNE(n_components=2,init='random').fit_transform(pca_reduction) # perform TSNE reduction in the first 100 dimensions
sc_tsne_pca=tsne_pca_reduction[0:sc.shape[1],:]
pseudobulk_tsne_pca_sum=tsne_pca_reduction[sc.shape[1]:sc.shape[1]+6,:]
pseudobulk_tsne_pca_mean=tsne_pca_reduction[sc.shape[1]+6:,:]

colormap=["#8d70c9","#6ca74d","#c8588c","#49adad","#cc5a43","#b69040"]

figura=plt.figure(figsize=(10,5))
((axs1,axs2))=figura.subplots(1,2)
i=0
for cell_type in cells:
   axs1.scatter(sc_tsne_pca[sc_labels==cell_type,0],sc_tsne_pca[sc_labels==cell_type,1],label=cell_type,c=colormap[i],s=8,marker=".",alpha=0.05)
   i+=1

axs1.scatter(pseudobulk_tsne_pca_mean[:,0],pseudobulk_tsne_pca_mean[:,1],marker='*',s=25,c=colormap,label='Pseudobulk')
axs1.set_title('Pseudobulk (mean) + sc',fontsize=12)
axs1.tick_params(axis='both', which='major', labelsize=4)
axs1.tick_params(axis='both', which='minor', labelsize=4)
axs1.set_xlabel('tSNE1',fontsize=8)
axs1.set_ylabel('tSNE2',fontsize=8)

i=0
for cell_type in cells:
   axs2.scatter(sc_tsne_pca[sc_labels==cell_type,0],sc_tsne_pca[sc_labels==cell_type,1],label=cell_type,c=colormap[i],s=8,marker=".",alpha=0.05)
   i+=1

axs2.scatter(pseudobulk_tsne_pca_sum[:,0],pseudobulk_tsne_pca_sum[:,1],marker='*',s=25,c=colormap,label='Pseudobulk')
axs2.set_title('Pseudobulk (sum) + sc',fontsize=12)
axs2.tick_params(axis='both', which='major', labelsize=4)
axs2.tick_params(axis='both', which='minor', labelsize=4)
axs2.set_xlabel('tSNE1',fontsize=8)
axs2.set_ylabel('tSNE2',fontsize=8)
axs2.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=8,markerscale=2)
plt.tight_layout()
plt.savefig(figures+"/tSNE+PCA_granja_p.pdf")
plt.clf()


##################################################################################################################

percentage=0.001

granjap_NK_rep_mean,granjap_NK_rep_sum=replicate3(granjap_NK,percentage)
granjap_pl_dendri_rep_mean,granjap_pl_dendri_rep_sum=replicate3(granjap_pl_dendri,percentage)
granjap_dendri_rep_mean,granjap_dendri_rep_sum=replicate3(granjap_dendri,percentage)
granjap_Bcell_rep_mean,granjap_Bcell_rep_sum=replicate3(granjap_Bcell,percentage)
granjap_Mono_rep_mean,granjap_Mono_rep_sum=replicate3(granjap_Mono,percentage)
granjap_Tcell_rep_mean,granjap_Tcell_rep_sum=replicate3(granjap_Tcell,percentage)

pseudomean_rep=pd.DataFrame(np.vstack((granjap_NK_rep_mean,granjap_pl_dendri_rep_mean,granjap_dendri_rep_mean,granjap_Bcell_rep_mean,granjap_Mono_rep_mean,granjap_Tcell_rep_mean)).T)
pseudomean_rep=(pseudomean_rep-pseudomean_rep.min())/(pseudomean_rep.max()-pseudomean_rep.min())
pseudosum_rep=pd.DataFrame(np.vstack((granjap_NK_rep_sum,granjap_pl_dendri_rep_sum,granjap_dendri_rep_sum,granjap_Bcell_rep_sum,granjap_Mono_rep_sum,granjap_Tcell_rep_sum)).T)
pseudosum_rep=(pseudosum_rep-pseudosum_rep.min())/(pseudosum_rep.max()-pseudosum_rep.min())
pseudomean_rep.index=sc.index
pseudosum_rep.index=sc.index

sc_all=pd.concat([sc,pseudosum,pseudomean,pseudomean_rep,pseudosum_rep],axis=1).fillna(0)
pca_reduction = PCA(n_components=100).fit_transform(sc_all.T) # perform PCA 
tsne_pca_reduction = TSNE(n_components=2,init='random').fit_transform(pca_reduction) # perform TSNE reduction in the first 100 dimensions
sc_tsne_pca=tsne_pca_reduction[0:sc.shape[1],:]
pseudobulk_tsne_pca_sum=tsne_pca_reduction[sc.shape[1]:sc.shape[1]+6,:]
pseudobulk_tsne_pca_mean=tsne_pca_reduction[sc.shape[1]+6:sc.shape[1]+12,:]
pseudobulk_tsne_pca_repsum=tsne_pca_reduction[sc.shape[1]+12:sc.shape[1]+30,:]
pseudobulk_tsne_pca_repmean=tsne_pca_reduction[sc.shape[1]+30:sc.shape[1]+48,:]

colormap=["#8d70c9","#6ca74d","#c8588c","#49adad","#cc5a43","#b69040"]
colormap_rep=["#8d70c9","#8d70c9","#8d70c9","#6ca74d","#6ca74d","#6ca74d","#c8588c","#c8588c","#c8588c","#49adad","#49adad","#49adad","#cc5a43","#cc5a43","#cc5a43","#b69040","#b69040","#b69040"]

figura=plt.figure(figsize=(10,5))
((axs1,axs2))=figura.subplots(1,2)
i=0
for cell_type in cells:
   axs1.scatter(sc_tsne_pca[sc_labels==cell_type,0],sc_tsne_pca[sc_labels==cell_type,1],label=cell_type,c=colormap[i],s=8,marker=".",alpha=0.05)
   i+=1

axs1.scatter(pseudobulk_tsne_pca_mean[:,0],pseudobulk_tsne_pca_mean[:,1],marker='*',s=25,c=colormap,label='Pseudobulk')
axs1.scatter(pseudobulk_tsne_pca_repmean[:,0],pseudobulk_tsne_pca_repmean[:,1],marker='^',s=25,c=colormap_rep,label='Pseudobulk replicas 0.001')
axs1.set_title('Pseudobulk (mean) + sc',fontsize=12)
axs1.tick_params(axis='both', which='major', labelsize=4)
axs1.tick_params(axis='both', which='minor', labelsize=4)
axs1.set_xlabel('tSNE1',fontsize=8)
axs1.set_ylabel('tSNE2',fontsize=8)

i=0
for cell_type in cells:
   axs2.scatter(sc_tsne_pca[sc_labels==cell_type,0],sc_tsne_pca[sc_labels==cell_type,1],label=cell_type,c=colormap[i],s=8,marker=".",alpha=0.05)
   i+=1

axs2.scatter(pseudobulk_tsne_pca_sum[:,0],pseudobulk_tsne_pca_sum[:,1],marker='*',s=25,c=colormap,label='Pseudobulk')
axs2.scatter(pseudobulk_tsne_pca_repsum[:,0],pseudobulk_tsne_pca_repsum[:,1],marker='^',s=25,c=colormap_rep,label='Pseudobulk replicas 0.001')
axs2.set_title('Pseudobulk (sum) + sc',fontsize=12)
axs2.tick_params(axis='both', which='major', labelsize=4)
axs2.tick_params(axis='both', which='minor', labelsize=4)
axs2.set_xlabel('tSNE1',fontsize=8)
axs2.set_ylabel('tSNE2',fontsize=8)
axs2.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=8,markerscale=2)
plt.tight_layout()
plt.savefig(figures+"/tSNE+PCA_granja_p_0.001.pdf")
plt.clf()


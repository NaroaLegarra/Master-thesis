################################################################################
## AUTHOR: NAROA LEGARRA MARCOS
## DATE: 10/03/22

## GENERATES PLOTS FOR MARKER GENE, ENRICHMENT ANALYSIS AND DEA COMPARISON
################################################################################

wd <- "~/Desktop/TFM" #set the working directory
pseudobulk <- "/DATA/pseudo/" #path to the folder to store the pseudobulks
figures <- "/FIGURES"

################################## READ DATA ###################################

#All samples
results_mean=readRDS(file=paste0(wd,"/DATA/DEA_results/results_mean_all_bulk.rds"))
results_sum=readRDS(file=paste0(wd,"/DATA/DEA_results/results_sum_all_bulk.rds"))

#Reduced bulk samples
results_mean=readRDS(file=paste0(wd,"/DATA/DEA_results/results_mean_4_bulk.rds"))
results_sum=readRDS(file=paste0(wd,"/DATA/DEA_results/results_sum_4_bulk.rds"))

results_mean=readRDS(file=paste0(wd,"/DATA/DEA_results/results_mean_4_rep_bulk.rds"))
results_sum=readRDS(file=paste0(wd,"/DATA/DEA_results/results_sum_4_rep_bulk.rds"))

results_mean=readRDS(file=paste0(wd,"/DATA/DEA_results/results_mean_4_rep2_bulk.rds"))
results_sum=readRDS(file=paste0(wd,"/DATA/DEA_results/results_sum_4_rep2_bulk.rds"))

########################### GENERATE MARKER GENE GRAPH #########################
# bind all results
results_all_markers <- rbind(process_results_markers(results_mean,"Pseudobulk (mean)"),
                             process_results_markers(results_sum,"Pseudobulk (sum)"))

results_all_markers$category <- str_replace_all(results_all_markers$category,"bulk","FACS-bulk RNAseq")
results_all_markers$category <- str_replace_all(results_all_markers$category,"sc","sc RNAseq")

# generate boxplot
PLOT <- results_all_markers %>% 
  ggplot(aes(x=category,y=value,fill=category))+
  geom_boxplot(alpha=0.5)+
  geom_jitter(aes(col=COMPARISON),size=0.5)+
  scale_fill_manual(values=c("#ff9dc5","#fae980"))+
  labs(x="",y="Percentage of marker genes detected",fill="Type of experiment")+
  facet_grid(cols = vars(EXPERIMENT))+
  theme_light()
PLOT

# save plot in pdf 
ggsave(plot = PLOT, width = 9, height = 6, dpi = 300, filename = paste0(wd,figures,"/MARKERS.pdf"))


############################## ENRICHMENT ANALYSIS #############################

HSC_MEP_mds_mean_GO <- GO_comparison("pseudobulk (mean)","HSC_MEP_mds",results_mean$HSC_MEP_bulkvssc_mds)
HSC_GMP_mds_mean_GO <- GO_comparison("pseudobulk (mean)","HSC_GMP_mds",results_mean$HSC_GMP_bulkvssc_mds)
MEP_GMP_mds_mean_GO <- GO_comparison("pseudobulk (mean)","MEP_GMP_mds",results_mean$MEP_GMP_bulkvssc_mds)
HSC_MEP_senior_mean_GO <- GO_comparison("pseudobulk (mean)","HSC_MEP_senior",results_mean$HSC_MEP_bulkvssc_senior)
HSC_GMP_senior_mean_GO <- GO_comparison("pseudobulk (mean)","HSC_GMP_senior",results_mean$HSC_GMP_bulkvssc_senior)
MEP_GMP_senior_mean_GO <- GO_comparison("pseudobulk (mean)","MEP_GMP_senior",results_mean$MEP_GMP_bulkvssc_senior)

HSC_MEP_mds_sum_GO <- GO_comparison("pseudobulk (sum)","HSC_MEP_mds",results_sum$HSC_MEP_bulkvssc_mds)
HSC_GMP_mds_sum_GO <- GO_comparison("pseudobulk (sum)","HSC_GMP_mds",results_sum$HSC_GMP_bulkvssc_mds)
MEP_GMP_mds_sum_GO <- GO_comparison("pseudobulk (sum)","MEP_GMP_mds",results_sum$MEP_GMP_bulkvssc_mds)
HSC_MEP_senior_sum_GO <- GO_comparison("pseudobulk (sum)","HSC_MEP_senior",results_sum$HSC_MEP_bulkvssc_senior)
HSC_GMP_senior_sum_GO <- GO_comparison("pseudobulk (sum)","HSC_GMP_senior",results_sum$HSC_GMP_bulkvssc_senior)
MEP_GMP_senior_sum_GO <- GO_comparison("pseudobulk (sum)","MEP_GMP_senior",results_sum$MEP_GMP_bulkvssc_senior)

GO_term_enrichment <- rbind(HSC_MEP_mds_mean_GO,HSC_GMP_mds_mean_GO,MEP_GMP_mds_mean_GO)
GOterms <- GO_term_enrichment %>% 
  mutate(sample=paste0(celltypes,"_",pseudo)) %>% 
  ggplot() + 
  geom_col(aes(x=experiment,y=total,fill="Unique GO terms"))+
  geom_col(aes(x=experiment,y=common,fill="common GO terms"))+
  theme_light()+
  theme(axis.text.x=element_text(angle = -70, hjust = 0))+
  labs(y="Significant GO term count",x=" ",fill="")+
  scale_fill_manual(values=c("#fae980","#ff9dc5","#82dde0"))+
  theme(plot.title = element_text(size=8))+
  facet_grid(sample~GO,scales='free')+
  theme(legend.position="top")
pdf(paste0(wd,figures,"/GO_terms_mean_mds.pdf"),height=11,width=10)
GOterms
dev.off()
GO_term_enrichment <- rbind(HSC_MEP_senior_mean_GO,HSC_GMP_senior_mean_GO,MEP_GMP_senior_mean_GO)
GOterms <- GO_term_enrichment %>% 
  mutate(sample=paste0(celltypes,"_",pseudo)) %>% 
  ggplot() + 
  geom_col(aes(x=experiment,y=total,fill="Unique GO terms"))+
  geom_col(aes(x=experiment,y=common,fill="common GO terms"))+
  theme_light()+
  theme(axis.text.x=element_text(angle = -70, hjust = 0))+
  labs(y="Significant GO term count",x=" ",fill="")+
  scale_fill_manual(values=c("#fae980","#ff9dc5","#82dde0"))+
  theme(plot.title = element_text(size=8))+
  facet_grid(sample~GO,scales='free')+
  theme(legend.position="top")
pdf(paste0(wd,figures,"/GO_terms_mean_senior.pdf"),height=11,width=10)
GOterms
dev.off()
GO_term_enrichment <- rbind(HSC_MEP_mds_sum_GO,HSC_GMP_mds_sum_GO,MEP_GMP_mds_sum_GO)
GOterms <- GO_term_enrichment %>% 
  mutate(sample=paste0(celltypes,"_",pseudo)) %>% 
  ggplot() + 
  geom_col(aes(x=experiment,y=total,fill="Unique GO terms"))+
  geom_col(aes(x=experiment,y=common,fill="common GO terms"))+
  theme_light()+
  theme(axis.text.x=element_text(angle = -70, hjust = 0))+
  labs(y="Significant GO term count",x=" ",fill="")+
  scale_fill_manual(values=c("#fae980","#ff9dc5","#82dde0"))+
  theme(plot.title = element_text(size=8))+
  facet_grid(sample~GO,scales='free')+
  theme(legend.position="top")
pdf(paste0(wd,figures,"/GO_terms_sum_mds.pdf"),height=11,width=10)
GOterms
dev.off()
GO_term_enrichment <- rbind(HSC_MEP_senior_sum_GO,HSC_GMP_senior_sum_GO,MEP_GMP_senior_sum_GO)
GOterms <- GO_term_enrichment %>% 
  mutate(sample=paste0(celltypes,"_",pseudo)) %>% 
  ggplot() + 
  geom_col(aes(x=experiment,y=total,fill="Unique GO terms"))+
  geom_col(aes(x=experiment,y=common,fill="common GO terms"))+
  theme_light()+
  theme(axis.text.x=element_text(angle = -70, hjust = 0))+
  labs(y="Significant GO term count",x=" ",fill="")+
  scale_fill_manual(values=c("#fae980","#ff9dc5","#82dde0"))+
  theme(plot.title = element_text(size=8))+
  facet_grid(sample~GO,scales='free')+
  theme(legend.position="top")
pdf(paste0(wd,figures,"/GO_terms_sum_senior.pdf"),height=11,width=10)
GOterms
dev.off()

GOterms <- GO_term_enrichment %>% 
  mutate(sample=paste0(celltypes,"_",pseudo)) %>% 
  ggplot() + 
  geom_col(aes(x=experiment,y=total,fill="Unique GO terms"))+
  geom_col(aes(x=experiment,y=common,fill="common GO terms"))+
  theme_light()+
  theme(axis.text.x=element_text(angle = -70, hjust = 0))+
  labs(y="Significant GO term count",x=" ",fill="")+
  scale_fill_manual(values=c("#fae980","#ff9dc5","#82dde0"))+
  theme(plot.title = element_text(size=8))+
  facet_grid(sample~GO,scales='free')+
  theme(legend.position="top")
pdf(paste0(wd,figures,"/GO_terms_sum_mds.pdf"),height=11,width=10)
GOterms
dev.off()

GO_term_enrichment <- rbind(HSC_MEP_mds_mean_GO,HSC_GMP_mds_mean_GO,MEP_GMP_mds_mean_GO,
                            HSC_MEP_senior_mean_GO,HSC_GMP_senior_mean_GO,MEP_GMP_senior_mean_GO,
                            HSC_MEP_mds_sum_GO,HSC_GMP_mds_sum_GO,MEP_GMP_mds_sum_GO,
                            HSC_MEP_senior_sum_GO,HSC_GMP_senior_sum_GO,MEP_GMP_senior_sum_GO)

GOterms <- GO_term_enrichment %>% 
  mutate(sample=paste0(celltypes,"_",pseudo)) %>% 
  ggplot() + 
  geom_col(aes(x=experiment,y=total,fill="Unique GO terms"))+
  geom_col(aes(x=experiment,y=common,fill="common GO terms"))+
  theme_light()+
  theme(axis.text.x=element_text(angle = -70, hjust = 0))+
  labs(y="Significant GO term count",x=" ",fill="")+
  scale_fill_manual(values=c("#fae980","#ff9dc5","#82dde0"))+
  theme(plot.title = element_text(size=8))+
  facet_grid(pseudo~GO,scales='free')+
  theme(legend.position="top")
pdf(paste0(wd,figures,"/GO_terms.pdf"),height=8,width=8)
GOterms
dev.off()

########################### GENERATE PLOTS WITH NO REPLICAS ####################
#save plots but without using the replicas
pdf(file=paste0(wd,figures,"/HSC_MEP_DEG_senior_NO_REP.pdf"))
compare_DEG_graph3(results_sum$HSC_MEP_bulkvssc_senior,"(sum)")
compare_DEG_graph3(results_mean$HSC_MEP_bulkvssc_senior,"(mean)")
results_mean$HSC_MEP_bulkvssc_senior$LFC_graph
results_sum$HSC_MEP_bulkvssc_senior$LFC_graph
cowplot::plot_grid(results_mean$HSC_MEP_bulkvssc_senior$sc_PCAplots,
                   results_mean$HSC_MEP_bulkvssc_senior$bulk_PCAplots,
                   ncol = 1)
cowplot::plot_grid(results_sum$HSC_MEP_bulkvssc_senior$sc_PCAplots,
                   results_sum$HSC_MEP_bulkvssc_senior$bulk_PCAplots,
                   ncol= 1)
dev.off()
#--------------------------------------------------------------------------------
pdf(file=paste0(wd,figures,"/HSC_GMP_DEG_senior_NO_REP.pdf"))
compare_DEG_graph3(results_sum$HSC_GMP_bulkvssc_senior,"(sum)")
compare_DEG_graph3(results_mean$HSC_GMP_bulkvssc_senior,"(mean)")
results_mean$HSC_GMP_bulkvssc_senior$LFC_graph
results_sum$HSC_GMP_bulkvssc_senior$LFC_graph
cowplot::plot_grid(results_mean$HSC_GMP_bulkvssc_senior$sc_PCAplots,
                   results_mean$HSC_GMP_bulkvssc_senior$bulk_PCAplots,
                   ncol = 1)
cowplot::plot_grid(results_sum$HSC_GMP_bulkvssc_senior$sc_PCAplots,
                   results_sum$HSC_GMP_bulkvssc_senior$bulk_PCAplots,
                   ncol= 1)
dev.off()
#--------------------------------------------------------------------------------
pdf(file=paste0(wd,figures,"/MEP_GMP_DEG_senior_NO_REP.pdf"))
compare_DEG_graph3(results_sum$MEP_GMP_bulkvssc_senior,"(sum)")
compare_DEG_graph3(results_mean$MEP_GMP_bulkvssc_senior,"(mean)")
results_mean$MEP_GMP_bulkvssc_senior$LFC_graph
results_sum$MEP_GMP_bulkvssc_senior$LFC_graph
cowplot::plot_grid(results_mean$MEP_GMP_bulkvssc_senior$sc_PCAplots,
                   results_mean$MEP_GMP_bulkvssc_senior$bulk_PCAplots,
                   ncol = 1)
cowplot::plot_grid(results_sum$MEP_GMP_bulkvssc_senior$sc_PCAplots,
                   results_sum$MEP_GMP_bulkvssc_senior$bulk_PCAplots,
                   ncol= 1)
dev.off()
#--------------------------------------------------------------------------------
pdf(file=paste0(wd,figures,"/HSC_MEP_DEG_mds_NO_REP.pdf"))
compare_DEG_graph3(results_sum$HSC_MEP_bulkvssc_mds,"(sum)")
compare_DEG_graph3(results_mean$HSC_MEP_bulkvssc_mds,"(mean)")
results_mean$HSC_MEP_bulkvssc_mds$LFC_graph
results_sum$HSC_MEP_bulkvssc_mds$LFC_graph
cowplot::plot_grid(results_mean$HSC_MEP_bulkvssc_mds$sc_PCAplots,
                   results_mean$HSC_MEP_bulkvssc_mds$bulk_PCAplots,
                   ncol = 1)
cowplot::plot_grid(results_sum$HSC_MEP_bulkvssc_mds$sc_PCAplots,
                   results_sum$HSC_MEP_bulkvssc_mds$bulk_PCAplots,
                   ncol= 1)
dev.off()
#--------------------------------------------------------------------------------
pdf(file=paste0(wd,figures,"/HSC_GMP_DEG_mds_NO_REP.pdf"))
compare_DEG_graph3(results_sum$HSC_GMP_bulkvssc_mds,"(sum)")
compare_DEG_graph3(results_mean$HSC_GMP_bulkvssc_mds,"(mean)")
results_mean$HSC_GMP_bulkvssc_mds$LFC_graph
results_sum$HSC_GMP_bulkvssc_mds$LFC_graph
cowplot::plot_grid(results_mean$HSC_GMP_bulkvssc_mds$sc_PCAplots,
                   results_mean$HSC_GMP_bulkvssc_mds$bulk_PCAplots,
                   ncol = 1)
cowplot::plot_grid(results_sum$HSC_GMP_bulkvssc_mds$sc_PCAplots,
                   results_sum$HSC_GMP_bulkvssc_mds$bulk_PCAplots,
                   ncol= 1)
dev.off()
#--------------------------------------------------------------------------------
pdf(file=paste0(wd,figures,"/MEP_GMP_DEG_mds_NO_REP.pdf"))
compare_DEG_graph3(results_sum$MEP_GMP_bulkvssc_mds,"(sum)")
compare_DEG_graph3(results_mean$MEP_GMP_bulkvssc_mds,"(mean)")
results_mean$MEP_GMP_bulkvssc_mds$LFC_graph
results_sum$MEP_GMP_bulkvssc_mds$LFC_graph
cowplot::plot_grid(results_mean$MEP_GMP_bulkvssc_mds$sc_PCAplots,
                   results_mean$MEP_GMP_bulkvssc_mds$bulk_PCAplots,
                   ncol = 1)
cowplot::plot_grid(results_sum$MEP_GMP_bulkvssc_mds$sc_PCAplots,
                   results_sum$MEP_GMP_bulkvssc_mds$bulk_PCAplots,
                   ncol= 1)
dev.off()


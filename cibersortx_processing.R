################################################################################
## AUTHOR: NAROA LEGARRA MARCOS
## DATE: 19/07/22

## Process the rds files to save the psuedobulk matrices in tsv format for cibersortx
################################################################################

wd <- "~/Desktop/TFM" #set the working directory
pseudobulk <- "/DATA/pseudo/" #path to the folder to store the pseudobulks

#################################################################################


save_tsv_pseudo_senior("pseudo_senior1_mean",
                "pseudo_senior2_mean",
                "pseudo_senior3_mean",
                "pseudo_senior_mean")

save_tsv_pseudo_mds("pseudo_mds1_mean",
                       "pseudo_mds3_mean",
                       "pseudo_mds5_mean",
                       "pseudo_mds10_mean",
                       "pseudo_mds_mean")

save_tsv_pseudo_senior("pseudo_senior1_sum",
                       "pseudo_senior2_sum",
                       "pseudo_senior3_sum",
                       "pseudo_senior_sum")

save_tsv_pseudo_mds("pseudo_mds1_sum",
                    "pseudo_mds3_sum",
                    "pseudo_mds5_sum",
                    "pseudo_mds10_sum",
                    "pseudo_mds_sum")

save_tsv_pseudo_mds("pseudo_mds1_rep_mean",
                       "pseudo_mds3_rep_mean",
                       "pseudo_mds5_rep_mean",
                       "pseudo_mds10_rep_mean",
                       "pseudo_mds_rep_mean")

save_tsv_pseudo_mds("pseudo_mds1_rep_sum",
                    "pseudo_mds3_rep_sum",
                    "pseudo_mds5_rep_sum",
                    "pseudo_mds10_rep_sum",
                    "pseudo_mds_rep_sum")

save_tsv_pseudo_senior("pseudo_senior1_rep_mean",
                       "pseudo_senior2_rep_mean",
                       "pseudo_senior3_rep_mean",
                       "pseudo_senior_rep_mean")

save_tsv_pseudo_senior("pseudo_senior1_rep_sum",
                       "pseudo_senior2_rep_sum",
                       "pseudo_senior3_rep_sum",
                       "pseudo_senior_rep_sum")

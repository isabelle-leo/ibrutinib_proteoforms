library(MultiAssayExperiment)
cll_metadata_shiny <- load("~/Downloads/multiomics_MAE.RData")
write.csv(colData(MAE_heidelberg), "CLL_clinical_sample_metadata.csv", row.names = FALSE)
MAE_heidelberg@metadata$treatment_num_legend
#small molecule = 3
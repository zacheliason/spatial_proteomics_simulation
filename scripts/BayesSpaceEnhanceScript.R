library(BayesSpace)

args <- commandArgs(trailingOnly = TRUE)

input_file <- file.path("./data_hallmarks", args[1], "sce.rds")
spatial_output_file <- file.path("./data_hallmarks", args[1], "sce_spatial_enhanced.rds")
feature_output_file <- file.path("./data_hallmarks", args[1], "sce_feature_enhanced.rds")
ST_sce <- readRDS(input_file)  

#check image orientation - swap imagerow and imagecol
if (cor(ST_sce$imagecol, ST_sce$col) < 0.9) {
  tmp <- ST_sce$imagecol
  ST_sce$imagecol <- ST_sce$imagerow
  ST_sce$imagerow <- tmp
}

#enhance single cell experiment object
ST_sce.enhanced <- spatialEnhance(ST_sce, q=length(unique(ST_sce$spatial.cluster)), d=15, platform="Visium", init=ST_sce$spatial.cluster, 
                                  nrep=200000, gamma=3, verbose = TRUE,  jitter_prior=0.3, save.chain=TRUE)


saveRDS(ST_sce.enhanced, spatial_output_file)


#load single cell experiment object at sub-spot level (without gene expression)
sce_enhanced <- readRDS(spatial_output_file)

#Impute genes at sub-spot resolution
sce_enhanced <- enhanceFeatures(sce_enhanced, readRDS(input_file),
                                model="xgboost",
                                nrounds=100,
                                assay.type = "SCT")
#save file
saveRDS(sce_enhanced, feature_output_file)

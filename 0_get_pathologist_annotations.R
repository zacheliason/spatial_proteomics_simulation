library(tidyverse)
library(Seurat)

root <- '/Users/zacheliason/Downloads/20304456'
full_sample_names <- list.files(root, full.names = TRUE, recursive = FALSE, pattern = "sample_\\d+.Robj")
sample_names <- list.files(root, full.names = FALSE, recursive = FALSE, pattern = "sample_\\d+.Robj") %>%
  str_remove("\\.Robj")

# Initialize an empty list to store the tibbles for each sample
tibble_list <- list()

# Loop through each sample
for (i in 1:length(sample_names)) {
  sample_name = sample_names[i]
  full_sample_name = full_sample_names[i]
  
  # Load the sample object using load
  load(full_sample_name)
  sample_obj <- get(sample_name)
  
  # Extract relevant information
  spots <- colnames(sample_obj@assays$SCT)
  pathologist_annotation <- sample_obj@meta.data$pathologist_anno.x
  cluster_annotation <- sample_obj@meta.data$cluster_annotations
  
  # Create a tibble
  tibble <- tibble(barcode = spots, pathologist_annotation = pathologist_annotation, cluster_annotation = cluster_annotation)
  
  write_csv(tibble, file.path(root, paste0(sample_name,"_pathologist_annotations", ".csv")))
}


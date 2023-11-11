library(tidyverse)
library(Seurat)


load(file = "./sample_1.Robj")
spots <- colnames(sample_1@assays$SCT)
pathologist_annotation <- sample_1@meta.data$pathologist_anno.x

df <- data.frame(Annotation = pathologist_annotation) 
rownames(df)<-spots

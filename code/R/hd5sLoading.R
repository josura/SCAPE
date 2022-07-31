remotes::install_github("mojaveazure/seurat-disk")

library(SeuratDisk)

Convert("/home/josura/Projects/thingsToRemove/trisicell/trisicell/datasets/real/acute_lymphocytic_leukemia1.h5ad",dest = "h5seurat")

b2905 <- LoadH5Seurat("/home/josura/Projects/thingsToRemove/trisicell/trisicell/datasets/real/acute_lymphocytic_leukemia1.h5seurat")


b2905Graph <- DimPlot(b2905)
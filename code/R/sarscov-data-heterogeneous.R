library(Matrix)
library(Seurat)

mtxFileLocation <- '/run/media/josura/Codigos/dati/sars-cov/Large-scale single-cell analysis reveals critical immune characteristics of COVID-19 patients/GSE158055_covid19_counts.mtx'
cellBarcodesFile <- '/run/media/josura/Codigos/dati/sars-cov/Large-scale single-cell analysis reveals critical immune characteristics of COVID-19 patients/GSE158055_covid19_barcodes.tsv'
featurefile <- '/run/media/josura/Codigos/dati/sars-cov/Large-scale single-cell analysis reveals critical immune characteristics of COVID-19 patients/GSE158055_covid19_features.tsv'

annotationFile <- '/run/media/josura/Codigos/dati/sars-cov/Large-scale single-cell analysis reveals critical immune characteristics of COVID-19 patients/GSE158055_cell_annotation.csv'


matrix <- ReadMtx(
  mtxFileLocation,
  cellBarcodesFile,
  featurefile,
  cell.column = 1,
  feature.column = 1,
  cell.sep = "\t",
  feature.sep = "\t",
  skip.cell = 0,
  skip.feature = 0,
  mtx.transpose = FALSE,
  unique.features = TRUE,
  strip.suffix = FALSE
)


firstPart <- Read10X("/run/media/josura/Codigos/dati/sars-cov/Large-scale single-cell analysis reveals critical immune characteristics of COVID-19 patients/part1dir/",gene.column = 1 )
Read10X("/run/media/josura/Codigos/dati/sars-cov/Large-scale single-cell analysis reveals critical immune characteristics of COVID-19 patients/part2dir/",gene.column = 1 )
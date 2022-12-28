library(Matrix)
library(Seurat)
library(dplyr)

# mtxFileLocation <- '/run/media/josura/Codigos/dati/sars-cov/Large-scale single-cell analysis reveals critical immune characteristics of COVID-19 patients/GSE158055_covid19_counts.mtx'
# cellBarcodesFile <- '/run/media/josura/Codigos/dati/sars-cov/Large-scale single-cell analysis reveals critical immune characteristics of COVID-19 patients/GSE158055_covid19_barcodes.tsv'
# featurefile <- '/run/media/josura/Codigos/dati/sars-cov/Large-scale single-cell analysis reveals critical immune characteristics of COVID-19 patients/GSE158055_covid19_features.tsv'
# 
# annotationFile <- '/run/media/josura/Codigos/dati/sars-cov/Large-scale single-cell analysis reveals critical immune characteristics of COVID-19 patients/GSE158055_cell_annotation.csv'
# 
# 
# matrix <- ReadMtx(
#   mtxFileLocation,
#   cellBarcodesFile,
#   featurefile,
#   cell.column = 1,
#   feature.column = 1,
#   cell.sep = "\t",
#   feature.sep = "\t",
#   skip.cell = 0,
#   skip.feature = 0,
#   mtx.transpose = FALSE,
#   unique.features = TRUE,
#   strip.suffix = FALSE
# )
# 
# 
# firstPart <- Read10X("/run/media/josura/Codigos/dati/sars-cov/Large-scale single-cell analysis reveals critical immune characteristics of COVID-19 patients/part1dir/",gene.column = 1 )
# Read10X("/run/media/josura/Codigos/dati/sars-cov/Large-scale single-cell analysis reveals critical immune characteristics of COVID-19 patients/part2dir/",gene.column = 1 )

postacute.dir <- "/run/media/josura/2217-021A/Data/post-acute COVID_longCOVID/files/"

postacute.metadata <- read.csv(paste(postacute.dir,"E-MTAB-10129.sdrf.txt",sep=""),sep = "\t")
#Derived.Array.Data.File contains all the genes
#Derived.Array.Data.File.1 contains the filtered genes that have enough counts
#Derived.Array.Data.File.2 contains the metadata for the cells isolated for the single patient
#Factor.Value.blood.draw contains the time of the sample , that is "healthy control","blood draw near diagnosis","blood draw a few days after diagnosis","blood draw at convalescence"
#all the values of the above matrices are real values (probably normalized)

test <- read.csv(paste(postacute.dir,"heathdc118_pbmc_pro_library_1053BW.txt",sep=""),sep = "\t",header = TRUE)
rownames(test) <- test[,1]
test[,1] <- NULL


# PREPARING THE MATRICES AND DIVIDE THE DATA BETWEEN TIME OF DIAGNOSIS, ALSO TAKING INTO ACCOUNT ONLY THE DERIVED ARRAY FILTERED
merge.scrnaseq.matrices <- function(dir,file.list){
  countmatrix.list <- list()
  for (file in file.list) {
    postacute.file <- read.csv(paste(dir,file,sep=""),sep = "\t",header = TRUE)
    rownames(postacute.file) <- postacute.file[,1]
    postacute.file[,1] <- NULL
    #phensim.embed$gene.name <- pathway.name
    countmatrix.list[[length(countmatrix.list)+1]] <- postacute.file
    print(paste("file read done for ",file,sep = ": "))
  }
  merged.countmatrices <- countmatrix.list[[1]]
  countmatrix.list <- countmatrix.list[-1]
  iteratori <- 1
  for (countmatrix in countmatrix.list) {
    merged.countmatrices <- bind_rows(merged.countmatrices,countmatrix)
    print(paste("merge done for pathway",file.list[iteratori],sep = ": "))
    iteratori <- iteratori+1
  }
  #merged.embeddings[is.na(merged.embeddings)] <- 0
  merged.countmatrices
}

merge.scrnaseq.matrices.filtering <- function(dir,file.list){
  countmatrix.list <- list()
  for (file in file.list) {
    postacute.file <- read.csv(paste(dir,file,sep=""),sep = "\t",header = TRUE)
    rownames(postacute.file) <- postacute.file[,1]
    postacute.file[,1] <- NULL
    #phensim.embed$gene.name <- pathway.name
    countmatrix.list[[length(countmatrix.list)+1]] <- postacute.file
    print(paste("file read done for ",file,sep = ": "))
  }
  merged.countmatrices <- countmatrix.list[[1]]
  countmatrix.list <- countmatrix.list[-1]
  iteratori <- 1
  for (countmatrix in countmatrix.list) {
    merged.countmatrices <- bind_rows(merged.countmatrices,countmatrix)
    print(paste("merge done for pathway",file.list[iteratori],sep = ": "))
    iteratori <- iteratori+1
  }
  #merged.embeddings[is.na(merged.embeddings)] <- 0
  merged.countmatrices
}


#92 samples without condition
#139 before diagnosis
#132 after diagnosis
#88 convalescence
#15 controls
#1 mix_donor1
beforediagnosis.filelist <- postacute.metadata %>%
  filter(Factor.Value.blood.draw. == "blood draw near diagnosis") %>%
  dplyr::select(Derived.Array.Data.File.1)
beforediagnosis.filelist <- beforediagnosis.filelist[,1]
afterdiagnosis.filelist <- postacute.metadata %>%
  filter(Factor.Value.blood.draw. == "blood draw a few days after diagnosis") %>%
  dplyr::select(Derived.Array.Data.File.1)
afterdiagnosis.filelist <- afterdiagnosis.filelist[,1]
convalescence.filelist <- postacute.metadata %>%
  filter(Factor.Value.blood.draw. == "blood draw at convalescence") %>%
  dplyr::select(Derived.Array.Data.File.1)
convalescence.filelist <- convalescence.filelist[,1]
controls.filelist <- postacute.metadata %>%
  filter(Factor.Value.blood.draw. == "healthy control") %>%   #also not considering the Mix_donor as control as well since the condition is normal
  dplyr::select(Derived.Array.Data.File.1)
controls.filelist <- controls.filelist[,1]

controls.matrix.filtered <- merge.scrnaseq.matrices(postacute.dir,controls.filelist)
#the genes are the same in all the filtered files, I think I will change the approach to take all the genes with at least 95% of the cells with a value or something along those lines

##full data
controls.filelist <- postacute.metadata %>%
  filter(Factor.Value.blood.draw. == "healthy control") %>%   #also not considering the Mix_donor as control as well since the condition is normal
  dplyr::select(Derived.Array.Data.File)
controls.filelist <- controls.filelist[,1]

controls.matrix.full <- merge.scrnaseq.matrices(postacute.dir,controls.filelist)
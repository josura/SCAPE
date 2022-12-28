
install.packages("tidyverse")
BiocManager::install("DESeq2")
BiocManager::install("tximport")

#after you've installed packages:
library("tidyverse")
library("tximport")
library("DESeq2")


fastqDirectories <- list.dirs("/home/josura/Projects/tesi/data/tesi/patient1ALL/SingleCell",recursive = FALSE, full=TRUE)
#fastqDirectories <- fastqDirectories[2:length(fastqDirectories)]

samples <- read.csv("/home/josura/Projects/tesi/data/tesi/patient1ALL/SingleCell/sample.csv")

# Transcript - Gene ID Data Frame
# read in the file from url, store in new object
tx2gene_map <- read_tsv("/home/josura/Projects/tesi/data/tesi/referenceGenomes/human/referenceTranscriptome//referenceTranscriptomeAnnotation.tsv/result.tsv")
# look at the first 6 lines
head(tx2gene_map)
# What is this file and how do I make one? 
# Need to associate Transcript IDs with GeneIDs - Salmon only provides Transcript IDs
# Making a tx2gene file is often a little different for each organism. 
# If your organism has a transcriptome (or *rna_from_genomic.fna.gz file) on RefSeq, you can often get the information from the *feature_table.txt.gz. 
# You might also be able to parse a gtf or gff file to produce the information you need.
# This information is sometimes also found in the fasta headers of the transcriptome file itself.
# https://bioconductor.org/packages/3.7/bioc/vignettes/tximport/inst/doc/tximport.html (under Import transcript-level estimates)
# https://support.bioconductor.org/p/106319/
Qfiles <- samples$quant_file


txi <- tximport(files = Qfiles, type = "salmon", tx2gene = tx2gene_map,ignoreTxVersion=TRUE)

names(txi)
head(txi$counts)
summary(txi)
#Change the column names to the names of the samples (experiments)
colnames(txi$counts) <- samples$sample

head(txi$counts)

#align to transcript 
for (directory in fastqDirectories){
  print(directory)
  #reads_1 <- readFastq(fastqFiles[1], withIds=TRUE)
  #reads_2 <- readFastq(fastqFiles[2], withIds=TRUE)
  quant <- file.path(directory,"salmon/quantification")
  print(quant)
  
}


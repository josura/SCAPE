BiocManager::install("Biostrings")
BiocManager::install("Rbowtie2")

BiocManager::install("GenomicFeatures")

library(Biostrings)
library(Rbowtie2)



fastqDirectories <- list.dirs("/home/josura/Projects/tesi/data/tesi/patient1ALL/SingleCell",recursive = FALSE, full=TRUE)
#fastqDirectories <- fastqDirectories[2:length(fastqDirectories)]

fastqDirectories <- list.dirs("/home/josura/Projects/tesi/data/tesi/patient2ALL/singleCell",recursive = FALSE, full=TRUE)
fastqDirectories <- list.dirs("/home/josura/Projects/tesi/data/tesi/patient3ALL/singleCell",recursive = FALSE, full=TRUE)
fastqDirectories <- list.dirs("/home/josura/Projects/tesi/data/tesi/patient4ALL/singleCell",recursive = FALSE, full=TRUE)

#removing library adapters
for (directory in fastqDirectories){
  print(directory)
  fastqFiles <- list.files(directory, full=TRUE, pattern = "\\.fastq$")
  print(fastqFiles[1])
  print(fastqFiles[2])
  #reads_1 <- readFastq(fastqFiles[1], withIds=TRUE)
  #reads_2 <- readFastq(fastqFiles[2], withIds=TRUE)
  reads_1 <- fastqFiles[1]
  reads_2 <- fastqFiles[2]

  #identify library adapters
  td <- tempdir()
  (adapters <-
      identify_adapters(file1=reads_1,file2=reads_2,
                        basename=file.path(directory,"trimmed/reads"),
                        "--threads 3",overwrite=TRUE))
  #remove library adapters
  (cmdout<-remove_adapters(file1=reads_1,file2=reads_2,adapter1 = adapters[1],
                            adapter2 = adapters[2],
                            output1=file.path(directory,"trimmed/reads_1.trimmed.fastq"),
                            output2=file.path(directory,"trimmed/reads_2.trimmed.fastq"),
                            basename=file.path(directory,"trimmed/reads.base"),overwrite=TRUE,"--threads 3"))
}

#prepare the index references for bowtie2 (i do not needed since I downloaded the index files directly)
td <- tempdir()
(cmdout<-bowtie2_build(references="/home/josura/Projects/tesi/data/tesi/referenceGenomes/human/GRCh38_fasta/chromosomes/noalt", 
                       bt2Index=file.path(td, "human_genome"), "--threads 4 --quiet",
                       overwrite=TRUE))

#align data and create bam files
for (directory in fastqDirectories){
  print(directory)
  #reads_1 <- readFastq(fastqFiles[1], withIds=TRUE)
  #reads_2 <- readFastq(fastqFiles[2], withIds=TRUE)
  reads_1 <- file.path(directory,"trimmed/reads_1.trimmed.fastq")
  reads_2 <- file.path(directory,"trimmed/reads_2.trimmed.fastq")
  print(reads_1)
  print(reads_2)
  
  if(file.exists("/home/josura/Projects/tesi/data/tesi/referenceGenomes/human/GRCh38_noalt_as")){
    (cmdout<-bowtie2_samtools(bt2Index = file.path("/home/josura/Projects/tesi/data/tesi/referenceGenomes/human/GRCh38_noalt_as/GRCh38_noalt_as"),
                              output = file.path(directory,"aligned/result"),
                              outputType = "sam",
                              seq1=reads_1,
                              seq2=reads_2,
                              overwrite=TRUE,
                              bamFile = NULL,
                              "--threads 6"))
    head(readLines(file.path(directory,"aligned/result.sam")))
  }
  
}



##  QUANTIFYING EXPRESSION ESTIMATES ##

# load gene bank
library("GenomicFeatures")

gtffile <- file.path("/home/josura/Projects/tesi/data/tesi/referenceGenomes/human/ensemblGene/gtf/Homo_sapiens.GRCh38.106.gtf")
(txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character()))
(ebg <- exonsBy(txdb, by="gene"))

# read counting step
library("GenomicAlignments")
library("BiocParallel")

bamfiles <- paste(fastqDirectories, "aligned/result.bam",sep = "/")

getNamesExperiments <- function(directories){
  firstPos <- unlist(gregexpr("SRR",directories[1]))
  return(substr(directories,firstPos,firstPos+10))
}

names(bamfiles) <- getNamesExperiments(fastqDirectories)


summarized.experiment <- summarizeOverlaps(features=ebg, reads=bamfiles,
                                           mode="Union",
                                           singleEnd=FALSE,
                                           ignore.strand=TRUE,
                                           fragments=TRUE )

countsPatient1 <- data.frame( summarized.experiment@assays@data@listData$counts )
rownames(countsPatient1) <-  ebg@partitioning@NAMES 
colnames(countsPatient1) <-  names(bamfiles)

#### PROBLEMS since the final results are 0 counts for everything
res <- colSums(countsPatient1==0)/nrow(countsPatient1)*100

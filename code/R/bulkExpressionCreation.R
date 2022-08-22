BiocManager::install("Biostrings")
BiocManager::install("Rbowtie2")

BiocManager::install("GenomicFeatures")

library(Biostrings)
library(Rbowtie2)



fastqDirectories <- list.dirs("/home/josura/Projects/tesi/data/tesi/patient1ALL/bulk",recursive = FALSE, full=TRUE)
#fastqDirectories <- fastqDirectories[2:length(fastqDirectories)]

fastqDirectories <- list.dirs("/home/josura/Projects/tesi/data/tesi/patient2ALL/bulk",recursive = FALSE, full=TRUE)
fastqDirectories <- list.dirs("/home/josura/Projects/tesi/data/tesi/patient3ALL/bulk",recursive = FALSE, full=TRUE)
fastqDirectories <- list.dirs("/home/josura/Projects/tesi/data/tesi/patient4ALL/bulk",recursive = FALSE, full=TRUE)

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

ctrlfastqDirectories <- list.dirs("/home/josura/Projects/tesi/data/tesi/controls/bulk",recursive = FALSE, full=TRUE)

#removing library adapters
for (directory in ctrlfastqDirectories){
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

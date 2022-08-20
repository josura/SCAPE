# download file: prefetch will download and save SRA file related to SRR accession in 
# the current directory under newly created SRA accession directory
$ prefetch  SRR5790106  # for a single file
$ prefetch  SRR5790106 SRR5790104  # multiple files

# convert to FASTQ: fastq-dump will convert SRR5790106.sra to SRR5790106.fastq
$ fastq-dump  SRR5790106  # single file
$ fastq-dump  SRR5790106  SRR5790104 # multiple files

# now you can also replace fastq-dump with fasterq-dump which is much faster 
# and efficient for large datasets
# by default it will use 6 threads (-e option)
$ fasterq-dump  SRR5790106  # single file
$ fasterq-dump  SRR5790106  SRR5790104 # multiple files

# for paired-end data use --split-files (fastq-dump) and -S or --split-files (fasterq-dump) option
# use --split-3 option, if the paired-end reads are not properly arranged (e.g. some reads has lack of mate pair)
$ fastq-dump --split-files SRR8296149
$ fasterq-dump -S SRR8296149

# download biological and technical reads (cell and sample barcodes) in case of single 
# cell RNA-seq (10x chromium)
$ fasterq-dump -S --include-technical SRR12564282

# download alignment files (SAM)
# make sure corresponding accession has alignment file at SRA database
$ sam-dump --output-file SRR1236468.sam SRR1236468

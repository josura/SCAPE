
library("tidyverse")
library("tximport")
library("DESeq2")

### PREPARATION OF DATA ###
### EXPRESSION SHOULD ONLY BE SEEN FOR GENES THAT ARE EXPRESSED AT LEAST IN ONE OF THE PATIENTS, THE UNION OF THESE GENES FOR THE PATIENTS WILL BE THE FEATURES TO CONSIDER  FIRST

normalizeMax <- function(x){(x-min(x))/(max(x)-min(x))}
normalizeTest <- function(x){(x)/(100000.0)}

## PATIENT 1

fastqDirectoriesPatient <- list.dirs("/home/josura/Projects/tesi/data/tesi/patient4ALL/bulk",recursive = FALSE, full=TRUE)
bulkTumor <- tail(fastqDirectoriesPatient,n=1)
#SRR1517768 is patient 4 tumor bone marrow RNA-seq

#getting all samples CTRL and ALL for patient 4 in a table
samples <- read.csv("/home/josura/Projects/tesi/data/tesi/patient4ALL/bulk/sample.csv")

# Transcript - Gene ID Data Frame
tx2gene_map <- read_tsv("/home/josura/Projects/tesi/data/tesi/referenceGenomes/human/referenceTranscriptome/referenceTranscriptomeAnnotation.tsv/result.tsv")
Qfiles <- samples$quant_file

txi <- tximport(files = Qfiles, type = "salmon", tx2gene = tx2gene_map,ignoreTxVersion=TRUE)


colnames(txi$counts) <- samples$sample

transcript.counts <-  t(txi$counts)

#count zeros as percentages in genes expressions
res <- colSums(transcript.counts==0)/nrow(transcript.counts)*100
#list genes with non-zero percentage of expression
res.nonzero <- res[which(res<100)]


### SCALE AND SELECT GENES ###

## PATIENT 4
#select only the genes that are expressed

transcript.counts.filtered <- transcript.counts[,names(res.nonzero)]

#normality tests on the columns
library(data.table)
resNormality <- data.table(transcript.counts.filtered) %>% mutate_all(function(x){shapiro.test(x)$p.value})

#normalization of counts
#transcript.counts.norm <-  data.table(data.table(transcript.counts.filtered) %>% mutate_all(scale))
transcript.counts.norm <-  data.table(data.table(transcript.counts.filtered) %>% mutate_all(normalizeMax))
transcript.counts.norm <-  data.table(data.table(transcript.counts.filtered) %>% mutate_all(normalizeTest))

write.csv(transcript.counts,"/home/josura/Projects/tesi/data/patient4ALL/allmerged-bulk.csv")
write.csv(transcript.counts.norm,"/home/josura/Projects/tesi/data/patient4ALL/allmergedNormalized-bulk.csv")
#results suggests that the lesser transcript gives more reads than the extended index for some reason




### differential analysis ###


## PATIENT 4
# Prepare this data for the differential expression algorithm
# Look through the documentation for this function to find out all the options you have here
# documentation in R is found by a command like "?DESeqDataSetFromTximport")

txi.augmented <- copy(txi)
txi.augmented$counts <- txi.augmented$counts +1

#the groups are the condition of the sample, in this case 1 patient from ALL and 3 CTRL patients


dds <- DESeqDataSetFromTximport(txi = txi.augmented,           #if the group in samples has only one value, there will be problems
                                design = ~ condition, #all variables in design should be the different groups in samples
                                colData = samples)


# Run DESeq on this dataset
dds <- DESeq(dds)
# Store DESeq results in a new object
res <- results(dds)
# Check out results
head(res)
# Store a subset of results in a new object, in this case, the ones with an adjusted p-value of < 0.05
res_sig <- subset(res, padj<.05)
# Out of the subset we created above, subset the results that changed between conditions
res_lfc <- subset(res_sig, abs(log2FoldChange) > 1) 
head(res_lfc)


# add a column of expression labelling
res$diffexpressed <- "NO"
# if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
res$diffexpressed[res$log2FoldChange > 0.5 & res$padj < 0.05] <- "UP"
# if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
res$diffexpressed[res$log2FoldChange < -0.5 & res$padj < 0.05] <- "DOWN"


### PHENSIM input preparation
library(biomaRt)

load("/home/josura/Projects/tesi/data/ensemblToGene4Bulk.RData")  #load mart object since the ensembl site goes down sometimes


mart <- useDataset("hsapiens_gene_ensembl",useMart("ensembl")) 


Glist4 <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene_id","description"),values = colnames(transcript.counts.filtered),mart = mart)
allGlist4 <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene_id","description"),values = colnames(transcript.counts),mart = mart)
save(mart,Glist4,allGlist4,file = "/home/josura/Projects/tesi/data/ensemblToGene4Bulk.RData")

res$ensembl_gene_id <- rownames(res)

joinedRes <- inner_join(as.data.frame(res), allGlist4)

allGenesTable <- joinedRes%>% 
  dplyr::select(hgnc_symbol,entrezgene_id,diffexpressed) %>%
  #distinct(entrezgene_id) %>% 
  na.omit()

nonExpressedGenesPatient4 <- allGenesTable[allGenesTable$diffexpressed=="NO",] %>% 
  distinct(entrezgene_id)


underExpressedGenesPatient4 <- allGenesTable[allGenesTable$diffexpressed=="UP",] %>% 
  distinct(entrezgene_id)
  
overExpressedGenesPatient4 <- allGenesTable[allGenesTable$diffexpressed=="DOWN",] %>% 
  distinct(entrezgene_id)

trueNonExpressed <- allGlist4 %>%
  filter(!(entrezgene_id %in% Glist4$entrezgene_id)) %>%
  dplyr::select(entrezgene_id)


write.table(joinedRes,
            file='/home/josura/Projects/tesi/data/patient4ALL/bulk-allGenesDE.tsv', quote=FALSE, sep='\t', col.names = TRUE,row.names = FALSE)

#only txt with entrezgene_id for phensim
write.table(nonExpressedGenesPatient4,
            file='/home/josura/Projects/tesi/data/patient4ALL/bulk-nonExpressed.txt', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)

write.table(overExpressedGenesPatient4,
            file='/home/josura/Projects/tesi/data/patient4ALL/bulk-overexpressed.txt', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)

write.table(underExpressedGenesPatient4,
            file='/home/josura/Projects/tesi/data/patient4ALL/bulk-downexpressed.txt', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)

write.table(trueNonExpressed,
            file='/home/josura/Projects/tesi/data/patient4ALL/bulk-trueNonExpressed.tsv', quote=FALSE, sep='\t', col.names = TRUE,row.names = FALSE)


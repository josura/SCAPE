install.packages("umap")
library("umap")

library("tidyverse")
library("tximport")
library("DESeq2")

### PREPARATION OF DATA ###
### EXPRESSION SHOULD ONLY BE SEEN FOR GENES THAT ARE EXPRESSED AT LEAST IN ONE OF THE PATIENTS, THE UNION OF THESE GENES FOR THE PATIENTS WILL BE THE FEATURES TO CONSIDER  FIRST

normalizeMax <- function(x){(x-min(x))/(max(x)-min(x))}
normalizeTest <- function(x){(x)/(100000.0)}

## PATIENT 1

fastqDirectories <- list.dirs("/home/josura/Projects/tesi/data/tesi/patient1ALL/SingleCell",recursive = FALSE, full=TRUE)


samples <- read.csv("/home/josura/Projects/tesi/data/tesi/patient1ALL/SingleCell/sample.csv")

#different counts with extended transcript from a different salmon  index, used to compare the two salmon outputs
samplesFull <- read.csv("/home/josura/Projects/tesi/data/tesi/patient1ALL/SingleCell/sampleFull.csv")

# Transcript - Gene ID Data Frame
tx2gene_map <- read_tsv("/home/josura/Projects/tesi/data/tesi/referenceGenomes/human/referenceTranscriptome/referenceTranscriptomeAnnotation.tsv/result.tsv")
Qfiles <- samples$quant_file
QfilesFull <- samplesFull$quant_file

txi <- tximport(files = Qfiles, type = "salmon", tx2gene = tx2gene_map,ignoreTxVersion=TRUE)
txiFull <- tximport(files = QfilesFull, type = "salmon", tx2gene = tx2gene_map,ignoreTxVersion=TRUE)

colnames(txi$counts) <- samples$sample

transcript.counts <-  t(txi$counts)
transcriptFull.counts <-  t(txiFull$counts)

#count zeros as percentages in genes expressions
res <- colSums(transcript.counts==0)/nrow(transcript.counts)*100
resFull <- colSums(transcriptFull.counts==0)/nrow(transcriptFull.counts)*100
#list genes with non-zero percentage of expression
res.nonzero <- res[which(res<100)]
resFull.nonzero <- resFull[which(resFull<100)]
#results suggests that the lesser transcript gives more reads than the extended index for some reason


## PATIENT 2


fastqDirectories2 <- list.dirs("/home/josura/Projects/tesi/data/tesi/patient2ALL/singleCell",recursive = FALSE, full=TRUE)
samples2 <- read.csv("/home/josura/Projects/tesi/data/tesi/patient2ALL/singleCell/sample.csv")

Qfiles2 <- samples2$quant_file

txi2 <- tximport(files = Qfiles2, type = "salmon", tx2gene = tx2gene_map,ignoreTxVersion=TRUE)

colnames(txi2$counts) <- samples2$sample

transcript2.counts <-  t(txi2$counts)

res2 <- colSums(transcript2.counts==0)/nrow(transcript2.counts)*100
res2.nonzero <- res2[which(res2<100)]



## PATIENT 3

fastqDirectories3 <- list.dirs("/home/josura/Projects/tesi/data/tesi/patient3ALL/singleCell",recursive = FALSE, full=TRUE)

samples3 <- read.csv("/home/josura/Projects/tesi/data/tesi/patient3ALL/singleCell/sample.csv")

Qfiles3 <- samples3$quant_file

txi3 <- tximport(files = Qfiles3, type = "salmon", tx2gene = tx2gene_map,ignoreTxVersion=TRUE)

colnames(txi3$counts) <- samples3$sample

transcript3.counts <-  t(txi3$counts)

res3 <- colSums(transcript3.counts==0)/nrow(transcript3.counts)*100
res3.nonzero <- res3[which(res3<100)]



## PATIENT 4

fastqDirectories4 <- list.dirs("/home/josura/Projects/tesi/data/tesi/patient4ALL/singleCell",recursive = FALSE, full=TRUE)

samples4 <- read.csv("/home/josura/Projects/tesi/data/tesi/patient4ALL/singleCell/sample.csv")

Qfiles4 <- samples4$quant_file

txi4 <- tximport(files = Qfiles4, type = "salmon", tx2gene = tx2gene_map,ignoreTxVersion=TRUE)

colnames(txi4$counts) <- samples4$sample

transcript4.counts <-  t(txi4$counts)

res4 <- colSums(transcript4.counts==0)/nrow(transcript4.counts)*100
res4.nonzero <- res4[which(res4<100)]


full.nonzero <- union(union(union(names(res.nonzero),names(res2.nonzero)),names(res3.nonzero)),names(res4.nonzero))



### SCALE AND SELECT GENES ###

## PATIENT 1
#select only the genes that are expressed

transcript.counts.filtered.singular <- transcript.counts[,names(res.nonzero)]
transcript.counts.filtered <- transcript.counts[,full.nonzero]


#normality tests on the columns
resNormality <- data.table(transcript.counts.filtered) %>% mutate_all(function(x){shapiro.test(x)$p.value})

#normalization of counts
#transcript.counts.norm <-  data.table(data.table(transcript.counts.filtered) %>% mutate_all(scale))
transcript.counts.norm <-  data.table(data.table(transcript.counts.filtered) %>% mutate_all(normalizeTest))
transcript.counts.norm.singular<-  data.table(data.table(transcript.counts.filtered.singular) %>% mutate_all(normalizeTest))

write.csv(transcript.counts,"/home/josura/Projects/tesi/data/tesi/patient1ALL/allmerged.csv")
write.csv(transcript.counts.filtered,"/home/josura/Projects/tesi/data/tesi/patient1ALL/allmerged_filtered.csv")

## PATIENT 2

#select only the genes that are expressed
transcript2.counts.filtered.singular <- transcript2.counts[,names(res2.nonzero)]
transcript2.counts.filtered <- transcript2.counts[,full.nonzero]

#normality tests on the columns
res2Normality <- data.table(transcript2.counts.filtered) %>% mutate_all(function(x){shapiro.test(x)$p.value})

#transcript2.counts.norm <-  data.table(data.table(transcript2.counts.filtered) %>% mutate_all(scale))
transcript2.counts.norm <-  data.table(data.table(transcript2.counts.filtered) %>% mutate_all(normalizeTest))
transcript2.counts.norm.singular <-  data.table(data.table(transcript2.counts.filtered.singular) %>% mutate_all(normalizeTest))


write.csv(transcript2.counts,"/home/josura/Projects/tesi/data/tesi/patient2ALL/allmerged.csv")
write.csv(transcript2.counts.filtered,"/home/josura/Projects/tesi/data/tesi/patient2ALL/allmerged_filtered.csv")

## PATIENT 3

#select only the genes that are expressed
transcript3.counts.filtered.singular <- transcript3.counts[,names(res3.nonzero)]
transcript3.counts.filtered <- transcript3.counts[,full.nonzero]

#normality tests on the columns
res3Normality <- data.table(transcript3.counts.filtered) %>% mutate_all(function(x){shapiro.test(x)$p.value})

#transcript3.counts.norm <-  data.table(data.table(transcript3.counts.filtered) %>% mutate_all(scale))
transcript3.counts.norm <-  data.table(data.table(transcript3.counts.filtered) %>% mutate_all(normalizeTest))
transcript3.counts.norm.singular <-  data.table(data.table(transcript3.counts.filtered.singular) %>% mutate_all(normalizeTest))


write.csv(transcript3.counts,"/home/josura/Projects/tesi/data/tesi/patient3ALL/allmerged.csv")
write.csv(transcript3.counts.filtered,"/home/josura/Projects/tesi/data/tesi/patient3ALL/allmerged_filtered.csv")

## PATIENT 4


#select only the genes that are expressed
transcript4.counts.filtered.singular <- transcript4.counts[,names(res4.nonzero)]
transcript4.counts.filtered <- transcript4.counts[,full.nonzero]

#normality tests on the columns
res4Normality <- data.table(transcript4.counts.filtered) %>% mutate_all(function(x){shapiro.test(x)$p.value})

#transcript4.counts.norm <-  data.table(data.table(transcript4.counts.filtered) %>% mutate_all(scale))
transcript4.counts.norm <-  data.table(data.table(transcript4.counts.filtered) %>% mutate_all(normalizeTest))
transcript4.counts.norm.singular <-  data.table(data.table(transcript4.counts.filtered.singular) %>% mutate_all(normalizeTest))


write.csv(transcript4.counts,"/home/josura/Projects/tesi/data/tesi/patient4ALL/allmerged.csv")
write.csv(transcript4.counts.filtered,"/home/josura/Projects/tesi/data/tesi/patient4ALL/allmerged_filtered.csv")




### UMAP AND CLUSTER ###

## PATIENT 1 ##


#custom umap configuration
custom.config = umap.defaults
custom.config$n_components = 3  #to change into the dimension we want

umapped.transcript <- umap(transcript.counts.norm,config = custom.config)
umapped.transcript.singular <- umap(transcript.counts.norm.singular,config = custom.config)

#clustering
distances <- dist(transcript.counts.norm, method = "euclidean")
distances.umap <- dist(umapped.transcript$layout, method = "euclidean")
distances.umap.singular <- dist(umapped.transcript.singular$layout, method = "euclidean")


res.clust <- hclust(distances, method = "centroid")
res.clust.umap <- hclust(distances.umap, method = "centroid")
res.clust.umap.singular <- hclust(distances.umap.singular, method = "centroid")

cuttree <- cutree(res.clust, k = 4)
cuttree.umap <- cutree(res.clust.umap, k = 4)
cuttree.umap.singular <- cutree(res.clust.umap.singular, k = 4)

#plot data
library(plotly)
plot_ly(x=umapped.transcript$layout[,1],y=umapped.transcript$layout[,2],z=umapped.transcript$layout[,3],type="scatter3d",color = cuttree.umap)



### PATIENT 2 ###

umapped.transcript2 <- umap(transcript2.counts.norm,config = custom.config)
umapped.transcript2.singular <- umap(transcript2.counts.norm.singular,config = custom.config)

#clustering
distances2 <- dist(transcript2.counts.norm, method = "euclidean")
distances.umap2 <- dist(umapped.transcript2$layout, method = "euclidean")
distances.umap2.singular <- dist(umapped.transcript2.singular$layout, method = "euclidean")


res.clust2 <- hclust(distances2, method = "centroid")
res.clust.umap2 <- hclust(distances.umap2, method = "centroid")
res.clust.umap2.singular <- hclust(distances.umap2.singular, method = "centroid")

cuttree2 <- cutree(res.clust2, k = 4)
cuttree.umap2 <- cutree(res.clust.umap2, k = 4)
cuttree.umap2.singular <- cutree(res.clust.umap2.singular, k = 4)

#plot data
plot_ly(x=umapped.transcript2$layout[,1],y=umapped.transcript2$layout[,2],z=umapped.transcript2$layout[,3],type="scatter3d",color = cuttree.umap2)


### PATIENT 3 ###


umapped.transcript3 <- umap(transcript3.counts.norm,config = custom.config)
umapped.transcript3.singular <- umap(transcript3.counts.norm.singular,config = custom.config)

#clustering
distances3 <- dist(transcript3.counts.norm, method = "euclidean")
distances.umap3 <- dist(umapped.transcript3$layout, method = "euclidean")
distances.umap3.singular <- dist(umapped.transcript3.singular$layout, method = "euclidean")


res.clust3 <- hclust(distances3, method = "centroid")
res.clust.umap3 <- hclust(distances.umap3, method = "centroid")
res.clust.umap3.singular <- hclust(distances.umap3.singular, method = "centroid")

cuttree3 <- cutree(res.clust3, k = 4)
cuttree.umap3 <- cutree(res.clust.umap3, k = 4)
cuttree.umap3.singular <- cutree(res.clust.umap3.singular, k = 4)

#plot data
plot_ly(x=umapped.transcript3$layout[,1],y=umapped.transcript3$layout[,2],z=umapped.transcript3$layout[,3],type="scatter3d",color = cuttree.umap3)

### PATIENT 4 ###

umapped.transcript4 <- umap(transcript4.counts.norm,config = custom.config)
umapped.transcript4.singular <- umap(transcript4.counts.norm.singular,config = custom.config)

#clustering
distances4 <- dist(transcript4.counts.norm, method = "euclidean")
distances.umap4 <- dist(umapped.transcript4$layout, method = "euclidean")
distances.umap4.singular <- dist(umapped.transcript4.singular$layout, method = "euclidean")


res.clust4 <- hclust(distances4, method = "centroid")
res.clust.umap4 <- hclust(distances.umap4, method = "centroid")
res.clust.umap4.singular <- hclust(distances.umap4.singular, method = "centroid")

cuttree4 <- cutree(res.clust4, k = 4)
cuttree.umap4 <- cutree(res.clust.umap4, k = 4)
cuttree.umap4.singular <- cutree(res.clust.umap4.singular, k = 4)

#plot data
plot_ly(x=umapped.transcript4$layout[,1],y=umapped.transcript4$layout[,2],z=umapped.transcript4$layout[,3],type="scatter3d",color = cuttree.umap4)


### differential analysis ###


## PATIENT 1 
# Prepare this data for the differential expression algorithm
# Look through the documentation for this function to find out all the options you have here
# documentation in R is found by a command like "?DESeqDataSetFromTximport")

txi.augmented <- copy(txi)
txi.augmented$counts <- txi.augmented$counts +1

samples["group"] <- factor(cuttree.umap)

dds <- DESeqDataSetFromTximport(txi = txi.augmented,           #if the group in samples has only one value, there will be problems
                                design = ~ group, #all variables in design should be the different groups in samples
                                colData = samples)
#transcript.counts.augmented <- transcript.counts.filtered + 1
#dds <- DESeqDataSetFromMatrix(t(transcript.counts.augmented),
#                              design = ~ group,
#                              colData = samples)

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

#We aren't scheduled to start visualizing yet, but it's fun to make graphs and you worked hard today so try this one: 
plotMA(res)

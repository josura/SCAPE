
library(Seurat)
library(Matrix)
library(biomaRt)


raw_counts.calu3 <- read.table(file=paste0("/home/josura/Projects/tesi/data/tesi/modernData/COVID/","GSE148729_Calu3_polyA_series2_readcounts_rev.tsv"),sep="\t")
raw_counts.calu3 <- read.csv("/home/josura/Projects/tesi/data/tesi/modernData/COVID/GSE148729_Calu3_polyA_series2_readcounts_rev.tsv",sep="\t")
#the matrix is formed by:
#- the three rows that are gene_id, gene_name, length
#- all the cell lines counts (21 for this example)

raw_counts.scrnaseq.calu3 <- read.csv("/home/josura/Projects/tesi/data/tesi/modernData/COVID/Calu3_scRNAseq_morethan1000genes_rawcounts_tr.txt",sep="\t")
# for Seurat, genes in the columns, samples-cells in the rows

metadata.scrnaseq.calu3 <- read.csv("/home/josura/Projects/tesi/data/tesi/modernData/COVID/Calu3_scRNAseq_morethan1000genes_metadata.txt",sep="\t")
# the important value to separate the various cell types into groups is the "orig.ident" that identifies the origin experiment, 
#the strain that identifies the patient condition {SARSCOV1,SARSCOV2,NA}, 
#and the "infect" feature used to separate controls from infected cells
#cell_id shouldn't be necessary since the rows/columns are already ordered and are identified

metadata.scrnaseq.calu3$strain == "SARSCoV1" && (metadata.scrnaseq.calu3$orig.ident == "Calu3-S1-4h-A" | metadata.scrnaseq.calu3$orig.ident == "Calu3-S2-4h-A")
sarscov.scrnaseq.calu3 <- raw_counts.scrnaseq.calu3[metadata.scrnaseq.calu3$strain == "SARSCoV1",]
sarscov.metadata.calu3 <- metadata.scrnaseq.calu3[metadata.scrnaseq.calu3$strain == "SARSCoV1",]

sarscov2.scrnaseq.calu3 <- raw_counts.scrnaseq.calu3[metadata.scrnaseq.calu3$strain == "SARSCoV2",]
sarscov2.metadata.calu3 <- metadata.scrnaseq.calu3[metadata.scrnaseq.calu3$strain == "SARSCoV2",]

uninf.scrnaseq.calu3 <- raw_counts.scrnaseq.calu3[metadata.scrnaseq.calu3$infect == "uninf",]
uninf.metadata.calu3 <- metadata.scrnaseq.calu3[metadata.scrnaseq.calu3$infect == "uninf",]


sarscov2.uninf.scrnaseq.calu3 <- raw_counts.scrnaseq.calu3[metadata.scrnaseq.calu3$strain == "SARSCoV2" | metadata.scrnaseq.calu3$infect == "uninf",]
sarscov2.uninf.metadata.calu3 <- metadata.scrnaseq.calu3[metadata.scrnaseq.calu3$strain == "SARSCoV2" | metadata.scrnaseq.calu3$infect == "uninf",]
write.csv(sarscov2.uninf.metadata.calu3,"/home/josura/Projects/tesi/data/tesi/modernData/COVID/sarscov-uninfmetadata.csv")
write.csv(sarscov2.uninf.scrnaseq.calu3,"/home/josura/Projects/tesi/data/tesi/modernData/COVID/sarscov-uninf.csv")


### Loading data directly
sarscov2.uninf.scrnaseq.calu3 <- read.csv("/home/josura/Projects/tesi/data/tesi/modernData/COVID/sarscov-uninf.csv")
sarscov2.uninf.metadata.calu3 <- read.csv("/home/josura/Projects/tesi/data/tesi/modernData/COVID/sarscov-uninfmetadata.csv")
sarscov2.uninf.metadataOth.calu3 <- read.csv("/home/josura/Projects/tesi/data/tesi/modernData/COVID/Calu3_Table_fullmetadata.csv")

### TODO see how the data is formed and reestablish the rownames and colnames in case

### PREPROCESSING OF THE METADATA SINCE NOT ALL THE SINGLE-CELLS ARE REPORTED WITH VIRAL LOAD(I MUST SAY THAT MOST OF THEM AREN'T EVEN INCLUDED)

#INNER JOIN
fulltest <- merge(x=sarscov2.uninf.metadata.calu3,y=sarscov2.uninf.metadataOth.calu3,by="cell_id")
fulltest <- fulltest[fulltest$type == "S2",]

library(dplyr)

test <- fulltest %>% select("cell_id","orig.ident","SCov2_Load")

sarscov2.uninf.scrnaseq.calu3$X <- NULL


filtered.points <- merge(x=test,y=sarscov2.uninf.scrnaseq.calu3,by="cell_id")

# I will be using the already available UMAP coordinates in the metadata dataset



library(ggplot2)
#not working
# p <- ggplot() +
#   geom_point(data=fulltest[fulltest$time=="4h",], aes(x=UMAPrna_1, y=UMAPrna_2, color=SCov2_Load), shape=21, size=3) +
#        scale_color_gradient(low="red", high="gray50") +
#        geom_point(data=fulltest[fulltest$time=="8h",], aes(x=UMAPrna_1, y=UMAPrna_2, color=SCov2_Load), shape=21, size=2) +
#        scale_fill_gradient(low="gray90", high="blue") +
#        geom_point(data=fulltest[fulltest$time=="12h",], aes(x=UMAPrna_1, y=UMAPrna_2, color=SCov2_Load), shape=21, size=2) +
#        scale_(low="blue", high="yellow")


p <- qplot(data=fulltest, x=UMAPrna_1, y=UMAPrna_2, colour=SCov2_Load, shape=time) +
     guides(name="viral load",colour = guide_colourbar(order=1),
            shape = guide_legend(order=2)
     )



p <- qplot(data=sarscov2.uninf.metadataOth.calu3[sarscov2.uninf.metadataOth.calu3$type=="S2",], x=UMAPrna_1, y=UMAPrna_2, colour=SCov2_Load, shape=time) +
  guides(fill=guide_legend(title = "SCov2 Load"),colour = guide_colourbar(order=1),
         shape = guide_legend(order=2)
  )

install.packages("ggnewscale")
library(ggnewscale)

p <- ggplot() +
     ggtitle("SCov2 load at different times") +
     geom_point(data=fulltest[fulltest$time=="4h",], aes(x=UMAPrna_1, y=UMAPrna_2, color=SCov2_Load)) +
          scale_color_gradient(low="red", high="gray50") +
          labs(colour="4h") +
          new_scale_color() +
          geom_point(data=fulltest[fulltest$time=="8h",], aes(x=UMAPrna_1, y=UMAPrna_2, color=SCov2_Load)) +
          scale_color_gradient(low="gray90", high="blue") +
          labs(colour="8h") +
          new_scale_color() +
          geom_point(data=fulltest[fulltest$time=="12h",], aes(x=UMAPrna_1, y=UMAPrna_2, color=SCov2_Load)) +
          scale_color_gradient(low="black", high="yellow") +
          labs(colour="12h")
  

# all the records

ggplot() +
  ggtitle("SCov2 load at 4h") +
  geom_point(data=sarscov2.uninf.metadataOth.calu3[sarscov2.uninf.metadataOth.calu3$type=="S2" & sarscov2.uninf.metadataOth.calu3$time=="4h",], aes(x=UMAPrna_1, y=UMAPrna_2, color=SCov2_Load)) +
  scale_color_gradient(low="red", high="gray50") +
  labs(colour="SCov2 load") 

ggplot() +
  ggtitle("SCov2 load at 8h") +
  geom_point(data=sarscov2.uninf.metadataOth.calu3[sarscov2.uninf.metadataOth.calu3$type=="S2" & sarscov2.uninf.metadataOth.calu3$time=="8h",], aes(x=UMAPrna_1, y=UMAPrna_2, color=SCov2_Load)) +
  scale_color_gradient(low="gray90", high="blue") +
  labs(colour="SCov2 load") 

ggplot() +
  ggtitle("SCov2 load at 12h") +
  geom_point(data=sarscov2.uninf.metadataOth.calu3[sarscov2.uninf.metadataOth.calu3$type=="S2" & sarscov2.uninf.metadataOth.calu3$time=="12h",], aes(x=UMAPrna_1, y=UMAPrna_2, color=SCov2_Load)) +
  scale_color_gradient(low="black", high="yellow") +
  labs(colour="SCov2 load")

p <- ggplot() +
  ggtitle("SCov2 load at different times") +
  geom_point(data=sarscov2.uninf.metadataOth.calu3[sarscov2.uninf.metadataOth.calu3$type=="S2" & sarscov2.uninf.metadataOth.calu3$time=="4h",], aes(x=UMAPrna_1, y=UMAPrna_2, color=SCov2_Load)) +
  scale_color_gradient(low="red", high="gray50") +
  labs(colour="4h") +
  new_scale_color() +
  geom_point(data=sarscov2.uninf.metadataOth.calu3[sarscov2.uninf.metadataOth.calu3$type=="S2" & sarscov2.uninf.metadataOth.calu3$time=="8h",], aes(x=UMAPrna_1, y=UMAPrna_2, color=SCov2_Load)) +
  scale_color_gradient(low="gray90", high="blue") +
  labs(colour="8h") +
  new_scale_color() +
  geom_point(data=sarscov2.uninf.metadataOth.calu3[sarscov2.uninf.metadataOth.calu3$type=="S2" & sarscov2.uninf.metadataOth.calu3$time=="12h",], aes(x=UMAPrna_1, y=UMAPrna_2, color=SCov2_Load)) +
  scale_color_gradient(low="black", high="yellow") +
  labs(colour="12h")

#mock data

p <- ggplot() +
  ggtitle("SCov2 load for mock data different times") +
  geom_point(data=sarscov2.uninf.metadataOth.calu3[sarscov2.uninf.metadataOth.calu3$type=="mock" & sarscov2.uninf.metadataOth.calu3$time=="4h",], aes(x=UMAPrna_1, y=UMAPrna_2, color=SCov2_Load)) +
  scale_color_gradient(low="red", high="gray50") +
  labs(colour="4h") +
  new_scale_color() +
  geom_point(data=sarscov2.uninf.metadataOth.calu3[sarscov2.uninf.metadataOth.calu3$type=="mock" & sarscov2.uninf.metadataOth.calu3$time=="12h",], aes(x=UMAPrna_1, y=UMAPrna_2, color=SCov2_Load)) +
  scale_color_gradient(low="black", high="yellow") +
  labs(colour="12h")




#LEFT JOIN to use for maybe regression for the viral load 
test <- merge(x=sarscov2.uninf.metadata.calu3,y=sarscov2.uninf.metadataOth.calu3,by="cell_id",all.x=TRUE)

###

### VIsualization of these points after filtering and preprocessing

sarscov2.scrnaseq.calu3 <- sarscov2.uninf.scrnaseq.calu3[sarscov2.uninf.metadata.calu3$strain == "SARSCoV2",]
sarscov2.metadata.calu3 <- sarscov2.uninf.metadata.calu3[sarscov2.uninf.metadata.calu3$strain == "SARSCoV2",]

sarscov2.uninf.scrnaseq.calu3$cell_id <- NULL
sarscov2.uninf.scrnaseq.calu3$X <- NULL
rownames(sarscov2.uninf.scrnaseq.calu3) <- sarscov2.uninf.metadata.calu3$cell_id

seurat.sarscov2.uninf.calu3 <- CreateSeuratObject(counts = sarscov2.uninf.scrnaseq.calu3.transpose, min.cells = 3, min.genes = 200, project = "calu3_sarscov2_scRNAseq")
#seurat.sarscov2.calu3 <- CreateSeuratObject(counts = sarscov2.scrnaseq.calu3, min.cells = 3, min.genes = 200, project = "calu3_sarscov2_scRNAseq")
seurat.sarscov2.calu3 <- CreateSeuratObject(counts = sarscov2.uninf.scrnaseq.calu3[sarscov2.uninf.metadata.calu3$strain == "SARSCoV2",], min.cells = 3, min.genes = 200, project = "calu3_sarscov2_scRNAseq")

rm(sarscov2.uninf.scrnaseq.calu3)

### LOADING THE SEURAT OBJECT DIRECTLY

saveRDS(seurat.sarscov2.uninf.calu3,file = "/home/josura/Projects/tesi/data/sarscov-seurat.rds")
seurat.sarscov2.uninf.calu3 <- readRDS("/home/josura/Projects/tesi/data/sarscov-seurat.rds")


### USELESS PART
#visualization of sarscov2 clusters based on "viral RNA quantities" or whatever since I do not see anything explicit in the metadata
seurat.sarscov2.calu3 <- ScaleData(seurat.sarscov2.calu3)
seurat.sarscov2.calu3 <- FindVariableFeatures(seurat.sarscov2.calu3, selection.method = "vst", nfeatures = 40)

seurat.sarscov2.calu3 <- RunPCA(seurat.sarscov2.calu3, npcs = 30, verbose = FALSE,approx=FALSE)
seurat.sarscov2.calu3 <- RunUMAP(seurat.sarscov2.calu3,reduction = "pca",dims=1:30)

seurat.sarscov2.calu3@reductions$umap@cell.embeddings #UMAP embeddings
### END USELESS PART



### DIFFERENTIAL ANALYSIS
seurat.sarscov2.uninf.calu3 <- NormalizeData(seurat.sarscov2.uninf.calu3)
seurat.sarscov2.uninf.calu3 <- FindVariableFeatures(seurat.sarscov2.uninf.calu3, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat.sarscov2.uninf.calu3)
all.cells <- colnames(seurat.sarscov2.uninf.calu3)
seurat.sarscov2.uninf.calu3 <- ScaleData(seurat.sarscov2.uninf.calu3, features = all.genes)

all.markers <- FindAllMarkers(seurat.sarscov2.uninf.calu3, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)


FeaturePlot(mydata4, features = all.markers$gene, min.cutoff = "q9")
#sarscov2.uninf.metadata.calu3$cell_id <- gsub("_","-",sarscov2.uninf.metadata.calu3$cell_id)
### adding the time metadata to the seurat object
filtered.metadata <- sarscov2.uninf.metadata.calu3[sarscov2.uninf.metadata.calu3$cell_id %in% all.cells,]
seurat.sarscov2.uninf.calu3@meta.data$origin.experiment <- filtered.metadata$orig.ident


seurat.sarscov2.uninf.calu3 <- SetIdent(seurat.sarscov2.uninf.calu3, value = "origin.experiment")
DimPlot(seurat.sarscov2.uninf.calu3, label = T , repel = T, label.size = 3) + NoLegend()

table(seurat.sarscov2.uninf.calu3@meta.data$origin.experiment)


## 4H
diffexpr.4h <- FindMarkers(seurat.sarscov2.uninf.calu3, ident.1 = "Calu3-S2-4h-A", ident.2 = "Calu3-mock-4h-A")
# view results
head(diffexpr.4h)

library(ggplot2)
# The basic scatter plot: x is "log2FoldChange", y is "pvalue"
ggplot(data=diffexpr.4h, aes(x=avg_log2FC, y=p_val)) + geom_point()
p <- ggplot(data=diffexpr.4h, aes(x=avg_log2FC, y=-log10(p_val))) + geom_point()
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +  geom_hline(yintercept=-log10(0.05), col="red")


###preparing phensim input
diffexpr.4h$diffexpressed <- "NO"
# if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
diffexpr.4h$diffexpressed[diffexpr.4h$avg_log2FC > 0.5 & diffexpr.4h$p_val < 0.05] <- "OVEREXPRESSION"
# if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
diffexpr.4h$diffexpressed[diffexpr.4h$avg_log2FC < -0.5 & diffexpr.4h$p_val < 0.05] <- "UNDEREXPRESSION"



mart <- useDataset("hsapiens_gene_ensembl",useMart("ensembl")) 


Glist.4h <- getBM(filters = "hgnc_symbol", attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene_id","description"),values = rownames(diffexpr.4h),mart = mart)

diffexpr.4h["hgnc_symbol"] <- rownames(diffexpr.4h)

diffexpr.4h.final <- merge(x=Glist.4h,y=diffexpr.4h)

diffexpr.4h.PHENSIM.input <- diffexpr.4h.final %>%
  filter(diffexpressed =="OVEREXPRESSION" | diffexpressed == "UNDEREXPRESSION") %>%
  dplyr::select(entrezgene_id,diffexpressed) %>%
  unique()



write.table(diffexpr.4h.PHENSIM.input,
            file='/home/josura/Projects/tesi/data/modernData/COVID/diffexpr4h-PHENSIMinput.txt', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)


## 8H
diffexpr.8h <- FindMarkers(seurat.sarscov2.uninf.calu3, ident.1 = "Calu3-S2-8h-A", ident.2 = "Calu3-mock-12h-A")
# view results
head(diffexpr.8h)

# The basic scatter plot: x is "log2FoldChange", y is "pvalue"
ggplot(data=diffexpr.8h, aes(x=avg_log2FC, y=p_val)) + geom_point()
p <- ggplot(data=diffexpr.8h, aes(x=avg_log2FC, y=-log10(p_val))) + geom_point()
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +  geom_hline(yintercept=-log10(0.05), col="red")


###preparing phensim input
diffexpr.8h$diffexpressed <- "NO"
# if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
diffexpr.8h$diffexpressed[diffexpr.8h$avg_log2FC > 0.5 & diffexpr.8h$p_val < 0.05] <- "OVEREXPRESSION"
# if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
diffexpr.8h$diffexpressed[diffexpr.8h$avg_log2FC < -0.5 & diffexpr.8h$p_val < 0.05] <- "UNDEREXPRESSION"


mart <- useDataset("hsapiens_gene_ensembl",useMart("ensembl")) 


Glist.8h <- getBM(filters = "hgnc_symbol", attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene_id","description"),values = rownames(diffexpr.8h),mart = mart)

diffexpr.8h["hgnc_symbol"] <- rownames(diffexpr.8h)

diffexpr.8h.final <- merge(x=Glist.8h,y=diffexpr.8h)

diffexpr.8h.PHENSIM.input <- diffexpr.8h.final %>%
  filter(diffexpressed =="OVEREXPRESSION" | diffexpressed == "UNDEREXPRESSION") %>%
  dplyr::select(entrezgene_id,diffexpressed)%>%
  unique()


write.table(diffexpr.8h.PHENSIM.input,
            file='/home/josura/Projects/tesi/data/modernData/COVID/diffexpr8h-PHENSIMinput.txt', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)




## 12H
diffexpr.12h <- FindMarkers(seurat.sarscov2.uninf.calu3, ident.1 = "Calu3-S2-12h-A", ident.2 = "Calu3-mock-12h-A")
# view results
head(diffexpr.12h)

# The basic scatter plot: x is "log2FoldChange", y is "pvalue"
ggplot(data=diffexpr.12h, aes(x=avg_log2FC, y=p_val)) + geom_point()
p <- ggplot(data=diffexpr.12h, aes(x=avg_log2FC, y=-log10(p_val))) + geom_point()
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +  geom_hline(yintercept=-log10(0.05), col="red")


###preparing phensim input
diffexpr.12h$diffexpressed <- "NO"
# if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
diffexpr.12h$diffexpressed[diffexpr.12h$avg_log2FC > 0.5 & diffexpr.12h$p_val < 0.05] <- "OVEREXPRESSION"
# if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
diffexpr.12h$diffexpressed[diffexpr.12h$avg_log2FC < -0.5 & diffexpr.12h$p_val < 0.05] <- "UNDEREXPRESSION"

Glist.12h <- getBM(filters = "hgnc_symbol", attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene_id","description"),values = rownames(diffexpr.12h),mart = mart)

diffexpr.12h["hgnc_symbol"] <- rownames(diffexpr.12h)

diffexpr.12h.final <- merge(x=Glist.12h,y=diffexpr.12h)

diffexpr.12h.PHENSIM.input <- diffexpr.12h.final %>%
  filter(diffexpressed =="OVEREXPRESSION" | diffexpressed == "UNDEREXPRESSION") %>%
  dplyr::select(entrezgene_id,diffexpressed) %>%
  unique()


write.table(diffexpr.12h.PHENSIM.input,
            file='/home/josura/Projects/tesi/data/modernData/COVID/diffexpr12h-PHENSIMinput.txt', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)


### BULK differential analysis



bulk.calu3.series1 <- read.csv("/home/josura/Projects/tesi/data/tesi/modernData/COVID/Calu3_polyA_series1_readcounts.txt",sep="\t")
bulk.calu3.series2 <- read.csv("/home/josura/Projects/tesi/data/tesi/modernData/COVID/GSE148729_Calu3_polyA_series2_readcounts_rev.tsv",sep="\t")

bulk.calu3.series1 <- bulk.calu3.series1[bulk.calu3.series1$gene_id %in% bulk.calu3.series2$gene_id,]
bulk.calu3.series2 <- bulk.calu3.series2[bulk.calu3.series2$gene_id %in% bulk.calu3.series1$gene_id,]

bulk.calu3 <- merge(x=bulk.calu3.series1,y=bulk.calu3.series2,by=c("gene_id","gene_name","length"))

bulk.calu3.metadata <- read.csv("/home/josura/Projects/tesi/data/tesi/modernData/COVID/Calu3_polyA_metadata.tsv",sep="\t")

colnames(bulk.calu3)

library("DESeq2")

#4h SarsCov2
bulk.4h.calu3.metadata <- bulk.calu3.metadata[bulk.calu3.metadata$time == "4h" & (bulk.calu3.metadata$condition == "S2" | bulk.calu3.metadata$condition == "mock"),]

bulk.4h.calu3 <- bulk.calu3 %>%
  dplyr::select(gene_name,bulk.4h.calu3.metadata$name) %>%
  group_by(gene_name) %>% 
  summarise(across(everything(),sum))

bulk.genes.4h.calu3 <- bulk.4h.calu3$gene_name
bulk.4h.calu3 <- as.data.frame(bulk.4h.calu3)
bulk.4h.calu3.genenames <- bulk.4h.calu3$gene_name 
bulk.4h.calu3$gene_name <- NULL

bulk.4h.calu3.metadata$condition <- as.factor(bulk.4h.calu3.metadata$condition) 

# transform to integers since deseq2 needs it
bulk.4h.calu3 <- sapply(bulk.4h.calu3,as.integer)
rownames(bulk.4h.calu3) <- bulk.4h.calu3.genenames


dds <- DESeqDataSetFromMatrix(countData = bulk.4h.calu3,
                              colData = bulk.4h.calu3.metadata,
                              design = ~ condition)

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
res$diffexpressed[res$log2FoldChange > 0.5 & res$padj < 0.05] <- "OVEREXPRESSION"
# if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
res$diffexpressed[res$log2FoldChange < -0.5 & res$padj < 0.05] <- "UNDEREXPRESSION"

mart <- useDataset("hsapiens_gene_ensembl",useMart("ensembl")) 

Glist.4h.bulk <- getBM(filters = "hgnc_symbol", attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene_id","description"),values = rownames(res),mart = mart)

res["hgnc_symbol"] <- rownames(res)

diffexpr.4h.bulk.final <- merge(x=Glist.4h.bulk,y=as.data.frame(res))

diffexpr.4h.bulk.PHENSIM.input <- diffexpr.4h.bulk.final %>%
  filter(diffexpressed =="OVEREXPRESSION" | diffexpressed == "UNDEREXPRESSION") %>%
  dplyr::select(entrezgene_id,diffexpressed) %>%
  unique()


write.table(diffexpr.4h.bulk.PHENSIM.input,
            file='/home/josura/Projects/tesi/data/modernData/COVID/diffexpr4h-bulk-PHENSIMinput.txt', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)

#8h SarsCov2
bulk.8h.calu3.metadata <- bulk.calu3.metadata[(bulk.calu3.metadata$time == "8h" | (bulk.calu3.metadata$time == "12h" & bulk.calu3.metadata$condition == "mock")) & (bulk.calu3.metadata$condition == "S2" | bulk.calu3.metadata$condition == "mock"),]

bulk.8h.calu3 <- bulk.calu3 %>%
  dplyr::select(gene_name,bulk.8h.calu3.metadata$name) %>%
  group_by(gene_name) %>% 
  summarise(across(everything(),sum))

bulk.genes.8h.calu3 <- bulk.8h.calu3$gene_name
bulk.8h.calu3 <- as.data.frame(bulk.8h.calu3)
bulk.8h.calu3.genenames <- bulk.8h.calu3$gene_name 
bulk.8h.calu3$gene_name <- NULL

bulk.8h.calu3.metadata$condition <- as.factor(bulk.8h.calu3.metadata$condition) 

# transform to integers since deseq2 needs it
bulk.8h.calu3 <- sapply(bulk.8h.calu3,as.integer)
rownames(bulk.8h.calu3) <- bulk.8h.calu3.genenames


dds <- DESeqDataSetFromMatrix(countData = bulk.8h.calu3,
                              colData = bulk.8h.calu3.metadata,
                              design = ~ condition)

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
res$diffexpressed[res$log2FoldChange > 0.5 & res$padj < 0.05] <- "OVEREXPRESSION"
# if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
res$diffexpressed[res$log2FoldChange < -0.5 & res$padj < 0.05] <- "UNDEREXPRESSION"

Glist.8h.bulk <- getBM(filters = "hgnc_symbol", attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene_id","description"),values = rownames(res),mart = mart)

res["hgnc_symbol"] <- rownames(res)

diffexpr.8h.bulk.final <- merge(x=Glist.8h.bulk,y=as.data.frame(res))

diffexpr.8h.bulk.PHENSIM.input <- diffexpr.8h.bulk.final %>%
  filter(diffexpressed =="OVEREXPRESSION" | diffexpressed == "UNDEREXPRESSION") %>%
  dplyr::select(entrezgene_id,diffexpressed) %>%
  unique()


write.table(diffexpr.8h.bulk.PHENSIM.input,
            file='/home/josura/Projects/tesi/data/modernData/COVID/diffexpr8h-bulk-PHENSIMinput.txt', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)


#12h SarsCov2
bulk.12h.calu3.metadata <- bulk.calu3.metadata[bulk.calu3.metadata$time == "12h" & (bulk.calu3.metadata$condition == "S2" | bulk.calu3.metadata$condition == "mock"),]

bulk.12h.calu3 <- bulk.calu3 %>%
  dplyr::select(gene_name,bulk.12h.calu3.metadata$name) %>%
  group_by(gene_name) %>% 
  summarise(across(everything(),sum))

bulk.genes.12h.calu3 <- bulk.12h.calu3$gene_name
bulk.12h.calu3 <- as.data.frame(bulk.12h.calu3)
bulk.12h.calu3.genenames <- bulk.12h.calu3$gene_name 
bulk.12h.calu3$gene_name <- NULL

bulk.12h.calu3.metadata$condition <- as.factor(bulk.12h.calu3.metadata$condition) 

# transform to integers since deseq2 needs it
bulk.12h.calu3 <- sapply(bulk.12h.calu3,as.integer)
rownames(bulk.12h.calu3) <- bulk.12h.calu3.genenames


dds <- DESeqDataSetFromMatrix(countData = bulk.12h.calu3,
                              colData = bulk.12h.calu3.metadata,
                              design = ~ condition)

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
res$diffexpressed[res$log2FoldChange > 0.5 & res$padj < 0.05] <- "OVEREXPRESSION"
# if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
res$diffexpressed[res$log2FoldChange < -0.5 & res$padj < 0.05] <- "UNDEREXPRESSION"

Glist.12h.bulk <- getBM(filters = "hgnc_symbol", attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene_id","description"),values = rownames(res),mart = mart)

res["hgnc_symbol"] <- rownames(res)

diffexpr.12h.bulk.final <- merge(x=Glist.12h.bulk,y=as.data.frame(res))

diffexpr.12h.bulk.PHENSIM.input <- diffexpr.12h.bulk.final %>%
  filter(diffexpressed =="OVEREXPRESSION" | diffexpressed == "UNDEREXPRESSION") %>%
  dplyr::select(entrezgene_id,diffexpressed) %>%
  unique()


write.table(diffexpr.12h.bulk.PHENSIM.input,
            file='/home/josura/Projects/tesi/data/modernData/COVID/diffexpr12h-bulk-PHENSIMinput.txt', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)

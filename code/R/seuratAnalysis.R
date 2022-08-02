
install.packages("rgeos")
### seurat and dependencies have some problems with linux
### download rgeos dependency with the package manager

install.packages("SeuratData")
library(Seurat)
library(Matrix)


raw_counts<-read.table(file=paste0("/home/josura/Projects/tesi/data/patient1ALL/","allmerged.csv"),sep=",")
samples.names <- raw_counts[,1][-1]

raw_counts$V1 <- NULL

genesNames <- as.vector((unlist(raw_counts[1,])))
all.genes <- genesNames
raw_counts <- raw_counts[-1,]


rownames(raw_counts) <- samples.names
colnames(raw_counts) <- genesNames

#count zeros as percentages in genes expressions
res <- colSums(raw_counts==0)/nrow(raw_counts)*100
#list genes with non-zero percentage of expression
res.nonzero <- res[which(res<100)]


raw_counts <- as.data.frame(t(raw_counts[,names(res.nonzero)]))

genesNames <- rownames(raw_counts)



## traduce ensembl gene ids to gene names

BiocManager::install("biomaRt")
library(biomaRt)

load("/home/josura/Projects/tesi/data/ensemblToGene.RData")  #load mart object since the ensembl site goes down sometimes

mart <- useDataset("hsapiens_gene_ensembl",useMart("ensembl")) 


Glist <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id","hgnc_symbol","description"),values = genesNames,mart = mart)
allGlist <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id","hgnc_symbol","description"),values = all.genes,mart = mart)
save(mart,Glist,allGlist,file = "/home/josura/Projects/tesi/data/ensemblToGene.RData")

library(dplyr) 
genesNames <- (Glist %>% distinct(ensembl_gene_id, .keep_all = TRUE))["hgnc_symbol"][,1]
all.genes <- (allGlist %>% distinct(ensembl_gene_id, .keep_all = TRUE))["hgnc_symbol"][,1]

raw_counts["genenames"] <- genesNames

raw_counts <- (raw_counts %>% distinct(genenames, .keep_all = TRUE))

genesNames.filtered <- raw_counts["genenames"]

raw_counts["genenames"] <- NULL


rownames(raw_counts) <- as.vector(unlist(genesNames.filtered))
colnames(raw_counts) <- samples.names


mydata <- CreateSeuratObject(counts = raw_counts, min.cells = 3, min.genes = 200, project = "mydata_scRNAseq")


mydata <- ScaleData(mydata)
mydata <- FindVariableFeatures(mydata, selection.method = "vst", nfeatures = 40)

mydata <- RunPCA(mydata, npcs = 30, verbose = FALSE,approx=FALSE)
mydata <- RunUMAP(mydata,reduction = "pca",dims=1:30)

mydata <- FindNeighbors(mydata, reduction = "pca", dims = 1:30)
mydata <- FindClusters(mydata, resolution = 0.5)


patient1Graph <- DimPlot(mydata, reduction = "umap", label = TRUE, repel = TRUE)

### installing metap package to find conserved markers

BiocManager::install('multtest')
install.packages('metap')

### doing differential analysis on clusters previously found with RNA counts and not with transformed data
DefaultAssay(mydata) <- "RNA"


mydata <- NormalizeData(mydata)
mydata <- FindVariableFeatures(mydata, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(mydata)
mydata <- ScaleData(mydata, features = all.genes)

all.markers <- FindAllMarkers(mydata, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)

# mydata[["groups"]] <- 
# patient1.markers <- FindConservedMarkers(mydata,grouping.var = "seurat_clusters", ident.1 = 0, verbose = FALSE)
# head(nk.markers)

FeaturePlot(mydata, features = all.markers$gene, min.cutoff = "q9")

### ANNOTATION USING THE MARKERS
BiocManager::install('SingleR')
BiocManager::install('celldex')

library(SingleR)
library(celldex)

# download the reference markers atlas for humans
hpca.ref <- celldex::HumanPrimaryCellAtlasData()


BiocManager::install('SingleCellExperiment')

# converting to DietSeurat for convenience
sce <- as.SingleCellExperiment(DietSeurat(mydata))

hpca.main <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.main)
hpca.fine <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.fine)

table(hpca.main$pruned.labels)
table(hpca.fine$pruned.labels)

mydata@meta.data$hpca.main <- hpca.main$pruned.labels
mydata@meta.data$hpca.fine <- hpca.fine$pruned.labels

mydata <- SetIdent(mydata, value = "hpca.fine")
DimPlot(mydata, label = T , repel = T, label.size = 3) + NoLegend()

### DIFFERENTIAL EXPRESSION TESTING

astrocyteVStcell.de.markers <- FindMarkers(mydata, ident.1 = "Astrocyte:Embryonic_stem_cell-derived", ident.2 = "T_cell:CCR10-CLA+1,25(OH)2_vit_D3/IL-12")
# view results
head(astrocyteVStcell.de.markers)

library(ggplot2)
# The basic scatter plot: x is "log2FoldChange", y is "pvalue"
ggplot(data=astrocyteVStcell.de.markers, aes(x=avg_log2FC, y=p_val)) + geom_point()
p <- ggplot(data=astrocyteVStcell.de.markers, aes(x=avg_log2FC, y=-log10(p_val))) + geom_point()
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +  geom_hline(yintercept=-log10(0.05), col="red")

# add a column of expression labelling
astrocyteVStcell.de.markers$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
astrocyteVStcell.de.markers$diffexpressed[astrocyteVStcell.de.markers$avg_log2FC > 0.6 & astrocyteVStcell.de.markers$p_val < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
astrocyteVStcell.de.markers$diffexpressed[astrocyteVStcell.de.markers$avg_log2FC < -0.6 & astrocyteVStcell.de.markers$p_val < 0.05] <- "DOWN"
#adding nemes as the rownames
p <- ggplot(data=astrocyteVStcell.de.markers, aes(x=avg_log2FC, y=-log10(p_val), col=diffexpressed, label=rownames(astrocyteVStcell.de.markers))) + geom_point() + theme_minimal() + geom_text()
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +  geom_hline(yintercept=-log10(0.05), col="red")



## preparing the data to use phensim


all.genes.filtered <- distinct(data.frame(all.genes[!all.genes %in% rownames(astrocyteVStcell.de.markers)]))

all.genes.filtered <- all.genes.filtered$all.genes..all.genes..in..rownames.astrocyteVStcell.de.markers..

diffExpressedSim <- data.frame(rep("NO",length(all.genes.filtered[!all.genes.filtered %in% rownames(astrocyteVStcell.de.markers)])),row.names = all.genes.filtered[!all.genes.filtered %in% rownames(astrocyteVStcell.de.markers)])
colnames(diffExpressedSim) <- "diffexpressed"

diffExpressedSim <- rbind(diffExpressedSim, select(astrocyteVStcell.de.markers,c("diffexpressed")))

diffExpressedSim["hgnc_symbol"] <- rownames(diffExpressedSim)
rownames(diffExpressedSim) <- NULL
# querying for entrez_id
entrez.list <- getBM(filters = "hgnc_symbol", attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene_id","description"),values = diffExpressedSim["hgnc_symbol"],mart = mart)
entrez.list <- (entrez.list %>% distinct(hgnc_symbol, .keep_all = TRUE))

full.phensim.input <- merge(diffExpressedSim,select(entrez.list,c("hgnc_symbol","entrezgene_id")))

full.phensim.input <- full.phensim.input[!is.na(full.phensim.input$entrezgene_id),]

full.phensim.input$diffexpressed[full.phensim.input$diffexpressed == "UP"] <- "OVEREXPRESSION"
full.phensim.input$diffexpressed[full.phensim.input$diffexpressed == "DOWN"] <- "UNDEREXPRESSION"

write.table((select(full.phensim.input,c("entrezgene_id","diffexpressed")))[!full.phensim.input$diffexpressed=="NO",],
            file='/home/josura/Projects/tesi/data/patient1ALL/full.phensim.input.tsv', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)


### patient 2 



raw_counts2<-read.table(file=paste0("/home/josura/Projects/tesi/data/patient2ALL/","allmerged.csv"),sep=",")
samplesNames2 <- raw_counts2[,1]

raw_counts2$V1 <- NULL
##useless since the gene names are the same for all the patients, not true if some genes were filtered before
genesNames2 <- as.vector((unlist(raw_counts2[1,])))
raw_counts2 <- raw_counts2[-1,]

raw_counts2 <- t(raw_counts2)


rownames(raw_counts2) <- genesNames2

mydata2 <- CreateSeuratObject(counts = raw_counts2, min.cells = 3, min.genes = 200, project = "mydata_scRNAseq",row.names = genesNames2)


mydata2 <- ScaleData(mydata2)
mydata2 <- FindVariableFeatures(mydata2, selection.method = "vst", nfeatures = 40)

mydata2 <- RunPCA(mydata2, npcs = 30, verbose = FALSE,approx=FALSE)
mydata2 <- RunUMAP(mydata2,reduction = "pca",dims=1:30)

mydata2 <- FindNeighbors(mydata2, reduction = "pca", dims = 1:30)
mydata2 <- FindClusters(mydata2, resolution = 0.5)


patient2Graph <- DimPlot(mydata2, reduction = "umap", label = TRUE, repel = TRUE)

### patient 3 


raw_counts3<-read.table(file=paste0("/home/josura/Projects/tesi/data/patient2ALL/","allmerged.csv"),sep=",")
samplesNames3 <- raw_counts3[,1]

raw_counts3$V1 <- NULL

genesNames3 <- as.vector((unlist(raw_counts3[1,])))
raw_counts3 <- raw_counts3[-1,]

raw_counts3 <- t(raw_counts3)


rownames(raw_counts3) <- genesNames3

mydata3 <- CreateSeuratObject(counts = raw_counts3, min.cells = 3, min.genes = 200, project = "mydata_scRNAseq",row.names = genesNames3)


mydata3 <- ScaleData(mydata3)
mydata3 <- FindVariableFeatures(mydata3, selection.method = "vst", nfeatures = 40)

mydata3 <- RunPCA(mydata3, npcs = 30, verbose = FALSE,approx=FALSE)
mydata3 <- RunUMAP(mydata3,reduction = "pca",dims=1:30)

mydata3 <- FindNeighbors(mydata3, reduction = "pca", dims = 1:30)
mydata3 <- FindClusters(mydata3, resolution = 0.5)


patient3Graph <- DimPlot(mydata3, reduction = "umap", label = TRUE, repel = TRUE)

### patient 4 



raw_counts4<-read.table(file=paste0("/home/josura/Projects/tesi/data/patient4ALL/","allmerged.csv"),sep=",")
samplesNames4 <- raw_counts4[,1]

raw_counts4$V1 <- NULL

genesNames4 <- as.vector((unlist(raw_counts4[1,])))
raw_counts4 <- raw_counts4[-1,]

raw_counts4 <- t(raw_counts4)


rownames(raw_counts4) <- genesNames4

mydata4 <- CreateSeuratObject(counts = raw_counts4, min.cells = 3, min.genes = 200, project = "mydata_scRNAseq",row.names = genesNames4)


mydata4 <- ScaleData(mydata4)
mydata4 <- FindVariableFeatures(mydata4, selection.method = "vst", nfeatures = 40)

mydata4 <- RunPCA(mydata4, npcs = 30, verbose = FALSE,approx=FALSE)
mydata4 <- RunUMAP(mydata4,reduction = "pca",dims=1:30)

mydata4 <- FindNeighbors(mydata4, reduction = "pca", dims = 1:30)
mydata4 <- FindClusters(mydata4, resolution = 0.5)


patient4Graph <- DimPlot(mydata4, reduction = "umap", label = TRUE, repel = TRUE)


DefaultAssay(mydata4) <- "RNA"


mydata4 <- NormalizeData(mydata4)
mydata4 <- FindVariableFeatures(mydata4, selection.method = "vst", nfeatures = 2000)
all.genes4 <- rownames(mydata4)
mydata4 <- ScaleData(mydata4, features = all.genes4)

all.markers4 <- FindAllMarkers(mydata4, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)

# mydata[["groups"]] <- 
# patient1.markers <- FindConservedMarkers(mydata,grouping.var = "seurat_clusters", ident.1 = 0, verbose = FALSE)
# head(nk.markers)




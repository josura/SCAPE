
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

write.table((select(full.phensim.input,c("entrezgene_id","diffexpressed")))[full.phensim.input$diffexpressed=="NO",],
            file='/home/josura/Projects/tesi/data/patient1ALL/fullNonExpressed.tsv', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)


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

samplesNames4 <- raw_counts4[,1][-1]

raw_counts4$V1 <- NULL

genesNames4 <- as.vector((unlist(raw_counts4[1,])))
all.genes4 <- genesNames4
raw_counts4 <- raw_counts4[-1,]

rownames(raw_counts4) <- samples.names4
colnames(raw_counts4) <- genesNames4

#count zeros as percentages in genes expressions
res4 <- colSums(raw_counts4==0)/nrow(raw_counts4)*100
#list genes with non-zero percentage of expression
res.nonzero4 <- res4[which(res4<100)]

raw_counts4 <- as.data.frame(t(raw_counts4[,names(res.nonzero4)]))  #t(raw_counts4)

genesNames4 <- rownames(raw_counts4)


## traduce ensembl gene ids to gene names

BiocManager::install("biomaRt")
library(biomaRt)

load("/home/josura/Projects/tesi/data/ensemblToGene4.RData")  #load mart object since the ensembl site goes down sometimes

mart <- useDataset("hsapiens_gene_ensembl",useMart("ensembl")) 


Glist4 <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene_id","description"),values = genesNames4,mart = mart)
allGlist4 <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene_id","description"),values = all.genes4,mart = mart)
save(mart,Glist4,allGlist4,file = "/home/josura/Projects/tesi/data/ensemblToGene4.RData")

library(dplyr) 
genesNames4 <- (Glist4 %>% distinct(ensembl_gene_id, .keep_all = TRUE))["hgnc_symbol"][,1]
all.genes4 <- (allGlist4 %>% distinct(ensembl_gene_id, .keep_all = TRUE))["hgnc_symbol"][,1]

raw_counts4["genenames"] <- genesNames4

raw_counts4 <- (raw_counts4 %>% distinct(genenames, .keep_all = TRUE))

genesNames.filtered4 <- raw_counts4["genenames"]

raw_counts4["genenames"] <- NULL


rownames(raw_counts4) <- as.vector(unlist(genesNames.filtered4))
colnames(raw_counts4) <- samplesNames4


#using the entrezgene_id since ensemble_gene_ids contain non-standard characters like -

entrez.list4 <- getBM(filters = "hgnc_symbol", attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene_id","description"),values = genesNames4,mart = mart)
entrez.list4 <- (entrez.list4 %>% distinct(entrezgene_id, .keep_all = TRUE)) %>% na.omit(entrezgene_id) %>% distinct(hgnc_symbol, .keep_all = TRUE)
  
  
raw_counts.entrez4 <- raw_counts4[genesNames.filtered4$genenames %in% entrez.list4$hgnc_symbol,]
rownames(raw_counts.entrez4) <- entrez.list4$entrezgene_id

raw_counts4 <- raw_counts4[genesNames.filtered4$genenames %in% entrez.list4$hgnc_symbol,]
rownames(raw_counts4) <- entrez.list4$hgnc_symbol

#mydata4 <- CreateSeuratObject(counts = raw_counts4, min.cells = 3, min.genes = 200, project = "mydata_scRNAseq",row.names = genesNames4)
mydata4 <- CreateSeuratObject(counts = raw_counts.entrez4, min.cells = 3, min.genes = 200, project = "mydata_scRNAseq")
mydata4 <- CreateSeuratObject(counts = raw_counts4, min.cells = 3, min.genes = 200, project = "mydata_scRNAseq")


mydata4 <- ScaleData(mydata4)
mydata4 <- FindVariableFeatures(mydata4, selection.method = "vst", nfeatures = 40)

mydata4 <- RunPCA(mydata4, npcs = 30, verbose = FALSE,approx=FALSE)
mydata4 <- RunUMAP(mydata4,reduction = "pca",dims=1:30)

mydata4 <- FindNeighbors(mydata4, reduction = "pca", dims = 1:30)
mydata4 <- FindClusters(mydata4, resolution = 0.5)


patient4Graph <- DimPlot(mydata4, reduction = "umap", label = TRUE, repel = TRUE)


#not using the reduced dimenzionality but the original counts
DefaultAssay(mydata4) <- "RNA"


mydata4 <- NormalizeData(mydata4)
mydata4 <- FindVariableFeatures(mydata4, selection.method = "vst", nfeatures = 2000)
all.genes4 <- rownames(mydata4)
mydata4 <- ScaleData(mydata4, features = all.genes4)

all.markers4 <- FindAllMarkers(mydata4, only.pos = T, min.pct = 0.5, logfc.threshold = 0.4)
#all.markers4 <- FindConservedMarkers(mydata4, only.pos = T, min.pct = 0.5, logfc.threshold = 0.4)

# mydata[["groups"]] <- 
# patient1.markers <- FindConservedMarkers(mydata,grouping.var = "seurat_clusters", ident.1 = 0, verbose = FALSE)
# head(nk.markers)


FeaturePlot(mydata4, features = all.markers4$gene, min.cutoff = "q9")

### ANNOTATION USING THE MARKERS

library(SingleR)
library(celldex)

# download the reference markers atlas for humans
hpca.ref <- celldex::HumanPrimaryCellAtlasData()

# converting to DietSeurat for convenience
sce4 <- as.SingleCellExperiment(DietSeurat(mydata4))

hpca.main4 <- SingleR(test = sce4,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.main)
hpca.fine4 <- SingleR(test = sce4,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.fine, aggr.ref = TRUE) #pseudobulk

clusters <- as.numeric(levels(mydata4$seurat_clusters))[mydata4$seurat_clusters]
hpca.fine4.clusters <- SingleR(test = sce4,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.fine, aggr.ref = TRUE,clusters = mydata4$seurat_clusters) #with clusters previously found

table(hpca.main4$pruned.labels)
table(hpca.fine4$pruned.labels)
table(hpca.fine4.clusters$pruned.labels)

mydata4@meta.data$hpca.main <- hpca.main4$pruned.labels
mydata4@meta.data$hpca.fine <- hpca.fine4$pruned.labels


clusterLabels <-  hpca.fine4.clusters$pruned.labels

mydata4@meta.data$hpca.fine.clusters[mydata4$seurat_clusters == 0] <- hpca.fine4.clusters$pruned.labels[1]
mydata4@meta.data$hpca.fine.clusters[mydata4$seurat_clusters == 1] <- hpca.fine4.clusters$pruned.labels[2]
mydata4@meta.data$hpca.fine.clusters[mydata4$seurat_clusters == 2] <- hpca.fine4.clusters$pruned.labels[3]




mydata4 <- SetIdent(mydata4, value = "hpca.fine.clusters")
DimPlot(mydata4, label = T , repel = T, label.size = 3) + NoLegend()

### DIFFERENTIAL EXPRESSION TESTING

DCmonocyteVSiPS.de.markers <- FindMarkers(mydata4, ident.1 = "DC:monocyte-derived:rosiglitazone", ident.2 = "iPS_cells:skin_fibroblast")
# view results
head(astrocyteVStcell.de.markers)

library(ggplot2)
# The basic scatter plot: x is "log2FoldChange", y is "pvalue"
ggplot(data=astrocyteVStcell.de.markers, aes(x=avg_log2FC, y=p_val)) + geom_point()
p <- ggplot(data=astrocyteVStcell.de.markers, aes(x=avg_log2FC, y=-log10(p_val))) + geom_point()
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +  geom_hline(yintercept=-log10(0.05), col="red")

# add a column of expression labelling
DCmonocyteVSiPS.de.markers$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
DCmonocyteVSiPS.de.markers$diffexpressed[DCmonocyteVSiPS.de.markers$avg_log2FC > 0.6 & DCmonocyteVSiPS.de.markers$p_val < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
DCmonocyteVSiPS.de.markers$diffexpressed[DCmonocyteVSiPS.de.markers$avg_log2FC < -0.6 & DCmonocyteVSiPS.de.markers$p_val < 0.05] <- "DOWN"
#adding nemes as the rownames
p <- ggplot(data=DCmonocyteVSiPS.de.markers, aes(x=avg_log2FC, y=-log10(p_val), col=diffexpressed, label=rownames(DCmonocyteVSiPS.de.markers))) + geom_point() + theme_minimal() + geom_text()
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +  geom_hline(yintercept=-log10(0.05), col="red")


##also for clusters
mydata4 <- SetIdent(mydata4, value = "seurat_clusters")
cluster1vsCluster2.de.markers4 <- FindMarkers(mydata4, ident.1 = 1, ident.2 = 2)
cluster1vsCluster2.de.markers4["hgnc_symbol"] <- rownames(cluster1vsCluster2.de.markers4) 
cluster1vsCluster0.de.markers4 <- FindMarkers(mydata4, ident.1 = 1, ident.2 = 0)
cluster1vsCluster0.de.markers4["hgnc_symbol"] <- rownames(cluster1vsCluster0.de.markers4) 
degenescluster1 <- union(rownames(cluster1vsCluster2.de.markers4),rownames(cluster1vsCluster0.de.markers4))
cluster1Markers <- bind_rows(cluster1vsCluster0.de.markers4,cluster1vsCluster2.de.markers4)
cluster1Markers <- cluster1Markers %>% group_by(hgnc_symbol) %>%
  arrange(avg_log2FC) %>%
  slice(1) %>% ungroup


# add a column of expression labelling
cluster1Markers$diffexpressed <- "NO"
# if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
cluster1Markers$diffexpressed[cluster1Markers$avg_log2FC > 0.5 & cluster1Markers$p_val < 0.05] <- "UP"
# if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
cluster1Markers$diffexpressed[cluster1Markers$avg_log2FC < -0.5 & cluster1Markers$p_val < 0.05] <- "DOWN"


#cluster2
cluster2vsCluster1.de.markers4 <- FindMarkers(mydata4, ident.1 = 2, ident.2 = 1)
cluster2vsCluster1.de.markers4["hgnc_symbol"] <- rownames(cluster2vsCluster1.de.markers4)
cluster2vsCluster0.de.markers4 <- FindMarkers(mydata4, ident.1 = 2, ident.2 = 0)
cluster2vsCluster0.de.markers4["hgnc_symbol"] <- rownames(cluster2vsCluster0.de.markers4) 
degenescluster2 <- union(rownames(cluster2vsCluster1.de.markers4),rownames(cluster2vsCluster0.de.markers4))
cluster2Markers <- bind_rows(cluster2vsCluster0.de.markers4,cluster2vsCluster1.de.markers4)
cluster2Markers <- cluster2Markers %>% group_by(hgnc_symbol) %>%
  arrange(avg_log2FC) %>%
  slice(1) %>% ungroup


# add a column of expression labelling
cluster2Markers$diffexpressed <- "NO"
# if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
cluster2Markers$diffexpressed[cluster2Markers$avg_log2FC > 0.5 & cluster2Markers$p_val < 0.05] <- "UP"
# if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
cluster2Markers$diffexpressed[cluster2Markers$avg_log2FC < -0.5 & cluster2Markers$p_val < 0.05] <- "DOWN"


#cluster0
cluster0vsCluster1.de.markers4 <- FindMarkers(mydata4, ident.1 = 0, ident.2 = 1)
cluster0vsCluster1.de.markers4["hgnc_symbol"] <- rownames(cluster0vsCluster1.de.markers4) 
cluster0vsCluster2.de.markers4 <- FindMarkers(mydata4, ident.1 = 0, ident.2 = 2)
cluster0vsCluster2.de.markers4["hgnc_symbol"] <- rownames(cluster0vsCluster2.de.markers4) 
degenescluster0 <- union(rownames(cluster0vsCluster1.de.markers4),rownames(cluster0vsCluster2.de.markers4))
cluster0Markers <- bind_rows(cluster0vsCluster1.de.markers4,cluster0vsCluster2.de.markers4)
cluster0Markers <- cluster0Markers %>% group_by(hgnc_symbol) %>%
  arrange(avg_log2FC) %>%
  slice(1) %>% ungroup


# add a column of expression labelling
cluster0Markers$diffexpressed <- "NO"
# if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
cluster0Markers$diffexpressed[cluster0Markers$avg_log2FC > 0.5 & cluster0Markers$p_val < 0.05] <- "UP"
# if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
cluster0Markers$diffexpressed[cluster0Markers$avg_log2FC < -0.5 & cluster0Markers$p_val < 0.05] <- "DOWN"


## preparing the data to use phensim


#cluster0

cluster0.all.genes.filtered4 <- distinct(data.frame(all.genes4[!all.genes4 %in% cluster0Markers$hgnc_symbol]))

cluster0.all.genes.filtered4 <- cluster0.all.genes.filtered4$all.genes4..all.genes4..in..cluster0Markers.hgnc_symbol.

diffExpressedSim4 <- data.frame(rep("NO",length(cluster0.all.genes.filtered4[!cluster0.all.genes.filtered4 %in% cluster0Markers$hgnc_symbol])),row.names = cluster0.all.genes.filtered4[!cluster0.all.genes.filtered4 %in% cluster0Markers$hgnc_symbol])
colnames(diffExpressedSim4) <- "diffexpressed"

diffExpressedSim <- rbind(diffExpressedSim, select(astrocyteVStcell.de.markers,c("diffexpressed")))

diffExpressedSim["hgnc_symbol"] <- rownames(diffExpressedSim)
rownames(diffExpressedSim) <- NULL

# symplyfing the list of entrez_id genes as input for phensim since it only needs a list
allNonexpressed <- data.frame(rep("NO",length(allGlist4$hgnc_symbol[!(allGlist4$hgnc_symbol %in% entrez.list4$hgnc_symbol)])),allGlist4$hgnc_symbol[!(allGlist4$hgnc_symbol %in% entrez.list4$hgnc_symbol)])
diffNonExpressedSim4 <- data.frame(rep("NO",length(entrez.list4$hgnc_symbol[!(entrez.list4$hgnc_symbol %in% cluster0Markers$hgnc_symbol)])),entrez.list4$hgnc_symbol[!(entrez.list4$hgnc_symbol %in% cluster0Markers$hgnc_symbol)])
colnames(diffNonExpressedSim4) <- c("diffexpressed","hgnc_symbol")
colnames(allNonexpressed) <- c("diffexpressed","hgnc_symbol")

#non considered are seen as non-expressed
nonConsideredEntries.cluster0 <- inner_join(allNonexpressed,select(allGlist4,c("hgnc_symbol","entrezgene_id")),by = "hgnc_symbol")
nonExpressedEntries.cluster0 <- inner_join(diffNonExpressedSim4,select(allGlist4,c("hgnc_symbol","entrezgene_id")),by = "hgnc_symbol")
expressedEntries.cluster0 <- inner_join(select(cluster0Markers,c("diffexpressed","hgnc_symbol")),select(allGlist4,c("hgnc_symbol","entrezgene_id")))

full.phensim.input.cluster0 <- inner_join(diffNonExpressedSim4,select(allGlist4,c("hgnc_symbol","entrezgene_id")),by = "hgnc_symbol")

full.phensim.input.cluster0 <- full.phensim.input[!is.na(full.phensim.input$entrezgene_id),]

#merging all non expressed genes
nonExpr.all.cluster0 <- union(nonConsideredEntries.cluster0$entrezgene_id,nonExpressedEntries.cluster0$entrezgene_id)

nonExpr.all.cluster0 <- union(nonExpr.all.cluster0,expressedEntries.cluster0[expressedEntries.cluster0$diffexpressed=="NO",]$entrezgene_id)

overexpressed.cluster0 <- expressedEntries.cluster0[expressedEntries.cluster0$diffexpressed=="UP",]$entrezgene_id
downexpressed.cluster0 <- expressedEntries.cluster0[expressedEntries.cluster0$diffexpressed=="DOWN",]$entrezgene_id

write.table(expressedEntries.cluster0,
            file='/home/josura/Projects/tesi/data/patient4ALL/cluster0-diffExpressedMarkers.tsv', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)

write.table(nonConsideredEntries.cluster0,
            file='/home/josura/Projects/tesi/data/patient4ALL/cluster0-nonConsidered.tsv', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)

write.table(nonExpressedEntries.cluster0,
            file='/home/josura/Projects/tesi/data/patient4ALL/cluster0-notdiffExpressed.tsv', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)

#only txt with entrezgene_id for phensim
write.table(nonExpr.all.cluster0,
            file='/home/josura/Projects/tesi/data/patient4ALL/cluster0-nonExpressed.txt', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)

write.table(overexpressed.cluster0,
            file='/home/josura/Projects/tesi/data/patient4ALL/cluster0-overexpressed.txt', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)

write.table(downexpressed.cluster0,
            file='/home/josura/Projects/tesi/data/patient4ALL/cluster0-downexpressed.txt', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)




#cluster1

cluster1.all.genes.filtered4 <- distinct(data.frame(all.genes4[!all.genes4 %in% cluster1Markers$hgnc_symbol]))

cluster1.all.genes.filtered4 <- cluster1.all.genes.filtered4$all.genes4..all.genes4..in..cluster1Markers.hgnc_symbol.

# symplyfing the list of entrez_id genes as input for phensim since it only needs a list
allNonexpressed <- data.frame(rep("NO",length(allGlist4$hgnc_symbol[!(allGlist4$hgnc_symbol %in% entrez.list4$hgnc_symbol)])),allGlist4$hgnc_symbol[!(allGlist4$hgnc_symbol %in% entrez.list4$hgnc_symbol)])
diffNonExpressedSim4 <- data.frame(rep("NO",length(entrez.list4$hgnc_symbol[!(entrez.list4$hgnc_symbol %in% cluster1Markers$hgnc_symbol)])),entrez.list4$hgnc_symbol[!(entrez.list4$hgnc_symbol %in% cluster1Markers$hgnc_symbol)])
colnames(diffNonExpressedSim4) <- c("diffexpressed","hgnc_symbol")
colnames(allNonexpressed) <- c("diffexpressed","hgnc_symbol")

#non considered are seen as non-expressed
nonConsideredEntries.cluster1 <- inner_join(allNonexpressed,select(allGlist4,c("hgnc_symbol","entrezgene_id")),by = "hgnc_symbol")
nonExpressedEntries.cluster1 <- inner_join(diffNonExpressedSim4,select(allGlist4,c("hgnc_symbol","entrezgene_id")),by = "hgnc_symbol")
expressedEntries.cluster1 <- inner_join(select(cluster1Markers,c("diffexpressed","hgnc_symbol")),select(allGlist4,c("hgnc_symbol","entrezgene_id")))

full.phensim.input.cluster1 <- inner_join(diffNonExpressedSim4,select(allGlist4,c("hgnc_symbol","entrezgene_id")),by = "hgnc_symbol")

full.phensim.input.cluster1 <- full.phensim.input.cluster1[!is.na(full.phensim.input.cluster1$entrezgene_id),]

#merging all non expressed genes
nonExpr.all.cluster1 <- union(nonConsideredEntries.cluster1$entrezgene_id,nonExpressedEntries.cluster1$entrezgene_id)

nonExpr.all.cluster1 <- union(nonExpr.all.cluster1,expressedEntries.cluster1[expressedEntries.cluster1$diffexpressed=="NO",]$entrezgene_id)

overexpressed.cluster1 <- expressedEntries.cluster1[expressedEntries.cluster1$diffexpressed=="UP",]$entrezgene_id
downexpressed.cluster1 <- expressedEntries.cluster1[expressedEntries.cluster1$diffexpressed=="DOWN",]$entrezgene_id

write.table(expressedEntries.cluster1,
            file='/home/josura/Projects/tesi/data/patient4ALL/cluster1-diffExpressedMarkers.tsv', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)

write.table(nonConsideredEntries.cluster1,
            file='/home/josura/Projects/tesi/data/patient4ALL/cluster1-nonConsidered.tsv', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)

write.table(nonExpressedEntries.cluster1,
            file='/home/josura/Projects/tesi/data/patient4ALL/cluster1-notdiffExpressed.tsv', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)

#only txt with entrezgene_id for phensim
write.table(nonExpr.all.cluster1,
            file='/home/josura/Projects/tesi/data/patient4ALL/cluster1-nonExpressed.txt', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)

write.table(overexpressed.cluster1,
            file='/home/josura/Projects/tesi/data/patient4ALL/cluster1-overexpressed.txt', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)

write.table(downexpressed.cluster1,
            file='/home/josura/Projects/tesi/data/patient4ALL/cluster1-downexpressed.txt', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)




#cluster2

cluster2.all.genes.filtered4 <- distinct(data.frame(all.genes4[!all.genes4 %in% cluster2Markers$hgnc_symbol]))

cluster2.all.genes.filtered4 <- cluster2.all.genes.filtered4$all.genes4..all.genes4..in..cluster2Markers.hgnc_symbol.

# symplyfing the list of entrez_id genes as input for phensim since it only needs a list
allNonexpressed <- data.frame(rep("NO",length(allGlist4$hgnc_symbol[!(allGlist4$hgnc_symbol %in% entrez.list4$hgnc_symbol)])),allGlist4$hgnc_symbol[!(allGlist4$hgnc_symbol %in% entrez.list4$hgnc_symbol)])
diffNonExpressedSim4 <- data.frame(rep("NO",length(entrez.list4$hgnc_symbol[!(entrez.list4$hgnc_symbol %in% cluster2Markers$hgnc_symbol)])),entrez.list4$hgnc_symbol[!(entrez.list4$hgnc_symbol %in% cluster2Markers$hgnc_symbol)])
colnames(diffNonExpressedSim4) <- c("diffexpressed","hgnc_symbol")
colnames(allNonexpressed) <- c("diffexpressed","hgnc_symbol")

#non considered are seen as non-expressed
nonConsideredEntries.cluster2 <- inner_join(allNonexpressed,select(allGlist4,c("hgnc_symbol","entrezgene_id")),by = "hgnc_symbol")
nonExpressedEntries.cluster2 <- inner_join(diffNonExpressedSim4,select(allGlist4,c("hgnc_symbol","entrezgene_id")),by = "hgnc_symbol")
expressedEntries.cluster2 <- inner_join(select(cluster2Markers,c("diffexpressed","hgnc_symbol")),select(allGlist4,c("hgnc_symbol","entrezgene_id")))

full.phensim.input.cluster2 <- inner_join(diffNonExpressedSim4,select(allGlist4,c("hgnc_symbol","entrezgene_id")),by = "hgnc_symbol")

full.phensim.input.cluster2 <- full.phensim.input.cluster2[!is.na(full.phensim.input.cluster2$entrezgene_id),]

#merging all non expressed genes
nonExpr.all.cluster2 <- union(nonConsideredEntries.cluster2$entrezgene_id,nonExpressedEntries.cluster2$entrezgene_id)

nonExpr.all.cluster2 <- union(nonExpr.all.cluster2,expressedEntries.cluster2[expressedEntries.cluster2$diffexpressed=="NO",]$entrezgene_id)

overexpressed.cluster2 <- expressedEntries.cluster2[expressedEntries.cluster2$diffexpressed=="UP",]$entrezgene_id
downexpressed.cluster2 <- expressedEntries.cluster2[expressedEntries.cluster2$diffexpressed=="DOWN",]$entrezgene_id

write.table(expressedEntries.cluster2,
            file='/home/josura/Projects/tesi/data/patient4ALL/cluster2-diffExpressedMarkers.tsv', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)

write.table(nonConsideredEntries.cluster2,
            file='/home/josura/Projects/tesi/data/patient4ALL/cluster2-nonConsidered.tsv', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)

write.table(nonExpressedEntries.cluster2,
            file='/home/josura/Projects/tesi/data/patient4ALL/cluster2-notdiffExpressed.tsv', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)

#only txt with entrezgene_id for phensim
write.table(nonExpr.all.cluster2,
            file='/home/josura/Projects/tesi/data/patient4ALL/cluster2-nonExpressed.txt', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)

write.table(overexpressed.cluster2,
            file='/home/josura/Projects/tesi/data/patient4ALL/cluster2-overexpressed.txt', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)

write.table(downexpressed.cluster2,
            file='/home/josura/Projects/tesi/data/patient4ALL/cluster2-downexpressed.txt', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)



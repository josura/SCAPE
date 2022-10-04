library(ggplot2)
library(Seurat)

metadata.scrnaseq.calu3 <- read.csv("/home/josura/Projects/tesi/data/tesi/modernData/COVID/Calu3_scRNAseq_metadata_full.txt")
### visualization by strain
ggplot() +
  ggtitle("Visualization for mock data at different times") +
  geom_point(data=metadata.scrnaseq.calu3[metadata.scrnaseq.calu3$strain=="nan",], aes(x=UMAP_1, y=UMAP_2, color=time),size=0.5)
ggplot() +
  ggtitle("Visualization for SarsCov1 data at different times") +
  geom_point(data=metadata.scrnaseq.calu3[metadata.scrnaseq.calu3$strain=="SARSCoV1",], aes(x=UMAP_1, y=UMAP_2, color=time),size=0.5)
ggplot() +
  ggtitle("Visualization for SarsCov2 data at different times") +
  geom_point(data=metadata.scrnaseq.calu3[metadata.scrnaseq.calu3$strain=="SARSCoV2",], aes(x=UMAP_1, y=UMAP_2, color=time),size=0.5)

### visualization by time
ggplot() +
  ggtitle("Visualization for 4h data") +
  geom_point(data=metadata.scrnaseq.calu3[metadata.scrnaseq.calu3$time=="4h",], aes(x=UMAP_1, y=UMAP_2, color=strain),size=0.5)
ggplot() +
  ggtitle("Visualization for 8h data") +
  geom_point(data=metadata.scrnaseq.calu3[metadata.scrnaseq.calu3$time=="8h",], aes(x=UMAP_1, y=UMAP_2, color=strain),size=0.5)
ggplot() +
  ggtitle("Visualization for 12h data") +
  geom_point(data=metadata.scrnaseq.calu3[metadata.scrnaseq.calu3$time=="12h",], aes(x=UMAP_1, y=UMAP_2, color=strain),size=0.5)

#visualization by origin
ggplot() +
  ggtitle("Visualization for mock data at different times") +
  geom_point(data=metadata.scrnaseq.calu3, aes(x=UMAP_1, y=UMAP_2, color=orig.ident),size=0.5)



###additional metadata
metadataOth.scrnaseq.calu3 <- read.csv("/home/josura/Projects/tesi/data/tesi/modernData/COVID/Calu3_Table_fullmetadata.csv")

library(ggnewscale)

p <- ggplot() +
  ggtitle("SCov2 load at different times") +
  geom_point(data=metadataOth.scrnaseq.calu3[metadataOth.scrnaseq.calu3$type=="S2" & metadataOth.scrnaseq.calu3$time=="4h",], aes(x=UMAPrna_1, y=UMAPrna_2, color=SCov2_Load)) +
  scale_color_gradient(low="red", high="gray50") +
  labs(colour="4h") +
  new_scale_color() +
  geom_point(data=metadataOth.scrnaseq.calu3[metadataOth.scrnaseq.calu3$type=="S2" & metadataOth.scrnaseq.calu3$time=="8h",], aes(x=UMAPrna_1, y=UMAPrna_2, color=SCov2_Load)) +
  scale_color_gradient(low="gray90", high="blue") +
  labs(colour="8h") +
  new_scale_color() +
  geom_point(data=metadataOth.scrnaseq.calu3[metadataOth.scrnaseq.calu3$type=="S2" & metadataOth.scrnaseq.calu3$time=="12h",], aes(x=UMAPrna_1, y=UMAPrna_2, color=SCov2_Load)) +
  scale_color_gradient(low="black", high="yellow") +
  labs(colour="12h")


p <- ggplot() +
  ggtitle("SCov1 load at different times") +
  geom_point(data=metadataOth.scrnaseq.calu3[metadataOth.scrnaseq.calu3$type=="S1" & metadataOth.scrnaseq.calu3$time=="4h",], aes(x=UMAPrna_1, y=UMAPrna_2, color=SCov1_Load)) +
  scale_color_gradient(low="red", high="gray50") +
  labs(colour="4h") +
  new_scale_color() +
  geom_point(data=metadataOth.scrnaseq.calu3[metadataOth.scrnaseq.calu3$type=="S1" & metadataOth.scrnaseq.calu3$time=="8h",], aes(x=UMAPrna_1, y=UMAPrna_2, color=SCov1_Load)) +
  scale_color_gradient(low="gray90", high="blue") +
  labs(colour="8h") +
  new_scale_color() +
  geom_point(data=metadataOth.scrnaseq.calu3[metadataOth.scrnaseq.calu3$type=="S1" & metadataOth.scrnaseq.calu3$time=="12h",], aes(x=UMAPrna_1, y=UMAPrna_2, color=SCov1_Load)) +
  scale_color_gradient(low="black", high="yellow") +
  labs(colour="12h")


### seurat visualization
seurat.all.calu3 <- readRDS("/home/josura/Projects/tesi/data/sarscov-seurat-all.rds")

seurat.all.calu3@meta.data$orig.ident <- metadata.scrnaseq.calu3$orig.ident

# normalize, always gc after every evaluation if memory is "low"
seurat.all.calu3 = NormalizeData(seurat.all.calu3)
# scale
seurat.all.calu3 = ScaleData(seurat.all.calu3)
# Find variable features
seurat.all.calu3 = FindVariableFeatures(seurat.all.calu3)
# Run PCA
seurat.all.calu3 = RunPCA(seurat.all.calu3, verbose=FALSE)
# Run UMAP on selected dims
seurat.all.calu3 = RunUMAP(seurat.all.calu3, dims = 1:30)

SetIdent(seurat.all.calu3,value="orig.ident")
graph.all <- DimPlot(seurat.all.calu3, reduction = "umap", label = TRUE, repel = TRUE)


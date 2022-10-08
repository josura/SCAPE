library(dplyr)

phensim_output.sarscov2.4h <- read.csv("/home/josura/Projects/tesi/data/modernData/COVID/sarscov2-4h/2640-sarscov2-4h-calu3-output.tsv",sep="\t")
phensim_output.sarscov2.8h <- read.csv("/home/josura/Projects/tesi/data/modernData/COVID/sarscov2-8h/2641-sarscov2-8h-calu3-output.tsv",sep="\t")
phensim_output.sarscov2.12h <- read.csv("/home/josura/Projects/tesi/data/modernData/COVID/sarscov2-12h/2642-sarscov2-12h-calu3-output.tsv",sep="\t")

phensim_output.sarscov2.4h$time <- "4h"
phensim_output.sarscov2.8h$time <- "8h"
phensim_output.sarscov2.12h$time <- "12h"


phensim_output.sarscov2.4h$strain <- "sarscov2"
phensim_output.sarscov2.8h$strain <- "sarscov2"
phensim_output.sarscov2.12h$strain <- "sarscov2"


phensim_output.sarscov1.4h <- read.csv("/home/josura/Projects/tesi/data/modernData/COVID/sarscov1-4h/2748-sarscov1-4h-calu3-output.tsv",sep="\t")
phensim_output.sarscov1.8h <- read.csv("/home/josura/Projects/tesi/data/modernData/COVID/sarscov1-8h/2749-sarscov1-8h-calu3-output.tsv",sep="\t")
phensim_output.sarscov1.12h <- read.csv("/home/josura/Projects/tesi/data/modernData/COVID/sarscov1-12h/2750-sarscov1-12h-calu3-output.tsv",sep="\t")

phensim_output.sarscov1.4h$strain <- "sarscov1"
phensim_output.sarscov1.8h$strain <- "sarscov1"
phensim_output.sarscov1.12h$strain <- "sarscov1"

phensim_output.sarscov2.4h.endpoints <- phensim_output.sarscov2.4h %>% group_by(X..Pathway.Id) %>%
  arrange(Is.Endpoint) %>%
  slice(1) %>% ungroup

phensim_output.sarscov2.4h.endpoints.ordered <- phensim_output.sarscov2.4h.endpoints %>%
  arrange(desc(Pathway.Activity.Score))


phensim_output.sarscov2.4h.endpoints.significant.ordered <- phensim_output.sarscov2.4h.endpoints %>%
  filter(Pathway.Adjusted.p.value < 0.05) %>% 
  arrange(Pathway.Activity.Score)


full.data <- bind_rows(phensim_output.sarscov1.4h, phensim_output.sarscov1.8h, phensim_output.sarscov1.12h,
                       phensim_output.sarscov2.4h, phensim_output.sarscov2.8h, phensim_output.sarscov2.12h)


### REACTOME


phensim_output.sarscov1.4h.reactome <- read.csv("/home/josura/Projects/tesi/data/modernData/COVID/sarscov1-4h/3591-sarscov1-4h-calu3-reactome-output.tsv",sep="\t")
phensim_output.sarscov1.8h.reactome <- read.csv("/home/josura/Projects/tesi/data/modernData/COVID/sarscov1-8h/3592-sarscov1-8h-calu3-reactome-output.tsv",sep="\t")
phensim_output.sarscov1.12h.reactome <- read.csv("/home/josura/Projects/tesi/data/modernData/COVID/sarscov1-12h/3593-sarscov1-12h-calu3-reactome-output.tsv",sep="\t")

phensim_output.sarscov2.4h.reactome <- read.csv("/home/josura/Projects/tesi/data/modernData/COVID/sarscov2-4h/3594-sarscov2-4h-calu3-reactome-output.tsv",sep="\t")
phensim_output.sarscov2.8h.reactome <- read.csv("/home/josura/Projects/tesi/data/modernData/COVID/sarscov2-8h/3595-sarscov2-8h-calu3-reactome-output.tsv",sep="\t")
phensim_output.sarscov2.12h.reactome <- read.csv("/home/josura/Projects/tesi/data/modernData/COVID/sarscov2-12h/3596-sarscov2-12h-calu3-reactome-output.tsv",sep="\t")

phensim_output.sarscov2.4h.reactome$time <- "4h"
phensim_output.sarscov2.8h.reactome$time <- "8h"
phensim_output.sarscov2.12h.reactome$time <- "12h"


phensim_output.sarscov2.4h.reactome$strain <- "sarscov2"
phensim_output.sarscov2.8h.reactome$strain <- "sarscov2"
phensim_output.sarscov2.12h.reactome$strain <- "sarscov2"


phensim_output.sarscov1.4h.reactome$time <- "4h"
phensim_output.sarscov1.8h.reactome$time <- "8h"
phensim_output.sarscov1.12h.reactome$time <- "12h"

phensim_output.sarscov1.4h.reactome$strain <- "sarscov1"
phensim_output.sarscov1.8h.reactome$strain <- "sarscov1"
phensim_output.sarscov1.12h.reactome$strain <- "sarscov1"

### interferons pathways

phensim_output.sarscov1.4h.reactome.interferon <- phensim_output.sarscov1.4h.reactome[grepl("interferon", phensim_output.sarscov1.4h.reactome$Pathway.Name,ignore.case = TRUE) | grepl("ifn", phensim_output.sarscov1.4h.reactome$Pathway.Name,ignore.case = TRUE) ,]

interferons.list <- unique(phensim_output.sarscov1.4h.reactome.interferon$X..Pathway.Id)
interferons.list <- gsub("-",".",interferons.list)

create.phensim.pathways.embeddings <- function(phensim_output){
  #TODO control input
  pathway.list <- unique(phensim_output$X..Pathway.Id)
  pathway.embedding.list <- list()
  for (pathway.name in pathway.list) {
    sub.phensim_output <- phensim_output %>%
      filter(X..Pathway.Id == pathway.name)
    phensim.embed <- data.frame(t(unlist(sub.phensim_output$Average.Node.Perturbation)),row.names = pathway.name)
    colnames(phensim.embed) <- sub.phensim_output$Node.Name
    #phensim.embed$gene.name <- pathway.name
    pathway.embedding.list[[length(pathway.embedding.list)+1]] <- phensim.embed
    print(paste("embed done for pathway",pathway.name,sep = ": "))
  }
  merged.embeddings <- pathway.embedding.list[[1]]
  pathway.embedding.list <- pathway.embedding.list[-1]
  for (path.embed in pathway.embedding.list) {
    merged.embeddings <- bind_rows(merged.embeddings,path.embed)
    print(paste("merge done for pathway",rownames(path.embed),sep = ": "))
  }
  #merged.embeddings[is.na(merged.embeddings)] <- 0
  merged.embeddings
}
## ALL THE PATHWAYS
phensim_output.sarscov1.4h.reactome.pathwayembed <-  create.phensim.pathways.embeddings(phensim_output.sarscov1.4h.reactome)
rownames(phensim_output.sarscov1.4h.reactome.pathwayembed) <- paste(rownames(phensim_output.sarscov1.4h.reactome.pathwayembed), "-4h-sarscov1", sep="")

phensim_output.sarscov1.8h.reactome.pathwayembed <-  create.phensim.pathways.embeddings(phensim_output.sarscov1.8h.reactome)
rownames(phensim_output.sarscov1.8h.reactome.pathwayembed) <- paste(rownames(phensim_output.sarscov1.8h.reactome.pathwayembed), "-8h-sarscov1", sep="")

phensim_output.sarscov1.12h.reactome.pathwayembed <- create.phensim.pathways.embeddings(phensim_output.sarscov1.12h.reactome)
rownames(phensim_output.sarscov1.12h.reactome.pathwayembed) <- paste(rownames(phensim_output.sarscov1.12h.reactome.pathwayembed), "-12h-sarscov1", sep="")

write.csv(phensim_output.sarscov1.4h.reactome.pathwayembed,file="/home/josura/Projects/tesi/data/modernData/COVID/sarscov1-4h-PHENSIM-embeddings.csv")
write.csv(phensim_output.sarscov1.8h.reactome.pathwayembed,file="/home/josura/Projects/tesi/data/modernData/COVID/sarscov1-8h-PHENSIM-embeddings.csv")
write.csv(phensim_output.sarscov1.12h.reactome.pathwayembed,file="/home/josura/Projects/tesi/data/modernData/COVID/sarscov1-12h-PHENSIM-embeddings.csv")


phensim_output.sarscov2.4h.reactome.pathwayembed <-  create.phensim.pathways.embeddings(phensim_output.sarscov2.4h.reactome)
rownames(phensim_output.sarscov2.4h.reactome.pathwayembed) <- paste(rownames(phensim_output.sarscov2.4h.reactome.pathwayembed), "-4h-sarscov2", sep="")

phensim_output.sarscov2.8h.reactome.pathwayembed <-  create.phensim.pathways.embeddings(phensim_output.sarscov2.8h.reactome)
rownames(phensim_output.sarscov2.8h.reactome.pathwayembed) <- paste(rownames(phensim_output.sarscov2.8h.reactome.pathwayembed), "-8h-sarscov2", sep="")

phensim_output.sarscov2.12h.reactome.pathwayembed <- create.phensim.pathways.embeddings(phensim_output.sarscov2.12h.reactome)
rownames(phensim_output.sarscov2.12h.reactome.pathwayembed) <- paste(rownames(phensim_output.sarscov2.12h.reactome.pathwayembed), "-12h-sarscov2", sep="")

## ONLY INTERFERONS
phensim_output.sarscov1.4h.reactome.pathwayembed.interferons <-  create.phensim.pathways.embeddings(phensim_output.sarscov1.4h.reactome[grepl("interferon", phensim_output.sarscov1.4h.reactome$Pathway.Name,ignore.case = TRUE) | grepl("ifn", phensim_output.sarscov1.4h.reactome$Pathway.Name,ignore.case = TRUE) ,])
rownames(phensim_output.sarscov1.4h.reactome.pathwayembed.interferons) <- paste(interferons.list, "-4h-sarscov1", sep="")

phensim_output.sarscov1.8h.reactome.pathwayembed.interferons <-  create.phensim.pathways.embeddings(phensim_output.sarscov1.8h.reactome[grepl("interferon", phensim_output.sarscov1.8h.reactome$Pathway.Name,ignore.case = TRUE) | grepl("ifn", phensim_output.sarscov1.8h.reactome$Pathway.Name,ignore.case = TRUE) ,])
rownames(phensim_output.sarscov1.8h.reactome.pathwayembed.interferons) <- paste(interferons.list, "-8h-sarscov1", sep="")

phensim_output.sarscov1.12h.reactome.pathwayembed.interferons <- create.phensim.pathways.embeddings(phensim_output.sarscov1.12h.reactome[grepl("interferon", phensim_output.sarscov1.12h.reactome$Pathway.Name,ignore.case = TRUE) | grepl("ifn", phensim_output.sarscov1.12h.reactome$Pathway.Name,ignore.case = TRUE) ,])
rownames(phensim_output.sarscov1.12h.reactome.pathwayembed.interferons) <- paste(interferons.list, "-12h-sarscov1", sep="")


phensim_output.sarscov2.4h.reactome.pathwayembed.interferons <-  create.phensim.pathways.embeddings(phensim_output.sarscov2.4h.reactome[grepl("interferon", phensim_output.sarscov2.4h.reactome$Pathway.Name,ignore.case = TRUE) | grepl("ifn", phensim_output.sarscov2.4h.reactome$Pathway.Name,ignore.case = TRUE) ,])
rownames(phensim_output.sarscov2.4h.reactome.pathwayembed.interferons) <- paste(interferons.list, "-4h-sarscov2", sep="")

phensim_output.sarscov2.8h.reactome.pathwayembed.interferons <-  create.phensim.pathways.embeddings(phensim_output.sarscov2.8h.reactome[grepl("interferon", phensim_output.sarscov2.8h.reactome$Pathway.Name,ignore.case = TRUE) | grepl("ifn", phensim_output.sarscov2.8h.reactome$Pathway.Name,ignore.case = TRUE) ,])
rownames(phensim_output.sarscov2.8h.reactome.pathwayembed.interferons) <- paste(interferons.list, "-8h-sarscov2", sep="")

phensim_output.sarscov2.12h.reactome.pathwayembed.interferons <- create.phensim.pathways.embeddings(phensim_output.sarscov2.12h.reactome[grepl("interferon", phensim_output.sarscov2.12h.reactome$Pathway.Name,ignore.case = TRUE) | grepl("ifn", phensim_output.sarscov2.12h.reactome$Pathway.Name,ignore.case = TRUE) ,])
rownames(phensim_output.sarscov2.12h.reactome.pathwayembed.interferons) <- paste(interferons.list, "-12h-sarscov2", sep="")

full.data.matrix <- bind_rows(phensim_output.sarscov1.4h.reactome.pathwayembed.interferons, phensim_output.sarscov1.8h.reactome.pathwayembed.interferons, phensim_output.sarscov1.12h.reactome.pathwayembed.interferons,
                              phensim_output.sarscov2.4h.reactome.pathwayembed.interferons, phensim_output.sarscov2.8h.reactome.pathwayembed.interferons, phensim_output.sarscov2.12h.reactome.pathwayembed.interferons)
full.data.metadata <- data.frame(rownames(full.data.matrix),row.names = rownames(full.data.matrix))
colnames(full.data.metadata) <- "names"

full.data.metadata$time <- c(rep("4h",times = length(interferons.list)),rep("8h",times = length(interferons.list)),rep("12h",times = length(interferons.list)),
                             rep("4h",times = length(interferons.list)),rep("8h",times = length(interferons.list)),rep("12h",times = length(interferons.list)))

full.data.metadata$strain <- c(rep("sarscov1",times = length(interferons.list)*3),
                               rep("sarscov2",times = length(interferons.list)*3))
full.data.matrix.nonzerocols <- full.data.matrix %>%
  #dplyr::select(where(~ any(. != 0 )))
  select_if(~!all(is.na(.) | . == 0))
write.csv(full.data.matrix.nonzerocols,file="/home/josura/Projects/tesi/data/modernData/COVID/calu3-sarscov-interferons-PHENSIM-avg_perturbation_genes-filtered.csv")
write.csv(full.data.matrix,file="/home/josura/Projects/tesi/data/modernData/COVID/calu3-sarscov-interferons-PHENSIM-avg_perturbation_genes.csv")
write.csv(full.data.metadata,file="/home/josura/Projects/tesi/data/modernData/COVID/interferons-PHENSIM-embeddings-metadata.csv")

### NULL models
phensim_output.sarscov1.4h.reactome.matrix <- t(read.csv("/home/josura/Projects/tesi/data/modernData/COVID/sarscov1-4h/3591-sarscov1-4h-calu3-reactome-pathway-matrix.tsv",sep="\t"))
phensim_output.sarscov1.8h.reactome.matrix <- t(read.csv("/home/josura/Projects/tesi/data/modernData/COVID/sarscov1-8h/3592-sarscov1-8h-calu3-reactome-pathway-matrix.tsv",sep="\t"))
phensim_output.sarscov1.12h.reactome.matrix <- t(read.csv("/home/josura/Projects/tesi/data/modernData/COVID/sarscov1-12h/3593-sarscov1-12h-calu3-reactome-pathway-matrix.tsv",sep="\t"))

phensim_output.sarscov2.4h.reactome.matrix <- t(read.csv("/home/josura/Projects/tesi/data/modernData/COVID/sarscov2-4h/3594-sarscov2-4h-calu3-reactome-pathway-matrix.tsv",sep="\t"))
phensim_output.sarscov2.8h.reactome.matrix <- t(read.csv("/home/josura/Projects/tesi/data/modernData/COVID/sarscov2-8h/3595-sarscov2-8h-calu3-reactome-pathway-matrix.tsv",sep="\t"))
phensim_output.sarscov2.12h.reactome.matrix <- t(read.csv("/home/josura/Projects/tesi/data/modernData/COVID/sarscov2-12h/3596-sarscov2-12h-calu3-reactome-pathway-matrix.tsv",sep="\t"))

phensim_output.sarscov2.12h.reactome.gene.matrix <- read.csv("/home/josura/Projects/tesi/data/modernData/COVID/sarscov2-12h/3596-sarscov2-12h-calu3-reactome-nodes-matrix.tsv",sep="\t")

phensim_output.sarscov1.4h.reactome.interferons.matrix <- phensim_output.sarscov1.4h.reactome.matrix[interferons.list,]
phensim_output.sarscov1.4h.reactome.interferons.matrix[is.na(phensim_output.sarscov1.4h.reactome.interferons.matrix)] <- 0
sarscov1.pathways.names.4h <- rownames(phensim_output.sarscov1.4h.reactome.matrix)
phensim_output.sarscov1.4h.reactome.interferons.matrix <- data.frame(phensim_output.sarscov1.4h.reactome.interferons.matrix)
rownames(phensim_output.sarscov1.4h.reactome.interferons.matrix) <- paste(interferons.list, "-4h-sarscov1", sep="")

phensim_output.sarscov1.8h.reactome.interferons.matrix <- phensim_output.sarscov1.8h.reactome.matrix[interferons.list,]
phensim_output.sarscov1.8h.reactome.interferons.matrix[is.na(phensim_output.sarscov1.8h.reactome.interferons.matrix)] <- 0
sarscov1.pathways.names.8h <- rownames(phensim_output.sarscov1.8h.reactome.matrix)
phensim_output.sarscov1.8h.reactome.interferons.matrix <- data.frame(phensim_output.sarscov1.8h.reactome.interferons.matrix)
rownames(phensim_output.sarscov1.8h.reactome.interferons.matrix) <- paste(interferons.list, "-8h-sarscov1", sep="")

phensim_output.sarscov1.12h.reactome.interferons.matrix <- phensim_output.sarscov1.12h.reactome.matrix[interferons.list,]
phensim_output.sarscov1.12h.reactome.interferons.matrix[is.na(phensim_output.sarscov1.12h.reactome.interferons.matrix)] <- 0
sarscov1.pathways.names.12h <- rownames(phensim_output.sarscov1.12h.reactome.matrix)
phensim_output.sarscov1.12h.reactome.interferons.matrix <- data.frame(phensim_output.sarscov1.12h.reactome.interferons.matrix)
rownames(phensim_output.sarscov1.12h.reactome.interferons.matrix) <- paste(interferons.list, "-12h-sarscov1", sep="")

phensim_output.sarscov2.4h.reactome.interferons.matrix <- phensim_output.sarscov2.4h.reactome.matrix[interferons.list,]
phensim_output.sarscov2.4h.reactome.interferons.matrix[is.na(phensim_output.sarscov2.4h.reactome.interferons.matrix)] <- 0
sarscov2.pathways.names.4h <- rownames(phensim_output.sarscov2.4h.reactome.matrix)
phensim_output.sarscov2.4h.reactome.interferons.matrix <- data.frame(phensim_output.sarscov2.4h.reactome.interferons.matrix)
rownames(phensim_output.sarscov2.4h.reactome.interferons.matrix) <- paste(interferons.list, "-4h-sarscov2", sep="")

phensim_output.sarscov2.8h.reactome.interferons.matrix <- phensim_output.sarscov2.8h.reactome.matrix[interferons.list,]
phensim_output.sarscov2.8h.reactome.interferons.matrix[is.na(phensim_output.sarscov2.8h.reactome.interferons.matrix)] <- 0
sarscov2.pathways.names.8h <- rownames(phensim_output.sarscov2.8h.reactome.matrix)
phensim_output.sarscov2.8h.reactome.interferons.matrix <- data.frame(phensim_output.sarscov2.8h.reactome.interferons.matrix)
rownames(phensim_output.sarscov2.8h.reactome.interferons.matrix) <- paste(interferons.list, "-8h-sarscov2", sep="")

phensim_output.sarscov2.12h.reactome.interferons.matrix <- phensim_output.sarscov2.12h.reactome.matrix[interferons.list,]
phensim_output.sarscov2.12h.reactome.interferons.matrix[is.na(phensim_output.sarscov2.12h.reactome.interferons.matrix)] <- 0
sarscov2.pathways.names.12h <- rownames(phensim_output.sarscov2.12h.reactome.matrix)
phensim_output.sarscov2.12h.reactome.interferons.matrix <- data.frame(phensim_output.sarscov2.12h.reactome.interferons.matrix)
rownames(phensim_output.sarscov2.12h.reactome.interferons.matrix) <- paste(interferons.list, "-12h-sarscov2", sep="")


full.data.matrix <- bind_rows(phensim_output.sarscov1.4h.reactome.interferons.matrix, phensim_output.sarscov1.8h.reactome.interferons.matrix, phensim_output.sarscov1.12h.reactome.interferons.matrix,
                              phensim_output.sarscov2.4h.reactome.interferons.matrix, phensim_output.sarscov2.8h.reactome.interferons.matrix, phensim_output.sarscov2.12h.reactome.interferons.matrix)

full.data.metadata <- data.frame(rownames(full.data.matrix),row.names = rownames(full.data.matrix))
colnames(full.data.metadata) <- "names"

full.data.metadata$time <- c(rep("4h",times = length(interferons.list)),rep("8h",times = length(interferons.list)),rep("12h",times = length(interferons.list)),
  rep("4h",times = length(interferons.list)),rep("8h",times = length(interferons.list)),rep("12h",times = length(interferons.list)))

full.data.metadata$strain <- c(rep("sarscov1",times = length(interferons.list)*3),
                             rep("sarscov2",times = length(interferons.list)*3))

library(umap)

custom.config <- umap.defaults
custom.config$random_state <- 123
full.data.matrix.distances <- as.matrix(dist(full.data.matrix))
full.data.matrix.umap.model.dist <- umap(full.data.matrix.distances, config=custom.config, input="dist")
#full.data.matrix.umap <- umap(full.data.matrix)

full.data.metadata$UMAPpathway_1 <-  full.data.matrix.umap.model.dist$layout[,1]
full.data.metadata$UMAPpathway_2 <-  full.data.matrix.umap.model.dist$layout[,2]

### VISUALIZATION

library(ggplot2)
library(ggnewscale)

#see UMAP values

p <- ggplot() +
  ggtitle("SCov strain at different times for PHENSIM output (perturbation for interferones genes)") +
  geom_point(data=full.data.metadata[ full.data.metadata$time=="4h",], aes(x=UMAPpathway_1, y=UMAPpathway_2, color=strain)) +
  scale_colour_manual(values = c("red", "blue")) +
  labs(colour="4h") +
  new_scale_color() +
  geom_point(data=full.data.metadata[ full.data.metadata$time=="8h",], aes(x=UMAPpathway_1, y=UMAPpathway_2, color=strain)) +
  scale_colour_manual(values = c("violet", "black")) +
  labs(colour="8h") +
  new_scale_color() +
  geom_point(data=full.data.metadata[ full.data.metadata$time=="12h",], aes(x=UMAPpathway_1, y=UMAPpathway_2, color=strain)) +
  scale_colour_manual(values = c("orange", "gray60")) +
  labs(colour="12h") 




### UMAP values for the whole samples 
sarscov1.4h.all.genes.PHENSIM <- data.frame(t(unlist(phensim_output.sarscov1.4h.reactome$Average.Node.Perturbation)),row.names = "sarscov1-4h")
colnames(sarscov1.4h.all.genes.PHENSIM) <- phensim_output.sarscov1.4h.reactome$X..Pathway.Id


sarscov1.8h.all.genes.PHENSIM <- data.frame(t(unlist(phensim_output.sarscov1.8h.reactome$Average.Node.Perturbation)),row.names = "sarscov1-8h")
colnames(sarscov1.4h.all.genes.PHENSIM) <- phensim_output.sarscov1.8h.reactome$X..Pathway.Id


sarscov1.12h.all.genes.PHENSIM <- data.frame(t(unlist(phensim_output.sarscov1.12h.reactome$Average.Node.Perturbation)),row.names = "sarscov1-12h")
colnames(sarscov1.4h.all.genes.PHENSIM) <- phensim_output.sarscov1.12h.reactome$X..Pathway.Id


sarscov2.4h.all.genes.PHENSIM <- data.frame(t(unlist(phensim_output.sarscov2.4h.reactome$Average.Node.Perturbation)),row.names = "sarscov2-4h")
colnames(sarscov2.4h.all.genes.PHENSIM) <- phensim_output.sarscov2.4h.reactome$X..Pathway.Id


sarscov2.8h.all.genes.PHENSIM <- data.frame(t(unlist(phensim_output.sarscov2.8h.reactome$Average.Node.Perturbation)),row.names = "sarscov2-8h")
colnames(sarscov2.4h.all.genes.PHENSIM) <- phensim_output.sarscov2.8h.reactome$X..Pathway.Id


sarscov2.12h.all.genes.PHENSIM <- data.frame(t(unlist(phensim_output.sarscov2.12h.reactome$Average.Node.Perturbation)),row.names = "sarscov2-12h")
colnames(sarscov2.4h.all.genes.PHENSIM) <- phensim_output.sarscov2.12h.reactome$X..Pathway.Id

bind_rows(sarscov1.4h.all.genes.PHENSIM,sarscov1.8h.all.genes.PHENSIM,sarscov1.12h.all.genes.PHENSIM,
          sarscov2.4h.all.genes.PHENSIM,sarscov2.8h.all.genes.PHENSIM,sarscov2.12h.all.genes.PHENSIM)
## INTERFERONS
sarscov1.4h.all.genes.PHENSIM.interferons <- data.frame(t(unlist(phensim_output.sarscov1.4h.reactome[grepl("interferon", phensim_output.sarscov1.4h.reactome$Pathway.Name,ignore.case = TRUE) | grepl("ifn", phensim_output.sarscov1.4h.reactome$Pathway.Name,ignore.case = TRUE) ,]$Average.Node.Perturbation)),row.names = "sarscov1-4h")
colnames(sarscov1.4h.all.genes.PHENSIM.interferons) <- phensim_output.sarscov1.4h.reactome[grepl("interferon", phensim_output.sarscov1.4h.reactome$Pathway.Name,ignore.case = TRUE) | grepl("ifn", phensim_output.sarscov1.4h.reactome$Pathway.Name,ignore.case = TRUE) ,]$Node.Name


sarscov1.8h.all.genes.PHENSIM.interferons <- data.frame(t(unlist(phensim_output.sarscov1.8h.reactome[grepl("interferon", phensim_output.sarscov1.4h.reactome$Pathway.Name,ignore.case = TRUE) | grepl("ifn", phensim_output.sarscov1.4h.reactome$Pathway.Name,ignore.case = TRUE) ,]$Average.Node.Perturbation)),row.names = "sarscov1-8h")
colnames(sarscov1.8h.all.genes.PHENSIM.interferons) <- phensim_output.sarscov1.8h.reactome[grepl("interferon", phensim_output.sarscov1.4h.reactome$Pathway.Name,ignore.case = TRUE) | grepl("ifn", phensim_output.sarscov1.4h.reactome$Pathway.Name,ignore.case = TRUE) ,]$Node.Name


sarscov1.12h.all.genes.PHENSIM.interferons <- data.frame(t(unlist(phensim_output.sarscov1.12h.reactome[grepl("interferon", phensim_output.sarscov1.4h.reactome$Pathway.Name,ignore.case = TRUE) | grepl("ifn", phensim_output.sarscov1.4h.reactome$Pathway.Name,ignore.case = TRUE) ,]$Average.Node.Perturbation)),row.names = "sarscov1-12h")
colnames(sarscov1.12h.all.genes.PHENSIM.interferons) <- phensim_output.sarscov1.12h.reactome[grepl("interferon", phensim_output.sarscov1.4h.reactome$Pathway.Name,ignore.case = TRUE) | grepl("ifn", phensim_output.sarscov1.4h.reactome$Pathway.Name,ignore.case = TRUE) ,]$Node.Name


sarscov2.4h.all.genes.PHENSIM.interferons <- data.frame(t(unlist(phensim_output.sarscov2.4h.reactome[grepl("interferon", phensim_output.sarscov1.4h.reactome$Pathway.Name,ignore.case = TRUE) | grepl("ifn", phensim_output.sarscov1.4h.reactome$Pathway.Name,ignore.case = TRUE) ,]$Average.Node.Perturbation)),row.names = "sarscov2-4h")
colnames(sarscov2.4h.all.genes.PHENSIM.interferons) <- phensim_output.sarscov2.4h.reactome[grepl("interferon", phensim_output.sarscov1.4h.reactome$Pathway.Name,ignore.case = TRUE) | grepl("ifn", phensim_output.sarscov1.4h.reactome$Pathway.Name,ignore.case = TRUE) ,]$Node.Name


sarscov2.8h.all.genes.PHENSIM.interferons <- data.frame(t(unlist(phensim_output.sarscov2.8h.reactome[grepl("interferon", phensim_output.sarscov1.4h.reactome$Pathway.Name,ignore.case = TRUE) | grepl("ifn", phensim_output.sarscov1.4h.reactome$Pathway.Name,ignore.case = TRUE) ,]$Average.Node.Perturbation)),row.names = "sarscov2-8h")
colnames(sarscov2.8h.all.genes.PHENSIM.interferons) <- phensim_output.sarscov2.8h.reactome[grepl("interferon", phensim_output.sarscov1.4h.reactome$Pathway.Name,ignore.case = TRUE) | grepl("ifn", phensim_output.sarscov1.4h.reactome$Pathway.Name,ignore.case = TRUE) ,]$Node.Name


sarscov2.12h.all.genes.PHENSIM.interferons <- data.frame(t(unlist(phensim_output.sarscov2.12h.reactome[grepl("interferon", phensim_output.sarscov1.4h.reactome$Pathway.Name,ignore.case = TRUE) | grepl("ifn", phensim_output.sarscov1.4h.reactome$Pathway.Name,ignore.case = TRUE) ,]$Average.Node.Perturbation)),row.names = "sarscov2-12h")
colnames(sarscov2.12h.all.genes.PHENSIM.interferons) <- phensim_output.sarscov2.12h.reactome[grepl("interferon", phensim_output.sarscov1.4h.reactome$Pathway.Name,ignore.case = TRUE) | grepl("ifn", phensim_output.sarscov1.4h.reactome$Pathway.Name,ignore.case = TRUE) ,]$Node.Name


all.interferons.activityscores <- bind_rows(sarscov1.4h.all.genes.PHENSIM.interferons,sarscov1.8h.all.genes.PHENSIM.interferons,sarscov1.12h.all.genes.PHENSIM.interferons,
          sarscov2.4h.all.genes.PHENSIM.interferons,sarscov2.8h.all.genes.PHENSIM.interferons,sarscov2.12h.all.genes.PHENSIM.interferons)

library(umap)

custom.config <- umap.defaults
custom.config$n_neighbors <- 3
#custom.config$random_state <- 123
full.data.interferons.distances <- as.matrix(dist(all.interferons.activityscores))
full.data.interferons.umap.model.dist <- umap(full.data.interferons.distances, config=custom.config, input="dist")
#full.data.matrix.umap <- umap(full.data.matrix)


full.data.metadata$UMAPpathway_1 <-  full.data.matrix.umap.model.dist$layout[,1]
full.data.metadata$UMAPpathway_2 <-  full.data.matrix.umap.model.dist$layout[,2]


#see pathways activity scores

phensim_output.sarscov2.4h.reactome.endpoints <- phensim_output.sarscov2.4h.reactome %>% group_by(X..Pathway.Id) %>%
  arrange(Is.Endpoint) %>%
  slice(1) %>% ungroup


phensim_output.sarscov2.4h.endpoints.reactome.ordered <- phensim_output.sarscov2.4h.reactome.endpoints %>%
  arrange(desc(Average.Pathway.Perturbation))


pathways.sarscov2.mostperturbed <- phensim_output.sarscov2.4h.endpoints.reactome.ordered[phensim_output.sarscov2.4h.endpoints.reactome.ordered$Average.Pathway.Perturbation]


phensim_output.sarscov2.4h.endpoints.significant.ordered <- phensim_output.sarscov2.4h.endpoints %>%
  filter(Pathway.Adjusted.p.value < 0.05) %>% 
  arrange(Pathway.Activity.Score)


# Change the colors manually
p <- ggplot(data=phensim_output.sarscov2.4h.endpoints.significant.ordered, aes(x=dose, y=len, fill=supp)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()
# Use custom colors
p + scale_fill_manual(values=c('#999999','#E69F00'))
# Use brewer color palettes
p + scale_fill_brewer(palette="Blues")
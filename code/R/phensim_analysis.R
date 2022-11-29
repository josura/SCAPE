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

### TGFb and autophagy pathways

phensim_output.sarscov1.4h.reactome.autophagy <- phensim_output.sarscov1.4h.reactome[grepl("tgfb", phensim_output.sarscov1.4h.reactome$Pathway.Name,ignore.case = TRUE) | grepl("autophagy", phensim_output.sarscov1.4h.reactome$Pathway.Name,ignore.case = TRUE) ,]

autophagy.list <- unique(phensim_output.sarscov1.4h.reactome.autophagy$X..Pathway.Id)
autophagy.names <- unique(phensim_output.sarscov1.4h.reactome.autophagy$Pathway.Name)
autophagy.list <- gsub("-",".",interferons.list)

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

full.data.metadata$names <- rep(autophagy.names,times=6)
full.data.matrix.nonzerocols <- full.data.matrix %>%
  #dplyr::select(where(~ any(. != 0 )))
  select_if(~!all(is.na(.) | . == 0))
write.csv(full.data.matrix.nonzerocols,file="/home/josura/Projects/tesi/data/modernData/COVID/calu3-sarscov-interferons-PHENSIM-avg_perturbation_genes-filtered.csv")
write.csv(full.data.matrix,file="/home/josura/Projects/tesi/data/modernData/COVID/calu3-sarscov-interferons-PHENSIM-avg_perturbation_genes.csv")
write.csv(full.data.metadata,file="/home/josura/Projects/tesi/data/modernData/COVID/interferons-PHENSIM-embeddings-metadata.csv")

## ONLY autophagyS
phensim_output.sarscov1.4h.reactome.pathwayembed.autophagys <-  create.phensim.pathways.embeddings(phensim_output.sarscov1.4h.reactome[grepl("autophagy", phensim_output.sarscov1.4h.reactome$Pathway.Name,ignore.case = TRUE) | grepl("tgfb", phensim_output.sarscov1.4h.reactome$Pathway.Name,ignore.case = TRUE) ,])
rownames(phensim_output.sarscov1.4h.reactome.pathwayembed.autophagys) <- paste(autophagy.list, "-4h-sarscov1", sep="")

phensim_output.sarscov1.8h.reactome.pathwayembed.autophagys <-  create.phensim.pathways.embeddings(phensim_output.sarscov1.8h.reactome[grepl("autophagy", phensim_output.sarscov1.8h.reactome$Pathway.Name,ignore.case = TRUE) | grepl("tgfb", phensim_output.sarscov1.8h.reactome$Pathway.Name,ignore.case = TRUE) ,])
rownames(phensim_output.sarscov1.8h.reactome.pathwayembed.autophagys) <- paste(autophagy.list, "-8h-sarscov1", sep="")

phensim_output.sarscov1.12h.reactome.pathwayembed.autophagys <- create.phensim.pathways.embeddings(phensim_output.sarscov1.12h.reactome[grepl("autophagy", phensim_output.sarscov1.12h.reactome$Pathway.Name,ignore.case = TRUE) | grepl("tgfb", phensim_output.sarscov1.12h.reactome$Pathway.Name,ignore.case = TRUE) ,])
rownames(phensim_output.sarscov1.12h.reactome.pathwayembed.autophagys) <- paste(autophagy.list, "-12h-sarscov1", sep="")


phensim_output.sarscov2.4h.reactome.pathwayembed.autophagys <-  create.phensim.pathways.embeddings(phensim_output.sarscov2.4h.reactome[grepl("autophagy", phensim_output.sarscov2.4h.reactome$Pathway.Name,ignore.case = TRUE) | grepl("tgfb", phensim_output.sarscov2.4h.reactome$Pathway.Name,ignore.case = TRUE) ,])
rownames(phensim_output.sarscov2.4h.reactome.pathwayembed.autophagys) <- paste(autophagy.list, "-4h-sarscov2", sep="")

phensim_output.sarscov2.8h.reactome.pathwayembed.autophagys <-  create.phensim.pathways.embeddings(phensim_output.sarscov2.8h.reactome[grepl("autophagy", phensim_output.sarscov2.8h.reactome$Pathway.Name,ignore.case = TRUE) | grepl("tgfb", phensim_output.sarscov2.8h.reactome$Pathway.Name,ignore.case = TRUE) ,])
rownames(phensim_output.sarscov2.8h.reactome.pathwayembed.autophagys) <- paste(autophagy.list, "-8h-sarscov2", sep="")

phensim_output.sarscov2.12h.reactome.pathwayembed.autophagys <- create.phensim.pathways.embeddings(phensim_output.sarscov2.12h.reactome[grepl("autophagy", phensim_output.sarscov2.12h.reactome$Pathway.Name,ignore.case = TRUE) | grepl("tgfb", phensim_output.sarscov2.12h.reactome$Pathway.Name,ignore.case = TRUE) ,])
rownames(phensim_output.sarscov2.12h.reactome.pathwayembed.autophagys) <- paste(autophagy.list, "-12h-sarscov2", sep="")

full.data.matrix <- bind_rows(phensim_output.sarscov1.4h.reactome.pathwayembed.autophagys, phensim_output.sarscov1.8h.reactome.pathwayembed.autophagys, phensim_output.sarscov1.12h.reactome.pathwayembed.autophagys,
                              phensim_output.sarscov2.4h.reactome.pathwayembed.autophagys, phensim_output.sarscov2.8h.reactome.pathwayembed.autophagys, phensim_output.sarscov2.12h.reactome.pathwayembed.autophagys)
full.data.metadata <- data.frame(rownames(full.data.matrix),row.names = rownames(full.data.matrix))
colnames(full.data.metadata) <- "names"

full.data.metadata$time <- c(rep("4h",times = length(autophagy.list)),rep("8h",times = length(autophagy.list)),rep("12h",times = length(autophagy.list)),
                             rep("4h",times = length(autophagy.list)),rep("8h",times = length(autophagy.list)),rep("12h",times = length(autophagy.list)))

full.data.metadata$strain <- c(rep("sarscov1",times = length(autophagy.list)*3),
                               rep("sarscov2",times = length(autophagy.list)*3))
full.data.matrix.nonzerocols <- full.data.matrix %>%
  #dplyr::select(where(~ any(. != 0 )))
  select_if(~!all(is.na(.) | . == 0))
write.csv(full.data.matrix.nonzerocols,file="/home/josura/Projects/tesi/data/modernData/COVID/calu3-sarscov-autophagys-PHENSIM-avg_perturbation_genes-filtered.csv")
write.csv(full.data.matrix,file="/home/josura/Projects/tesi/data/modernData/COVID/calu3-sarscov-autophagys-PHENSIM-avg_perturbation_genes.csv")
write.csv(full.data.metadata,file="/home/josura/Projects/tesi/data/modernData/COVID/autophagys-PHENSIM-embeddings-metadata.csv")


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


#see pathways and genes activity scores

#4h
##sarscov2
phensim_output.sarscov2.4h.reactome.endpoints <- phensim_output.sarscov2.4h.reactome %>% group_by(X..Pathway.Id) %>%
  arrange(Is.Endpoint) %>%
  slice(1) %>% ungroup

phensim_output.sarscov2.4h.endpoints.reactome.ordered <- phensim_output.sarscov2.4h.reactome.endpoints %>%
  arrange(desc(Average.Pathway.Perturbation))


phensim_output.sarscov2.4h.endpoints.significant.ordered <- phensim_output.sarscov2.4h.endpoints %>%
  filter(Pathway.Adjusted.p.value < 0.05) %>% 
  arrange(Pathway.Activity.Score)
##sarscov1
phensim_output.sarscov1.4h.reactome.endpoints <- phensim_output.sarscov1.4h.reactome %>% group_by(X..Pathway.Id) %>%
  arrange(Is.Endpoint) %>%
  slice(1) %>% ungroup

phensim_output.sarscov1.4h.endpoints.reactome.ordered <- phensim_output.sarscov1.4h.reactome.endpoints %>%
  arrange(desc(Average.Pathway.Perturbation))


pathways.sarscov1.mostperturbed <- phensim_output.sarscov1.4h.endpoints.reactome.ordered[phensim_output.sarscov1.4h.endpoints.reactome.ordered$Average.Pathway.Perturbation]


phensim_output.sarscov1.4h.endpoints.significant.ordered <- phensim_output.sarscov1.4h.endpoints %>%
  filter(Pathway.Adjusted.p.value < 0.05) %>% 
  arrange(Pathway.Activity.Score)

## data binded
pathways.mostperturbed.4h <- bind_rows(phensim_output.sarscov2.4h.reactome.endpoints,phensim_output.sarscov1.4h.reactome.endpoints)

pathways.mostperturbed.4h.ordered <- pathways.mostperturbed.4h %>%
  arrange(desc(Average.Pathway.Perturbation))

###pathways
pathways.mostperturbed.4h.lesserfirstquartile <- pathways.mostperturbed.4h.ordered[pathways.mostperturbed.4h.ordered$Average.Pathway.Perturbation < quantile(pathways.mostperturbed.4h.ordered$Average.Pathway.Perturbation,probs = c(0.25)),]
pathways.mostperturbed.4h.higherthirdquartile <- pathways.mostperturbed.4h.ordered[pathways.mostperturbed.4h.ordered$Average.Pathway.Perturbation > quantile(pathways.mostperturbed.4h.ordered$Average.Pathway.Perturbation,probs = c(0.75)),]
### genes
genes.mostperturbed.4h.all <- bind_rows(phensim_output.sarscov1.4h.reactome,phensim_output.sarscov2.4h.reactome)
genes.mostperturbed.4h.lesserfirstquartile <- genes.mostperturbed.4h.all[genes.mostperturbed.4h.all$Average.Node.Perturbation < quantile(genes.mostperturbed.4h.all$Average.Node.Perturbation,probs = c(0.25)),]
genes.mostperturbed.4h.higherthirdquartile <- genes.mostperturbed.4h.all[genes.mostperturbed.4h.all$Average.Node.Perturbation > quantile(genes.mostperturbed.4h.all$Average.Node.Perturbation,probs = c(0.75)),]

#8h
##sarscov2
phensim_output.sarscov2.8h.reactome.endpoints <- phensim_output.sarscov2.8h.reactome %>% group_by(X..Pathway.Id) %>%
  arrange(Is.Endpoint) %>%
  slice(1) %>% ungroup

phensim_output.sarscov2.8h.endpoints.reactome.ordered <- phensim_output.sarscov2.8h.reactome.endpoints %>%
  arrange(desc(Average.Pathway.Perturbation))


phensim_output.sarscov2.8h.endpoints.significant.ordered <- phensim_output.sarscov2.8h.endpoints %>%
  filter(Pathway.Adjusted.p.value < 0.05) %>% 
  arrange(Pathway.Activity.Score)
##sarscov1
phensim_output.sarscov1.8h.reactome.endpoints <- phensim_output.sarscov1.8h.reactome %>% group_by(X..Pathway.Id) %>%
  arrange(Is.Endpoint) %>%
  slice(1) %>% ungroup

phensim_output.sarscov1.8h.endpoints.reactome.ordered <- phensim_output.sarscov1.8h.reactome.endpoints %>%
  arrange(desc(Average.Pathway.Perturbation))


pathways.sarscov1.mostperturbed <- phensim_output.sarscov1.8h.endpoints.reactome.ordered[phensim_output.sarscov1.8h.endpoints.reactome.ordered$Average.Pathway.Perturbation]


phensim_output.sarscov1.8h.endpoints.significant.ordered <- phensim_output.sarscov1.8h.endpoints %>%
  filter(Pathway.Adjusted.p.value < 0.05) %>% 
  arrange(Pathway.Activity.Score)

## data binded
pathways.mostperturbed.8h <- bind_rows(phensim_output.sarscov2.8h.reactome.endpoints,phensim_output.sarscov1.8h.reactome.endpoints)

pathways.mostperturbed.8h.ordered <- pathways.mostperturbed.8h %>%
  arrange(desc(Average.Pathway.Perturbation))

###pathways
pathways.mostperturbed.8h.lesserfirstquartile <- pathways.mostperturbed.8h.ordered[pathways.mostperturbed.8h.ordered$Average.Pathway.Perturbation < quantile(pathways.mostperturbed.8h.ordered$Average.Pathway.Perturbation,probs = c(0.25)),]
pathways.mostperturbed.8h.higherthirdquartile <- pathways.mostperturbed.8h.ordered[pathways.mostperturbed.8h.ordered$Average.Pathway.Perturbation > quantile(pathways.mostperturbed.8h.ordered$Average.Pathway.Perturbation,probs = c(0.75)),]
### genes
genes.mostperturbed.8h.all <- bind_rows(phensim_output.sarscov1.8h.reactome,phensim_output.sarscov2.8h.reactome)
genes.mostperturbed.8h.lesserfirstquartile <- genes.mostperturbed.8h.all[genes.mostperturbed.8h.all$Average.Node.Perturbation < quantile(genes.mostperturbed.8h.all$Average.Node.Perturbation,probs = c(0.25)),]
genes.mostperturbed.8h.higherthirdquartile <- genes.mostperturbed.8h.all[genes.mostperturbed.8h.all$Average.Node.Perturbation > quantile(genes.mostperturbed.8h.all$Average.Node.Perturbation,probs = c(0.75)),]



#12h
##sarscov2
phensim_output.sarscov2.12h.reactome.endpoints <- phensim_output.sarscov2.12h.reactome %>% group_by(X..Pathway.Id) %>%
  arrange(Is.Endpoint) %>%
  slice(1) %>% ungroup

phensim_output.sarscov2.12h.endpoints.reactome.ordered <- phensim_output.sarscov2.12h.reactome.endpoints %>%
  arrange(desc(Average.Pathway.Perturbation))


phensim_output.sarscov2.12h.endpoints.significant.ordered <- phensim_output.sarscov2.12h.endpoints %>%
  filter(Pathway.Adjusted.p.value < 0.05) %>% 
  arrange(Pathway.Activity.Score)
##sarscov1
phensim_output.sarscov1.12h.reactome.endpoints <- phensim_output.sarscov1.12h.reactome %>% group_by(X..Pathway.Id) %>%
  arrange(Is.Endpoint) %>%
  slice(1) %>% ungroup

phensim_output.sarscov1.12h.endpoints.reactome.ordered <- phensim_output.sarscov1.12h.reactome.endpoints %>%
  arrange(desc(Average.Pathway.Perturbation))


pathways.sarscov1.mostperturbed <- phensim_output.sarscov1.12h.endpoints.reactome.ordered[phensim_output.sarscov1.12h.endpoints.reactome.ordered$Average.Pathway.Perturbation]


phensim_output.sarscov1.12h.endpoints.significant.ordered <- phensim_output.sarscov1.12h.endpoints %>%
  filter(Pathway.Adjusted.p.value < 0.05) %>% 
  arrange(Pathway.Activity.Score)

## data binded
pathways.mostperturbed.12h <- bind_rows(phensim_output.sarscov2.12h.reactome.endpoints,phensim_output.sarscov1.12h.reactome.endpoints)

pathways.mostperturbed.12h.ordered <- pathways.mostperturbed.12h %>%
  arrange(desc(Average.Pathway.Perturbation))

###pathways
pathways.mostperturbed.12h.lesserfirstquartile <- pathways.mostperturbed.12h.ordered[pathways.mostperturbed.12h.ordered$Average.Pathway.Perturbation < quantile(pathways.mostperturbed.12h.ordered$Average.Pathway.Perturbation,probs = c(0.25)),]
pathways.mostperturbed.12h.higherthirdquartile <- pathways.mostperturbed.12h.ordered[pathways.mostperturbed.12h.ordered$Average.Pathway.Perturbation > quantile(pathways.mostperturbed.12h.ordered$Average.Pathway.Perturbation,probs = c(0.75)),]
### genes
genes.mostperturbed.12h.all <- bind_rows(phensim_output.sarscov1.12h.reactome,phensim_output.sarscov2.12h.reactome)
genes.mostperturbed.12h.lesserfirstquartile <- genes.mostperturbed.12h.all[genes.mostperturbed.12h.all$Average.Node.Perturbation < quantile(genes.mostperturbed.12h.all$Average.Node.Perturbation,probs = c(0.25)),]
genes.mostperturbed.12h.higherthirdquartile <- genes.mostperturbed.12h.all[genes.mostperturbed.12h.all$Average.Node.Perturbation > quantile(genes.mostperturbed.12h.all$Average.Node.Perturbation,probs = c(0.75)),]


# load Venn diagram package
install.packages("VennDiagram")
library("VennDiagram")



# create pairwise Venn diagram

create.venndiagram.pairwise <- function(pathways){
  area.sarscov1 <- dim(pathways[pathways$strain == "sarscov1",])[1]
  area.sarscov2 <- dim(pathways[pathways$strain == "sarscov2",])[1]
  grid.newpage()
  draw.pairwise.venn(area1=area.sarscov1, 
                     area2=area.sarscov2,
                     cross.area=length(intersect(pathways[pathways$strain == "sarscov1",]$X..Pathway.Id,
                                                 pathways[pathways$strain == "sarscov2",]$X..Pathway.Id)),
                     category=c("sarscov1","sarscov2"),fill=c("Red","Yellow"))
}

## 4h

create.venndiagram.pairwise(pathways.mostperturbed.4h.lesserfirstquartile)
create.venndiagram.pairwise(pathways.mostperturbed.4h.higherthirdquartile)

create.venndiagram.pairwise(genes.mostperturbed.4h.lesserfirstquartile)
create.venndiagram.pairwise(genes.mostperturbed.4h.higherthirdquartile)


## 8h
create.venndiagram.pairwise(pathways.mostperturbed.8h.lesserfirstquartile)
create.venndiagram.pairwise(pathways.mostperturbed.8h.higherthirdquartile)

create.venndiagram.pairwise(genes.mostperturbed.8h.lesserfirstquartile)
create.venndiagram.pairwise(genes.mostperturbed.8h.higherthirdquartile)

## 12h
create.venndiagram.pairwise(pathways.mostperturbed.12h.lesserfirstquartile)
create.venndiagram.pairwise(pathways.mostperturbed.12h.higherthirdquartile)


create.venndiagram.pairwise(genes.mostperturbed.12h.lesserfirstquartile)
create.venndiagram.pairwise(genes.mostperturbed.12h.higherthirdquartile)


#box plots
library(ggplot2)

# create a data frame
all.pathways.avgperturbation <- bind_rows(pathways.mostperturbed.4h,
                                          pathways.mostperturbed.8h,
                                          pathways.mostperturbed.12h)
variety=all.pathways.avgperturbation$time
strain=all.pathways.avgperturbation$strain
Average.Pathway.Perturbation=all.pathways.avgperturbation$Average.Pathway.Perturbation
data=data.frame(variety, strain ,  Average.Pathway.Perturbation)



# grouped boxplot
times <- factor(all.pathways.avgperturbation$time, level = c("4h", "8h", "12h"))

ggplot(data, aes(x=times, y=Average.Pathway.Perturbation, fill=strain)) + 
  geom_boxplot()


# most 20 perturbed pathways for sarscov2
#12h
##sarscov2
most20perturbed.sarscov2.12h <- phensim_output.sarscov2.12h.endpoints.reactome.ordered[1:20,]

##sarscov1
most20perturbed.sarscov2.12h.all <-bind_rows(most20perturbed.sarscov2.12h, phensim_output.sarscov1.12h.endpoints.reactome.ordered[phensim_output.sarscov1.12h.endpoints.reactome.ordered$X..Pathway.Id %in% most20perturbed.sarscov2.12h$X..Pathway.Id,])
#8h
##sarscov2
most20perturbed.sarscov2.12h.all <- bind_rows(most20perturbed.sarscov2.12h.all, phensim_output.sarscov2.8h.endpoints.reactome.ordered[phensim_output.sarscov2.8h.endpoints.reactome.ordered$X..Pathway.Id %in% most20perturbed.sarscov2.12h.all$X..Pathway.Id,])

##sarscov1
most20perturbed.sarscov2.12h.all <-bind_rows(most20perturbed.sarscov2.12h.all, phensim_output.sarscov1.8h.endpoints.reactome.ordered[phensim_output.sarscov1.8h.endpoints.reactome.ordered$X..Pathway.Id %in% most20perturbed.sarscov2.12h.all$X..Pathway.Id,])
#4h
##sarscov2
most20perturbed.sarscov2.12h.all <- bind_rows(most20perturbed.sarscov2.12h.all, phensim_output.sarscov2.4h.endpoints.reactome.ordered[phensim_output.sarscov2.4h.endpoints.reactome.ordered$X..Pathway.Id %in% most20perturbed.sarscov2.12h.all$X..Pathway.Id,])

##sarscov1
most20perturbed.sarscov2.12h.all <-bind_rows(most20perturbed.sarscov2.12h.all, phensim_output.sarscov1.4h.endpoints.reactome.ordered[phensim_output.sarscov1.4h.endpoints.reactome.ordered$X..Pathway.Id %in% most20perturbed.sarscov2.12h.all$X..Pathway.Id,])



variety=most20perturbed.sarscov2.12h.all$time
strain=most20perturbed.sarscov2.12h.all$strain
pathway.name <- most20perturbed.sarscov2.12h.all$Pathway.Name
Average.Pathway.Perturbation <- most20perturbed.sarscov2.12h.all$Average.Pathway.Perturbation
data=data.frame(variety, strain ,  pathway.name,Average.Pathway.Perturbation)



# grouped boxplot
times <- factor(variety, level = c("4h", "8h", "12h"))
pathways <- factor(pathway.name, level=most20perturbed.sarscov2.12h$Pathway.Name)

ggplot(data, aes(x=times, y=Average.Pathway.Perturbation, color=strain)) + 
  geom_point() +
  facet_wrap(facets = vars(pathway.name))

  # Change the colors manually
p <- ggplot(data=phensim_output.sarscov2.4h.endpoints.significant.ordered, aes(x=dose, y=len, fill=supp)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()
# Use custom colors
p + scale_fill_manual(values=c('#999999','#E69F00'))
# Use brewer color palettes
p + scale_fill_brewer(palette="Blues")
library(dplyr)

cluster0.nonconsidered <- read.csv("/home/josura/Projects/tesi/data/patient4ALL/cluster0-nonConsidered.tsv",header = FALSE, sep = "\t")

Phensim.input.cluster0.nonconsidered <- cluster0.nonconsidered %>%
  select(V3,V1)

write.table(Phensim.input.cluster0.nonconsidered,
            file='/home/josura/Projects/tesi/data/patient4ALL/cluster0-nonConsideredPHENSIMinput.tsv', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)


#cluster0
cluster0.PHENSIM.input <- read.csv("/home/josura/Projects/tesi/data/patient4ALL/cluster0-diffExpressedMarkers.tsv",header = FALSE, sep = "\t")

cluster0.PHENSIM.input <- cluster0.PHENSIM.input %>%
  filter(V1 =="UP" | V1 == "DOWN") %>%
  dplyr::select(V3,V1)

# change UP to OVEREXPRESSED
cluster0.PHENSIM.input$V1[cluster0.PHENSIM.input$V1 == "UP"] <- "OVEREXPRESSION"
# change DOWN to UNDEREXPRESSED
cluster0.PHENSIM.input$V1[cluster0.PHENSIM.input$V1 == "DOWN"] <- "UNDEREXPRESSION"
write.table(cluster0.PHENSIM.input,
            file='/home/josura/Projects/tesi/data/patient4ALL/cluster0-PHENSIMinput.txt', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)


#cluster1
cluster1.PHENSIM.input <- read.csv("/home/josura/Projects/tesi/data/patient4ALL/cluster1-diffExpressedMarkers.tsv",header = FALSE, sep = "\t")

cluster1.PHENSIM.input <- cluster1.PHENSIM.input %>%
  filter(V1 =="UP" | V1 == "DOWN") %>%
  dplyr::select(V3,V1)

# change UP to OVEREXPRESSED
cluster1.PHENSIM.input$V1[cluster1.PHENSIM.input$V1 == "UP"] <- "OVEREXPRESSION"
# change DOWN to UNDEREXPRESSED
cluster1.PHENSIM.input$V1[cluster1.PHENSIM.input$V1 == "DOWN"] <- "UNDEREXPRESSION"
write.table(cluster1.PHENSIM.input,
            file='/home/josura/Projects/tesi/data/patient4ALL/cluster1-PHENSIMinput.txt', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)


#cluster2
cluster2.PHENSIM.input <- read.csv("/home/josura/Projects/tesi/data/patient4ALL/cluster2-diffExpressedMarkers.tsv",header = FALSE, sep = "\t")

cluster2.PHENSIM.input <- cluster2.PHENSIM.input %>%
  filter(V1 =="UP" | V1 == "DOWN") %>%
  dplyr::select(V3,V1)

# change UP to OVEREXPRESSED
cluster2.PHENSIM.input$V1[cluster2.PHENSIM.input$V1 == "UP"] <- "OVEREXPRESSION"
# change DOWN to UNDEREXPRESSED
cluster2.PHENSIM.input$V1[cluster2.PHENSIM.input$V1 == "DOWN"] <- "UNDEREXPRESSION"
write.table(cluster2.PHENSIM.input,
            file='/home/josura/Projects/tesi/data/patient4ALL/cluster2-PHENSIMinput.txt', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)



#bulk
bulk.PHENSIM.input <- read.csv("/home/josura/Projects/tesi/data/patient4ALL/bulk-allGenesDE.tsv",header = TRUE, sep = "\t")

bulk.PHENSIM.input <- bulk.PHENSIM.input %>%
  filter(diffexpressed =="UP" | diffexpressed == "DOWN") %>%
  dplyr::select(entrezgene_id,diffexpressed) %>%
  distinct()

# change UP to OVEREXPRESSED
bulk.PHENSIM.input$diffexpressed[bulk.PHENSIM.input$diffexpressed == "UP"] <- "OVEREXPRESSION"
# change DOWN to UNDEREXPRESSED
bulk.PHENSIM.input$diffexpressed[bulk.PHENSIM.input$diffexpressed == "DOWN"] <- "UNDEREXPRESSION"
write.table(bulk.PHENSIM.input,
            file='/home/josura/Projects/tesi/data/patient4ALL/bulk-PHENSIMinput.txt', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)

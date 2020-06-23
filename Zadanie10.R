# biblioteka ggplot2 do wykresów
library(ggplot2)
# biblioteki GO
library(goseq)
library(topGO)
# biblioteki operacji na danych
library(tidyr)
library(stringr)

# --- KOD Z ĆWICZEŃ --- START
dane <- read.csv("/Dydaktyka/counts.txt", skip = 1, sep = "\t")
gene.bad.id <- dane$Geneid
geneID2GO <- readMappings(file = '/Dydaktyka/geneid2go.map2')
r <- read.csv('/Dydaktyka/DE_wyniki.csv', row.names = 1)
gene.bad.id <- dane$Geneid
popraw.nazwy <- function(gen.id){
  dobry.gen.id <- strsplit(as.character(gen.id), split = ".", fixed = T)[[1]][1]
  return(dobry.gen.id)
}
r1 <- as.data.frame(r)
de.gene <- row.names(r1)
de.gene <- sapply(de.gene, popraw.nazwy)
assigned.gene <- sapply(gene.bad.id, popraw.nazwy)
gene.vector <- as.integer(assigned.gene%in%de.gene)
names(gene.vector) <- assigned.gene
length.data <- dane$Length
pwf <- nullp(gene.vector,"unnOrg","unnGen", length.data)
# --- KOD Z ĆWICZEŃ --- KONIEC

GO.wall <- goseq(pwf,"unnOrg","unnGen", gene2cat = geneID2GO, test.cats = c("GO:CC", "GO:BP", "GO:MF"), method = "Wallenius")
GO.wall.CC <- GO.wall[GO.wall$ontology == 'CC',]
enriched.GO <- GO.wall.CC$category[p.adjust(GO.wall.CC$over_represented_pvalue, method = "BH") < .05]
GO.wall1.CC <- GO.wall.CC[GO.wall.CC$category%in%enriched.GO,]
GO.wall1.CC$FDR <- p.adjust(GO.wall.CC$over_represented_pvalue,method = "BH")[p.adjust(GO.wall.CC$over_represented_pvalue,method = "BH") < .05]
GO.wall1.CC <- GO.wall1.CC[!is.na(GO.wall1.CC$FDR),]
GO.wall.MF <- GO.wall[GO.wall$ontology == 'MF',]
enriched.GO <- GO.wall.MF$category[p.adjust(GO.wall.MF$over_represented_pvalue,method = "BH")<.05]
GO.wall1.MF <- GO.wall.MF[GO.wall.MF$category%in%enriched.GO,]
GO.wall1.MF$FDR <- p.adjust(GO.wall.MF$over_represented_pvalue,method = "BH")[p.adjust(GO.wall.MF$over_represented_pvalue,method = "BH") < .05]
GO.wall1.MF <- GO.wall1.MF[!is.na(GO.wall1.MF$FDR),]
GO.wall.BP <- GO.wall[GO.wall$ontology == 'BP',]
enriched.GO <- GO.wall.BP$category[p.adjust(GO.wall.BP$over_represented_pvalue, method = "BH")<.05]
GO.wall1.BP <- GO.wall.BP[GO.wall.BP$category%in%enriched.GO,]
GO.wall1.BP$FDR <- p.adjust(GO.wall.BP$over_represented_pvalue,method="BH")[p.adjust(GO.wall.BP$over_represented_pvalue,method = "BH") < .05]
GO.wall1.BP <- GO.wall1.BP[!is.na(GO.wall1.BP$FDR),]
polaczone <- rbind(GO.wall1.BP, GO.wall1.MF, GO.wall1.CC)

row.names(r) <- de.gene
r <- cbind(rownames(r), data.frame(r, row.names = NULL))
colnames(r) <- c('loc','logFC','logCPM','PValue','FDR')
genes_ <- data.frame(unlist(geneID2GO, use.names = TRUE))
geneid <- cbind(rownames(genes_), data.frame(genes_, row.names = NULL))
colnames(geneid) <- c('loc','GO')
geneid$loc <- substr(geneid$loc,1,11)
geneMerged <- merge(r,geneid)
geneMerged <- na.omit(geneMerged)

polaczone$plus <- 0
polaczone$minus <- 0
for (m in 1:length(geneMerged$logFC)){
  if (geneMerged$logFC[m] > 0) {
    for (n in 1:length(polaczone$category)){
      if (polaczone$category[n] == geneMerged$GO[m]){
        polaczone$plus[n] <- polaczone$plus[n] + 1
      }
    }
  } else if(geneMerged$logFC[m] <= 0){
    for (n in 1:length(polaczone$category)){
      if (polaczone$category[n] == geneMerged$GO[m]){
        polaczone$minus[n] <- polaczone$minus[n] + 1
      }
    }
  }
}

plotData <- data.frame(polaczone$category, polaczone$ontology, polaczone$plus, polaczone$minus)
colnames(plotData) <- c("category","ontology","Up","Down")
plotData$DGEs_Number <- 0
plotData <- rbind(plotData,plotData)
plotData$DGE <- ''
for (s in 1:(length(plotData$category)/2)){
  plotData[s,5] <- plotData[s,3]
  plotData[s,6] <- "Up"
  s <- s+length(plotData$category)/2
  plotData[s,5] <- plotData[s,4]
  plotData[s,6] <- "Down"
}

ggplot(data = plotData, aes(x = category, y = DGEs_Number, fill = DGE)) +
  geom_bar(stat="identity", position=position_dodge()) +
  coord_flip() +
  theme_minimal() +
  scale_fill_manual("", values = c("Up" = "red", "Down" = "blue")) +
  ylab("DGEs Number") +
  xlab("Category") +
  geom_text(aes(label = DGEs_Number), position=position_dodge(width=0.1), size = 2) +
  facet_grid(ontology ~ ., scales = "free")
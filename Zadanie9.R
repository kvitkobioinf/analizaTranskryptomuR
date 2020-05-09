# Zadanie 8 -- start

BiocManager::install(c("DESeq2"))
library(DESeq2)

data22 <- read.csv('source_files/counts_22.txt', sep = '\t', skip = 1)
dataERCC92 <- read.csv('source_files/counts_ERCC92.txt', sep = '\t', skip = 1)

normalizeDDS <- function(data){
  countData <- data[,7:12]
  rownames(countData) = data$Geneid
  samples <- names(countData)
  cond_1 <- rep("cond1", 3)
  cond_2 <- rep("cond2", 3)
  condition <- factor(c(cond_1, cond_2))
  colData <- data.frame(samples = samples, condition = condition)
  dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~condition)
  
  return(dds)
}

dds_22 <- normalizeDDS(data22)
dds_ERCC99 <- normalizeDDS(dataERCC92)

# Zadanie 8 -- koniec



# Zadanie 9 -- start

summarizeDDS <- function(ddsData){
  dds <- ddsData
  
  dds <- DESeq(dds)
  res <- results(dds)
  
  r <- res[res$baseMean != 0,] # pominięcie genów bez ekspresji (zerowej)
  r <- r[r$log2FoldChange > 1 | r$log2FoldChange < -1,] # założenie, że istotne są geny z co najmniej dwukrotną zmianą ekspresji 2^1 lub 2^-1
  x <- !is.na(r$padj) # usunięcie danych brakujących (NA)
  r <- r[x,]
  r <- r[r$padj < 0.05,] # wybór genów o istotności statystycznej < 5%
  print(paste("Geny o istotnie statystycznie zmienionym poziomie ekspresji pomiędzy badanymi warunkami:", nrow(r), sep = ""))
  print(r)
}

summarizeDDS(dds_22)
summarizeDDS(dds_ERCC99)

# Zadanie 9 -- koniec
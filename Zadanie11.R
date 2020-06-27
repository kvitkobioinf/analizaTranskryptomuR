library(DESeq2)

dataERCC92 <- read.delim("~/cwiczenia/zad11/counts_ERCC92.txt", comment.char = "#")
dataChrom22 <- read.delim("~/cwiczenia/zad11/counts_22.txt", comment.char = "#")
dataProvided <- read.delim("~/cwiczenia/zad11/dane11.txt", comment.char = "#")

normalizeDDS <- function(data){
  countData <- data[,7:12]
  rownames(countData) <- data$Geneid
  samples <- names(countData)
  cond_1 <- rep("cond1", 3)
  cond_2 <- rep("cond2", 3)
  condition <- factor(c(cond_1, cond_2))
  colData <- data.frame(samples = samples, condition = condition)
  dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~condition)
  
  return(dds)
}

analyseDESeq <- function(Data, ymin, ymax){
  dds <- DESeq(Data)
  res <- results(dds)
  
  layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
  plotMA(res, ylim=c(ymin,ymax))
  plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
  plotCounts(dds, gene=which.max(res$padj), intgroup="condition")
  
  resOrdered <- res[order(res$padj),]
  resSig <- subset(resOrdered, padj<0.05)
  return(resSig)
}

normalizedERCC92 <- normalizeDDS(dataERCC92)
normalizedChrom22 <- normalizeDDS(dataChrom22)

dddERCC92 <- analyseDESeq(normalizedERCC92, -4, 4)
dddChrom22 <- analyseDESeq(normalizedChrom22, -10, 10)

dddERCC92 <- cbind(rownames(dddERCC92), data.frame(dddERCC92, row.names=NULL))
colnames(dddERCC92) <- c("ERCC.ID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")

Test.corr <- data.frame(dddERCC92$ERCC.ID, dddERCC92$log2FoldChange)
colnames(Test.corr) <- c("ERCC.ID","log2FoldChange")
rownames(Test.corr) <- Test.corr$ERCC.ID
compare <- data.frame(dataProvided[,2],dataProvided[,7])
colnames(compare) <- c("ERCC.ID","log2")
rownames(compare) <- compare$ERCC.ID

Test.corr <- merge(Test.corr, compare, by=0, all=TRUE)
Test.corr <- na.omit(Test.corr)
rownames(Test.corr) <- Test.corr$Row.names

corRes <- cor.test(Test.corr$log2, Test.corr$log2FoldChange,method = "pearson", )
print(corRes)
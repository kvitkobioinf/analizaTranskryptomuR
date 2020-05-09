#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.10")

#BiocManager::install(c("Rsubread"))
#library(Rsubread)
#HBR_UHR <- c("../../Zadanie 8/_pliki/BAM/HBR_1.bam",
#           "../../Zadanie 8/_pliki/BAM/HBR_2.bam",
#           "../../Zadanie 8/_pliki/BAM/HBR_3.bam",
#           "../../Zadanie 8/_pliki/BAM/UHR_1.bam",
#           "../../Zadanie 8/_pliki/BAM/UHR_2.bam",
#           "../../Zadanie 8/_pliki/BAM/UHR_3.bam")
#ERCC92 <- "../../Zadanie 8/_pliki/ERCC92.gtf"
#HBR_UHR_ERCC92_fe <- featureCounts(files = HBR_UHR, annot.ext = ERCC92, isGTFAnnotationFile = T)
#head(HBR_UHR_ERCC92_fe)
#ostatecznie pliki przetworzono z wersji linuksowej programu featureCounts

BiocManager::install(c("DESeq2"))
library(DESeq2)

data22 <- read.csv('source_files/counts_22.txt', sep = '\t', skip = 1)
dataERCC92 <- read.csv('source_files/counts_ERCC92.txt', sep = '\t', skip = 1)

normalize <- function(data){
  countData <- data[,7:12]
  rownames(countData) = data$Geneid
  samples <- names(countData)
  cond_1 <- rep("cond1", 3)
  cond_2 <- rep("cond2", 3)
  condition <- factor(c(cond_1, cond_2))
  colData <- data.frame(samples = samples, condition = condition)
  dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~condition)
  dss1 <- estimateSizeFactors(dds)
  
  log_data <- rlog(dds)
  normalized_data<- assay(log_data)
  normalized_data <- as.data.frame(normalized_data)
  
  return(normalized_data)
}

normalized_data22 <- normalize(data22)
normalized_ERCC99 <- normalize(dataERCC92)




#BiocManager::install(c("ComplexHeatmap"))
library(ComplexHeatmap)

drawHeatmap <- function(data){
  threshold <- 10
  data <- data[
                data$X.home.bioinformatyka.s119494.BAM_indeks.BAM_FILES.HBR_1.bam > threshold | 
                data$X.home.bioinformatyka.s119494.BAM_indeks.BAM_FILES.HBR_2.bam > threshold | 
                data$X.home.bioinformatyka.s119494.BAM_indeks.BAM_FILES.HBR_3.bam > threshold | 
                data$X.home.bioinformatyka.s119494.BAM_indeks.BAM_FILES.UHR_1.bam > threshold | 
                data$X.home.bioinformatyka.s119494.BAM_indeks.BAM_FILES.UHR_2.bam > threshold | 
                data$X.home.bioinformatyka.s119494.BAM_indeks.BAM_FILES.UHR_3.bam > threshold
                , 
              ]
  renamed_data <- data
  colnames(renamed_data) = c("HBR1","HBR2","HBR3", "UHR1", "UHR2", "UHR3")
  
  Heatmap(renamed_data , cluster_columns = FALSE,
          row_names_side = "left",
          row_dend_sid = "left",
          row_names_gp=gpar(cex = 0.8))
}

drawHeatmap(normalized_data22)
drawHeatmap(normalized_ERCC99)

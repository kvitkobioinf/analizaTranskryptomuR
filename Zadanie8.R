#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.10")

#BiocManager::install(c("Rsubread"))
library(Rsubread)

HBR_UHR <- c("../../Zadanie 8/_pliki/BAM/HBR_1.bam",
           "../../Zadanie 8/_pliki/BAM/HBR_2.bam",
           "../../Zadanie 8/_pliki/BAM/HBR_3.bam",
           "../../Zadanie 8/_pliki/BAM/UHR_1.bam",
           "../../Zadanie 8/_pliki/BAM/UHR_2.bam",
           "../../Zadanie 8/_pliki/BAM/UHR_3.bam")
ERCC92 <- "../../Zadanie 8/_pliki/ERCC92.gtf"
HBR_UHR_ERCC92_fe <- featureCounts(files = HBR_UHR, annot.ext = ERCC92, isGTFAnnotationFile = T)

# Zadanie 1 -- start
raw_data <- read.csv("source_files\\counts.txt", sep = "\t", skip = 1)
geneLengths <- raw_data[, c(1, 6:9)]
TPM_step1 <- geneLengths
TPM_step1$bam.flower.bam <- TPM_step1$bam.flower.bam / TPM_step1$Length
TPM_step1$bam.stem.bam <- TPM_step1$bam.stem.bam / TPM_step1$Length
TPM_step1$bam.leaf.bam <- TPM_step1$bam.leaf.bam / TPM_step1$Length
TPM_step2 <- TPM_step1
TPM_step2$bam.flower.bam <- TPM_step2$bam.flower.bam / (sum(TPM_step2$bam.flower.bam) / 1000000)
TPM_step2$bam.stem.bam <- TPM_step2$bam.stem.bam / (sum(TPM_step2$bam.stem.bam) / 1000000)
TPM_step2$bam.leaf.bam <- TPM_step2$bam.leaf.bam / (sum(TPM_step2$bam.leaf.bam) / 1000000)
dane_TPM <- TPM_step2[,c(1,3:5)]
colnames(dane_TPM)[2:4] <- c("liść_TPM", "pęd_TPM", "kwiat_TPM")
# Zadanie 1 -- koniec

# Zadanie 2 -- start
MYB <- read.csv("source_files\\MYB.txt", sep = "\t")
library(dplyr)
dane_myb <- dane_TPM %>% filter(Geneid %in% MYB$Locus)
# Zadanie 2 -- koniec

MYB_Zad3 <- read.csv("source_files\\MYB_Zad3.txt", sep = "\t")
CEN <- read.csv("source_files\\CEN.txt", sep = "\t")

dane_myb <- dane_TPM %>% filter(Geneid %in% MYB_Zad3$Locus)
dane_cen <- dane_TPM %>% filter(Geneid %in% CEN$Locus)

dane_myb_AND_cen <- dane_myb %>% filter(Geneid %in% CEN$Locus)

dane_myb_OR_cen <- rbind(dane_myb, dane_cen, )
sum(duplicated(dane_myb_OR_cen$Geneid))
dim(dane_myb_OR_cen)


# Weryfikacja poprawności szukania duplikatów wybierając losowy wiersz z tabeli MYC
Geneid <- c("LOC109339520")
liść_TPM <- c(0.12450058)
pęd_TPM <- c(0.12183480)
kwiat_TPM <- (0.52207619)
dummyDataFrame <- data.frame(Geneid,	liść_TPM,	pęd_TPM,	kwiat_TPM)
dataFrameWithDuplicatesTest <- rbind(dane_myb_OR_cen, dummyDataFrame)
sum(duplicated(dataFrameWithDuplicatesTest$Geneid))

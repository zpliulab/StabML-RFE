## clear 
rm(list=ls())

## packages
library(dplyr)       
library(tidyr)
library(tidyverse)    


## load data
setwd('D:\\E\\²©Ê¿\\R_³ÌÐò\\HGSOC')

# feature gene  data
Data1 = read.table("Data\\GSE69428_scale.txt", header = T, check.names = FALSE)
gene = read.csv("datanew\\bio.csv", header=TRUE, sep = ',')
dim(Data1)

colnames(gene) <- c('gene')
Data2 <- cbind(rownames(Data1), Data1)
colnames(Data2) <- c('gene', colnames(Data1))

genedata <- merge(gene, Data2, by = "gene")#[,-2]
genedata1 <- genedata %>% tibble::column_to_rownames(colnames(.)[1])
genedata2 <- rbind(Data1[1,],genedata1)
rownames(genedata2) <- c('Lable', rownames(genedata1))
dim(genedata2)    #  19 20
# write.table(genedata2,"datanew\\GSE69428_bio_forDEGs.txt",quote=F,sep="\t") 



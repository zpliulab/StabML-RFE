
##
rm(list = ls())


setwd("D:\\E\\博士\\R_程序\\HGSOC\\DataGEOTogather")

## 读入数据
data1 = read.table("GSE69428_outcome.txt", header = T, sep='\t', fill=TRUE, strip.white = T, check.names = F)
data2 = read.table("GSE27651_outcome.txt", header = T, sep='\t', fill=TRUE, strip.white = T, check.names = F)


## merge
data22 <- rbind(data2[1,],data2[rownames(data1)[-1],])
data <- cbind(data1,data22)
dim(data)

# my_scale -------------------------------------------------------------------

my_scale <- function(x){
  x1 <- cbind(t(x[1,]), scale(t(x[-1,])))
  x2 <- t(x1)
  return(x2)
}

# 重组 ----------------------------------------------------------------------

data1 <- my_scale(data)
dim(data1)    # 17690    48
# write.table(data1,"GSE69428_27651scale.txt",quote=F,sep="\t")


# DE  ---------------------------------------------------------------------

rm(list = ls())
setwd("D:\\E\\博士\\R_程序\\HGSOC\\DataGEOtogather")

Data = read.table("GSE69428_27651scale.txt", header = T, check.names = FALSE)


# phenotype label ----------------------------------------------------------------------

sample <- as.matrix(colnames(Data))
label <- t(Data[1,])
class <- cbind(sample, label)
colnames(class) <- c("sample","class")
# write.csv(class, file = "phenotype.csv", row.names = F)


# Rank and newdata ----------------------------------------------------------------------

Data1 <- Data[,which(Data[1,] == 1)]
Data0 <- Data[,which(Data[1,] == 0)]
data <- cbind(Data0,Data1)[-1,]    # 注意！！-- 要与数据对应好，哪儿是0/1


phe = read.csv("phenotype.csv", header=TRUE, sep = ',')
phe1 <- phe[which(phe[,2] == 1),]
phe1[,2] <- c("Normal")
phe0 <- phe[which(phe[,2] == 0),]
phe0[,2] <- c("HGSOC")
phenotype <- rbind(phe0,phe1)
# write.table(phenotype, file = "phe_zh.txt", quote=F, sep="\t", row.names = F)

# DE analysis ----------------------------------------------------------------
disease = read.table("phe_zh.txt",header = TRUE, sep = "\t")
disease <- factor(disease[,"class"])


library(simpleaffy)
library(affyPLM)
library(RColorBrewer)
library(graph)
library(genefilter)
library(affycoretools)
library(affy)
library(limma)



design <- model.matrix(~-1+disease)
contrast.matrix <- makeContrasts(contrasts = "diseaseNormal - diseaseHGSOC", levels = design)
fit <- lmFit(data,design)
fit1 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit1)
dif <- topTable(fit2,coef = "diseaseNormal - diseaseHGSOC",n = nrow(fit2),lfc = log2(2))
dif0.05 <- dif[dif[,"adj.P.Val"] < 0.05,]    #  1768    6
dim(dif0.05)
dif_FC <- dif0.05[abs(dif0.05[,"logFC"]) > 1.3,]     # 533   6
dim(dif_FC)
# write.csv(dif, file = "dif.csv")  # 
# write.csv(dif_FC, file = "dif_FC.csv")

gene_adjp <- rownames(dif_FC)
# write.csv(gene_adjp, file = "gene_adjp.csv", row.names = F)



# 差异gene的表达数据 -------------------------------------------------------------
rm(list = ls())


# merge
setwd("D:\\E\\博士\\R_程序\\HGSOC\\DataGEOtogather")
data <- read.table("GSE69428_27651scale.txt", header = T, check.names = FALSE)
data1 <- cbind(rownames(data),data)
colnames(data1) <- c("gene_symbol", colnames(data))
data2 <- data1[-1,]    # 删除第一行label

gene_adjp = read.csv("gene_adjp.csv", header=TRUE, sep = ',')
colnames(gene_adjp) <- c("gene_adjp")

data3 <- merge(gene_adjp, data2, by.x="gene_adjp",by.y = "gene_symbol",all=FALSE) 
data4 <- rbind(data1[1,-1],data3[,-1])
row.names(data4) <- c("Label",as.character(data3[,1]))
# View(data4[1:10,1:10])
dim(data4)    # 518  20
# write.table(data4, file = "GSE69428_27651_outcome.txt",quote = F, sep = "\t")



# Data split --------------------------------------------------------------


datalabt <- data.frame(t(data4))
datalabt[2,2]

set.seed(2022)

library(caret)
library(dplyr)
trainingsamples <- datalabt$Label %>% createDataPartition(p = 0.7, list = FALSE)
## Delete label, obtain the matrix that RFE needs
datanum <- data4[-1,]
traindata  <- datanum[, trainingsamples]
testdata <- datanum[, -trainingsamples]
traindata[100,3]

# write.table(traindata, file = "matrix_DEtrain.txt",quote = F, sep = "\t")
# write.table(testdata, file = "matrix_DEtest.txt",quote = F, sep = "\t")


label <- data4[1,]
colnames(label) <- c("outcome")
trainLab <- as.numeric(label[trainingsamples])
testLab <- as.numeric(label[-trainingsamples])
# View(trainLab)
sum(trainLab)

# write.table(trainLab, file = "trainLab.txt",quote = F, sep = "\t")
# write.table(testLab, file = "testLab.txt",quote = F, sep = "\t")



## clear 
rm(list=ls())

## packages
library(caret)
library(lattice)
library(ggplot2)

setwd("D:\\E\\²©Ê¿\\R_³ÌÐò\\HGSOC\\DataGEOtogather")

eset <- read.table("matrix_DEtrain.txt",header = TRUE,sep = "\t")
View(eset[,1:10])
colnames(eset)

##################  NNet-rfe  ##################

Control <- rfeControl(functions = caretFuncs, method = "cv",
                      verbose = FALSE , returnResamp = "final")

trControl1 <- trainControl( method = "cv",
                            classProbs=TRUE,
                            summaryFunction = twoClassSummary)

# disease <- as.factor(c(rep("SPTB",69), rep("Term",160)))
lab <- read.table("trainLab.txt",header = TRUE,sep = "\t")
lab[which(lab[,1] == "1"),1] <- c("Normal")
lab[which(lab[,1] == "0"),1] <- c("HGSOC")
disease <- as.factor(c(rep("HGSOC",8), rep("Normal",6),
                           rep("HGSOC",14), rep("Normal",5),
                           rep("HGSOC",1)))


rf2 <- rfe(
  t(eset), disease, sizes = c(533),
  rfeControl = Control, trControl = trControl1, method = "nnet",
  tuneGrid = expand.grid(size = c(8), decay = c(0.1)),
  maxit = 30, MaxNWts = 100000
)
feature_sele <- rf2$optVariables
# write.table(feature_sele, file = "DataPython\\ranklist_NNETrfe.txt", quote=F, sep="\t")

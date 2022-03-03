## clear 
rm(list=ls())

## packages
library(caret)
library(lattice)
library(ggplot2)

setwd("D:\\E\\²©Ê¿\\R_³ÌÐò\\HGSOC\\Data")

eset <- read.table("matrix_DE.txt",header = TRUE,sep = "\t")

##################  NNet-rfe  ##################

Control <- rfeControl(functions = caretFuncs, method = "cv",
                      verbose = FALSE , returnResamp = "final")

trControl1 <- trainControl( method = "cv",
                            classProbs=TRUE,
                            summaryFunction = twoClassSummary)

disease <- as.factor(c(rep("HGSOC",10), rep("Normal",10)))

rf2 <- rfe(
  t(eset), disease, sizes = c(517),
  rfeControl = Control, trControl = trControl1, method = "nnet",
  tuneGrid = expand.grid(size = c(8), decay = c(0.1)),
  maxit = 30, MaxNWts = 100000
)
feature_sele <- rf2$optVariables
# write.table(feature_sele, file = "datanew\\ranklist_NNETrfe.txt", quote=F, sep="\t")

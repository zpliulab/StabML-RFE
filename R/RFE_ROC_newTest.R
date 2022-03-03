
## clear
rm(list = ls())

## packages
library(e1071)
library(limma)
library(ggplot2)
library(reshape2)
library(gmodels)


## load data
setwd("D:\\E\\≤© ø\\R_≥Ã–Ú\\HGSOC")
eset <- as.matrix(read.table("Data\\matrix_DE.txt",header = TRUE,sep = "\t"))


# ∂¡»Î≈≈–ÚŒƒº˛ ------------------------------------------------------------------
rfe_svm <- read.table("DataPython\\ranklist_svmrfe.txt",stringsAsFactors = FALSE)
colnames(rfe_svm) <- c("rank","chipID")
rfe_svm_dif<- head(rfe_svm,n = 40L)
eset_svm <- eset[rfe_svm_dif$chipID,]


# DE gene ----------------------------------------------------------------

library(limma)
data1 <- eset_svm
disease = read.table("data\\phe_zh.txt",header = TRUE, sep = "\t")
disease <- factor(disease[,"class"])
design <- model.matrix(~-1+disease)
contrast.matrix <- makeContrasts(contrasts = "diseaseNormal - diseaseHGSOC", levels = design)
fit_feature <- lmFit(data1,design)
fit1_feature <- contrasts.fit(fit_feature,contrast.matrix)
fit2_feature <- eBayes(fit1_feature)
dif <- topTable(fit2_feature,coef = "diseaseNormal - diseaseHGSOC",n = nrow(fit2_feature),lfc = log2(2))
dif0.05 <- dif[dif[,"adj.P.Val"] < 0.05,]    #  2089    6
dif_feature <- dif0.05[abs(dif0.05[,"logFC"]) > 1.5,]
# write.csv(dif_feature, file = "Datanew\\dif_feature_nb.csv")

# —µ¡∑---≤‚ ‘ -----------------------------------------------------------------
x <- t(eset[rownames(dif_feature)[1:20],])
gene_ref <- colnames(x)

intersect(gene_ref, rfe_svm_dif[1:40,2])

y1 <- as.matrix(c(rep(0,10), rep(1,10)))
disease <- as.factor(c(rep("HGSOC",10), rep("Normal",10)))
y <- disease


# Acc ---------------------------------------------------------------------
svm_acc <- list()
cost <- c(rep(10,4),rep(100,4))
gamma <- rep(c(0.1,0.01,0.001,0.0001),2)
parameters <- data.frame(cost,gamma)
for(k in 1:dim(parameters)[1]){
  svm_acc_tmp <- c()
  costtmp <- parameters[k,"cost"]
  gammatmp <- parameters[k,"gamma"]
  all_test_label <- c()
  all_pred_label <- c()
  
  set.seed(666) # for reproducing results
  rowIndices <- 1 : nrow(x) # prepare row indices
  sampleSize <- 0.70 * length(rowIndices) # training sample size
  trainingRows <- sample (rowIndices, sampleSize) # random sampling
  trainingData <- x[trainingRows, ] # training data
  testData <- x[-trainingRows, ] # test data
  trainingLabel <- y[trainingRows]
  testLabel <- y[-trainingRows]
  
  # trainingLabel <- y1[trainingRows]
  # testLabel <- y1[-trainingRows]

  svmfit <- svm (trainingData,trainingLabel, kernel = "radial", cost = costtmp, gamma=gammatmp, scale = FALSE) # radial svm, scaling turned OFF
  eset_pred<- predict(svmfit, testData)
  all_test_label <- c(all_test_label,as.vector(testLabel))
  all_pred_label <- c(all_pred_label,as.vector(eset_pred))

  svm_acc <- c(svm_acc,mean(all_test_label == all_pred_label))
}

library(gmodels) # CrossTable
CrossTable(x=all_test_label,y=all_pred_label, prop.chisq=FALSE)

parametertypes <- c()
for(k in 1:dim(parameters)[1]){
  costtmp <- parameters[k,"cost"]
  costtmp <- paste("cost:",costtmp,sep = "")
  gammatmp <- parameters[k,"gamma"]
  gammatmp <- paste("gamma:",gammatmp,sep = "")
  parametertmp <- paste(costtmp,gammatmp)
  parametertypes <- c(parametertypes,parametertmp)
}
# View(parametertypes)
names(svm_acc) <- parametertypes
svm_acc<- data.frame(svm_acc)
# View(svm_acc)
library(ggplot2)
library(reshape2)
svm_melt<- melt(svm_acc)
colnames(svm_melt) <- c("Parameter","Accuracy")
svm_melt$Accuracy <- round(svm_melt$Accuracy,3)
p <- ggplot(data = svm_melt,aes(x = Parameter,y = Accuracy,fill = Parameter))+
  geom_bar(stat = 'identity', width = 0.6, show.legend = F) +  # show.legend = NA, #Õº¿˝
  geom_text(aes(label = Accuracy),vjust=-0.3) +
  labs(title = "Accuracy of XGB-RFE in different parameter") +
  theme_bw() +
  # theme(axis.text.x = element_text(angle=30,size=10)) +
  theme(legend.text = element_text(color = 'black',size = 12,),
        axis.text = element_text(color = 'black',size = 15),
        axis.text.x = element_blank(),
        axis.title = element_text(color = 'black',size = 15,),
        axis.ticks = element_line(color = 'black'),
        axis.ticks.x = element_blank())
p

pdf(file = "datanew\\RFExgbSVMAcc.pdf",width = 5,height = 5)  # 7.5
p
dev.off()






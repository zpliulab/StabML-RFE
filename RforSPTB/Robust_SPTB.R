## 2022.1.23 计算稳定性的指标

rm(list = ls())

# install.packages("stabm")
# install.packages("ggdendro")
library(stabm)
library(ggdendro)

# my data -----------------------------------------------------------------

# load Rdata --------------------------------------------------------------
setwd('D:\\E\\博士\\R_程序\\HGSOC\\DataTcgaGtex') 
myfile = list.files("Robust")    
dir = paste("./Robust/", myfile, sep="")     #用paste命令构建路径变量dir
n = length(dir) 



library(stringr)
myname <- function(x){
  # x <- 1
  name <- dir[x]
  name1 <- str_split_fixed(name, "./", 2);
  name2 <- str_split_fixed(name1[2], "/", 2)
  name3 <- str_split_fixed(name2[2], "[_]", 2);
  name4 <- str_split_fixed(name3[2], "[_]", 2);
  name5 <- str_split_fixed(name4[2], "[.]", 2);
  name6 <- name5[1]
  return(name6)
}


methodlist <- c(myname(1), myname(2), myname(3), myname(4))

Biorfe <- as.matrix(read.csv(file = dir[1],header=TRUE, sep = ','))
ABrfe <- as.matrix(read.csv(file = dir[2],header=TRUE, sep = ','))
GBMrfe <- as.matrix(read.csv(file = dir[3],header=TRUE, sep = ','))
NNETrfe <- as.matrix(read.csv(file = dir[4],header=TRUE, sep = ','))


allmethod <- cbind(Biorfe, ABrfe, GBMrfe, NNETrfe)
colnames(allmethod) <- c("Bio-RFE", "AB-RFE", "GBDT-RFE", "NNET-RFE")




features <- list(Biorfe, ABrfe, GBMrfe, NNETrfe)
listname <- alist(Biorfe, ABrfe, GBMrfe, NNETrfe)
names(features) <- c(myname(1), myname(2), myname(3), myname(4))

listwithname <- function(features) {
  names(features) <- eval(substitute(listname))
  return(features)
}

listwithname(features)

all.feats = unique(unlist(features, use.names = FALSE))
all.feats
p <- length(all.feats)


# Selected in each repetition:
Reduce(intersect, features)


# Sorted selection frequency across all 30 repetitions:
sort(table(unlist(features)), decreasing = TRUE)






################# all possiable Combination including 4 subsets  ###############
n <- 4
k2 <- 2



## extract all combinations
A2 <- as.matrix(combn(n,k2))
colnames(A2) <- rep(1:choose(n, k2), 1)


## C^2_8 + C^3_8 + C^4_8 + ... + C^7_8： stability
var <- NULL
stability <- NULL
mycomb <- function(n,k,A){
  # k <- 8
  # A <- A8
  for (i in 1:choose(n, k)) {
    p <- length(unique(unlist(features[A[,i]], use.names = FALSE)))
    stab <- stabilityHamming(features[A[,i]], p)
    var <- c(var, i)
    stability <- c(stability, stab) 
  }
  return(stability)
}
Stability <- c(mycomb(n,k2,A2))
ll <- length(Stability)
Combination <- rep(1:ll,1)
plot(Combination, Stability, type='l')
l <- which(Stability == max(Stability))


# plot --------------------------------------------------------------------

# pdf(file = "Stability_GEOtogather.pdf",width = 15,height = 5)  # 7.5
par(cex=1.2)     # 坐标大小
plot(Combination, Stability, type='l', col=rgb(0,0,1), xaxt="n")  
axis(side=1, at=seq(1,ll,1), tck = -0.02, labels = Combination)
points(l, Stability[l], col='red',cex = 1.0, pch=18)
axis(3)
# dev.off()


A <- list(A2)
l <- which(Stability == max(Stability))
l
A[[6]][,7]
order(Stability, decreasing = T)
sort(Stability, decreasing = T)
Stability[l]


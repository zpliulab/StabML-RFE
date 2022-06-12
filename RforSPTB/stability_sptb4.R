## 2022.1.23 计算稳定性的指标

rm(list = ls())

# install.packages("stabm")
# install.packages("ggdendro")
library(stabm)
library(ggdendro)

# my data -----------------------------------------------------------------

# load Rdata --------------------------------------------------------------
setwd('D:\\E\\博士\\R_程序\\HGSOC\\DataTcgaGtex') 
myfile = list.files("DataPythonKauc")    
dir = paste("./DataPythonKauc/", myfile, sep="")     #用paste命令构建路径变量dir
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

ABrfe <- as.matrix(read.csv(file = dir[1],header=TRUE, sep = ','))
NBrfe <- as.matrix(read.csv(file = dir[2],header=TRUE, sep = ','))
NNETrfe <- as.matrix(read.csv(file = dir[3],header=TRUE, sep = ','))
SVMrfe <- as.matrix(read.csv(file = dir[4],header=TRUE, sep = ','))

allmethod <- cbind(ABrfe,  NBrfe, NNETrfe,SVMrfe)
colnames(allmethod) <- c("AB-RFE", "NB-RFE", "NNET-RFE", "XGB-RFE")
# write.csv(allmethod, file = "datanew\\AllMethodUpset.csv", row.names = F)
  

features <- list(ABrfe, NNETrfe, NBrfe, SVMrfe)
listname <- alist(ABrfe, NNETrfe, NBrfe, SVMrfe)
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
k3 <- k2+1
k4 <- k3+1


## extract all combinations
A2 <- as.matrix(combn(n,k2))
colnames(A2) <- rep(1:choose(n, k2), 1)
A3 <- as.matrix(combn(n,k3))
colnames(A3) <- rep((choose(n, k2)+1):(choose(n, k2) + choose(n, k3)), 1)
A4 <- as.matrix(combn(n,k4))
colnames(A4) <- rep((choose(n, k2)+choose(n, k3)+1):(choose(n, k2)+choose(n, k3) + choose(n, k4)), 1)

## save
getwd()
# write.csv(A2, "Datanew\\stabnum\\A2.csv", row.names=F)
# write.csv(A3, "Datanew\\stabnum\\A3.csv", row.names=F)
# write.csv(A4, "Datanew\\stabnum\\A4.csv", row.names=F)
# write.csv(A5, "Datanew\\stabnum\\A5.csv", row.names=F)
# write.csv(A6, "Datanew\\stabnum\\A6.csv", row.names=F)
# write.csv(A7, "Datanew\\stabnum\\A7.csv", row.names=F)
# write.csv(A8, "Datanew\\stabnum\\A8.csv", row.names=F)


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
# Stability <- c(mycomb(n,k2,A2), mycomb(n,k3,A3), mycomb(n,k4,A4), mycomb(n,k5,A5))
# Combination <- rep(1:56,1)
# Stability <- c(mycomb(n,k2,A2), mycomb(n,k3,A3), mycomb(n,k4,A4))
Stability <- c(mycomb(n,k2,A2), mycomb(n,k3,A3))
ll <- length(Stability)
Combination <- rep(1:ll,1)
plot(Combination, Stability, type='l')
l <- which(Stability == max(Stability))


# plot --------------------------------------------------------------------

# pdf(file = "Stability_SPTB4.pdf",width = 5,height = 5)  # 7.5
par(cex=1.2)     # 坐标大小
plot(Combination, Stability, type='l', col=rgb(0,0,1), xaxt="n")
axis(side=1, at=seq(1,ll,1), tck = -0.02, labels = Combination)
points(l, Stability[l], col='red',cex = 1.5, pch=18)
axis(3)
# dev.off()




# pdf(file = "datanew\\Stability_28.pdf",width = 6,height = 5)  # 7.5
lll <- choose(8,2)
par(cex=1.2)     # 坐标大小
plot(Combination[1:28], Stability[1:28], type='l', col=rgb(0,0,1), xaxt="n")  
axis(side=1, at=seq(1,lll,1), tck = -0.02, labels = Combination[1:lll])
points(5, Stability[5], col='red',cex = 1.5,pch=18)
axis(3)
# dev.off()

A <- list(A2,A3,A4)
l <- which(Stability == max(Stability))
l
A[[2]][,1]
order(Stability, decreasing = T)
sort(Stability, decreasing = T)
Stability[l]


# lab <- A[[1]][,5]
# lab <- c(1,2,6)
lab <- A[[2]][,1]

# Selected in each method
Reduce(intersect, features[lab])
bio <- Reduce(intersect, features[lab])
# write.csv(bio, file = "datanew\\bio.csv", row.names = F)
# Sorted selection frequency 
sort(table(unlist(features[lab])), decreasing = TRUE)
View(sort(table(unlist(features[lab])), decreasing = TRUE))
## The selection frequency can be visualized with the plotFeatures() function:
plotFeatures(features[lab])



# Sorted selection frequency across all 30 repetitions:
sort(table(unlist(features[lab])), decreasing = TRUE)
fren <- sort(table(unlist(features[lab])), decreasing = TRUE)
View(fren)
# write.csv(fren, file = "fren4.csv", row.names = F)



plotFeaturesmy = function(features, listname, sim.mat = NULL) {
  
  ## my add begin
  listwithname <- function(features){
    names(features) <- eval(substitute(listname))
    return(features)
  }
  ## my add end
  
  
  packages = c("ggplot2", "cowplot", "ggdendro")
  rn = lapply(packages, requireNamespace)
  
  # Checks
  checkmate::assertList(features, any.missing = FALSE, min.len = 2L,
                        types = c("integerish", "character"))
  type.character = sapply(features, is.character)
  if (any(type.character) && !all(type.character)) {
    stop("All features must numeric or all features must be character")
  }
  
  if (!is.null(sim.mat)) {
    pck = attr(class(sim.mat), "package")
    
    if (is.null(pck) || pck != "Matrix") {
      checkmate::assertMatrix(sim.mat, any.missing = FALSE, min.rows = 1L, min.cols = 1L, null.ok = FALSE)
      checkmate::assertTRUE(isSymmetric(unname(sim.mat)))
      checkmate::assertNumeric(sim.mat, lower = 0, upper = 1)
    } else {
      checkmate::assertTRUE(Matrix::isSymmetric(sim.mat))
      checkmate::assertNumeric(sim.mat@x, lower = 0, upper = 1)
    }
    
    if (any(type.character)) {
      checkmate::assertNames(colnames(sim.mat))
      rownames(sim.mat) = colnames(sim.mat)
      F.all = colnames(sim.mat)
    } else {
      F.all = seq_len(ncol(sim.mat))
    }
    
    lapply(features, function(f) {
      checkmate::assertVector(f, any.missing = FALSE, unique = TRUE, max.len = length(F.all))
      checkmate::assertSubset(f, F.all)
    })
    
  } else {
    lapply(features, function(f) {
      checkmate::assertVector(f, any.missing = FALSE, unique = TRUE)
    })
  }
  
  all.feats = unique(unlist(features, use.names = FALSE))
  
  if (length(all.feats) == 0) {
    stop("No feature selected in any set!")
  }
  
  mat = do.call(rbind, lapply(features, function(f) all.feats %in% f))
  colnames(mat) = NULL
  rownames(mat) = listname
  
  d.repls = dist(mat, method = "manhattan")
  hc.repls = hclust(d.repls, method = "average")
  o.repls = hc.repls$order
  dd.repls = as.dendrogram(hc.repls)
  
  if (length(all.feats) > 1) {
    if (is.null(sim.mat)) {
      d.feats = dist(t(mat), method = "manhattan")
    } else {
      d.feats = as.dist(1 - sim.mat[all.feats, all.feats])
    }
    hc.feats = hclust(d.feats, method = "average")
    o.feats = hc.feats$order
    dd.feats = as.dendrogram(hc.feats)
  } else {
    o.feats = 1L
  }
  
  mat = mat[o.repls, o.feats, drop = FALSE]
  colnames(mat) = paste0("V", all.feats[o.feats])
  ## the following: my add annotation
  # rownames(mat) = paste0("S", o.repls)      
  
  # this is a poor man's melt of mat
  mat.data = data.frame(
    repl = factor(rep(rownames(mat), ncol(mat)), levels = rownames(mat)),
    feature = factor(rep(colnames(mat), each = nrow(mat)), levels = colnames(mat)),
    selected = factor(ifelse(as.logical(mat), "Yes", "No"), levels = c("No", "Yes"))
  )
  
  # nchar
  max.char = max(nchar(all.feats))
  if (max.char > 2) {
    angle.feats = 90
  } else {
    angle.feats = 0
  }
  
  heat.plot = ggplot2::ggplot(mat.data, ggplot2::aes_string(x = "feature", y = "repl")) +
    ggplot2::geom_tile(ggplot2::aes_string(fill = "selected"), colour = "white") +
    # ggplot2::scale_fill_grey(name = "Selected", start = 0.9, end = 0.2, drop = FALSE) +
    ggplot2::scale_fill_hue(c=50, l= 80, name = "Selected", drop = FALSE) +
    ggplot2::theme_void() +
    ggplot2::labs(y = "Method", x = "Feature", title = "") +
    ggplot2::scale_y_discrete(expand = c(0, 0), labels = listname) +
    ggplot2::scale_x_discrete(expand = c(0, 0), labels = all.feats[o.feats]) +
    ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                   title = ggplot2::element_text(size = 1),
                   legend.position = "right",
                   legend.title = ggplot2::element_text(size = 12, color = 'black'),  # legend 标题
                   legend.text = ggplot2::element_text(size = 10),    # legend 内容
                   axis.title = ggplot2::element_text(size = 12, color = 'black'), #横纵标签
                   axis.title.y = ggplot2::element_text(angle = 90),
                   axis.text = ggplot2::element_text(size = 10),
                   axis.text.x = ggplot2::element_text(angle = angle.feats, 
                                                       hjust = 1, vjust = 0.5),
                   axis.text.y = ggplot2::element_text(angle = angle.feats))
  
  
  final.plot = heat.plot
  
  dendro.data.repls = ggdendro::dendro_data(dd.repls, type = "rectangle")
  dendro.repls = cowplot::axis_canvas(heat.plot, axis = "y", coord_flip = TRUE) +
    ggplot2::geom_segment(data = ggdendro::segment(dendro.data.repls),
                          ggplot2::aes_string(y = "y", x = "x", xend = "xend", yend = "yend"), size = 0.5) +
    ggplot2::coord_flip() +
    ggplot2::theme(plot.margin = ggplot2::unit(c(0, 1, 0, 0), "lines"))
  final.plot = cowplot::insert_yaxis_grob(final.plot, dendro.repls,
                                          grid::unit(0.2, "null"), position = "right")
  
  if (length(all.feats) > 1) {
    dendro.data.feats = ggdendro::dendro_data(dd.feats, type = "rectangle")
    dendro.feats = cowplot::axis_canvas(heat.plot, axis = "x") +
      ggplot2::geom_segment(data = ggdendro::segment(dendro.data.feats),
                            ggplot2::aes_string(x = "x", y = "y", xend = "xend", yend = "yend"), size = 0.5) +
      ggplot2::theme(plot.margin = ggplot2::unit(c(1, 0, 0, 0), "lines"))
    final.plot = cowplot::insert_xaxis_grob(final.plot, dendro.feats,
                                            grid::unit(0.2, "null"), position = "top")
  }
  cowplot::ggdraw(final.plot)
}


pdf(file = "stability_overlap3.pdf",width = 11, height = 5)  # 7.5  13
# source("D:\\E\\博士\\R_程序\\HGSOC\\R\\plotFeaturesmy.R")
fea <- features[lab]
listname <- c("AB-RFE", "DBGT-RFE", "NNET-RF")
plotFeaturesmy(fea, listname)
dev.off()


## 2022.2.14 以下没用










# similarity matrix
# data1 = read.table("GSE69428_scale.txt", header = T, check.names = FALSE)
data1 = read.table("GSE69428_outcome.txt", header = T, check.names = FALSE)
# View(data1[,1:10])
data2 <- t(data1[all.feats,])
# View(data2[,1:10])
simmat = abs(cor(data2))
# View(simmat[,1:10])


## each of the 30 feature sets includes either medv or cmedv:
plotFeatures(features, sim.mat = simmat)

## use one of the stability measures
stabilityIntersectionCount(features, sim.mat = simmat, N = 1000)
stabilityZucknick(features = features, sim.mat = simmat)


## considering identity matrix as similarity matrix and neglecting all similarities.
## Without taking into account the feature similarities, the stability value is:
nosimmat = diag(nrow(simmat))
colnames(nosimmat) = row.names(nosimmat) = colnames(simmat)
stabilityIntersectionCount(features, sim.mat = nosimmat)
stabilityZucknick(features = features, sim.mat = nosimmat)






## 图片的纵坐标有问题，修正了一下，并且调整了颜色，调整了字号
source("D:\\E\\博士\\R_程序\\HGSOC\\R\\plotFeaturesmy.R")
listStabilityMeasures()
plotFeaturesmy(features, listname, simmat)







# code in function --------------------------------------------------------

listwithname <- function(features) {
  names(features) <- eval(substitute(listname))
  return(features)
}

listwithname(features)

all.feats = unique(unlist(features, use.names = FALSE))
p <- length(all.feats)


# mat matrix --------------------------------------------------------------
mat = do.call(rbind, lapply(features, function(f) all.feats %in% f))
colnames(mat) = all.feats
rownames(mat) = listname
dim(mat)    # 5 107
dim(t(mat))


# repls -------------------------------------------------------------------
d.repls = dist(mat, method = "manhattan")
hc.repls = hclust(d.repls, method = "average")
o.repls = hc.repls$order
dd.repls = as.dendrogram(hc.repls)


# feats -------------------------------------------------------------------
d.feats = dist(t(mat), method = "manhattan")
View(d.feats)
hc.feats = hclust(d.feats, method = "average")
o.feats = hc.feats$order
dd.feats = as.dendrogram(hc.feats)


# mat cluster -------------------------------------------------------------
mat = mat[o.repls, o.feats, drop = FALSE]
# colnames(mat) = paste0("V", all.feats[o.feats])
# rownames(mat0) = paste0("S", o.repls)


# a poor man's melt of mat ------------------------------------------------
mat.data = data.frame(
  repl = factor(rep(rownames(mat), ncol(mat)), levels = rownames(mat)),
  feature = factor(rep(colnames(mat), each = nrow(mat)), levels = colnames(mat)),
  selected = factor(ifelse(as.logical(mat), "Yes", "No"), levels = c("No", "Yes"))
)

# nchar
max.char = max(nchar(all.feats))
if (max.char > 2) {
  angle.feats = 90
} else {
  angle.feats = 0
}


# 画热图 ---------------------------------------------------------------------
heat.plot = ggplot2::ggplot(mat.data, ggplot2::aes_string(x = "feature", y = "repl")) +
  ggplot2::geom_tile(ggplot2::aes_string(fill = "selected"), colour = "white") +
  # ggplot2::scale_fill_grey(name = "Selected", start = 0.9, end = 0.2, drop = FALSE) +
  ggplot2::scale_fill_hue(c=50, l= 80, name = "Selected", drop = FALSE) +
  ggplot2::theme_void() +
  ggplot2::labs(y = "Method", x = "Feature", title = "") +
  ggplot2::scale_y_discrete(expand = c(0, 0), labels = listname) +
  ggplot2::scale_x_discrete(expand = c(0, 0), labels = all.feats[o.feats]) +
  ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                 title = ggplot2::element_text(size = 1),
                 legend.position = "right",
                 legend.title = ggplot2::element_text(size = 14, color = 'black',
                                                      family = 'Times', face = 'plain'),  # legend 标题
                 legend.text = ggplot2::element_text(size = 10, family = 'Times'),    # legend 内容
                 axis.title = ggplot2::element_text(size = 14, color = 'black',
                                                    family = 'Times', face = 'plain'), #横纵标签
                 axis.title.y = ggplot2::element_text(angle = 90),
                 axis.text = ggplot2::element_text(size = 7, family = 'Times'),
                 axis.text.x = ggplot2::element_text(angle = angle.feats, 
                                                     hjust = 1, vjust = 0.5))


heat.plot

# pdf(file = "stability.pdf",width = 10, height = 5)
# heat.plot
# dev.off()




# 添加聚类线条 ------------------------------------------------------------------
final.plot = heat.plot

dendro.data.repls = ggdendro::dendro_data(dd.repls, type = "rectangle")
dendro.repls = cowplot::axis_canvas(heat.plot, axis = "y", coord_flip = TRUE) +
  ggplot2::geom_segment(data = ggdendro::segment(dendro.data.repls),
                        ggplot2::aes_string(y = "y", x = "x", xend = "xend", yend = "yend"), size = 0.5) +
  ggplot2::coord_flip() +
  ggplot2::theme(plot.margin = ggplot2::unit(c(0, 1, 0, 0), "lines"))
final.plot = cowplot::insert_yaxis_grob(final.plot, dendro.repls,
                                        grid::unit(0.2, "null"), position = "right")


dendro.data.feats = ggdendro::dendro_data(dd.feats, type = "rectangle")
dendro.feats = cowplot::axis_canvas(heat.plot, axis = "x") +
  ggplot2::geom_segment(data = ggdendro::segment(dendro.data.feats),
                        ggplot2::aes_string(x = "x", y = "y", xend = "xend", yend = "yend"), size = 0.5) +
  ggplot2::theme(plot.margin = ggplot2::unit(c(1, 0, 0, 0), "lines"))
final.plot = cowplot::insert_xaxis_grob(final.plot, dendro.feats,
                                        grid::unit(0.2, "null"), position = "top")

cowplot::ggdraw(final.plot)

# code in function end ------------------------------------------------------------


listStabilityMeasures()
aa <- as.matrix(listStabilityMeasures())
# write.csv(aa, file = "StabilityResult.csv", row.names = F)


feats <- features
# p: Total number of features in the datasets.
mat_matrix = 0.92 ^ abs(outer(1:p, 1:p, "-"))
colnames(mat_matrix) <- c(all.feats[o.feats])
rownames(mat_marix) <- colnames(mat_matrix)
# View(mat_marix[1:p,80:p])




feats <- features
mat_matrix <- simmat

# index only features -----------------------------------------------------
## Stability Measure Dice
stabilityDice(features = feats)

## Stability Measure Jaccard
stabilityJaccard(features = feats)

## Stability Measure Novoviˇcová
stabilityNovovicova(features = feats)

## Stability Measure Ochiai
stabilityOchiai(features = feats)




# index feature-num -------------------------------------------------------

## Stability Measure Davis
stabilityDavis(features = feats, p)

## Stability Measure Hamming
stabilityHamming(features = feats, p)

## Stability Measure Kappa
stabilityKappa(features = feats, p)

## Stability Measure Lustgarten
stabilityLustgarten(features = feats, p)

## Stability Measure Nogueira
stabilityNogueira(features = feats, p)

## Stability Measure Phi
stabilityPhi(features = feats, p)

## Stability Measure Somol
stabilitySomol(features = feats, p)

## Stability Measure Unadjusted
stabilityUnadjusted(features = feats, p)

## Stability Measure Wald
stabilityWald(features = feats, p)



# index feature-matrix -------------------------------------------------------
## Stability Measure Sechidis
stabilitySechidis(features = feats, sim.mat = mat_matrix)

## Stability Measure Zucknick
stabilityZucknick(features = feats, sim.mat = mat_matrix)
stabilityZucknick(features = feats, sim.mat = nosimmat)

## Stability Measure Adjusted Intersection Count
stabilityIntersectionCount(features = feats, sim.mat = mat_matrix, N = 1000)

## Stability Measure Adjusted Intersection Greedy
stabilityIntersectionGreedy(features = feats, sim.mat = mat_matrix, N = 1000)

## Stability Measure Adjusted Intersection MBM
stabilityIntersectionMBM(features = feats, sim.mat = mat_matrix, N = 1000)

## Stability Measure Adjusted Intersection Mean
stabilityIntersectionMean(features = feats, sim.mat = mat_matrix, N = 1000)

## Stability Measure Yu
stabilityYu(features = feats, sim.mat = mat_matrix, N = 1000)












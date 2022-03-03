## 2022.2.14 corrected

## clear
rm(list = ls())

## package
library(clusterProfiler)
library(topGO)
library(Rgraphviz)
library(carData)
library(org.Hs.eg.db)
library(GOplot)
library(stringr)

# load data --------------------------------------------------------------------
setwd("D:\\E\\博士\\R_程序\\HGSOC")
gene <- as.matrix(read.csv("Datanew\\bio.csv", header=TRUE, sep = ','))
colnames(gene) <- c('gene')


# genes symbol ----------------------------------------------------------
genelist <- as.character(gene)
eg <- bitr(genelist, 
           fromType="SYMBOL", 
           toType=c("ENTREZID","GENENAME"), 
           OrgDb="org.Hs.eg.db"); 
head(eg)

library(stargazer)
stargazer(as.matrix(cbind(eg$SYMBOL, eg$ENTREZID, eg$GENENAME)))

# go ----------------------------------------------------------------------
geneList <- eg$ENTREZID
go <- enrichGO(gene = geneList,
               OrgDb = org.Hs.eg.db,
               ont='ALL',
               pAdjustMethod = 'BH',
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               keyType = 'ENTREZID')
head(go)
# write.csv(go@result, file = "datanew\\go_BP.csv", row.names = F)   ## use

# 进行简单的可视化 ----------------------------------------------------------------
# barplot(go_BP, showCategory=10,drop=T, title=NULL)    
dotplot(go, showCategory=10) #泡泡图
# plotGOgraph(go_BP) 	#GO图，看不清楚可以尝试左上角另存为pdf
library(ggnewscale)

## 富集到的GO terms之间的基因重叠关系
library(enrichplot)
x2 <- pairwise_termsim(go)
emapplot(x2, showCategory = 10)


enrichGO = DOSE::setReadable(go_BP, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
enrichGO
library(ggnewscale)
cnetplot(enrichGO)



# KEGG --------------------------------------------------------------------

enrichKK <- enrichKEGG(gene         =  geneList,
                       organism     = 'hsa',
                       #universe     = gene_all,
                       pvalueCutoff = 0.2,
                       qvalueCutoff = 0.2)
head(enrichKK)[,1:6] 
dotplot(enrichKK)

barplot(enrichKK,showCategory=20)
dotplot(enrichKK)

#(3)展示top5通路的共同基因，要放大看。
#Gene-Concept Network 
# install.packages("ggnewscale")
library(ggnewscale)
cnetplot(enrichKK)


## 2022.2.17 add
kegg <- enrichKEGG(gene = geneList, 
                   organism = 'hsa', 
                   keyType = 'kegg', 
                   pvalueCutoff = 0.1, 
                   pAdjustMethod = 'BH', 
                   minGSSize = 3, 
                   maxGSSize = 500, 
                   qvalueCutoff = 0.2, 
                   use_internal_data = FALSE)
head(kegg)
barplot(kegg,showCategory=10,drop=T)

#pathway映射
browseKEGG(kegg, "hsa00770") 

## gene GO relationship
library(enrichplot)
library(DOSE)
x2 <- go

cnetplot(x2)
# use `layout` to change the layout of map
cnetplot(x2, layout = "star")
# use `showCategory` to select the displayed terms. It can be a number of a vector of terms.
cnetplot(x2, showCategory = 7)

cnetplot(x2, categorySize="pvalue", foldChange=genelist)

categorys <- c( "nuclear division",
                "sister chromatid segregation",
                "organelle fission",
                "water-soluble vitamin metabolic process",
                "nuclear chromosome segregation",
                "chromosome separation",
                "mitotic nuclear division",
                "chromosome segregation",
                "vitamin metabolic process",
                "meiotic chromosome separation")

x <- read.csv("Datanew\\biologFC.csv", header=TRUE, sep = ',')
foldchange <- x$logFC
names(foldchange) <- x$bio

plot <-  cnetplot(go_BP,
                  showCategory = categorys,
                  # foldChange = NULL,
                  foldChange =foldchange, 
                  layout = "kk",    # kk, gem
                  colorEdge = T,
                  circular = F,
                  node_label = "all", # "gene",
                  cex_category = 1.0,
                  cex_gene = 1.0,
                  cex_label_category = 0.6,
                  cex_label_gene = 0.8,
                  shadowtext = "none")
plot
# pdf(file = "datapython20\\cnetplot_cluster.pdf", width = 6.5, height = 4.5)
# plot
# dev.off()





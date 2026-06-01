rm(list = ls())
library(TCGAbiolinks)
library(readr)
library(readxl)
library(tidyverse)
library(lance)
setwd("/data/nas1/chenpeiru/44_GYZK-30212-7/09_GSEA")
load("/data/nas1/chenpeiru/44_GYZK-30212-7/00_raw_data/01_TCGA/dat.tcga_count425.rda")
dat <- dat.final
group <- read.delim2('/data/nas1/chenpeiru/44_GYZK-30212-7/01_TCGA_DEG/group.xls')
#group$sample<-gsub('.','-',group$sample,fixed = T)
table(group$group)
tumor.sample <- group[which(group$group=='Tumor'),]

dat <- dat[,tumor.sample$sample]
dat<-round(dat,digits = 0)

#deseq2-------
colData <- read.table('/data/nas1/chenpeiru/44_GYZK-30212-7/04_lasso/risk.train.group.txt',
                      sep = "\t",header = T)
head(colData)
colData <- colData[colData$id %in% tumor.sample$sample,] 

library(DESeq2)
#colData$group<-ifelse(colData$group=='High CCNB1','High_CCNB1','Low_CCNB1')
#table(colData$group)
head(colData)
head(dat)

colData$group<-factor(colData$group,levels = c('high','low'))
dat <- dat[,colData$id]
head(dat)
dds<-DESeqDataSetFromMatrix(countData = dat,colData=colData,design = ~group)
dds = dds[rownames(counts(dds)) > 1,]
dds<-DESeq(dds)
# dds
res =results(dds, contrast = c("group","high","low"))
res =res[order(res$padj),]
head(res)
summary(res)
table(res$padj<0.05)
allGeneSets<-as.data.frame(res)
allGeneSets<-na.omit(allGeneSets)
logFCcutoff <- 0.5
allGeneSets$change = as.factor(
  ifelse(allGeneSets$pvalue < 0.05 & abs(allGeneSets$log2FoldChange) > logFCcutoff,
         ifelse(allGeneSets$log2FoldChange > logFCcutoff,'UP','DOWN'),'NOT')
)

DEGeneSets <- subset(allGeneSets,
                     allGeneSets$pvalue < 0.05 & abs(allGeneSets$log2FoldChange) > logFCcutoff)#733
genelist <- allGeneSets$log2FoldChange
names(genelist) <- rownames(allGeneSets)
geneList <- sort(genelist, decreasing = T)
DEGeneSets <- DEGeneSets[order(DEGeneSets$pvalue),]
dim(DEGeneSets)# 5824    7
## GSEA KEGG
library(clusterProfiler)
library(enrichplot)
#setwd("/data/nas1/chenpeiru/41_YQNC-10615-6/10_GSEA")
#KEGG--------
kegg_set<- read.gmt("c2.cp.kegg.v7.5.1.symbols.gmt")
set.seed(1)
kegg_gsea <- GSEA(geneList, TERM2GENE = kegg_set, pvalueCutoff = 1)
kegg_result <- kegg_gsea@result
dim(kegg_result)#184，11
kegg_result <- subset(kegg_result,pvalue<0.05)
write.table(kegg_result,file = 'GSEA_KEGG_pvalue0.05.xls',sep = '\t',quote = F,row.names = F)
pdf(file = '01.KEGG_GSEA.pdf',w=12,h=8)
par(family="serif")
gseaplot2(kegg_gsea,c(1:5),color = c('#7B68EE','#CD3333','#20B2AA','#FF8C00','#FF6666'),
          title = 'KEGG GSEA',
          base_size = 11,
          rel_heights = c(1.5, 0.3, 0.5))
dev.off()
#GO--------
kegg_set<- read.gmt("c5.go.v7.5.1.symbols.gmt")
set.seed(1)
kegg_gsea <- GSEA(geneList, TERM2GENE = kegg_set, pvalueCutoff = 1)
kegg_result <- kegg_gsea@result
dim(kegg_result)#184，11
kegg_result <- subset(kegg_result,pvalue<0.05)
write.table(kegg_result,file = 'GSEA_GO_pvalue0.05.xls',sep = '\t',quote = F,row.names = F)
pdf(file = '02.GO_GSEA.pdf',w=12,h=8)
par(family="serif")
gseaplot2(kegg_gsea,c(1:5),color = c('#7B68EE','#CD3333','#20B2AA','#FF8C00','#FF6666'),
          title = 'GO GSEA',
          base_size = 11,
          rel_heights = c(1.5, 0.3, 0.5))
dev.off()
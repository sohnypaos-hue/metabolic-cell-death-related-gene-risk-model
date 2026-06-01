#venn------
library(dplyr)
rm(list = ls())

setwd("/data/nas1/chenpeiru/44_GYZK-30212-7/02_venn")
##VENN-------
library(readxl)
library(VennDiagram)

#when logfc=1, get venn--------- 
#rm(list = ls())
DEGs1 <- read.csv('/data/nas1/chenpeiru/44_GYZK-30212-7/01_TCGA_DEG/DEG_logfc1_padj0.05.csv',
                  sep = ',',header = T,row.names = 1)
head(DEGs1)
# "HIST1H3J" %in% rownames(DEGs1)
MCD_RGs <- read.csv("/data/nas1/chenpeiru/44_GYZK-30212-7/00_raw_data/03_third_genes/MCD_RGs.csv",
                    sep=',',header=T) %>% as.data.frame()
head(MCD_RGs)
#DEGs2
# DEGs2 <- read.table('/data/nas1/chenpeiru/43_JNZK-20704-5/03_train_DEGs/diffSig_mRNA_logFC1_adj.P.Val0.05.txt',
#                     sep = '\t',header = T,row.names = 1)

s3 <- list(
  DEGs1 = rownames(DEGs1),
  MCD_RGs = MCD_RGs$MCD_RGs
  
)

v2 <- venn.diagram(x = s3, filename = NULL,
                   scaled = F, # 根据比例显示大小
                   alpha= 0.4, #透明度
                   lwd=1,lty=1,col=c('#F4A460','#6A5ACD'), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
                   label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
                   cex = 2, # 数字大小
                   fontface = "bold",  # 字体粗细；加粗bold
                   fill=c('#F4A460','#6A5ACD'), # 填充色 配色https://www.58pic.com/
                   category.names = c("DEGs1", "MCD_RGs") , #标签名
                   cat.dist = 0.05, # 标签距离圆圈的远近
                   print.mode = c("percent", "raw"),#显示百分比
                   cat.pos = c(-40, -320), # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
                   cat.cex = 2, #标签字体大小
                   cat.fontface = "bold",  # 标签字体加粗
                   cat.col='black' ,   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
                   cat.default.pos = "outer",  # 标签位置, outer内;text 外
                   output=TRUE
)
# 
cowplot::plot_grid(v2)
setwd("/data/nas1/chenpeiru/44_GYZK-30212-7/02_venn")
pdf("DE-MCD-50.pdf",height =  9,width = 9)
par(family="serif")
cowplot::plot_grid(v2)
dev.off()

inter <- get.venn.partitions(s3)
interset_mRNA<-as.data.frame(inter$..values..[1])

#for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.csv(interset_mRNA,'DE-MCD-50.csv', row.names = FALSE, quote = FALSE)
#go kegg--------
rm(list=ls())
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(stringr)
library(tidyr)
##GO term enrichment
data <- read.csv("DE-MCD-50.csv",
                 header = T,sep = ',')
Gene = bitr(data$X1,fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")  #将SYMBOL转换为ENTREZID
Gene <- na.omit(Gene)  #去除未能成功转换的
head(Gene)
go.enrich = enrichGO(gene =Gene$ENTREZID,
                     OrgDb = org.Hs.eg.db,
                     keyType = "ENTREZID",
                     ont = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.1,
                     qvalueCutoff = 1,
                     readable = TRUE)

go_result = data.frame(go.enrich)
#setwd("/data/nas1/chenpeiru/43_JNZK-20704-5/06_GO_kegg")
head(go_result)
write.csv(go_result, '01_GO_enrichment.csv', quote = F, row.names = F)

#筛选Top10
go_result$term = paste(go_result$ID, go_result$Description, sep = ': ') #将ID与Description合并成新的一列
go_result$term = factor(go_result$term, levels = go_result$term,ordered = T) #转成因子，防止重新排列

go_result_BP = go_result[go_result$ONTOLOGY == 'BP', ]
go_result_MF = go_result[go_result$ONTOLOGY == 'MF', ]
go_result_CC = go_result[go_result$ONTOLOGY == 'CC', ]
table(go_result$ONTOLOGY)
go_result_Top10 = rbind(go_result_BP[1:10,], go_result_CC[1:10,],go_result_MF[1:10,])

#纵向柱状图——-根据pvalue值绘制#
library(ggplot2)
p1 = ggplot(go_result_Top10,aes(y=term,x=Count,fill=pvalue))+  #x、y轴定义；根据pvalue填充颜色
  geom_bar(stat = "identity",width=0.8)+ #柱状图宽度设置
  scale_fill_gradient(low = "#FDB338",high ="#1065AB" )+
  labs(title = "GO Term Enrich (Top10)",  #设置标题、x轴和Y轴名称
       x = "Gene number", 
       y = "GO Terms")+
  theme(axis.title.x = element_text(face = "bold",size = 16),
        axis.title.y = element_text(face = "bold",size = 16),
        legend.title = element_text(face = "bold",size = 16))+
  theme_bw()
head(go_result_Top10)
pdf('02_GO_enrichment_Top10.pdf', 15, 8)
#根据ONTOLOGY分类信息添加分组框
p1+facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')
dev.off()
# KEGG ANALYSIS--------
##1.2 KEGG富集分析############
kegg <- enrichKEGG(gene=Gene$ENTREZID,
                   organism = 'hsa',
                   keyType = 'kegg',
                   pvalueCutoff =0.05,
                   qvalueCutoff=1,
                   pAdjustMethod = 'BH'
)

kegg <- setReadable(kegg, 'org.Hs.eg.db', 'ENTREZID')
head(kegg,10)
library(GOplot)
library(dplyr)
library(stringr)
kegg_all <- kegg@result
head(kegg_all)
#setwd("/data/nas1/chenpeiru/43_JNZK-20704-5/06_GO_kegg")
write.csv(as.data.frame(kegg_all),"02_kegg-enrich.csv",row.names =F)
head(kegg_all,10)

#kegg analysis-------

kegg_all$Description = factor(kegg_all$Description, levels = kegg_all$Description,ordered = T) #转成因子，防止重新排列
head(kegg_all)
kegg.enrich.Top10 = kegg_all[1:10, ]

#纵向柱状图——-根据pvalue值绘制#
p2 = ggplot(kegg.enrich.Top10,aes(y=Description,x=Count,fill=pvalue))+  #x、y轴定义；根据pvalue填充颜色
  geom_bar(stat = "identity",width=0.8)+ #柱状图宽度设置
  scale_fill_gradient(high = "#FDB338",low ="#1065AB" )+
  labs(title = "KEGG Pathway Enrich (Top10)",  #设置标题、x轴和Y轴名称
       x = "Gene number",
       y = "KEGG pathways")+
  theme(axis.title.x = element_text(face = "bold",size = 16),
        axis.title.y = element_text(face = "bold",size = 16),
        legend.title = element_text(face = "bold",size = 16))+
  theme_bw()
p2
#setwd("/data/nas1/chenpeiru/43_JNZK-20704-5/06_GO_kegg")
pdf('04_KEGG_enrichment_Top10.pdf', 8, 8)
p2
dev.off()

png('04_KEGG_enrichment_Top10.png', 800, 800)
p2
dev.off()
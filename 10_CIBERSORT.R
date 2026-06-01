#免疫浸润分析 ssGSEA####
rm(list = ls())
#source('Cibersort_source.R')
library(RColorBrewer)
library(tidyr)
#library(radiant.data)
library(aplot)
library(dplyr)
library(e1071)
library(tidyverse)
library(reshape2)
library(psych)
library(corrplot)
library(cowplot)
library(ggcorrplot)
library(msigdbr)
library(GSEABase)
library(GSVA)
setwd("/data/nas1/chenpeiru/44_GYZK-30212-7/10_CIBERSORT")
load("/data/nas1/chenpeiru/44_GYZK-30212-7/00_raw_data/01_TCGA/dat_fpkm425.rda")
dat <- dat_fpkm

group <- read.delim2('/data/nas1/chenpeiru/44_GYZK-30212-7/01_TCGA_DEG/group.xls')
tumor.sample <- group[which(group$group=='Tumor'),]

head(group)
#group$sample<-gsub('.','-',group$sample,fixed = T)
dat <- dat[,tumor.sample$sample]#399
# dat<-read.csv('TCGA_fpkm_symbol.csv',row.names = 1,sep = ",",header = T)%>%lc.tableToNum()

ssGSEA_data <- dat
gmtFile="mmc3.gmt"
#分析免疫细胞含量
#载入背景基因
geneSet=getGmt(gmtFile,
               geneIdType=SymbolIdentifier())
#开始ssGSEA分析,数据一定要是矩阵，不是数据框
ssgseaScore=gsva(as.matrix(ssGSEA_data), geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)

range(ssgseaScore)
ssgseaOut <- as.data.frame(t(ssgseaScore))

#整理数据格式
ssgsea_result1 = ssgseaOut

write.csv(ssgsea_result1,"ssGSEA_result.csv",row.names = T,quote = F)

ssgsea_result1 <- read.csv("/data/nas1/chenpeiru/44_GYZK-30212-7/10_CIBERSORT/ssGSEA_result.csv",
                           header = T,row.names =1)
group10588 <- read.table("/data/nas1/chenpeiru/44_GYZK-30212-7/04_lasso/risk.train.group.txt",
                         sep = "\t",header = T) %>% lc.tableToNum()
head(group10588)
group10588 <- group10588[,c(1,13)]
colnames(group10588) <- c("Sample","group")
head(group10588)
head(ssgsea_result1)
ssgsea_result1 <- data.frame(rownames(ssgsea_result1),ssgsea_result1)
head(ssgsea_result1)
colnames(ssgsea_result1)[1] <- "Sample"
ssgsea_result1 <- merge(group10588,ssgsea_result1,by="Sample")
head(ssgsea_result1)
rownames(ssgsea_result1) <- ssgsea_result1[,1]
head(ssgsea_result1)

ssgsea_result1 <- ssgsea_result1[,-1]
head(ssgsea_result1)
ssgsea_result1 <- ssgsea_result1[order(ssgsea_result1$group, decreasing = TRUE), ]
re1 <- ssgsea_result1
head(re1)
#去除分组信息
colnames(re1)[1] <- "Type"
re2 = re1[,-1]       #去除第一列分组信息
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
#提取数据，多行变成多列，要多学习‘tidyr’里面的三个函数
dat_cell <- re2 %>% as.data.frame() %>%rownames_to_column("Sample") %>%gather(key = Cell_type,value = Proportion,-Sample)
#提取数据
dat_group = gather(re1,Cell_type,Proportion,-Type )
#合并分组
dat = cbind(dat_cell,dat_group$Type)
colnames(dat)[4] <- "Type"
##2.3柱状图############
# p1 <- ggplot(dat,aes(Sample,Proportion,fill = Cell_type)) +
#   geom_bar(stat = "identity") +
#   labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") +
#   theme_bw() +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         legend.position = "bottom") +
#   scale_y_continuous(expand = c(0.01,0)) +
#   scale_fill_manual(values = mypalette(28))
# p1
# dev.off()
# ggsave('barchart_ImmuneGeme_ssGSEA.pdf', plot = p1,width=15,height = 6)
#热图####
annotation_col <- data.frame(sample = c(rep("low", 158), rep("high",241)))
rownames(annotation_col) <- rownames(re2)
pdf('heatmap_ImmuneGeme.pdf',height = 9,width = 12)
pheatmap(t(re2),
         cluster_rows = T,
         cluster_cols = F,
         annotation_col = annotation_col,
         show_colnames = F,
         name = " ")
dev.off()

#筛选差异免疫浸润细胞 ####
box=dat
theme_zg2 <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          panel.background=element_rect(fill='transparent',color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),#去网格线
          axis.line = element_line(colour = "black"),
          #axis.title.x = element_blank(),#去x轴标签
          axis.title.y=element_text(face = "bold",size = 14),#y轴标签加粗及字体大小
          axis.title.x=element_text(face = "bold",size = 14),#X轴标签加粗及字体大小
          axis.text.y = element_text(face = "bold",size = 12),#y坐标轴刻度标签加粗
          axis.text.x = element_text(face = "bold",size = 10, vjust = 1, hjust = 1, angle = 45),#x坐标轴刻度标签加粗
          axis.ticks = element_line(color='black'),
          # axis.ticks.margin = unit(0.8,"lines"),
          legend.title=element_blank(),
          legend.position=c(0.9, 0.85),#图例在绘图区域的位置
          # legend.position="top",
          legend.direction = "horizontal",
          legend.text = element_text(face = "bold",size = 12),
          # legend.background = element_rect( linetype="solid",colour ="black")
    )
}
e1 <- ggplot(box,aes(x=Cell_type,y=Proportion),palette = "jco", add = "jitter")+
  geom_boxplot(aes(fill=Type),width=0.6)+
  labs(x = "Cell Type", y = "Estimated Proportion")+
  scale_fill_manual(values = c("#1874CD","#FED82F")) +
  theme_zg2() + stat_compare_means(aes(group = Type),label = "p.signif",method = 'wilcox.test')
e1

#小提琴图
e1 <- ggplot(box,aes(x=Cell_type,y=Proportion),palette = "jco", add = "jitter")+
  geom_violin(aes(fill=Type),position=position_dodge(0.5),width=0.6)+
  labs(x = "Cell Type", y = "Immune cell content")+
  scale_fill_manual(values = c("red","blue")) +
  theme_zg2()
e1 = e1 + stat_compare_means(aes(group = Type),label = "p.signif")
e1
#保存图片
ggsave('box_ImmuneGeme_ssGSEA.pdf', plot = e1,width=15,height = 6)

#免疫细胞相关性分析 ####
head(ssgsea_result1)
DEcells_spearman <- ssgsea_result1[,-1]
colnames(DEcells_spearman)
#删除没有丰度的细胞#1,2,3,5,10,11,12,13,15,16,18,19,20,23,24,25,28
DEcells_spearman <- DEcells_spearman[,-c(6,27)]
colnames(DEcells_spearman)
#DEcells_cor <- cor(DEcells_spearman,method = "spearman") 
#class(DEcells_spearman$Activated.CD8.T.cell)
#class(DEcells_spearman$Activated.CD4.T.cell)

ciber.res<-DEcells_spearman#[,colSums(DEcells_spearman)>0]
#ciber.res<-t(ciber.res)
cor_ciber<-cor(ciber.res,method = 'spearman')
cor_ciber<-round(cor(ciber.res),2)
pdf("03_cor.plot_DEcells_spearman.pdf",width=13,height=9) #保存图为pdf
col1=colorRampPalette(colors =c("blue","white","red"),space="Lab")
testRes = cor.mtest(cor_ciber, method="spearman",conf.level = 0.95)
corrplot(cor_ciber, 
         p.mat = testRes$p, diag = T, type = 'upper',col = col1(10),
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.2,
         insig = 'label_sig', pch.col = 'grey20', order = 'AOE')
#corrplot(cor_ciber, method = "number", type = "lower",col = addcol(100), 
#tl.col = "n", tl.cex = 0.8, tl.pos = "n",order = 'AOE',
#add = T)
dev.off()
write.csv(cor_ciber, "DEcells-DEcells_correlation.csv", quote = F, row.names = T)
testrespvalue <- testRes$p
write.csv(testrespvalue, "DEcells-DEcells_pvalue_testRes.csv", quote = F, row.names = T)





# env.p <- cor_pmat(DEcells_spearman,method = "spearman")
# 
# decol <- colorRampPalette(c("#1874CD", "white", "#FED82F"))(100)
# 
# pdf("heatmap_Spearman_ssGSEA.pdf",height = 12,width = 18)
# corrplot(corr =DEcells_cor,type="upper",tl.pos="tp",tl.col="black",p.mat = env.p,tl.cex = 1.6 ,
#          insig = "label_sig", sig.level = c(.01, .05),pch.cex=2,pch.col = "black",order = "AOE")
# corrplot(corr = DEcells_cor,type="lower",add=TRUE,method="circle",tl.pos="n",tl.col="black",
#          diag=FALSE, cl.pos="n",order = "AOE")
# dev.off()
#关键基因和差异免疫细胞的相关性####
# which(colnames(ssgsea_result1) == 'Type.2.T.helper.cell')
# DEcells <- data.frame(T.cells.follicular.helper = cibersort_results1[,9])
# rownames(DEcells) <- rownames(cibersort_results1)
keygenes <- c(
  'WWTR1',
  'SLC7A11',
  'G6PD',
  'ZEB1',
  'EGR1',
  'FLNA',
  'FADS2',
  'SLC39A14'
  
)
exprsTCGA<-read.table('/data/nas1/chenpeiru/44_GYZK-30212-7/04_lasso/risk.train.group.txt',
                      row.names = 1,sep = "\t",header = T)%>%lc.tableToNum()
head(exprsTCGA)
keyGs_exp <- exprsTCGA[,c(3:10)]
keyGs_exp <- t(keyGs_exp)
head(keyGs_exp)
#keyGs_exp <- keyGs_exp[,rownames(DEcells_spearman)]
genelist <- rownames(keyGs_exp)
keyGs_exp <- apply(keyGs_exp,2,function(x){as.numeric(x)})
rownames(keyGs_exp) <- genelist

#批量计算相关性
gene <- genelist
immuscore <- function(gene){
  y <- as.numeric(keyGs_exp[gene,])
  colnames <- colnames(DEcells_spearman)
  do.call(rbind,lapply(colnames, function(x){
    dd  <- corr.test(as.numeric(DEcells_spearman[,x]), y , method="spearman",adjust = "fdr")
    data.frame(gene=gene,immune_cells=x,cor=dd$r,p.value=dd$p )
  }))
}

#批量计算关键跟差异免疫浸润细胞相关性的结果
KGs_DEcells_spearman_data <- do.call(rbind,lapply(genelist,immuscore))
head(KGs_DEcells_spearman_data)
write.csv(KGs_DEcells_spearman_data, "KGs_DEcells_correlation.csv", quote = F, row.names = F)


#热图
KGs_DEcells_data <- KGs_DEcells_spearman_data
KGs_DEcells_data$sig <- ifelse(KGs_DEcells_data$p.value < 0.05,
                               ifelse(KGs_DEcells_data$p.value < 0.01,
                                      ifelse(KGs_DEcells_data$p.value < 0.001,'***','**'),'*'),' ')
pdf('KGs_DEcells_correlation_ssGSEA.pdf',width = 12,height = 10)
ggplot(KGs_DEcells_data,aes(x = immune_cells,y = gene, fill=cor))+ 
  geom_tile()+  
  scale_fill_gradient2(low = '#1874CD',mid = 'white',high ='#FED82F',
                       limits=c(-1,1),breaks=c(-1,-0.5,0,0.5,1))+
  labs(x=NULL,y=NULL)+theme_bw(base_size = 15)+
  geom_text(aes(label = sig), color = 'black', size = 8)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"))+
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","*** p < 0.001","\n\n","Correlation"))
dev.off()
#ESTIMATE算法计算高低风险组的基质评分、免疫评分、ESTIMATE评分####
# library(estimate)
# load("/data/nas1/chenpeiru/41_YQNC-10615-6/00_raw_data/01_TCGA/dat_fpkm425.rda")
# dat <- dat_fpkm
# #group$sample<-gsub('.','-',group$sample,fixed = T)
# 
# group <- read.delim2('/data/nas1/chenpeiru/41_YQNC-10615-6/01_TCGA_DEG/time30/group.xls')
# tumor.sample <- group[which(group$group=='Tumor'),]
# 
# head(group)
# dat <- dat[,tumor.sample$sample]#399
# 
# #dat <- data.frame(rownames(dat),dat)
# #colnames(dat)[1] <- "Sample"
# 
# group10588 <- read.table("/data/nas1/chenpeiru/41_YQNC-10615-6/06_km_roc/risk.train.group.txt",
#                          sep = "\t",header = T,row.names = 1) %>% lc.tableToNum()
# head(group10588)
# risk.train.group <- group10588[,9]%>%as.data.frame()
# rownames(risk.train.group) <- rownames(group10588)
# colnames(risk.train.group)[1] <- "group"
# #colnames(risk.train.group) <- c("Sample","group")
# dat_t <- t(dat)
# TCGA_fpkm_843 <- merge(risk.train.group,dat_t,by="row.names")
# TCGA_fpkm_843t <- t(TCGA_fpkm_843)
# head(TCGA_fpkm_843t)
# TCGA_fpkm_843t <- TCGA_fpkm_843t[-2,]
# colnames(TCGA_fpkm_843t) <- TCGA_fpkm_843t[1,]
# TCGA_fpkm_843t <- TCGA_fpkm_843t[-1,]
# head(TCGA_fpkm_843t)
# 
# write.table(TCGA_fpkm_843t,"TCGA_fpkm_estimate.txt",quote = F,row.names = T,sep = "\t")
# 
# TCGA_fpkm_843t_p <- read.table("TCGA_fpkm_estimate.txt",header = T,row.names = 1,sep = "\t")
# 
# 
# filterCommonGenes(input.f = 'TCGA_fpkm_estimate.txt',
#                   output.f = 'riskscore_immune.gct',
#                   id = "GeneSymbol")
# estimateScore(input.ds = 'riskscore_immune.gct',
#               output.ds= 'estimate_score.gct', 
#               platform="affymetrix")#还可以选择其他测序平台
# 
# estimate_score <- read.table("estimate_score.gct", skip = 2, header = TRUE,check.names = F)
# ##写出csv
# write.csv(estimate_score,"estimate_score.csv",row.names = FALSE,quote = F)
# # estimate_score <- read.csv("../data/10_immune/estimate_score.csv",header = T,check.names = F,row.names = 1)
# estimate_score <- estimate_score[,-1]
# estimate_score <- as.data.frame(t(estimate_score))
# colnames(estimate_score) <- estimate_score[1,]
# estimate_score <-estimate_score[-1,]
# head(estimate_score)
# head(risk.train.group)
# rownames(estimate_score)<-gsub('.','-',rownames(estimate_score),fixed = T)
# estimate_score <- data.frame(rownames(estimate_score),estimate_score)
# head(estimate_score)
# 
# colnames(estimate_score)[1] <- "Sample"
# head(risk.train.group)
# estimate_score <- estimate_score[-1,]
# estimate_score1 <- merge(risk.train.group,estimate_score,by="row.names")
# head(estimate_score1)
# rownames(estimate_score1) <- estimate_score1$Row.names
# estimate_score1 <- estimate_score1[,-1]
# head(estimate_score1)
# #estimate_score$group <- risk.train.group$group
# 
# #小提琴图  ####
# head(estimate_score1)
# estimate_score1 <- estimate_score1[order(estimate_score1$group),]
# compaired <- list(c( "low","high"))
# e_im <- ggplot(estimate_score1,aes(x=group,y=as.numeric(ImmuneScore)),palette = "jco", add = "jitter")+
#   geom_violin(aes(fill=group))+
#   geom_boxplot(width=0.1,cex=1.2,color = 'grey')+
#   labs(x = "Group", y = "ImmuneScore")+
#   scale_fill_manual(values = c("#104E8B","#56B4E9")) + theme_zg2()+
#   geom_signif(comparisons = compaired,step_increase = 0.3,
#               map_signif_level = F,test = wilcox.test,textsize = 6)
# 
# e_st <- ggplot(estimate_score1,aes(x=group,y=as.numeric(StromalScore)),palette = "jco", add = "jitter")+
#   geom_violin(aes(fill=group))+
#   geom_boxplot(width=0.1,cex=1.2,color = 'grey')+
#   labs(x = "Group", y = "StromalScore")+
#   scale_fill_manual(values = c("#104E8B","#56B4E9")) + theme_zg2()+
#   geom_signif(comparisons = compaired,step_increase = 0.3,
#               map_signif_level = F,test = wilcox.test,textsize = 6)
# 
# e_es <- ggplot(estimate_score1,aes(x=group,y=as.numeric(ESTIMATEScore)),palette = "jco", add = "jitter")+
#   geom_violin(aes(fill=group))+
#   geom_boxplot(width=0.1,cex=1.2,color = 'grey')+
#   labs(x = "Group", y = "ESTIMATEScore")+
#   scale_fill_manual(values = c("#104E8B","#56B4E9")) + theme_zg2()+
#   geom_signif(comparisons = compaired,step_increase = 0.3,map_signif_level = F,
#               test = wilcox.test,textsize = 6)
# 
# pdf('ESTIMATE_violin.pdf',width = 30,height = 8)
# e_im + e_st + e_es
# dev.off()
#tide score----

setwd("/data/nas1/chenpeiru/44_GYZK-30212-7/11_checkpoints")
rm(list = ls())
load("/data/nas1/chenpeiru/44_GYZK-30212-7/00_raw_data/01_TCGA/dat_fpkm425.rda")
dat <- dat_fpkm

group <- read.delim2('/data/nas1/chenpeiru/44_GYZK-30212-7/01_TCGA_DEG/group.xls')
tumor.sample <- group[which(group$group=='Tumor'),]

head(group)
#group$sample<-gsub('.','-',group$sample,fixed = T)
dat <- dat[,tumor.sample$sample]#399

rt <- dat
risk <- read.table('/data/nas1/chenpeiru/44_GYZK-30212-7/04_lasso/risk.train.group.txt',
                   sep = '\t',header = T,check.names = F)
head(risk)
check<-read_xlsx("checkpoints.xlsx",sheet=1) %>% as.data.frame()

int <- intersect(check$checkpoints,rownames(rt))
#write.table(int,'E:\\project\\24_GY0316\\GY0316\\GY0316_Result\\13.checkpoint\\intersect47.txt',sep = '\t',quote = F,row.names = T)



rt1<-rt[rownames(rt)%in%check$checkpoints,]
head(rt1)
write.table(rt1,'checkpoint exp_new.txt',sep = '\t',quote = F,row.names = T)

data<-read.table('checkpoint exp_new.txt',sep = '\t',header = T,check.names = F,row.names = 1)
head(data)

library(ggpubr)
library(ggplot2)
library(reshape2)
data<-data.frame(t(data))
data$id<-rownames(data)
head(data)
#colnames(data)<-gsub('[.]', '-', colnames(data))
#head(data)
head(risk)
cc<-merge(data,risk[,c(1,13)],by='id')
head(cc)
rownames(cc)<-cc$id
cc<-cc[,-1]
colnames(cc)
colnames(cc)[40] <- "risk"
box2 <- melt(cc, id.vars = "risk", variable.name = "Checkpoint",
             value.name = "Expression")
head(box2)
box2[,3] <- as.numeric(box2$Expression)
head(box2)
#box2[,3] <- log2(box2[,3] +1)

box <- ggplot(box2, aes(x = Checkpoint, y = Expression))+ 
  geom_boxplot(aes(fill = risk),
               position = position_dodge(0.5),
               width = 0.4)+
  scale_fill_manual(values = c("#FED82F","#1874CD")) +
  stat_compare_means(aes(group = risk),label = "p.signif")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 10),axis.text.y = element_text(size = 12))
box
ggsave('CheckpointBoxExp-new.pdf', box, width = 12,height = 7)

#高低风险组HLA表达
setwd("/data/nas1/chenpeiru/44_GYZK-30212-7/11_checkpoints")
rm(list = ls())
load("/data/nas1/chenpeiru/44_GYZK-30212-7/00_raw_data/01_TCGA/dat_fpkm425.rda")
dat <- dat_fpkm

group <- read.delim2('/data/nas1/chenpeiru/44_GYZK-30212-7/01_TCGA_DEG/group.xls')
tumor.sample <- group[which(group$group=='Tumor'),]

head(group)
#group$sample<-gsub('.','-',group$sample,fixed = T)
dat <- dat[,tumor.sample$sample]#399

rt <- dat
risk <- read.table('/data/nas1/chenpeiru/44_GYZK-30212-7/04_lasso/risk.train.group.txt',
                   sep = '\t',header = T,check.names = F)
head(risk)
check<-read_xlsx("checkpoints.xlsx",sheet=2) %>% as.data.frame()

int <- intersect(check$checkpoints,rownames(rt))
#write.table(int,'E:\\project\\24_GY0316\\GY0316\\GY0316_Result\\13.checkpoint\\intersect47.txt',sep = '\t',quote = F,row.names = T)



rt1<-rt[rownames(rt)%in%check$checkpoints,]
head(rt1)
write.table(rt1,'02_checkpoint exp_new_HLA.txt',sep = '\t',quote = F,row.names = T)

data<-read.table('02_checkpoint exp_new_HLA.txt',sep = '\t',header = T,check.names = F,row.names = 1)
head(data)

library(ggpubr)
library(ggplot2)
library(reshape2)
data<-data.frame(t(data))
data$id<-rownames(data)
head(data)
#colnames(data)<-gsub('[.]', '-', colnames(data))
#head(data)
head(risk)
cc<-merge(data,risk[,c(1,13)],by='id')
head(cc)
rownames(cc)<-cc$id
cc<-cc[,-1]
colnames(cc)
colnames(cc)[20] <- "risk"
box2 <- melt(cc, id.vars = "risk", variable.name = "Checkpoint",
             value.name = "Expression")
head(box2)
box2[,3] <- as.numeric(box2$Expression)
head(box2)
#box2[,3] <- log2(box2[,3] +1)

box <- ggplot(box2, aes(x = Checkpoint, y = Expression))+ 
  geom_boxplot(aes(fill = risk),
               position = position_dodge(0.5),
               width = 0.4)+
  scale_fill_manual(values = c("#FED82F","#1874CD")) +
  stat_compare_means(aes(group = risk),label = "p.signif")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 10),axis.text.y = element_text(size = 12))
box
ggsave('02_CheckpointBoxExp-new_HLA.pdf', box, width = 12,height = 7)
# cc$risk=ifelse(cc$risk=="high", "High-risk", "Low-risk")
# group=levels(factor(cc$risk))
# cc$risk=factor(cc$risk, levels=c("Low-risk", "High-risk"))
# group=levels(factor(cc$risk))
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# head(cc)
# gg1=ggviolin(cc, x="risk", y="PD-L1", fill = "risk",
#              xlab="",
#              ylab="PD-L1",
#              legend.title="",
#              palette = c("#00AFBB","#FC4E07"),
#              add = "boxplot", add.params = list(fill = "white"))+
#   stat_compare_means(comparisons = my_comparisons)
# pdf(file="PD-L1.pdf", width=6, height=5)
# print(gg1)
# dev.off()
#免疫浸润分析 ssGSEA####
# rm(list = ls())
# setwd("/data/nas1/chenpeiru/41_YQNC-10615-6/13_checkpoints")
# #source('/data/nas1/chenpeiru/39_YQWLMQ-10122-3/09_CIBERSORT/Cibersort_source.R')
# library(RColorBrewer)
# library(tidyr)
# #library(radiant.data)
# library(aplot)
# library(dplyr)
# library(e1071)
# library(tidyverse)
# library(reshape2)
# library(psych)
# library(corrplot)
# library(cowplot)
# library(ggcorrplot)
# library(msigdbr)
# library(GSEABase)
# library(GSVA)
# #setwd("/data/nas1/chenpeiru/39_YQWLMQ-10122-3/09_CIBERSORT")\
# load("/data/nas1/chenpeiru/41_YQNC-10615-6/00_raw_data/01_TCGA/dat_fpkm425.rda")
# dat <- dat_fpkm
# 
# group <- read.delim2('/data/nas1/chenpeiru/41_YQNC-10615-6/01_TCGA_DEG/time30/group.xls')
# tumor.sample <- group[which(group$group=='Tumor'),]
# 
# head(group)
# #group$sample<-gsub('.','-',group$sample,fixed = T)
# dat <- dat[,tumor.sample$sample]#399
# 
# #dat <- mRNA_FPKM_log
# 
# ssGSEA_data <- dat
# gmtFile="immune_pathway.gmt"
# #分析免疫细胞含量
# #载入背景基因
# geneSet=getGmt(gmtFile,
#                geneIdType=SymbolIdentifier())
# #开始ssGSEA分析,数据一定要是矩阵，不是数据框
# ssgseaScore=gsva(as.matrix(ssGSEA_data), geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
# 
# range(ssgseaScore)
# ssgseaOut <- as.data.frame(t(ssgseaScore))
# 
# #整理数据格式
# ssgsea_result1 = ssgseaOut
# 
# write.csv(ssgsea_result1,"ssGSEA_result.csv",row.names = T,quote = F)
# 
# library(lance)
# ssgsea_result1 <- read.csv("ssGSEA_result.csv",header = T,row.names =1)
# group10588 <-  read.table("/data/nas1/chenpeiru/41_YQNC-10615-6/06_km_roc/risk.train.group.txt",
#                           sep = "\t",header = T) %>% lc.tableToNum()
# head(group10588)
# group10588 <- group10588[,c(1,10)]
# colnames(group10588) <- c("Sample","group")
# head(group10588)
# head(ssgsea_result1)
# ssgsea_result1 <- data.frame(rownames(ssgsea_result1),ssgsea_result1)
# head(ssgsea_result1)
# colnames(ssgsea_result1)[1] <- "Sample"
# ssgsea_result1 <- merge(group10588,ssgsea_result1,by="Sample")
# head(ssgsea_result1)
# rownames(ssgsea_result1) <- ssgsea_result1[,1]
# head(ssgsea_result1)
# 
# ssgsea_result1 <- ssgsea_result1[,-1]
# head(ssgsea_result1)
# ssgsea_result1 <- ssgsea_result1[order(ssgsea_result1$group, decreasing = TRUE), ]
# re1 <- ssgsea_result1
# 
# head(re1)
# #去除分组信息
# colnames(re1)[1] <- "Type"
# re2 = re1[,-1]       #去除第一列分组信息
# mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
# #提取数据，多行变成多列，要多学习‘tidyr’里面的三个函数
# dat_cell <- re2 %>% as.data.frame() %>%rownames_to_column("Sample") %>%gather(key = Cell_type,value = Proportion,-Sample)
# #提取数据
# dat_group = gather(re1,Cell_type,Proportion,-Type )
# #合并分组
# dat = cbind(dat_cell,dat_group$Type)
# head(dat)
# colnames(dat)[4] <- "Type"
# ##2.3柱状图############
# 
# #热图####
# annotation_col <-  data.frame(sample = c(rep("low", 247), rep("high", 152)))
# rownames(annotation_col) <- rownames(re2)
# pdf('heatmap_ImmuneGeme.pdf',height = 9,width = 12)
# pheatmap(t(re2),
#          cluster_rows = T,
#          cluster_cols = F,
#          annotation_col = annotation_col,
#          show_colnames = F,
#          name = " ")
# dev.off()
# 
# #筛选差异免疫浸润细胞 ####
# box=dat
# theme_zg2 <- function(..., bg='white'){
#   require(grid)
#   theme_classic(...) +
#     theme(rect=element_rect(fill=bg),
#           plot.margin=unit(rep(0.5,4), 'lines'),
#           panel.background=element_rect(fill='transparent',color='black'),
#           panel.border=element_rect(fill='transparent', color='transparent'),
#           panel.grid=element_blank(),#去网格线
#           axis.line = element_line(colour = "black"),
#           #axis.title.x = element_blank(),#去x轴标签
#           axis.title.y=element_text(face = "bold",size = 14),#y轴标签加粗及字体大小
#           axis.title.x=element_text(face = "bold",size = 14),#X轴标签加粗及字体大小
#           axis.text.y = element_text(face = "bold",size = 12),#y坐标轴刻度标签加粗
#           axis.text.x = element_text(face = "bold",size = 10, vjust = 1, hjust = 1, angle = 45),#x坐标轴刻度标签加粗
#           axis.ticks = element_line(color='black'),
#           # axis.ticks.margin = unit(0.8,"lines"),
#           legend.title=element_blank(),
#           legend.position=c(0.9, 0.85),#图例在绘图区域的位置
#           # legend.position="top",
#           legend.direction = "horizontal",
#           legend.text = element_text(face = "bold",size = 12),
#           # legend.background = element_rect( linetype="solid",colour ="black")
#     )
# }
# e1 <- ggplot(box,aes(x=Cell_type,y=Proportion),palette = "jco", add = "jitter")+
#   geom_boxplot(aes(fill=Type),width=0.6)+
#   labs(x = "Cell Type", y = "Estimated Proportion")+
#   scale_fill_manual(values = c("#1874CD","#FED82F")) +
#   theme_zg2() + stat_compare_means(aes(group = Type),label = "p.signif",method = 'wilcox.test')
# e1
# ggsave('box_ImmuneGeme_ssGSEA.pdf', plot = e1,width=15,height = 6,family="serif")
# 
# #免疫细胞相关性分析 ####
# head(ssgsea_result1)
# DEcells_spearman <- ssgsea_result1[,-1]
# colnames(DEcells_spearman)
# #删除没有丰度的细胞#1,2,3,5,10,11,12,13,15,16,18,19,20,23,24,25,28
# DEcells_spearman <- DEcells_spearman[,c(4,5,6,7,11,16)]
# colnames(DEcells_spearman)
# #DEcells_cor <- cor(DEcells_spearman,method = "spearman") 
# #class(DEcells_spearman$Activated.CD8.T.cell)
# #class(DEcells_spearman$Activated.CD4.T.cell)
# 
# ciber.res<-DEcells_spearman#[,colSums(DEcells_spearman)>0]
# #ciber.res<-t(ciber.res)
# cor_ciber<-cor(ciber.res,method = 'spearman')
# cor_ciber<-round(cor(ciber.res),2)
# # pdf("cor.plot_DEcells_spearman.pdf",width=13,height=9) #保存图为pdf
# # col1=colorRampPalette(colors =c("blue","white","red"),space="Lab")
# # testRes = cor.mtest(cor_ciber, method="spearman",conf.level = 0.95)
# # corrplot(cor_ciber, 
# #          p.mat = testRes$p, diag = T, type = 'upper',col = col1(10),
# #          sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.2,
# #          insig = 'label_sig', pch.col = 'grey20', order = 'AOE')
# # #corrplot(cor_ciber, method = "number", type = "lower",col = addcol(100), 
# # #tl.col = "n", tl.cex = 0.8, tl.pos = "n",order = 'AOE',
# # #add = T)
# # dev.off()
# # write.csv(cor_ciber, "DEcells-DEcells_correlation.csv", quote = F, row.names = T)
# # testrespvalue <- testRes$p
# # write.csv(testrespvalue, "DEcells-DEcells_pvalue_testRes.csv", quote = F, row.names = T)
# 
# 
# 
# 
# 
# # env.p <- cor_pmat(DEcells_spearman,method = "spearman")
# # 
# # decol <- colorRampPalette(c("#1874CD", "white", "#FED82F"))(100)
# # 
# # pdf("heatmap_Spearman_ssGSEA.pdf",height = 12,width = 18)
# # corrplot(corr =DEcells_cor,type="upper",tl.pos="tp",tl.col="black",p.mat = env.p,tl.cex = 1.6 ,
# #          insig = "label_sig", sig.level = c(.01, .05),pch.cex=2,pch.col = "black",order = "AOE")
# # corrplot(corr = DEcells_cor,type="lower",add=TRUE,method="circle",tl.pos="n",tl.col="black",
# #          diag=FALSE, cl.pos="n",order = "AOE")
# # dev.off()
# #关键基因和差异免疫细胞的相关性####
# rm(list = ls())
# load("/data/nas1/chenpeiru/41_YQNC-10615-6/00_raw_data/01_TCGA/dat_fpkm425.rda")
# dat <- dat_fpkm
# 
# group <- read.delim2('/data/nas1/chenpeiru/41_YQNC-10615-6/01_TCGA_DEG/time30/group.xls')
# tumor.sample <- group[which(group$group=='Tumor'),]
# 
# head(group)
# #group$sample<-gsub('.','-',group$sample,fixed = T)
# dat <- dat[,tumor.sample$sample]#399
# 
# keygenes <- c(
#   'RAD9A',
#   'PTGIS',
#   'P4HB',
#   'MYC',
#   'G6PD'
# )
# exprsTCGA <- dat
# head(exprsTCGA)
# #keyGs_exp <- exprsTCGA[,c(3:7)]
# #keyGs_exp <- t(keyGs_exp)
# #head(keyGs_exp)
# #exprsTCGA <- t(exprsTCGA)
# keyGs_exp <- exprsTCGA[keygenes,]
# keyGs_exp <- keyGs_exp[,rownames(DEcells_spearman)]
# genelist <- rownames(keyGs_exp)
# keyGs_exp <- apply(keyGs_exp,2,function(x){as.numeric(x)})
# rownames(keyGs_exp) <- genelist
# 
# #批量计算相关性
# gene <- genelist
# immuscore <- function(gene){
#   y <- as.numeric(keyGs_exp[gene,])
#   colnames <- colnames(DEcells_spearman)
#   do.call(rbind,lapply(colnames, function(x){
#     dd  <- corr.test(as.numeric(DEcells_spearman[,x]), y , method="spearman",adjust = "fdr")
#     data.frame(gene=gene,immune_cells=x,cor=dd$r,p.value=dd$p )
#   }))
# }
# 
# #批量计算关键跟差异免疫浸润细胞相关性的结果
# KGs_DEcells_spearman_data <- do.call(rbind,lapply(genelist,immuscore))
# head(KGs_DEcells_spearman_data)
# #KGs_DEcells_spearman_data <- subset(KGs_DEcells_spearman_data,immune_cells==)
# write.csv(KGs_DEcells_spearman_data, "KGs_DEcells_correlation.csv", quote = F, row.names = F)
# 
# 
# #热图
# KGs_DEcells_data <- KGs_DEcells_spearman_data
# KGs_DEcells_data$sig <- ifelse(KGs_DEcells_data$p.value < 0.05,
#                                ifelse(KGs_DEcells_data$p.value < 0.01,
#                                       ifelse(KGs_DEcells_data$p.value < 0.001,'***','**'),'*'),' ')
# pdf('KGs_DEcells_correlation_ssGSEA.pdf',width = 6,height = 10)
# ggplot(KGs_DEcells_data,aes(x = immune_cells,y = gene, fill=cor))+ 
#   geom_tile()+  
#   scale_fill_gradient2(low = '#1874CD',mid = 'white',high ='#FED82F',
#                        limits=c(-1,1),breaks=c(-1,-0.5,0,0.5,1))+
#   labs(x=NULL,y=NULL)+theme_bw(base_size = 15)+
#   geom_text(aes(label = sig), color = 'black', size = 8)+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),
#         axis.text.y = element_text(size = 14, face = "bold"))+
#   labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","*** p < 0.001","\n\n","Correlation"))
# dev.off()
# 
# #ips------
# rm(list = ls())
# setwd("/data/nas1/chenpeiru/41_YQNC-10615-6/13_checkpoints")
# LL<-read.csv('IPS.csv',
#              sep = ',',header = T,check.names = F)
# 
# risk<-read.table('/data/nas1/chenpeiru/41_YQNC-10615-6/06_km_roc/risk.train.group.txt',
#                  sep = '\t',header = T,check.names = F)
# head(risk)
# risk$id<-substr(risk$id,1,12)
# head(risk)
# 
# LL1<-merge(LL,risk[,c(1,10)],by.y='id')
# head(LL1)
# LL1=LL1[!duplicated(LL1),]
# 
# rownames(LL1)<-LL1$id
# # write.csv(LL1,"LL1.csv")
# LL1<-LL1[,-1]
# data<-LL1
# head(data)
# colnames(data)[5] <- "risk"
# 
# data$risk=ifelse(data$risk=="high", "High-risk", "Low-risk")
# group=levels(factor(data$risk))
# data$risk=factor(data$risk, levels=c("Low-risk", "High-risk"))
# group=levels(factor(data$risk))
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# gg1=ggviolin(data, x="risk", y="ips_ctla4_neg_pd1_neg", fill = "risk",
#              xlab="",
#              ylab="ips_ctla4_neg_pd1_neg",
#              legend.title="",
#              palette = c("#1874CD","#FED82F"),
#              add = c("jitter", "mean_sd")#,
#              #add.params = list(fill = "white")
# )+
#   stat_compare_means(comparisons = my_comparisons)
# pdf(file="ips_ctla4_neg_pd1_neg.pdf", width=6, height=5)
# print(gg1)
# dev.off()
# 
# gg2=ggviolin(data, x="risk", y="ips_ctla4_neg_pd1_pos", fill = "risk",
#              xlab="",
#              ylab="ips_ctla4_neg_pd1_pos",
#              legend.title="",
#              palette = c("#1874CD","#FED82F"),
#              add = c("jitter", "mean_sd"))+
#   stat_compare_means(comparisons = my_comparisons)
# pdf(file="02-ips_ctla4_neg_pd1_pos.pdf", width=6, height=5)
# print(gg2)
# dev.off()
# 
# gg3=ggviolin(data, x="risk", y="ips_ctla4_pos_pd1_neg", fill = "risk",
#              xlab="",
#              ylab="ips_ctla4_pos_pd1_neg",
#              legend.title="",
#              palette = c("#1874CD","#FED82F"),
#              add = c("jitter", "mean_sd"))+
#   stat_compare_means(comparisons = my_comparisons)
# pdf(file="03-ips_ctla4_pos_pd1_neg.pdf", width=6, height=5)
# print(gg3)
# dev.off()
# 
# gg4=ggviolin(data, x="risk", y="ips_ctla4_pos_pd1_pos", fill = "risk",
#              xlab="",
#              ylab="ips_ctla4_pos_pd1_pos",
#              legend.title="",
#              palette = c("#1874CD","#FED82F"),
#              add = c("jitter", "mean_sd"))+
#   stat_compare_means(comparisons = my_comparisons)
# pdf(file="04-ips_ctla4_pos_pd1_pos.pdf", width=6, height=5)
# print(gg4)
# dev.off()
# 
# library("IMvigor210CoreBiologies")





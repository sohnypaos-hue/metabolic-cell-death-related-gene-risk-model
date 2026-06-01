# #train boxplot----
# rm(list = ls())
# setwd("/data/nas1/chenpeiru/44_GYZK-30212-7/00_raw_data/01_TCGA")
# keygenes <- c('WWTR1',
#               'SLC7A11',
#               'G6PD',
#               'ZEB1',
#               'EGR1',
#               'FLNA',
#               'FADS2',
#               'SLC39A14'
#               
#               
# )
# load("/data/nas1/chenpeiru/44_GYZK-30212-7/00_raw_data/01_TCGA/dat_fpkm425.rda")
# exprs55457 <- dat_fpkm
# group <- read.delim2("/data/nas1/chenpeiru/44_GYZK-30212-7/01_TCGA_DEG/group.xls")
# range(exprs55457)
# head(group)
# gse55457_rocgenes_exp <- exprs55457[keygenes,group$sample]
# gse55457_df <- as.data.frame(t(gse55457_rocgenes_exp)) 
# 
# # group55457 <-  read_xlsx("/data/nas1/chenpeiru/43_JNZK-20704-5/00_raw_data/01_train/group.xlsx",
# #                          sheet = 1) %>% as.data.frame()
# group55457 <- group
# colnames(group55457) <- c("Sample","group")
# 
# gse55457_df$group <- group55457$group
# head(gse55457_df)
# gse55457_df_melt <- reshape2::melt(gse55457_df , value.name = "value")
# colnames(gse55457_df_melt) <- c("group","Gene","value")
# #gse55457_df_melt$value <- log2(gse55457_df_melt$value+1)
# 
# #zd_genes_sig <- gse55457_topGenes[keygenes,]
# #zd_genes_sig$label <- ifelse(zd_genes_sig$P.Value < 0.05,'*',
# #                             ifelse(zd_genes_sig$P.Value < 0.01,'**',
# #                                   ifelse(zd_genes_sig$P.Value < 0.001,'***','ns')))
# #zd_genes_sig$gene <- rownames(zd_genes_sig)
# #pdf('E:/project/11_LZZK-10411-1/05_machine_learning/GSE55457_BOXPLOT.pdf',width = 12,height = 8)
# theme_zg <- function(..., bg='white'){
#   require(grid)
#   theme_classic(...) +
#     theme(rect=element_rect(fill=bg),
#           plot.margin=unit(rep(0.5,4), 'lines'),
#           panel.background=element_rect(fill='transparent',color='black'),
#           panel.border=element_rect(fill='transparent', color='transparent'),
#           panel.grid=element_blank(),#去网格线
#           axis.line = element_line(colour = "black"),
#           #axis.title.x = element_blank(),#去x轴标签
#           axis.title.y=element_text(face = "bold",size = 18),#y轴标签加粗及字体大小
#           axis.title.x=element_text(face = "bold",size = 18),#X轴标签加粗及字体大小
#           axis.text.y = element_text(face = "bold",size = 14),#y坐标轴刻度标签加粗
#           axis.text.x = element_text(face = "bold",size = 12, vjust = 1, hjust = 1, angle = 45),#x坐标轴刻度标签加粗
#           axis.ticks = element_line(color='black'),
#           # axis.ticks.margin = unit(0.8,"lines"),
#           legend.title=element_blank(),
#           # legend.position=c(0.9, 0.8),#图例在绘图区域的位置
#           legend.position="right",
#           legend.box = "vertical",
#           legend.direction = "vertical",
#           legend.text = element_text(face = "bold",size = 12,margin = ggplot2::margin(r=8))
#           # legend.background = element_rect( linetype="solid",colour ="black")
#     )
# }
# #setwd("/data/nas1/chenpeiru/41_YQNC-10615-6/16_GENE_EXPRESSION")
# pdf("/data/nas1/chenpeiru/44_GYZK-30212-7/06_independent_prognosis/01_train_boxplot.pdf",height =  12,width = 15)
# head(gse55457_df_melt)
# p25 <- ggplot(gse55457_df_melt,aes(x=factor(Gene),y=value),palette = "jco", add = "jitter")+
#   geom_boxplot(aes(fill=group),position=position_dodge(0.3),width=0.3)+
#   labs(x = "Gene", y = "Expression (log2)")+
#   scale_fill_manual(values = c("blue","#DC0000FF")) + theme_zg()+
#   stat_compare_means(label = "p.signif", aes(group = group),size = 6#,
#                      #method.args = list(alternative = "less")
#   )
# # geom_text(data = zd_genes_sig, aes(x = gene, y = pvalue,label = label),nudge_y = 4)
# p25
# dev.off()

#new
library(limma)
library(pheatmap)
library(ggplot2)
library(data.table)
library(dplyr)

#####训练集GSE94019训练集####
rm(list=ls())
keygenes <- c('WWTR1',
              'SLC7A11',
              'G6PD',
              'ZEB1',
              'EGR1',
              'FLNA',
              'FADS2',
              'SLC39A14'
              
              
)
load("/data/nas1/chenpeiru/44_GYZK-30212-7/00_raw_data/01_TCGA/dat_fpkm425.rda")
exprs55457 <- dat_fpkm

group <- read.delim2("/data/nas1/chenpeiru/44_GYZK-30212-7/01_TCGA_DEG/group.xls")
range(exprs55457)
head(group)
data_candi <- exprs55457[keygenes,group$sample]
group_N <- group$sample[group$group == "Normal"]
##候选基因
#gene <- read.csv("../03_Lasso/lasso_genes.csv")
#data_candi <- data[rownames(data)%in%gene$x,]
####把长列变短gather()
data_candi$gene <- rownames(data_candi)
violin_dat <- gather(data_candi, sample, value='expr', -c("gene"))
violin_dat$group <- ifelse(violin_dat$sample %in% group_N, "Control","Tumor") 

library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)
colnames(violin_dat)

library(rstatix)

stat.test<-violin_dat%>%
  group_by(gene)%>%
  wilcox_test(expr ~ group)%>%
  adjust_pvalue(method = 'BH')
#write.csv(stat.test,file = 'train_wilcox_result.csv')
# stat.test$p.signif <- ifelse(stat.test$p<0.0001,"****",ifelse(stat.test$p<0.001,"***",ifelse(stat.test$p<0.01,"**",ifelse(stat.test$p<0.05,"*","ns"))))
# stat.test <- stat.test %>% add_xy_position(x="gene")

#差异结果---
pval <- read.csv("/data/nas1/chenpeiru/44_GYZK-30212-7/01_TCGA_DEG/DEG_logfc1_padj0.05.csv",
                 header = T,row.names = 1)
pval <- pval[stat.test$gene,]
colnames(pval)
stat.test$padj <-pval$padj
stat.test$p.signif <- ifelse(stat.test$padj<0.0001,"****",ifelse(stat.test$padj<0.001,"***",ifelse(stat.test$padj<0.01,"**",ifelse(stat.test$padj<0.05,"*","ns"))))
#write.csv(stat.test,file = 'train_DESeq2_result.csv')

stat.test <- stat.test %>% add_xy_position(x="gene")
head(violin_dat)
violin_plot <- ggplot(violin_dat, aes(x=gene,  y=expr))+ 
  #  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               aes(fill=group),
               width=0.1,
               position = position_dodge(0.9)) +
  geom_violin(aes(fill=group))+
  #geom_boxplot(width=0.7,
  #             aes(fill=group),
  #             position=position_dodge(0.9),
   #            outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#0074b3","#982b2b"), name = "Group")+#"#0074b3","#d8d8d8","#982b2b"
  #labs(title="Immune Checkpoint", x="", y = "log2(expr+1)",size=20) 
  labs(x="", y = "log2(expr+1)",size=20)+ 
  stat_pvalue_manual(stat.test,
                     # group = gene,
                     #y.position=14.3,
                     label = "p.signif",
                     tip.length = 0.01,
                     color='black',
                     family = "Times",
                     #parse = T,
                     face = "bold")+
  theme_bw()+
theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
      axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=12), 
      axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
      axis.title.x=element_text(size=16,face="bold"),
      axis.title.y=element_text(size=16,face="bold"),
      legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
      legend.title = element_text(face = "bold", size = 14),
      legend.position = "top",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())

# facet_wrap(~gene,scales = "free",nrow = 1) 
violin_plot
ggsave('/data/nas1/chenpeiru/44_GYZK-30212-7/06_independent_prognosis/01_train_boxplot_new.pdf',
       violin_plot,w=12,h=5)
ggsave('/data/nas1/chenpeiru/44_GYZK-30212-7/06_independent_prognosis/01_train_boxplot_new.png',
       violin_plot,w=12,h=5)
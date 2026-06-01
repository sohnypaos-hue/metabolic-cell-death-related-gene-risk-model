
rm(list=ls())
setwd("/data/nas1/chenpeiru/44_GYZK-30212-7/01_TCGA_DEG")
library(tidyverse)
library(magrittr)
library(stringr)
library(limma)
library(readxl)
library(DESeq2)
library(lance)
##DESEQ2--------
load("/data/nas1/chenpeiru/44_GYZK-30212-7/00_raw_data/01_TCGA/dat_fpkm_tumor406_survival.rda")
colnames(fpkm_tumor_p_t_survival)
colnames(fpkm_tumor_p_t_survival)[1] <- "fustat"
colnames(fpkm_tumor_p_t_survival)[2] <- "futime"
fpkm_tumor_p_t_survival$futime <- as.numeric(fpkm_tumor_p_t_survival$futime)
fpkm_tumor_p_t_survival <- subset(fpkm_tumor_p_t_survival,futime>30)#399

setwd("/data/nas1/chenpeiru/44_GYZK-30212-7/00_raw_data/01_TCGA")
dat <- read.delim2('dat.tcga_count.xls',row.names = 1)%>%lc.tableToNum()
head(dat)

colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)
dat_tumor <- dat[,-c(1:19)]
dat_tumor_30 <- dat_tumor[,colnames(dat_tumor)%in%rownames(fpkm_tumor_p_t_survival)]
dat_30 <- cbind(dat[,1:19],dat_tumor_30)
colnames(dat_30)
dat_30<-round(dat_30,digits = 0)

# mete=data.frame(colnames(dat))  # 取第一行样本id
# for (i in 1:length(mete[,1])) {
#   num=as.numeric(as.character(substring(mete[i,1],14,15)))
#   if(num %in% seq(1,9)){mete[i,2]="T"}
#   if(num %in% seq(10,29)){mete[i,2]="N"}
# }
# names(mete)=c("id","group")
# table(mete$group)#  N 19 T 406
# mete$group=as.factor(mete$group)
# mete=subset(mete,mete$group=="T")
# exp_tumor<-dat[,which(colnames(dat)%in%mete$id)]
# exp_tumor<-as.data.frame(exp_tumor)
# survival<-read.delim2('TCGA-BLCA.survival.tsv')#448,4
# head(survival)
#survival <- subset(survival,survival$OS.time >90)
# exp_tumor<-exp_tumor[,colnames(exp_tumor)%in%survival$sample]#375
# 
# exp_control<-dat[,which(!colnames(dat)%in%mete$id)]#19
# exp_control<-as.data.frame(exp_control)
# # 19
# dat.final<-cbind(exp_control,exp_tumor)
# head(dat.final)#
# dat <- dat.final
# save(dat,file = "data_survival_3month.rda")
library(DESeq2)

colData <- data.frame(sample=colnames(dat_30),
                      group=c(rep('Normal',19),rep('Tumor',399)))
colData$group<-factor(colData$group,levels = c('Normal','Tumor'))
setwd("/data/nas1/chenpeiru/44_GYZK-30212-7/01_TCGA_DEG")
write.table(colData,file = 'group.xls',sep = '\t',quote = F,row.names = F)
row_sums <- rowSums(dat_30 == 0)
summary(row_sums)
# 创建一个新的数据框，只包括那些不超过3个0的行
#14647-13903=744
# filtered_data <- dat[row_sums <= 100, ]#只包括那些不超过100个0的行



dds <- DESeqDataSetFromMatrix(countData =dat_30,colData=colData,design = ~group)#16013
#keep <- rowSums(counts(dds) >= 10) >= 3   
#dds <- dds[keep, ] 
dds <- dds[rownames(counts(dds)) > 1,]
dds <- estimateSizeFactors(dds)
##提取标准化后的数据
#normalized_counts <- counts(dds,normalized=T)
#write.table(normalized_counts,file = 'normalized.counts.xls',sep = '\t',row.names = T,quote = F)
dds<-DESeq(dds)
## 提取差异结果
res =results(dds, contrast = c("group","Tumor","Normal"))
res =res[order(res$padj),]
head(res)
write.csv(res, "TCGA-blca_diffall.csv", quote = F, row.names = T)

summary(res)
table(res$padj<0.05)
#logfc2
rm(list=ls())
res <- read.csv("/data/nas1/chenpeiru/44_GYZK-30212-7/01_TCGA_DEG/TCGA-blca_diffall.csv",
                header = T,sep = ",",row.names = 1)
colnames(res)
#logfc2
DEG <- subset(res, padj < 0.05 & abs(log2FoldChange) >2 )#1935
DEG <- as.data.frame(DEG)
setwd("/data/nas1/chenpeiru/44_GYZK-30212-7/01_TCGA_DEG")
save(DEG,file = "DEG_logfc2_padj0.05.rda")
write.csv(DEG,"DEG_logfc2_padj0.05.csv",row.names = T,quote = F)
#logfc1
DEG1 <- subset(res, padj < 0.05 & abs(log2FoldChange) >1 )#4661
DEG1 <- as.data.frame(DEG1)
save(DEG1,file = "DEG_logfc1_padj0.05.rda")
write.csv(DEG1,"DEG_logfc1_padj0.05.csv",row.names = T,quote = F)
#logfc0.5

DEG2 <- subset(res, padj < 0.05 & abs(log2FoldChange) > 0.5 )#7742
DEG2 <- as.data.frame(DEG2)#7742
save(DEG2,file = "DEG_logfc0.5_padj0.05.rda")
write.csv(DEG2,"DEG_logfc0.5_padj0.05.csv",row.names = T,quote = F)

#logfc0.5 pvalue 0.05
# rm(list=ls())
# res <- read.csv("/data/nas1/chenpeiru/41_YQNC-10615-6/01_TCGA_DEG/TCGA-blca_diffall.csv",
#                 header = T,sep = ",",row.names = 1)
# colnames(res)
# DEG2 <- subset(res, pvalue < 0.05 & abs(log2FoldChange) > 0.5 )#8229
# DEG2 <- as.data.frame(DEG2)#8229
# setwd("/data/nas1/chenpeiru/41_YQNC-10615-6/01_TCGA_DEG")
# save(DEG2,file = "DEG_logfc0.5_pvalue0.05.rda")
# write.csv(DEG2,"DEG_logfc0.5_pvalue0.05.csv",row.names = T,quote = F)
# # DEG<-na.omit(DEG)
# dim(DEG)
# head(DEG)
# volcano plot------
rm(list = ls())

load("/data/nas1/chenpeiru/44_GYZK-30212-7/01_TCGA_DEG/DEG_logfc1_padj0.05.rda")#4661

GSE_diff <- DEG1
GSE_diff$symbol<-rownames(GSE_diff)
DEG <- read.csv('/data/nas1/chenpeiru/44_GYZK-30212-7/01_TCGA_DEG/TCGA-blca_diffall.csv',
                header = T,check.names = F,row.names = 1)

head(DEG)
DEG$change <- as.factor(ifelse(DEG$padj<0.05&abs(DEG$log2FoldChange)>1,
                               ifelse(DEG$log2FoldChange>1,"up","down"),
                               "no_diff"))
table(DEG$change)
# 对差异表达基因的FDR值进行从小到大排序
DEG <- DEG[order(DEG$padj), ]
head(DEG)
# 高表达的基因中，选择adj.P.Val最小的10个
up.genes <- head(rownames(DEG)[which(DEG$change == "up")], 10)

# 低表达的adj.P.Val最小的10个
down.genes <- head(rownames(DEG)[which(DEG$change == "down")], 10)


# 将up.genes和down.genes合并，并加入到Label中
deg.top10.genes <- c(as.character(up.genes), as.character(down.genes))
#DEGs_file$Label[match(deg.top10.genes, rownames(DEGs_file))] <- deg.top10.genes
library(dplyr)
best_in_class <- DEG[rownames(DEG)%in%deg.top10.genes,]



df <- data.frame(DEG$log2FoldChange, DEG$padj, DEG$change)
head(df)
colnames(df) <- c("logFC", "adj.P.Val", "group")
head(best_in_class)
colnames(best_in_class)[2] <- "logFC"
colnames(best_in_class)[6] <- "adj.P.Val"
colnames(best_in_class)[7] <- "group"

head(df)
head(DEG)
head(best_in_class)
#-log10(3.683861e-179)
options(ggrepel.max.overlaps = Inf)
p <- ggplot(df, aes(x = logFC, y = -log10(adj.P.Val), color = group)) +
  geom_point(size = 1.5, alpha = 0.7) +
  labs(y = "-log10(adj.P.Val)", x = "logFC") +
  scale_color_manual(values=c("red", "blue","black"), limits = c("up", "down", "no_diff")) +
  geom_vline(xintercept = c(-1, 1), lty = 4, col = "black", lwd = 0.5) +
  geom_hline(yintercept = -log10(0.05), lty = 2, col = "black", lwd = 0.5) +
  scale_x_continuous(breaks = c(-20,-10,0,10,20),limits = c(-25,25)) +
  scale_y_continuous(limits = c(0,60))+
  theme(legend.title = element_blank(), 
        panel.background = element_rect(color = "black", fill = "transparent")) +
  ggrepel::geom_label_repel(
    aes(label = rownames(best_in_class)),
    data = best_in_class
  )
p
setwd("/data/nas1/chenpeiru/44_GYZK-30212-7/01_TCGA_DEG")
ggsave(p,filename = paste0("TCGA DEG volcano.pdf"),height = 7)
ggsave(p,filename = paste0("TCGA DEG volcano.png"),height = 7)

#pheatmap
#热图
##全部差异基因表达矩阵提取 
#rm(list = ls())
setwd("/data/nas1/chenpeiru/44_GYZK-30212-7/01_TCGA_DEG")
load("/data/nas1/chenpeiru/44_GYZK-30212-7/00_raw_data/01_TCGA/dat_fpkm425.rda")
data <- dat_fpkm %>% as.data.frame()
group <- read.delim2("group.xls")
head(group)
red_de_expr_418 <- data[,colnames(data)%in%group$sample] #注意这里是data_GEO,不是data
write.csv(red_de_expr_418,file="/data/nas1/chenpeiru/44_GYZK-30212-7/01_TCGA_DEG/dat_fpkm418.txt",quote=F)


red_de_expr <- data[rownames(data) %in% rownames(GSE_diff),colnames(data)%in%group$sample] #注意这里是data_GEO,不是data
range(red_de_expr)
#red_de_expr <- log2(red_de_expr+1)
write.table(red_de_expr,file="diffSigExp.txt",sep="\t",quote=F)

#heatmap ####
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
## DEGs expression
#top_up <- DE_reslut %>% filter(change == "Up") 
#top_down <- DE_reslut %>% filter(change == "Down") 
#top <- c(rownames(top_up),rownames(top_down))

de_expr <- red_de_expr
range(de_expr)
head(de_expr)
de_expr <- de_expr[rownames(de_expr)%in%deg.top10.genes,]
madt<-as.matrix(de_expr)
madt2<-t(scale(t(madt)))
range(madt2)
Group <- read.delim2("/data/nas1/chenpeiru/44_GYZK-30212-7/01_TCGA_DEG/group.xls") %>% as.data.frame()
head(Group)
colnames(Group) <- c("Sample","Type")
head(Group)

col_fun <- colorRamp2(c(min(madt2), 0, max(madt2)),c("#4F94CD", "white", "#8B1A1A"))
ha1 = HeatmapAnnotation(group = as.factor(c(rep("Normal", table(Group$Type)[1][[1]]),rep("Tumor", table(Group$Type)[2][[1]]))),
                        col =list(group = c("Tumor" = "#E9967A", "Normal" = "#4682B4")))
pdf("heatmap_top10.pdf", width=10, height=8)
densityHeatmap(madt2,col=colorRampPalette(c("#4F94CD", "white", "#8B1A1A"))(50),quantile_gp = gpar(fontsize = 9),title = " ",
               ylab = "Expression",top_annotation = ha1) %v% 
  Heatmap(madt2, name = "expression",col =col_fun,show_row_names = FALSE,show_column_names = FALSE, cluster_rows = TRUE,height = unit(11, "cm"))
dev.off()


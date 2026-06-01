rm(list = ls())
library(lance)
library(tidyverse)
load("/data/nas1/chenpeiru/44_GYZK-30212-7/00_raw_data/01_TCGA/dat_fpkm425.rda")
dat <- dat_fpkm

group <- read.delim2('/data/nas1/chenpeiru/44_GYZK-30212-7/01_TCGA_DEG/group.xls')
tumor.sample <- group[which(group$group=='Tumor'),]

head(group)
#group$sample<-gsub('.','-',group$sample,fixed = T)
dat <- dat[,tumor.sample$sample]#399

# colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)
group <- read.table("/data/nas1/chenpeiru/44_GYZK-30212-7/04_lasso/risk.train.group.txt",
                    sep = "\t",header = T) %>% lc.tableToNum()

table(group$group)
high.sample <- group$id[which(group$group=='high')]
low.sample <- group$id[which(group$group=='low')]
head(group)
dat<-dat[,group$id]
#BiocManager::install('pRRophetic',type = "source")
#install.packages("pRRophetic")
library(pRRophetic)
library(ggplot2)
set.seed(12345)
colnames(group)

hubgene <- colnames(group)[4:11]
model_expr<-dat[hubgene,]
setwd("/data/nas1/chenpeiru/44_GYZK-30212-7/12_IC50")
drug<-read.table(file = 'drugs.txt',sep='\t',header=F)
ic50<-data.frame(group$id)
a<-data.frame(row.names=group$id,group=group$group)
colnames(a)<-'group'

cnt<-1

while (cnt < 139) {
  
  predictedPtype <- pRRopheticPredict(as.matrix(model_expr), drug[cnt,],selection=1)
  
  Tipifarnib<-data.frame(predictedPtype)
  
  colnames(Tipifarnib)<-drug[cnt,]
  
  a<-cbind(a,Tipifarnib)
  
  cnt = cnt + 1
}

write.table(a,'IC50.xls',sep='\t',quote=F)
write.table(a,'IC50.txt',sep='\t',quote=F)

b<-a
b[b<0]<-NA
# 先写成函数的形式，方便调用
removeRowsAllNa  <- function(x){x[apply(x, 1, function(y) any(!is.na(y))),]}
removeColsAllNa  <- function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}
c<-removeColsAllNa(b)
na_flag <- apply(is.na(c), 2, sum)
x <- c[, which(na_flag == 0)]
View(x)
dim(x)
# [1] 513,108
medicinal_result <- t(subset(x, select = -group)) 
high_group <- high.sample
low_group <- low.sample
pvalue = padj = log2FoldChange <- matrix(0, nrow(medicinal_result), 1)
for (i in 1:nrow(medicinal_result)){
  pvalue[i, 1] = p.value = wilcox.test(medicinal_result[i, high_group],
                                       medicinal_result[i, low_group])$p.value
  log2FoldChange[i, 1] = mean(medicinal_result[i, high_group]) - 
    mean(medicinal_result[i, low_group])
}
padj <- p.adjust(as.vector(pvalue), "fdr", n = length(pvalue))
rTable <- data.frame(log2FoldChange, 
                     pvalue, 
                     padj,
                     row.names = rownames(medicinal_result))
high_group_res <- signif(apply(medicinal_result[rownames(rTable), high_group], 
                               1,
                               median), 4)
low_group_res <- signif(apply(medicinal_result[rownames(rTable), low_group], 
                              1, 
                              median), 4)
rTable <- data.frame(high_group_res, 
                     low_group_res,
                     rTable[, c("padj", "pvalue", "log2FoldChange")])
rTable$drugs <- rownames(rTable)
rTable$sig <- ifelse(rTable$padj < 0.05,
                     ifelse(rTable$padj < 0.01, 
                            ifelse(rTable$padj < 0.001,
                                   ifelse(rTable$padj < 0.0001,
                                          paste(rTable$drugs, "****",  sep = ""),
                                          paste(rTable$drugs, "***", sep = "")),
                                   paste(rTable$drugs, "**", sep = "")),
                            paste(rTable$drugs, "*",  sep = "")), 
                     rTable$drugs)

write.table(rTable,
            file = "drugs_wilcox_test.xls",
            quote = F,
            row.names = F)
write.table(rTable,
            file = "drugs_wilcox_test.txt",
            quote = F,
            row.names = F,sep = "\t")
### 发散条形图绘制

#install.packages('ggprism')
library(ggprism)
## 横坐标药物，纵坐标：IC50(H)/IC50(L)-1
dat_plot<-data.frame(drug=rownames(rTable),
                     'IC50(H)/IC50(L)-1'=(rTable$high_group_res/rTable$low_group_res-1),
                     pvalue=rTable$pvalue,
                     padj=rTable$padj)
dat_plot$threshold=factor(ifelse(dat_plot$padj<0.05&dat_plot$pvalue<0.05,'P<0.05 & FDR<0.05',ifelse(dat_plot$pvalue<0.05&dat_plot$padj>0.05,'P<0.05 & FDR>0.05','P>0.05 & FDR>0.05')))
dat_plot<-dat_plot%>%arrange(desc(dat_plot$IC50.H..IC50.L..1))

dat_plot$drug<-factor(dat_plot$drug,levels=dat_plot$drug)

p <- ggplot(data = dat_plot,aes(x = drug,y = IC50.H..IC50.L..1,fill = threshold)) +
  geom_col()+
  scale_fill_manual(values = c('P<0.05 & FDR<0.05'= '#5F9EA0','P>0.05 & FDR>0.05'='#cccccc','P<0.05 & FDR>0.05'='#FFD700')) +
  xlab('') +
  ylab('Exp(Median IC50(H))/Exp(Median IC50(L))-1') +
  theme_prism(border = T) +
  theme(
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 65,size = 7,
                               hjust = 1,vjust = 1),
    axis.text.y = element_text(size = 13),
    legend.position = c(0.85,0.85),
    legend.text = element_text(size = 8,face = 'bold')
  )
p
ggsave(p,filename = '01.allplot.pdf',family="serif",w=10,h=6)
ggsave(p,filename = '01.allplot.png',w=10,h=6)


### 差异最显著的6个-----
head(dat_plot)
dat_plot <- dat_plot[order(dat_plot$padj),]

top5 <- c('Vorinostat','ABT.263','Gefitinib','BIBW2992','X681640')#药物最敏感
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)
#head(drugs_res)
all(rownames(rTable) == rownames(medicinal_result))
drugs_res <- data.frame(drugs=rownames(medicinal_result), medicinal_result, pvalue=rTable$pvalue)
drugs_res <- drugs_res[top5,]
violin_dat <- gather(drugs_res, key=indivs, value=score, -c("drugs","pvalue"))
table(group$group)
violin_dat$indivs <- ifelse(gsub("\\.","-",violin_dat$indivs) %in% high_group,
                            "High", "Low")
violin_dat$indivs <- factor(violin_dat$indivs, levels = c("High", "Low"))
head(violin_dat)
drugs_hub_boxplot1 <- ggboxplot(violin_dat, x = "indivs", y = "score",
                                color = "indivs", palette = c("#A73030FF", "#0073C2FF"),
                                add = "jitter",
                                short.panel.labs = T,
                                ggtheme = theme_bw()) +
  stat_compare_means(label = "p.signif", label.x = 1.4, vjust = 0.5)
drugs_hub_boxplot <- facet(drugs_hub_boxplot1,
                           facet.by = "drugs",
                           short.panel.labs = T,
                           panel.labs.background = list(fill = "white"),
                           ncol = 3,
                           scales = "free_y") + xlab("") + ylab("IC(50)") +
  # geom_text(data=data_text,
  #           mapping=aes(x=x,y=y,label=label),nudge_x=0.1,nudge_y=0.1)+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        # strip.background = element_blank(),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold",angle = 45,hjust = 1),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        text = element_text(size = 13, face = "bold"))
drugs_hub_boxplot
ggsave(filename = "IC50_drugs.plot-top5.pdf", family="serif",height = 6, width = 7)
ggsave(filename = "IC50_drugs.plot-top5.png", height = 6, width = 7)
#correlation-riskscore-IC50------
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
head(keyGs_exp)
keyGs_exp <- t(keyGs_exp)
head(keyGs_exp)
#keyGs_exp <- keyGs_exp[,rownames(DEcells_spearman)]
genelist <- rownames(keyGs_exp)
keyGs_exp <- apply(keyGs_exp,2,function(x){as.numeric(x)})
rownames(keyGs_exp) <- genelist

#批量计算相关性
gene <- genelist
head(medicinal_result)
medicinal_result <- medicinal_result[top5,]
DEcells_spearman <- t(medicinal_result)
head(DEcells_spearman)

immuscore <- function(gene){
  y <- as.numeric(keyGs_exp[gene,])
  colnames <- colnames(DEcells_spearman)
  do.call(rbind,lapply(colnames, function(x){
    dd  <- corr.test(as.numeric(DEcells_spearman[,x]), y , method="spearman",adjust = "fdr")
    data.frame(gene=gene,immune_cells=x,cor=dd$r,p.value=dd$p )
  }))
}

#批量计算
KGs_DEcells_spearman_data <- do.call(rbind,lapply(genelist,immuscore))
head(KGs_DEcells_spearman_data)
riskscore_ic50 <- KGs_DEcells_spearman_data#[c(),]#riskscore
write.csv(riskscore_ic50, "riskscore_ic50.csv", quote = F, row.names = F)


#热图
#setwd("/data/nas1/chenpeiru/28_YQ874-4/12_IC50")
riskscore_ic50 <- read.csv("riskscore_ic50.csv",header=T,sep=",")
head(riskscore_ic50)#excel-----
# riskscore_ic50 <- subset(riskscore_ic50,p.value<0.05)
# riskscore_ic50 <- riskscore_ic50[order(riskscore_ic50$cor,decreasing = T),]
# head(riskscore_ic50)
# riskscore_ic50_high10 <- riskscore_ic50[c(1:5),]
# head(riskscore_ic50_high10)
# # write_csv(riskscore_ic50_high10,"riskscore_ic50_corelation_high10.csv",quote = F)
# riskscore_ic50 <- riskscore_ic50[order(riskscore_ic50$cor,decreasing = F),]
# head(riskscore_ic50)
# riskscore_ic50_low10 <- riskscore_ic50[c(1:5),]
riskscore_ic50_10 <- riskscore_ic50
#riskscore_ic50_10 <- rbind(riskscore_ic50_high10,riskscore_ic50_low10)
head(riskscore_ic50_10)
KGs_DEcells_data <- riskscore_ic50_10
KGs_DEcells_data$sig <- ifelse(KGs_DEcells_data$p.value < 0.05,
                               ifelse(KGs_DEcells_data$p.value < 0.01,
                                      ifelse(KGs_DEcells_data$p.value < 0.001,'***','**'),'*'),' ')
pdf('riskscore_ic50_top10.pdf',width = 8,height = 10)
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
head(KGs_DEcells_data)
table(KGs_DEcells_data$immune_cells)
# ### 差异最显著的5个-----
# top5 <- riskscore_ic50_high10$immune_cells#药物最敏感
# library(tidyr)
# library(ggplot2)
# library(ggpubr)
# library(Ipaper)
# head(drugs_res)
# all(rownames(rTable) == rownames(medicinal_result))
# drugs_res <- data.frame(drugs=rownames(medicinal_result), medicinal_result, pvalue=rTable$pvalue)
# drugs_res <- drugs_res[top5,]
# violin_dat <- gather(drugs_res, key=indivs, value=score, -c("drugs","pvalue"))
# table(group$group)
# violin_dat$indivs <- ifelse(gsub("\\.","-",violin_dat$indivs) %in% high_group,
#                             "High", "Low") 
# violin_dat$indivs <- factor(violin_dat$indivs, levels = c("High", "Low"))
# head(violin_dat)
# drugs_hub_boxplot1 <- ggboxplot(violin_dat, x = "indivs", y = "score",
#                                 color = "indivs", palette = c("#A73030FF", "#0073C2FF"),
#                                 add = "jitter",
#                                 short.panel.labs = T,
#                                 ggtheme = theme_bw()) +
#   stat_compare_means(label = "p.signif", label.x = 1.4, vjust = 0.5)
# drugs_hub_boxplot <- facet(drugs_hub_boxplot1,
#                            facet.by = "drugs",
#                            short.panel.labs = T,
#                            panel.labs.background = list(fill = "white"),
#                            ncol = 3,
#                            scales = "free_y") + xlab("") + ylab("IC(50)") +
#   # geom_text(data=data_text,
#   #           mapping=aes(x=x,y=y,label=label),nudge_x=0.1,nudge_y=0.1)+
#   theme(panel.grid = element_blank(),
#         legend.position = "none",
#         # strip.background = element_blank(),
#         axis.title.y = element_text(size = 18, face = "bold"),
#         axis.title.x = element_blank(),
#         axis.text.y = element_text(size = 10, face = "bold"),
#         axis.text.x = element_text(size = 10, face = "bold",angle = 45,hjust = 1),
#         plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
#         text = element_text(size = 13, face = "bold"))
# drugs_hub_boxplot
# ggsave(filename = "IC50_drugs.plot-high-top5.pdf", height = 6, width = 7)
# ggsave(filename = "IC50_drugs.plot-high-top5.png", height = 6, width = 7)
# 
# 
# top_low <-riskscore_ic50_low10$immune_cells #药物最敏感
# library(tidyr)
# library(ggplot2)
# library(ggpubr)
# library(Ipaper)
# head(drugs_res)
# all(rownames(rTable) == rownames(medicinal_result))
# drugs_res <- data.frame(drugs=rownames(medicinal_result), medicinal_result, pvalue=rTable$pvalue)
# drugs_res <- drugs_res[top_low,]
# violin_dat <- gather(drugs_res, key=indivs, value=score, -c("drugs","pvalue"))
# table(group$group)
# violin_dat$indivs <- ifelse(gsub("\\.","-",violin_dat$indivs) %in% high_group,
#                             "High", "Low") 
# violin_dat$indivs <- factor(violin_dat$indivs, levels = c("High", "Low"))
# head(violin_dat)
# drugs_hub_boxplot1 <- ggboxplot(violin_dat, x = "indivs", y = "score",
#                                 color = "indivs", palette = c("#A73030FF", "#0073C2FF"),
#                                 add = "jitter",
#                                 short.panel.labs = T,
#                                 ggtheme = theme_bw()) +
#   stat_compare_means(label = "p.signif", label.x = 1.4, vjust = 0.5)
# drugs_hub_boxplot <- facet(drugs_hub_boxplot1,
#                            facet.by = "drugs",
#                            short.panel.labs = T,
#                            panel.labs.background = list(fill = "white"),
#                            ncol = 3,
#                            scales = "free_y") + xlab("") + ylab("IC(50)") +
#   # geom_text(data=data_text,
#   #           mapping=aes(x=x,y=y,label=label),nudge_x=0.1,nudge_y=0.1)+
#   theme(panel.grid = element_blank(),
#         legend.position = "none",
#         # strip.background = element_blank(),
#         axis.title.y = element_text(size = 18, face = "bold"),
#         axis.title.x = element_blank(),
#         axis.text.y = element_text(size = 10, face = "bold"),
#         axis.text.x = element_text(size = 10, face = "bold",angle = 45,hjust = 1),
#         plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
#         text = element_text(size = 13, face = "bold"))
# drugs_hub_boxplot
# ggsave(filename = "IC50_drugs.plot-low-top5.pdf", height = 6, width = 7)
# ggsave(filename = "IC50_drugs.plot-low-top5.png", height = 6, width = 7)

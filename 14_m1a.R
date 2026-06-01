rm(list = ls())
setwd("/data/nas1/chenpeiru/44_GYZK-30212-7/14_m1a")
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
#m1a-----
m1a <- read_xlsx("/data/nas1/chenpeiru/44_GYZK-30212-7/14_m1a/m1a.xlsx",sheet=1)%>% as.data.frame()
load("/data/nas1/chenpeiru/44_GYZK-30212-7/00_raw_data/01_TCGA/dat_fpkm425.rda")
dat <- dat_fpkm

group <- read.delim2('/data/nas1/chenpeiru/44_GYZK-30212-7/01_TCGA_DEG/group.xls')
tumor.sample <- group[which(group$group=='Tumor'),]

head(group)
#group$sample<-gsub('.','-',group$sample,fixed = T)
dat <- dat[,tumor.sample$sample]#399
dat <- t(dat)
head(dat)
m1a_exp <- dat[,m1a$m1a]
DEcells_spearman <- m1a_exp

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
write.csv(KGs_DEcells_spearman_data, "01_m1a_KGs_correlation.csv", quote = F, row.names = F)


#热图
KGs_DEcells_data <- KGs_DEcells_spearman_data
KGs_DEcells_data$sig <- ifelse(KGs_DEcells_data$p.value < 0.05,
                               ifelse(KGs_DEcells_data$p.value < 0.01,
                                      ifelse(KGs_DEcells_data$p.value < 0.001,'***','**'),'*'),' ')
pdf('01_m1a_KGs_correlation_ssGSEA.pdf',width = 12,height = 10)
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
#m5c---
m5c <- read_xlsx("/data/nas1/chenpeiru/44_GYZK-30212-7/14_m1a/m1a.xlsx",sheet=2)%>% as.data.frame()

m5c_exp <- dat[,m5c$m5c]
DEcells_spearman <- m5c_exp

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
write.csv(KGs_DEcells_spearman_data, "01_m5c_KGs_correlation.csv", quote = F, row.names = F)


#热图
KGs_DEcells_data <- KGs_DEcells_spearman_data
KGs_DEcells_data$sig <- ifelse(KGs_DEcells_data$p.value < 0.05,
                               ifelse(KGs_DEcells_data$p.value < 0.01,
                                      ifelse(KGs_DEcells_data$p.value < 0.001,'***','**'),'*'),' ')
pdf('02_m5c_KGs_correlation_ssGSEA.pdf',width = 12,height = 10)
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
#m6a----
m6a <- read_xlsx("/data/nas1/chenpeiru/44_GYZK-30212-7/14_m1a/m1a.xlsx",sheet=3)%>% as.data.frame()

m6a_exp <- dat[,m6a$m6a]
DEcells_spearman <- m6a_exp

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
write.csv(KGs_DEcells_spearman_data, "03_m6a_KGs_correlation.csv", quote = F, row.names = F)


#热图
KGs_DEcells_data <- KGs_DEcells_spearman_data
KGs_DEcells_data$sig <- ifelse(KGs_DEcells_data$p.value < 0.05,
                               ifelse(KGs_DEcells_data$p.value < 0.01,
                                      ifelse(KGs_DEcells_data$p.value < 0.001,'***','**'),'*'),' ')
pdf('03_m6a_KGs_correlation_ssGSEA.pdf',width = 12,height = 10)
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
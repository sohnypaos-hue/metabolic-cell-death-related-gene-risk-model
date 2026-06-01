# 01 获取数据集----------
rm(list=ls())
library(TCGAbiolinks)
library(readr)
library(readxl)
library(tidyverse)
library(lance)
setwd("/data/nas1/chenpeiru/44_GYZK-30212-7/00_raw_data/01_TCGA")
#library(lance)
## 读取从xena下载的数据
tcga.expr1 <- read_tsv(file = 'TCGA-BLCA.htseq_counts.tsv')%>%lc.tableToNum()#Rows: 60488 Columns: 431 
head(tcga.expr1)
# tcga.expr <- tcga.expr1 %>%lc.tableToNum()
# class(tcga.expr$`GTEX-S4Q7-0003-SM-3NM8M`)
rownames(tcga.expr1) <- tcga.expr1$Ensembl_ID
tcga.expr <- tcga.expr1[,-1]
rownames(tcga.expr) <- rownames(tcga.expr1)
head(tcga.expr)

## xena下载的数据经过了log2+1转化，需要将其还原
tcga.expr<-2^tcga.expr-1
## 对数据进行id转化
genecode<-read.table(file = 'gencode.v22.annotation.gene.probeMap')
probe2symbol<-genecode[,(1:2)]
colnames(probe2symbol)<-c('ID','symbol')
probe2symbol<-probe2symbol[-1,]
dat.tcga<-tcga.expr
dat.tcga$ID <- rownames(dat.tcga)
dat.tcga$ID<-as.character(dat.tcga$ID)
probe2symbol$ID<-as.character(probe2symbol$ID)
dat.tcga<-dat.tcga %>%
  inner_join(probe2symbol,by='ID')%>% 
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('TCGA',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
dim(dat.tcga)#[1] 58387   430

## 筛选癌症组织，去掉癌旁组织。01-09为肿瘤，10-19为正常对照
mete=data.frame(colnames(dat.tcga))  # 取第一行样本id
for (i in 1:length(mete[,1])) {
  num=as.numeric(as.character(substring(mete[i,1],14,15)))
  if(num %in% seq(1,9)){mete[i,2]="T"}
  if(num %in% seq(10,29)){mete[i,2]="N"}
}
names(mete)=c("id","group")
table(mete$group)#  N 19 T 411 
mete$group=as.factor(mete$group)
mete=subset(mete,mete$group=="T")
exp_tumor<-dat.tcga[,which(colnames(dat.tcga)%in%mete$id)]
exp_tumor<-as.data.frame(exp_tumor)
# 411
## 保留有生存数据的
survival<-read.delim2('TCGA-BLCA.survival.tsv')#448,4

write.table(survival,file = 'survival.xls',sep = '\t',row.names = F,quote = F)
head(survival)
exp_tumor<-exp_tumor[,colnames(exp_tumor)%in%survival$sample]
#406
head(exp_tumor)
pcg <- read.delim2('/data/nas1/chenpeiru/28_YQ874-4/00_raw_data/PCG.xls(v22)')
exp_tumor_p <- exp_tumor[pcg$gene_name,]
exp_tumor_p <- na.omit(exp_tumor_p)#19814
head(exp_tumor_p)
exp_tumor_p <- na.omit(exp_tumor_p)
save(exp_tumor_p,file = "TCGA_Tumor406_count.rda")#row is gene, col is sample

exp_tumor_p_t <- t(exp_tumor_p)
# head(exp_tumor_p_t)
library(tibble)

exp_tumor_p_t <- data.frame(names = row.names(exp_tumor_p_t), exp_tumor_p_t)
colnames(exp_tumor_p_t)[1] <- "sample"

exp_tumor_p_t_survival <- merge(survival,exp_tumor_p_t,by.x="sample",all=FALSE)#406
rownames(exp_tumor_p_t_survival) <- exp_tumor_p_t_survival$sample
exp_tumor_p_t_survival <- exp_tumor_p_t_survival[,-c(1,3)]
save(exp_tumor_p_t_survival,file = "TCGA_Tumor406_count_survival.rda")


exp_control<-dat.tcga[,which(!colnames(dat.tcga)%in%mete$id)]#19
exp_control<-as.data.frame(exp_control)
# 19
dat.final<-cbind(exp_control,exp_tumor)
head(dat.final)#425
##425
pcg <- read.delim2('/data/nas1/chenpeiru/28_YQ874-4/00_raw_data/PCG.xls(v22)')
dat.final <- dat.final[pcg$gene_name,]
dat.final <- na.omit(dat.final)#19814
dat.final <- na.omit(dat.final)
write.table(dat.final,file = 'dat.tcga_count.xls',sep = '\t',quote = F,row.names = T)
save(dat.final,file = "dat.tcga_count425.rda")#row is gene

colnames(dat.final)

mete=data.frame(colnames(dat.final))  # 取第一行样本id
for (i in 1:length(mete[,1])) {
  num=as.numeric(as.character(substring(mete[i,1],14,15)))
  if(num %in% seq(1,9)){mete[i,2]="T"}
  if(num %in% seq(10,29)){mete[i,2]="N"}
}
names(mete)=c("id","group")
table(mete$group)#N 19, T 406
##fpkm
setwd("/data/nas1/chenpeiru/44_GYZK-30212-7/00_raw_data/01_TCGA")
expr_fpkm1 <- read_tsv(file = 'TCGA-BLCA.htseq_fpkm.tsv') %>%lc.tableToNum()#Rows: 60483 Columns: 431   
rownames(expr_fpkm1) <- expr_fpkm1$Ensembl_ID
expr_fpkm <- expr_fpkm1[,-1]
rownames(expr_fpkm) <- rownames(expr_fpkm1)
head(expr_fpkm)

## 
## 对数据进行id转化
dat_fpkm<-expr_fpkm
dat_fpkm$ID <- rownames(dat_fpkm)
dat_fpkm$ID<-as.character(dat_fpkm$ID)
probe2symbol$ID<-as.character(probe2symbol$ID)
dat_fpkm<-dat_fpkm %>%
  inner_join(probe2symbol,by='ID')%>% 
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('TCGA',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
dim(dat_fpkm)# 58387   430
# dat_fpkm<-dat_fpkm[mRNA$gene_name,]
dat_fpkm<-dat_fpkm[,colnames(dat.final)]#418 sample
dat_fpkm<-dat_fpkm[pcg$gene_name,]
dat_fpkm<-na.omit(dat_fpkm)#19814,425
save(dat_fpkm,file = "dat_fpkm425.rda")
write.table(dat_fpkm,file = 'dat.fpkm.xls',sep = '\t',row.names = T,quote = F)
#select tumor and os time-----
dat_fpkm_tumor <-dat_fpkm[,colnames(exp_tumor)]#406 sample
dat_fpkm_tumor<-dat_fpkm_tumor[pcg$gene_name,]
dat_fpkm_tumor<-na.omit(dat_fpkm_tumor)#19814,406
dat_fpkm_tumor_t <- t(dat_fpkm_tumor)
range(dat_fpkm_tumor_t)
dat_fpkm_tumor_t <- data.frame(names = row.names(dat_fpkm_tumor_t), dat_fpkm_tumor_t)
colnames(dat_fpkm_tumor_t)[1] <- "sample"

fpkm_tumor_p_t_survival <- merge(survival,dat_fpkm_tumor_t,by.x="sample",all=FALSE)#406
rownames(fpkm_tumor_p_t_survival) <- fpkm_tumor_p_t_survival$sample
fpkm_tumor_p_t_survival <- fpkm_tumor_p_t_survival[,-c(1,3)]
save(fpkm_tumor_p_t_survival,file = "dat_fpkm_tumor406_survival.rda")

# write.table(dat_fpkm,file = 'dat.fpkm.xls',sep = '\t',row.names = T,quote = F)
# fpkm转TPM
# FPKM2TPM <- function(fpkm){
#   exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
# }
# 
# dat_tpm <- apply(dat_fpkm,2,FPKM2TPM)
# setwd("/data/nas1/chenpeiru/28_YQ874-4/00_raw_data")
# 
# write.table(dat_tpm,file = 'dat.tpm.xls',sep = '\t',row.names = T,quote = F)


## 整理表型数据的
rm(list=ls())
setwd("/data/nas1/chenpeiru/44_GYZK-30212-7/00_raw_data/01_TCGA")

phenotype1 <- read.delim2('TCGA-BLCA.GDC_phenotype.tsv')#
colnames(phenotype1)
select_list <- c("submitter_id.samples","age_at_initial_pathologic_diagnosis",
                 "gender.demographic",
                 'pathologic_M','pathologic_N','pathologic_T',
                 'tumor_stage.diagnoses',
                 'number_of_lymphnodes_positive_by_he')
phenotype <- phenotype1[,select_list]
head(phenotype)
colnames(phenotype) <- c("id.samples","Age","Gender","M_stage",
                         "N_stage","T_stage","Stage",'Lymphnodes_positive')
head(phenotype)
phenotype <- na.omit(phenotype)#454
head(phenotype)
write.csv(phenotype,file = 'phenotype.csv',row.names = F)

library(readxl)
phenotype_label <- read_xlsx("phenotype_label.xlsx.xlsx",
                             sheet=1) %>% as.data.frame()# 这个就是在上面csv文件基础上去除了NA，excel里面操作的
head(phenotype_label)
phenotype_label$Gender <- factor(phenotype_label$Gender,
                                 levels = c("male","female"),labels = c(1,2))
head(phenotype_label)
#M
M=phenotype_label$M_stage

M=gsub("M1.*","M1",M)
M=gsub("MX.*",NA,M)

phenotype_label$M_stage <- M
phenotype_label$M_stage <- factor(phenotype_label$M_stage,
                                  levels = c("M0","M1"),
                                  labels = c(1,2))
head(phenotype_label)
#N
phenotype_label$N_stage <- gsub("^NX$",NA,phenotype_label$N_stage)

phenotype_label$N_stage <- factor(phenotype_label$N_stage,
                                  levels = c("N0","N1","N2","N3"),
                                  labels = c(1,2,3,4))
head(phenotype_label,10)
#T
T=phenotype_label$T_stage

T=gsub("T0.*","T0",T)

T=gsub("T1.*","T1",T)
T=gsub("T2.*","T2",T)
T=gsub("T3.*","T3",T)
T=gsub("T4.*","T4",T)
T=gsub("TX.*",NA,T)

phenotype_label$T_stage <- T

phenotype_label$T_stage<- factor(phenotype_label$T_stage,
                                 levels = c("T0","T1","T2","T3","T4"),
                                 labels = c(0,1,2,3,4))
head(phenotype_label,10)
#STAGE
Stage=phenotype_label$Stage

# Stage=gsub("stage ia","stage i",Stage)
# Stage=gsub("stage ib","stage i",Stage)
# 
# Stage=gsub("stage iia","stage ii",Stage)
# Stage=gsub("stage iib","stage ii",Stage)
# 
# Stage=gsub("stage iiia","stage iii",Stage)
# Stage=gsub("stage iiib","stage iii",Stage)

# Stage=gsub("stage iva","stage iv",Stage)
# Stage=gsub("stage ivb","stage iv",Stage)

phenotype_label$Stage <- Stage

phenotype_label$Stage <- factor(phenotype_label$Stage,
                                levels = c("stage i","stage ii","stage iii","stage iv"),
                                labels = c(1,2,3,4))
head(phenotype_label,10)
phenotype_label <- na.omit(phenotype_label)
head(phenotype_label)#144
write.csv(phenotype_label,
          file = 'phenotype_change_lable.csv',row.names = F)


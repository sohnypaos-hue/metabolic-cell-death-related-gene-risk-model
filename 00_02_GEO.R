#0.数据整理####
rm(list=ls())
library(GEOquery)
library(dplyr)
library(tidyverse)
library(limma)
library(pheatmap)
library(readxl)
setwd("/data/nas1/chenpeiru/44_GYZK-30212-7/00_raw_data/02_geo")
data_GEO <- read.table("GEO_clean.txt",header = T,sep = "\t")
##读取GPL注释文件
GPL <- read_xlsx("GPL6102.xlsx",sheet=1) %>% as.data.frame()
##读取分组信息
Group <- read_xlsx("group_survival.xlsx",sheet = 1) %>% as.data.frame()
head(Group)
# select 165 cancer---
data_GEO <- data_GEO[,c("ID_REF", Group$`!Sample_geo_accession`)]
#探针ID转换
colnames(GPL)[1] <- "ID_REF"  #需要讲GPL的探针列名和表达矩阵的对应
data_GEO <- merge(GPL,data_GEO,by="ID_REF")  #根据ID_REF进行合并
data_GEO <- data_GEO[,-1]  #去除探针列的信息
exp1 <- data_GEO %>% as_tibble() %>% 
  separate_rows("Symbol", sep = " /// ")
colnames(exp1)[1] <- "Symbol"
# class(exp1$GSM773540)
# datexpr2=as.data.frame(lapply(data_GEO[,2:118],as.numeric))
datexpr2 <- exp1
data <- aggregate(.~Symbol,datexpr2,mean)#相同基因取平均
#data <- data[,colSums(data>0)>nrow(data)*0.8]
dim(data)# 24357   166
#data <- data[complete.cases(data[,1]),]
data <- subset(data,nchar(data$Symbol)>0)
rownames(data) <- data[,1]
data <- data[,-1]
range(data)
write.table(data,file="matrix.txt",sep="\t",row.names=T,quote=F)
save(data,file = "matrix.rda")
# data <- data[,Group$`!Sample_geo_accession`]
# colnames(Group)
# table(Group$`!Sample_characteristics_ch1`)
# head(data)
# data <- log2(data+1)
# head(data)
#提取出癌症组----
rm(list=ls())
setwd("/data/nas1/chenpeiru/44_GYZK-30212-7/00_raw_data/02_geo")
load("matrix.rda")
data_GEO <- data
GSE <- t(data_GEO)
GSE <- data.frame(id.small = rownames(GSE),GSE) 
head(GSE)

survial <- read_xlsx('group_survival.xlsx',sheet=1) %>% as.data.frame()
head(survial)
colnames(survial) <- c("id.small","patientid","fustat","futime")


# combined expression and survival--------
expr_OS_geo <- merge(survial,GSE,by = "id.small")#
head(expr_OS_geo)
rownames(expr_OS_geo) <- expr_OS_geo$id.small
head(expr_OS_geo)
expr_OS_geo <- expr_OS_geo[,-c(1,2)]
head(expr_OS_geo)
#expr_OS_geo$fustat <- factor(expr_OS_geo$fustat,levels = c("dead","alive"),labels = c(1,0))
# survial <- data.frame(
#   New_Column_Name = rownames(survial),
#   survial,
#   row.names = NULL
# )
# head(survial)
# EXP_survival <- merge(data_t,survial,by = "New_Column_Name")
# EXP_survival$CDH17
# EXP_survival$FSCN1
# EXP_survival <- EXP_survival[,c(1,25443,25444,2:25442)]
# class(EXP_survival$status)
# EXP_survival$CDH17
# EXP_survival$FSCN1
# EXP_survival$status <- ifelse(EXP_survival$status == "D",1,0)

write.table(expr_OS_geo,file="GSE-Exp_Survival.txt",sep="\t",row.names=F,quote=F)
save(expr_OS_geo,file = "expr_OS_geo.rda")
##



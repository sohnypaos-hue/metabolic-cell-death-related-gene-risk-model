
#boxplot-----
rm(list = ls())
setwd("/data/nas1/chenpeiru/44_GYZK-30212-7/00_raw_data/01_TCGA")
phenotype_TCGA <- read.csv('phenotype_change_lable.csv',
                           header = T,sep = ',',row.names = 1)

head(phenotype_TCGA)
riskscore <- read.table("/data/nas1/chenpeiru/44_GYZK-30212-7/04_lasso/risk.train.group.txt",
                        header = T,sep = "\t",row.names = 1)
phenotype_TCGA_riskscore <- merge(phenotype_TCGA,riskscore,by="row.names")
colnames(phenotype_TCGA_riskscore)
phenotype_TCGA_riskscore <- phenotype_TCGA_riskscore[,c(1:10,19,20)]
head(phenotype_TCGA_riskscore)
rownames(phenotype_TCGA_riskscore) <- phenotype_TCGA_riskscore$Row.names
phenotype_TCGA_riskscore <- phenotype_TCGA_riskscore[,-1]
head(phenotype_TCGA_riskscore)
setwd("/data/nas1/chenpeiru/44_GYZK-30212-7/08_cli_boxplot")
save(phenotype_TCGA_riskscore,file = "phenotype_TCGA_riskscore.rda")
write.csv(phenotype_TCGA_riskscore,"phenotype_TCGA_riskscore.csv",quote = F,row.names = T)


#1. 将临床特征改为两组 ####
rm(list = ls())
load("/data/nas1/chenpeiru/44_GYZK-30212-7/08_cli_boxplot/phenotype_TCGA_riskscore.rda")
head(phenotype_TCGA_riskscore)
#stage:Ⅰ、Ⅱ、Ⅲ、Ⅳ
#M： M0， M1 12
#N：N0 N1 N2 N3 1234
#T: T0 T1 T2 T3 T4 01234
#MALE 1 FEMALE 2
#AGE <=65, >65
phenotype_TCGA_riskscore$Gender <- ifelse(phenotype_TCGA_riskscore$Gender%in%1,"Male","Female")
phenotype_TCGA_riskscore$M_stage <- ifelse(phenotype_TCGA_riskscore$M_stage %in% 1,"M0","M1")
phenotype_TCGA_riskscore$N_stage <- ifelse(phenotype_TCGA_riskscore$N_stage %in% c(1,2),"N0-N1","N2-N3")
phenotype_TCGA_riskscore$T_stage <- ifelse(phenotype_TCGA_riskscore$T_stage %in% c(0,1,2),"T0-T1-T2","T3-T4")
phenotype_TCGA_riskscore$Stage <- ifelse(phenotype_TCGA_riskscore$Stage%in%c(1,2),"Stage I-II","Stage III-IV")
phenotype_TCGA_riskscore$Age <- ifelse(phenotype_TCGA_riskscore$Age > 65,">65","<=65")
phenotype_TCGA_riskscore$Lymphnodes_positive <- ifelse(phenotype_TCGA_riskscore$Lymphnodes_positive > 0,
                                                       "Lymphnodes_Positive","Lymphnodes_Negative")

head(phenotype_TCGA_riskscore)

#2.1 KM-age-------
library(survival)
library("survminer")
library(tidyr)

#3.1 age----
library(tidyverse)
library(ggpubr)
library(ggsignif)
library(rstatix)
# load("E:/project/18_BJTC-572-2/07_independent_prognosis/TCGA_phenotype_survival.rda")
head(phenotype_TCGA_riskscore)
data_long <- phenotype_TCGA_riskscore
data_long$Age <- factor(data_long$Age)

violin_data_clinical <- data_long
head(violin_data_clinical)
violin_data_age <- violin_data_clinical[,c(1,10)]
#violin_data_gender$gender <- ifelse(violin_data_gender$gender == 1,'male','female')
violin_data_age <- violin_data_age[order(violin_data_age$Age),]
head(violin_data_age,20)
compaired <- list(c("<=65",">65"))
#setwd("/data/nas1/chenpeiru/41_YQNC-10615-6/09_cli_feature")
library(showtext)
#font_add('Arial','/Library/Fonts/Arial.ttf') #加载字体，MAC 中字体库在 /Library/Fonts
#showtext_auto() #自动调用showtext，否则无法在ggsave()中使用，因为ggsave会自动打开和关闭图形设备。

theme_zg <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          panel.background=element_rect(fill='transparent',color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),#去网格线
          axis.line = element_line(colour = "black"),
          #axis.title.x = element_blank(),#去x轴标签
          axis.title.y=element_text(face = "bold",size = 20),#y轴标签加粗及字体大小
          axis.title.x=element_text(face = "bold",size = 20),#X轴标签加粗及字体大小
          axis.text.y = element_text(face = "bold",size = 16),#y坐标轴刻度标签加粗
          axis.text.x = element_text(face = "bold",size = 16, vjust = 1, hjust = 1, angle = 45),#x坐标轴刻度标签加粗
          axis.ticks = element_line(color='black'),
          # axis.ticks.margin = unit(0.8,"lines"),
          title = element_text(size=16),
          legend.title=element_blank(),
          # legend.position=c(0.9, 0.8),#图例在绘图区域的位置
          legend.position="right",
          legend.box = "vertical",
          legend.direction = "vertical",
          legend.text = element_text(face = "bold",size = 18,margin = margin(r=8))
          # legend.background = element_rect( linetype="solid",colour ="black")
    )
}
pdf('riskscore_age.pdf',width = 7,height = 7)
par(family="serif")
ggplot(violin_data_age,aes(x=Age,y=riskScore),palette = "jco", add = "jitter")+
  geom_violin(aes(fill=Age))+
  geom_boxplot(width=0.1,cex=1.2)+
  labs(x = "Age", y = "riskscore")+theme_zg()+
  scale_fill_manual(values = c("#104E8B","#E85600FF")) +
  geom_signif(comparisons = compaired,step_increase = 0.3,map_signif_level = F,test = wilcox.test,textsize = 5)
dev.off()
#3.2gender--------
head(data_long)
#violin_data_N <- data_long[,c(6,4)]
violin_data_clinical <- data_long
violin_data_gender <- violin_data_clinical[,c(2,10)]
#violin_data_gender$gender <- ifelse(violin_data_gender$gender == 1,'male','female')
violin_data_gender <- violin_data_gender[order(violin_data_gender$Gender),]
head(violin_data_gender)
compaired <- list(c("Female","Male"))
pdf('riskscore_gender.pdf',width = 7,height = 7)
par(family="serif")
ggplot(violin_data_gender,aes(x=Gender,y=riskScore),palette = "jco", add = "jitter")+
  geom_violin(aes(fill=Gender))+
  geom_boxplot(width=0.1,cex=1.2)+
  labs(x = "gender", y = "riskscore")+
  scale_fill_manual(values = c("#104E8B","#E85600FF")) + theme_zg()+
  geom_signif(comparisons = compaired,step_increase = 0.3,map_signif_level = F,test = wilcox.test,textsize = 5)
dev.off()
#3.3 M--------
#pathologic_N
head(data_long,20)
violin_data_M <- data_long[,c(3,10)]
# violin_data_N$pathologic_N <- ifelse(violin_data_N$pathologic_N == 0,'N0',
#                                      ifelse(violin_data_N$pathologic_N == 1,'N1',
#                                             ifelse(violin_data_N$pathologic_N == 2,'N2','N3')))

violin_data_M <- violin_data_M[order(violin_data_M$M_stage),]

compaired <- list(c("M0","M1"))
pdf('riskscore_M-stage.pdf',width = 7,height = 7)
par(family="serif")
ggplot(violin_data_M,aes(x=M_stage,y=riskScore),palette = "jco", add = "jitter")+
  geom_violin(aes(fill=M_stage))+
  geom_boxplot(width=0.1,cex=1.2)+
  labs(x = "M_stage", y = "riskscore")+
  scale_fill_manual(values = c("#104E8B","#E85600FF")) + theme_zg()+
  geom_signif(comparisons = compaired,step_increase = 0.3,map_signif_level = F,
              test = wilcox.test,textsize = 5)
dev.off()
#3.4 N--------
#pathologic_N
head(data_long,10)
violin_data_N <- data_long[,c(4,10)]
# violin_data_N$pathologic_N <- ifelse(violin_data_N$pathologic_N == 0,'N0',
#                                      ifelse(violin_data_N$pathologic_N == 1,'N1',
#                                             ifelse(violin_data_N$pathologic_N == 2,'N2','N3')))

violin_data_N <- violin_data_N[order(violin_data_N$N_stage),]

compaired <- list(c("N0-N1","N2-N3"))
pdf('riskscore_N_stage.pdf',width = 7,height = 7)
par(family="serif")

ggplot(violin_data_N,aes(x=N_stage,y=riskScore),palette = "jco", add = "jitter")+
  geom_violin(aes(fill=N_stage))+
  geom_boxplot(width=0.1,cex=1.2)+
  labs(x = "N_stage", y = "riskscore")+
  scale_fill_manual(values = c("#104E8B","#E85600FF")) + theme_zg()+
  geom_signif(comparisons = compaired,step_increase = 0.3,map_signif_level = F,
              test = wilcox.test,textsize = 5)
dev.off()
#3.4 T--------
#T----
head(data_long,10)
violin_data_T <- data_long[,c(5,10)]
# violin_data_N$pathologic_N <- ifelse(violin_data_N$pathologic_N == 0,'N0',
#                                      ifelse(violin_data_N$pathologic_N == 1,'N1',
#                                             ifelse(violin_data_N$pathologic_N == 2,'N2','N3')))

violin_data_T <- violin_data_T[order(violin_data_T$T_stage),]

compaired <- list(c("T0-T1-T2","T3-T4"))
pdf('riskscore_T_stage.pdf',width = 7,height = 7)
par(family="serif")
ggplot(violin_data_T,aes(x=T_stage,y=riskScore),palette = "jco", add = "jitter")+
  geom_violin(aes(fill=T_stage))+
  geom_boxplot(width=0.1,cex=1.2)+
  labs(x = "T_stage", y = "riskscore")+
  scale_fill_manual(values = c("#104E8B","#E85600FF")) + theme_zg()+
  geom_signif(comparisons = compaired,step_increase = 0.3,map_signif_level = F,
              test = wilcox.test,textsize = 5)
dev.off()
#3.3stage--------
head(data_long,10)
violin_data_stage <- data_long[,c(6,10)]
# violin_data_N$pathologic_N <- ifelse(violin_data_N$pathologic_N == 0,'N0',
#                                      ifelse(violin_data_N$pathologic_N == 1,'N1',
#                                             ifelse(violin_data_N$pathologic_N == 2,'N2','N3')))

violin_data_stage <- violin_data_stage[order(violin_data_stage$Stage),]
# theme_v <- function(..., bg='white'){
#   require(grid)
#   theme_classic(...) +
#     theme(rect=element_rect(fill=bg),
#           plot.margin=unit(rep(0.5,4), 'lines'),
#           panel.background=element_rect(fill='transparent',color='black'),
#           panel.border=element_rect(fill='transparent', color='transparent'),
#           panel.grid=element_blank(),
#           axis.line = element_line(colour = "black"),
#           #axis.title.x = element_blank(),
#           axis.title.y=element_text(face = "bold",size = 18),
#           axis.title.x=element_text(face = "bold",size = 18),
#           axis.text.y = element_text(face = "bold",size = 16),
#           axis.text.x = element_text(face = "bold",size = 16, vjust = 1, hjust = 1, angle = 45),
#           axis.ticks = element_line(color='black'),
#           # axis.ticks.margin = unit(0.8,"lines"),
#           legend.title=element_blank(),
#           # legend.position=c(0.9, 0.8),
#           legend.position="top",
#           legend.box = "vertical",
#           # legend.direction = "vertical",
#           legend.text = element_text(face = "bold",size = 18,margin = margin(r=8))
#           # legend.background = element_rect( linetype="solid",colour ="black")
#     )
# }
compaired <- list(c("Stage I-II","Stage III-IV"))
pdf('riskscore_stage.pdf',width = 7,height = 7)
par(family="serif")
ggplot(violin_data_stage,aes(x=Stage,y=riskScore),palette = "jco", add = "jitter")+
  geom_violin(aes(fill=Stage))+
  geom_boxplot(width=0.1,cex=1.2)+
  labs(x = "Stage", y = "riskscore")+
  scale_fill_manual(values = c("#104E8B","#E85600FF","#6D2E71FF")) +theme_zg()+
  geom_signif(comparisons = compaired,step_increase = 0.3,map_signif_level = F,
              test = wilcox.test,textsize = 5)
dev.off()
#Lymphnodes_positive------
head(data_long,10)
violin_data_N <- data_long[,c(7,10)]
# violin_data_N$pathologic_N <- ifelse(violin_data_N$pathologic_N == 0,'N0',
#                                      ifelse(violin_data_N$pathologic_N == 1,'N1',
#                                             ifelse(violin_data_N$pathologic_N == 2,'N2','N3')))

violin_data_N <- violin_data_N[order(violin_data_N$Lymphnodes_positive),]

compaired <- list(c("Lymphnodes_Negative","Lymphnodes_Positive"))
pdf('riskscore_Lymphnodes_positive.pdf',width = 7,height = 7)
par(family="serif")

ggplot(violin_data_N,aes(x=Lymphnodes_positive,y=riskScore),palette = "jco", add = "jitter")+
  geom_violin(aes(fill=Lymphnodes_positive))+
  geom_boxplot(width=0.1,cex=1.2)+
  labs(x = "Lymphnodes_positive", y = "riskscore")+
  scale_fill_manual(values = c("#104E8B","#E85600FF")) + theme_zg()+
  geom_signif(comparisons = compaired,step_increase = 0.3,map_signif_level = F,
              test = wilcox.test,textsize = 5)
dev.off()
#heatmap------
rm(list = ls())
#heatmap---
library(tidyverse)
library(ComplexHeatmap)
library(ggsci) #颜色
library(circlize) #
#setwd("/data/nas1/chenpeiru/41_YQNC-10615-6/09_cli_feature")
load("phenotype_TCGA_riskscore.rda")

phenotype_TCGA_riskscore$Gender <- ifelse(phenotype_TCGA_riskscore$Gender%in%1,"Male","Female")
phenotype_TCGA_riskscore$M_stage <- ifelse(phenotype_TCGA_riskscore$M_stage %in% 1,"M0","M1")
phenotype_TCGA_riskscore$N_stage <- ifelse(phenotype_TCGA_riskscore$N_stage %in% c(1,2),"N0-N1","N2-N3")
phenotype_TCGA_riskscore$T_stage <- ifelse(phenotype_TCGA_riskscore$T_stage %in% c(0,1,2),"T0-T1-T2","T3-T4")
phenotype_TCGA_riskscore$Stage <- ifelse(phenotype_TCGA_riskscore$Stage%in%c(1,2),"Stage I-II","Stage III-IV")
phenotype_TCGA_riskscore$Age <- ifelse(phenotype_TCGA_riskscore$Age > 65,">65","<=65")
phenotype_TCGA_riskscore$Lymphnodes_positive <- ifelse(phenotype_TCGA_riskscore$Lymphnodes_positive > 0,
                                                       "Lymphnodes_Positive","Lymphnodes_Negative")
phenotype_TCGA_riskscore$fustat <- ifelse(phenotype_TCGA_riskscore$fustat==1,
                                          "Dead","Alive")
head(phenotype_TCGA_riskscore)


riskScore_cli4 <- phenotype_TCGA_riskscore[order(phenotype_TCGA_riskscore$group), ]
colnames(riskScore_cli4)
riskScore_cli4 <- riskScore_cli4[,-10]#delete riskscore
head(riskScore_cli4)
# 提取想展示的临床数据
# riskScore_cli2 <- riskScore_cli2 %>% 
#   select(riskScore:tumor_stage,Age) %>% 
#   select(- "age")
# 构建列注释块

#
col_fun_time <- colorRamp2(
  c(0, 10, 20),  #根据值的范围设置
  c("#DC0000FF", "grey", "#1f78b4")
)
#
riskScore_cli4$group <- factor(riskScore_cli4$group,)
head(riskScore_cli4,20)
# colnames(riskScore_cli4)[6] <- "Stage"
# colnames(riskScore_cli4)[4] <- "N_stage"
# colnames(riskScore_cli4)[5] <- "T_stage"
# colnames(riskScore_cli4)[7] <- "OS.Time"
# colnames(riskScore_cli4)[10] <- "riskgroup"
# head(riskScore_cli4,10)

#order(riskScore_cli4$riskScore)
# head(riskScore_cli4)
ha=HeatmapAnnotation(df=riskScore_cli4)
# 构建zero矩阵
zero_row_mat=matrix(nrow=0, ncol=nrow(riskScore_cli4))
colnames(riskScore_cli4)
# ha <- HeatmapAnnotation(
#   Risk = riskScore_cli4$group,
#   Stage = riskScore_cli4$Stage,
#   OS.Status = riskScore_cli4$fustat,
#   OS.Time = riskScore_cli4$futime,
#   Gender = riskScore_cli4$Gender ,
#   Age = riskScore_cli4$Age,
#   M_stage <- riskScore_cli4$M_stage,
#   N_stage <- riskScore_cli4$N_stage,
#   T_stage <- riskScore_cli4$T_stage,
#   Lymphnodes_positive <- riskScore_cli4$Lymphnodes_positive,
#   col = list( 
#     Risk = c("high" = "#BC3C29FF", "low" = "#0072B5FF"),
#     OS.Status = c("Alive" = "#E18727FF", "Dead" = "#20854EFF"), #分类
#     OS.Time = col_fun_time , #连续
#     Gender = c("Female" = "#AB3282", "Male" = "#3A6963"),
#     Age = c("<=65" = "#712820", ">65" = "#E4C755"),
#     Stage = c("Stage I-II" = "#E64B35FF", "Stage III-IV" = "#4DBBD5FF"),
#     M_stage = c("M0" = "#FFB6C1", "M1" = "#87CEFA"),
#     N_stage = c("N0-N1" = "#DB7093", "N1-N2" = "#7B68EE"),
#     T_stage = c("T0-T1-T2" = "#DA70D6", "T3-T4" = "#4DBBD5FF"),
#     Lymphnodes_positive = c("Lymphnodes_Positive" = "#DA70D6", "Lymphnodes_Negative" = "#E4C755")
#      )
# )
head(ha)
head(zero_row_mat)
pdf('heatmap.pdf',width = 8,height = 6)
Hm <- Heatmap(zero_row_mat, top_annotation=ha)
draw(Hm, merge_legend = TRUE, 
     heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom",
     width = unit(16, "cm"), height = unit(1, "cm")
)
dev.off()



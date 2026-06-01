rm(list = ls())

setwd("/data/nas1/chenpeiru/44_GYZK-30212-7/05_test")
# if (! dir.exists("./05_Verif")){
#   dir.create("./05_Verif")
# }
# setwd("./05_Verif")


library(lance)
library(readxl)
library(readr)
library(tidyverse)
###读取数据——————
# 01 获取数据集----------
#disease <- "BLCA"
load("/data/nas1/chenpeiru/44_GYZK-30212-7/00_raw_data/02_geo/expr_OS_geo.rda")
expr_OS_geo$futime <- expr_OS_geo$futime*30
expr_OS_geo <- subset(expr_OS_geo,futime>30)
# survival<-read.csv(file = '../00_Rawdata/survival(GSE32894).csv',header=TRUE, row.names=1,)
# data <- read.csv("../00_Rawdata/dat(GSE32894).csv",header=TRUE, row.names=1, check.names=FALSE)
# survival_dat<-t(data)
# survival_dat<-survival_dat[rownames(survival_dat)%in%survival$sample,]
###候选基因
lasso_geneids <- readRDS("../04_lasso/lasso_geneids.rds")

## 合并生存数据
train_dat <- expr_OS_geo[,c("fustat","futime",lasso_geneids)]
# train_dat<-as.data.frame(train_dat)
# train_dat$sample<-rownames(train_dat)
# train_dat<-merge(survival,train_dat,by='sample')
# rownames(train_dat)<-train_dat$sample
# train_data<-train_dat[,-1]
# train_data <- na.omit(train_data)
# colnames(train_data)


###风险系数
outCol <- readRDS("../04_lasso/outCol.rds")
coef.min <- readRDS("../04_lasso/coef.min.rds")
coef.min <- coef.min[lasso_geneids,]

#风险模型的构建与验证------

train_dat <- train_dat[,-c(1:2)]
risk_out <- train_dat
#risk_out$risk_out <- NA
risk_out$riskScore_out <- NA

#计算riskscore
cnt <- 1
while (cnt < length(rownames(train_dat))+1) {
  risk_out$riskScore_out[cnt] <- sum(coef.min*train_dat[cnt,])
  cnt = cnt + 1
}
class(risk_out$riskScore_out)

head(risk_out)


##合并生存数据
survival <- expr_OS_geo[,c(1:2)]
head(survival)
#risk_out$sample <- rownames(risk_out)
risk_out2<-merge(survival,risk_out,by='row.names')
#rownames(risk_out2) <- risk_out2$sample
head(risk_out2)
rownames(risk_out2) <- risk_out2$Row.names
risk_out2 <- risk_out2[,-1]
colnames(risk_out2)
colnames(risk_out2)[11] <- "riskScore"
write.table(risk_out2, file = "risk.test.txt", sep = "\t", quote = F, row.names = T,
            col.names = T)
# ############    ROC曲线####
setwd("/data/nas1/chenpeiru/44_GYZK-30212-7/05_test")
library(survivalROC)
heatmap_train = read.table("risk.test.txt",header=T,sep="\t",
                           check.names=F,row.names=1)

colnames(heatmap_train)

heatmap_train_tumor <- heatmap_train
head(heatmap_train_tumor)
pdf(file = "04_ROC1-2-3.pdf",family="serif",height = 10,width = 10)
par(oma = c(2, 2.5, 2, 2.5), mar = c(3, 2, 2, 0.95),
    fig = c(0.1, 0.95, 0.1, 0.9), xpd = NA, family="serif")
roc=survivalROC(Stime=heatmap_train_tumor$futime, status=heatmap_train_tumor$fustat,
                marker = heatmap_train_tumor$riskScore,predict.time =1*365, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='#FA8072',
     xlab="False positive rate", ylab="True positive rate",
     #main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     main="ROC curve",
     lwd = 2, cex.main=3, cex.lab=3, cex.axis=2, font.lab=1)
abline(0,1)
aucText=c()
rocCol <- c('#FA8072','#63B8FF','#FFC1C1','#ADFF2F','#FFFF00')
aucText=c(aucText,paste0("1 year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
j =0
for (i in c(2,3)){
  roc1=survivalROC(Stime=heatmap_train_tumor$futime, status=heatmap_train_tumor$fustat,
                   marker = heatmap_train_tumor$riskScore,predict.time =i*365, method="KM")
  j=j+1
  aucText=c(aucText,paste0(i," year"," (AUC=",sprintf("%.3f",roc1$AUC),")"))
  lines(roc1$FP, roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j+1],lwd = 3)
}
legend("bottomright",aucText, lwd=2,bty="n",col=rocCol,cex = 2)
abline(0,1)
dev.off()


########KM####better
library(survival)
library("survminer")
library(tidyr)
risk.train = read.table("/data/nas1/chenpeiru/44_GYZK-30212-7/05_test/risk.test.txt",
                        header=T,sep="\t",row.names=1)
head(risk.train)
risk.train$futime=risk.train$futime/365
#median----
# group <- as.vector(ifelse(risk.train$riskScore> median(risk.train$riskScore),"high","low"))
#bestthreshold.surv---
bestthreshold.surv <- surv_cutpoint(risk.train,time = "futime",
                                    event = "fustat",
                                    'riskScore',
                                    minprop = 0.35,
                                    progressbar = TRUE)
bestthreshold.train <- bestthreshold.surv$cutpoint$cutpoint#4.433
bestthreshold.train
group <- ifelse(risk.train$riskScore > bestthreshold.train, 'high','low')
#here
diff_train=survdiff(Surv(futime, fustat) ~group,data = risk.train)
pValue=1-pchisq(diff_train$chisq,df=1)
pValue=round(pValue,5)
fit_train <- survfit(Surv(futime, fustat) ~ group, data = risk.train)
summary(fit_train)
risk.train.group <- cbind(risk.train,group)
risk.train.group <- data.frame(id=row.names(risk.train.group),risk.train.group)
write.table(risk.train.group,file="risk.test.group.txt",
            sep="\t",quote=F,row.names=F)

####        KM曲线
pdf(file="survival.coefficient.test.pdf",family="serif",height = 9,width = 9)
ggsurvplot(fit_train,
           conf.int=TRUE,
           pval=TRUE,
           risk.table=TRUE,
           xlab="time (years)",
           #legend.labs=legend.labs,
           legend.title="Risk",
           #linetype = lty,
           palette = c( "red", "blue"),
           title="Kaplan-Meier Curve for Survival",
           risk.table.height=.15) %>% print()
dev.off()


###############  风险曲线####
loc_train = match(c('futime','fustat','riskScore'),colnames(heatmap_train))
datalast_train = heatmap_train[,loc_train]
datalast_train = data.frame(datalast_train)
datalast_train$futime=datalast_train$futime/365
#datalast_train$riskScore[which(datalast_train$riskScore >50)] <- 50

pdf(file="risk curve_TCGA_test.pdf",height = 12,width = 12)
par(mfrow=c(2,1),family="serif")
par(mar=c(2,4.5,1,1))
col=c()
col[sort(datalast_train$riskScore) <= bestthreshold.train]="blue"
col[sort(datalast_train$riskScore) > bestthreshold.train]="red"
plot(sort(datalast_train$riskScore),axes=F,xlab = NA,ylab = "Risk Score",col=col,
     mgp=c(2.5,1,0),cex.lab=2)
box(lwd=2)
abline(v = length(which(sort(datalast_train$riskScore) <=  bestthreshold.train))+0.5,lty="dashed")
abline(h =  bestthreshold.train,lty="dashed")
text(100, bestthreshold.train+0.5,paste0('threshold = ',round( bestthreshold.train,4)),cex = 1.5,font=2)
axis(2,seq(0,50,1),cex.axis=2)
axis(1,seq(0,nrow(datalast_train),50),cex.axis=2)
legend("topleft", c("High risk","Low risk"), pch=16:16, col=c("red","blue"),cex = 1.5)
#图2
datalastSORT = datalast_train[order(datalast_train[,"riskScore"]),]
col=c()
col[datalastSORT[,"fustat"]==1]= "red"
col[datalastSORT[,"fustat"]==0]= "blue"
par(mar=c(2,4.5,1,1))
plot(datalastSORT[,"futime"],col=col,pch=16,axes=F,xlab = NA,mgp=c(2.5,1,0),cex.lab=1.8,
     ylab = "Following up (years)")
box(lwd=2)
abline(v = length(which(sort(datalast_train$riskScore) <= bestthreshold.train))+0.5,lty="dashed")
axis(2,seq(0,max(datalastSORT[,"futime"]),2),cex.axis=2)
axis(1,seq(0,nrow(datalast_train),10),cex.axis=2)
legend("topleft", c("Dead","Alive"), pch=16:16, col=c("red","blue"),cex = 1.5)

dev.off()
###########热图####
library(pheatmap)

risk.train = read.table("risk.test.group.txt",header=T,sep="\t",row.names=1)
head(risk.train)
risk.train <- risk.train[order(risk.train$group),]
head(risk.train)

#risk.train$futime=risk.train$futime/365
#risk.train$group <- as.vector(ifelse(risk.train$riskScore> median(risk.train$riskScore),"high","low"))
head(risk.train)
exp.gene <- t(risk.train[,c(3:10)])
head(exp.gene)
group <- data.frame(
  New_Column_Name = rownames(risk.train),
  risk.train,
  row.names = NULL
)
head(group)
group <- group[,c(1,13)]
head(group)
group <- group[order(group$group),]
head(group)
annotation_col = data.frame(
  Group = factor(group$group)
)
rownames(annotation_col) = group$New_Column_Name
color.key <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")
ann_colors = list(
  Group = c(low="royalblue1", high="#BC3C29FF"))
head(exp.gene)
pdf(file = '07_riskscore_heatmap.pdf',width = 12,height = 9)

pheatmap(exp.gene,cluster_cols=F,show_colnames=F,scale="row",annotation_col = annotation_col,
         color = colorRampPalette(color.key)(50),
         annotation_colors = ann_colors,fontsize = 20)


dev.off()#manual save 
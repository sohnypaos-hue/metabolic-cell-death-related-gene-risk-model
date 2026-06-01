rm(list = ls())
setwd("/data/nas1/chenpeiru/44_GYZK-30212-7/04_lasso")
library(lance)
library(readxl)
library(readr)
library(tidyverse)
###读取数据——————
# 01 获取数据集----------
cox_train <- read.table("/data/nas1/chenpeiru/44_GYZK-30212-7/03_cox/hub_expr.txt",
                        header=T,row.names=1)
gene <- read.table("/data/nas1/chenpeiru/44_GYZK-30212-7/03_cox/uniCoxResult.Sig.txt",
                   header=T,sep="\t",row.names=1,check.names=F)
train_dat <- cox_train[,c("fustat","futime",rownames(gene))]


### LASSO
library(survival)
library(survminer)
library(glmnet)
train_data<-train_dat
x_all <- as.matrix(train_data[,colnames(train_data)[-(1:2)]])
y_all <- data.matrix(Surv(train_data$futime,train_data$fustat))

set.seed(5566)#123456
cvfit <- cv.glmnet(x_all, y_all, family="cox", maxit = 1000)
fit.train <- cvfit$glmnet.fit
coef <- coef(cvfit, s = cvfit$lambda.min)
fit <- glmnet(x_all,y_all,family = "cox")
cvfit$lambda.min# 0.0329713
x <- coef(fit)  
tmp <- as.data.frame(as.matrix(x))
tmp$coef <- row.names(tmp)
tmp <- reshape::melt(tmp, id = "coef")
tmp$variable <- as.numeric(gsub("s", "", tmp$variable))
tmp$coef <- gsub('_','-',tmp$coef)
tmp$lambda <- fit$lambda[tmp$variable+1] # extract the lambda values
#tmp$norm <- apply(abs(x[-1,]), 2, sum)[tmp$variable+1] # compute L1 norm 
head(tmp)

log(cvfit$lambda.min) #-
log(cvfit$lambda.1se) #-
library(ggsci)
p.coef<-ggplot(tmp,aes(log(lambda),value,color = coef)) + 
  geom_vline(xintercept = log(cvfit$lambda.min),size=0.8,color='grey60',alpha=0.8,linetype=2)+
  geom_vline(xintercept = log(cvfit$lambda.1se),size=0.8,color='grey60',alpha=0.8,linetype=2)+
  geom_line(size=1) + 
  xlab("Lambda (log scale)") + 
  #xlab("L1 norm")+
  ylab('Coefficients')+
  theme_bw(base_rect_size = 2)+ 
  scale_color_manual(values = c(pal_npg()(10),pal_d3()(10),pal_jco()(7)))+
  scale_x_continuous(expand = c(0.01,0.01))+
  scale_y_continuous(expand = c(0.01,0.01))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=15,color='black'),
        axis.text = element_text(size=12,color='black'),
        legend.title = element_blank(),
        legend.text = element_text(size=12,color='black'),
        legend.position = 'right')+
  annotate('text',x = -4,y=0.3,label='Lambda.min = -3.998',color='black',size=5)+
  annotate('text',x = -3,y=0.2,label='Lambda.lse = -2.509',color='black',size=5)+
  guides(col=guide_legend(ncol = 1))
p.coef
ggsave('lasso.coefficients.venalty.pdf',p.coef,family="serif",w=8,h=7)
ggsave('lasso.coefficients.venalty.png',p.coef,w=8,h=7)


xx <- data.frame(lambda=cvfit[["lambda"]],cvm=cvfit[["cvm"]],cvsd=cvfit[["cvsd"]],
                 cvup=cvfit[["cvup"]],cvlo=cvfit[["cvlo"]],nozezo=cvfit[["nzero"]])
xx$ll <- log(xx$lambda)
xx$NZERO <- paste0(xx$nozezo,' vars')

a<-ggplot(xx,aes(ll,cvm,color=NZERO))+
  geom_errorbar(aes(x=ll,ymin=cvlo,ymax=cvup),width=0.05,size=1)+
  geom_vline(xintercept = log(cvfit$lambda.min),size=0.8,color='grey60',alpha=0.8,linetype=2)+
  geom_vline(xintercept = log(cvfit$lambda.1se),size=0.8,color='grey60',alpha=0.8,linetype=2)+
  geom_point(size=2)+
  xlab("Log Lambda")+ylab('Partial Likelihood Deviance')+
  theme_bw(base_rect_size = 1.5)+ 
  scale_color_manual(values = c(pal_npg()(10),pal_jco()(8)))+
  scale_x_continuous(expand = c(0.02,0.02))+
  scale_y_continuous(expand = c(0.02,0.02))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=15,color='black'),
        axis.text = element_text(size=12,color='black'),
        legend.title = element_blank(),
        legend.text = element_text(size=12,color='black'),
        legend.position = 'bottom')+
  annotate('text',x = -4,y=12.3,label='Lambda.min = -3.998',color='black',size=5)+
  annotate('text',x = -3,y=12.5,label='Lambda.lse = -2.509',color='black',size=5)+
  
  guides(col=guide_legend(ncol = 3))

a
ggsave('lasso_verify.pdf',a,family="serif",w=8,h=7)
ggsave('lasso_verify.png',a,w=8,h=7)

index <- which(as.numeric(coef) != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
coef.min = coef(cvfit, s = "lambda.min") 

write.csv(lassoGene, "lasso_genes.csv")
saveRDS(lassoGene, file="lasso_geneids.rds")
saveRDS(coef.min, file="coef.min.rds")

riskScore=predict(cvfit,newx = as.matrix(x_all),s=cvfit$lambda.min) #predict(模型对象，newx数据框拟合值)
riskScore<-as.numeric(riskScore) #转换为数值类型
class(riskScore) #获取R语言对象的类别信息
coxGene=lassoGene
outCol=c("fustat","futime",coxGene)
outCol
saveRDS(outCol, file="outCol.rds")
cox.data.step <- train_data[,outCol]
cox.data.plot <- cbind(cox.data.step,riskScore)
cox.data.plot$futime <- as.numeric(as.character(cox.data.plot$futime))
cox.data.plot$riskScore <- as.numeric(cox.data.plot$riskScore)
cox.data.plot <- cox.data.plot[order(cox.data.plot$riskScore),]
write.table(cox.data.plot, file = "risk.train.txt", sep = "\t", quote = F, row.names = T,
            col.names = T)
##
# ############    ROC曲线####
rm(list=ls())
setwd("/data/nas1/chenpeiru/44_GYZK-30212-7/04_lasso")
library(survivalROC)
heatmap_train = read.table("risk.train.txt",header=T,sep="\t",
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
risk.train = read.table("/data/nas1/chenpeiru/44_GYZK-30212-7/04_lasso/risk.train.txt",
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
bestthreshold.train <- bestthreshold.surv$cutpoint$cutpoint#1.131085
bestthreshold.train#1.498
group <- ifelse(risk.train$riskScore > bestthreshold.train, 'high','low')
#here
diff_train=survdiff(Surv(futime, fustat) ~group,data = risk.train)
pValue=1-pchisq(diff_train$chisq,df=1)
pValue=round(pValue,5)
fit_train <- survfit(Surv(futime, fustat) ~ group, data = risk.train)
summary(fit_train)
risk.train.group <- cbind(risk.train,group)
risk.train.group <- data.frame(id=row.names(risk.train.group),risk.train.group)
write.table(risk.train.group,file="risk.train.group.txt",
            sep="\t",quote=F,row.names=F)

####        KM曲线
pdf(file="survival.coefficient.train.pdf",family="serif",height = 9,width = 9)
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

pdf(file="risk curve_TCGA_tarin.pdf",height = 12,width = 12)
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
text(250, bestthreshold.train+0.5,paste0('threshold = ',round( bestthreshold.train,4)),cex = 1.5,font=2)
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

risk.train = read.table("risk.train.group.txt",header=T,sep="\t",row.names=1)
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

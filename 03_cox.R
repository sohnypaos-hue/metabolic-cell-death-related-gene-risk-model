setwd("/data/nas1/chenpeiru/44_GYZK-30212-7/03_cox")
# if (! dir.exists("./04_cox")){
#   dir.create("./04_cox")
# }
# setwd("./04_cox")
rm(list = ls())
load("/data/nas1/chenpeiru/44_GYZK-30212-7/00_raw_data/01_TCGA/dat_fpkm_tumor406_survival.rda")

colnames(fpkm_tumor_p_t_survival)
colnames(fpkm_tumor_p_t_survival)[1] <- "fustat"
colnames(fpkm_tumor_p_t_survival)[2] <- "futime"
fpkm_tumor_p_t_survival <- subset(fpkm_tumor_p_t_survival,futime>30)#399

hub <- read.csv("/data/nas1/chenpeiru/44_GYZK-30212-7/02_venn/DE-MCD-50.csv",
                  sep = ",",header = T)
head(hub)
hub_expr <- fpkm_tumor_p_t_survival[,colnames(fpkm_tumor_p_t_survival)%in%c(hub$X1,"futime","fustat")]

write.table(hub_expr,"hub_expr.txt",sep = "\t",row.names = T,quote = F)
save(hub_expr,file = "hub_expr.rda")
#test---
# load("/data/nas1/chenpeiru/41_YQNC-10615-6/00_raw_data/02_geo/expr_OS_geo.rda")
# validSet <- expr_OS_geo
# rm(expr_OS_geo)
# validSet$futime <- validSet$futime*30
# #validSet <- subset(validSet,futime>90)
# hub <- read.table("/data/nas1/chenpeiru/41_YQNC-10615-6/02_DE_SC_IS_Venn/scRNA_TCGA_IS-44.txt",
#                   sep = "\t",header = T)
# head(hub)
# #validSet$OBFC2B
# hub_expr <- validSet[,colnames(validSet)%in%c(hub$X1,"futime","fustat")]
# #colnames(hub_expr)[8] <- "NABP2"
# survival_GEO <- hub_expr
# save(survival_GEO,file = "survival_GEO.rda")


#hub_expr <- t(hub_expr) 
#hub_expr <- as.data.frame(hub_expr)
#hub_expr1 <- data.frame(id.small = rownames(hub_expr),hub_expr)
#######################    单因素####----
rm(list = ls())
library(survival)
setwd("/data/nas1/chenpeiru/44_GYZK-30212-7/03_cox")
set.seed(386)
cox_train <- read.table("hub_expr.txt",header=T,row.names=1)

colnames(cox_train)
#cox_train[,3:ncol(cox_train)] <- apply(cox_train[,3:ncol(cox_train)],2,function(x){log2(x+1)})
#head(cox_train)
#pFilter=0.05 #定义单因素显著性
pFilter=0.01#0.001
#colnames(cox_train)[1:2] <- c("fustat","futime") 
head(cox_train)
outTab=data.frame()
for(i in colnames(cox_train[,3:ncol(cox_train)])){
  cox <- coxph(Surv(futime, fustat) ~ cox_train[,i], data = cox_train)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=i,
                     z=coxSummary$coefficients[,"z"],
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}

outTab = outTab[is.na(outTab$pvalue)==FALSE,]
outTab=outTab[order(as.numeric(as.vector(outTab$pvalue))),]
#write.table(outTab,file="uniCoxResult.txt",sep="\t",row.names=F,quote=F)
#输出单因素显著的结果
sigTab=outTab[as.numeric(as.vector(outTab$pvalue))<pFilter,]

write.table(sigTab,file="01_single_uniCoxResult.Sig.txt",sep="\t",row.names=F,quote=F)
#输出单因素显著AS的PSI值，用于后续建模
sigGenes=c("futime","fustat")
sigGenes=c(sigGenes,as.vector(sigTab[,1]))
uniSigExp=cox_train[,sigGenes]# genes
#PH
# ph_res <- data.frame()
# for(i in 3:ncol(uniSigExp)){
#   fml <- as.formula(paste0('Surv(futime,fustat)~',colnames(uniSigExp)[i]))
#   cox.zph.fit <- cox.zph(coxph(fml, data=uniSigExp[,c(1:2,i)],id = rownames(uniSigExp)))$table
#   ph_res <- rbind(ph_res,cox.zph.fit[1,])
# }
# colnames(ph_res) <- c('chisq','df','p')
# rownames(ph_res) <- colnames(uniSigExp)[-(1:2)]
# ph_res
# write.csv(ph_res,"PH_results.csv",quote = F)
# # fml <- as.formula(paste0('Surv(futime,fustat)~',paste0(colnames(uniSigExp)[-c(1:2)],collapse = '+')))
# # f <- coxph(fml, data=uniSigExp,id = rownames(uniSigExp))
# # # summary(cox)
# # cox.zph.fit <- cox.zph(f)
# # cox.zph.fit
# # pdf('/data/nas1/chenpeiru/29_YQNJ-10132-5/05_cox/PH_analysis.pdf',width = 20,height = 16)
# # ggcoxzph(cox.zph.fit, font.caption=16,font.main = 16,font.submain=14, font.x = 16, font.y = 16,font.tickslab=14)
# # dev.off()
# # ph_res <- as.data.frame(cox.zph.fit$table)
# ph_genes <- subset(rownames(ph_res),ph_res$p > 0.05)
# ph_genes
# if('GLOBAL' %in% ph_genes){
#   ph_genes <- ph_genes[-length(ph_genes)]
# }else{
#   ph_genes <- ph_genes
# }
# ph_genes2 <- c(c("futime","fustat"),as.vector(ph_genes))
# if(length(ph_genes2) >= 4){
#   
#   #去除p<0.05的结果
#   uniSigExp = uniSigExp[,ph_genes2]
# }
#ph 检验 new---------
###-- 存Schoenfeld图 （ph假定检验的结果）------------
# setwd(path)
# if (! dir.exists("./02_Schoenfeld")){
#   dir.create("./02_Schoenfeld")
# }
# setwd("./02_Schoenfeld")
#PH
####
library(survminer)
ph_res <- data.frame()

for(i in 3:ncol(uniSigExp)){
  fml <- as.formula(paste0('Surv(futime,fustat)~',colnames(uniSigExp)[i]))
  coxph_model <- coxph(fml, data=uniSigExp[,c(1:2,i)],id = rownames(uniSigExp))
  cox.zph.fit <- cox.zph(coxph_model)$table
  
  
  cox_zph <- cox.zph(coxph_model)
  
  plot_cox <- ggcoxzph(cox_zph,var = 1)
  
  pdf_filename <- paste0(i-2,".Schoenfeld_", colnames(uniSigExp)[i], ".pdf")
  png_filename <- paste0(i-2,".Schoenfeld_", colnames(uniSigExp)[i], ".png")
  
  # 打开 PDF 设备
  pdf(file = pdf_filename, height = 4, width = 8, onefile = FALSE)
  a <- dev.cur()
  
  # 打开 PNG 设备
  png(filename = png_filename, height = 4, width = 8, units = 'in', res = 600, family = 'Times')
  dev.control("enable")
  
  # 绘制图表
  print(plot_cox)
  
  # 将图表复制到 PDF 设备
  dev.copy(which = a)
  
  # 关闭设备
  dev.off()
  dev.off()
  
  
  
  ph_res <- rbind(ph_res,cox.zph.fit[1,])
}
colnames(ph_res) <- c('chisq','df','p')
rownames(ph_res) <- colnames(uniSigExp)[-(1:2)]
ph_res
# fml <- as.formula(paste0('Surv(futime,fustat)~',paste0(colnames(uniSigExp)[-c(1:2)],collapse = '+')))
# f <- coxph(fml, data=uniSigExp,id = rownames(uniSigExp))
# for(j in 3:40){
#     gene <- colnames(uniSigExp)[j]
#     cox.zph(coxph(Surv(futime,fustat)~ NR2F1, data=uniSigExp[,c(1,2,39)],id = rownames(uniSigExp)))##PH假定
# }
#   cox.zph.fit <- cox.zph(f)
# ph_res <- as.data.frame(cox.zph.fit$table)
ph_genes <- subset(rownames(ph_res),ph_res$p > 0.05)
if('GLOBAL' %in% ph_genes){
  ph_genes <- ph_genes[-length(ph_genes)]
}else{
  ph_genes <- ph_genes
}
ph_genes2 <- c(c("futime","fustat"),as.vector(ph_genes))
# if(length(ph_genes2) >= 4){

#去除p<0.05的结果
uniSigExp = uniSigExp[,ph_genes2]
head(uniSigExp)
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
head(uniSigExp)
write.table(uniSigExp,file="uniSigExp.txt",sep="\t",row.names=F,quote=F)
head(sigTab)
sigTab_filter <- sigTab[sigTab$id%in%ph_genes2,]
head(sigTab_filter)
write.table(sigTab_filter,file="uniCoxResult.Sig.txt",sep="\t",row.names=F,quote=F)


setwd("/data/nas1/chenpeiru/44_GYZK-30212-7/03_cox")
univariateCox_Result=read.table("uniCoxResult.Sig.txt",header=T,sep="\t",row.names=1,check.names=F)
head(univariateCox_Result)
p.value<-signif(univariateCox_Result[,5], digits=3)
#wald.test<-signif(x$wald["test"], digits=2)
Coefficient<-signif(univariateCox_Result[,1], digits=3);#coeficient beta
HR <-signif(univariateCox_Result[,2], digits=3);#exp(beta)
HR.confint.lower <- signif(univariateCox_Result[,3], digits=3)
HR.confint.upper <- signif(univariateCox_Result[,4], digits=3)
HR.combine <- paste0(HR, " (",
                     HR.confint.lower, "-", HR.confint.upper, ")")
rescox.temp.1 <-cbind(HR, HR.confint.lower,HR.confint.upper,HR.combine,p.value)
names(rescox.temp.1)<-c("HR", "HR.confint.lower", "HR.confint.upper",'HR.combine',
                        "p.value")
rownames(rescox.temp.1) <- rownames(univariateCox_Result)

univOut_sig.plot.univar <- rescox.temp.1
gene <- rownames(univOut_sig.plot.univar)
hr <- univOut_sig.plot.univar[,1]
hrLow <- univOut_sig.plot.univar[,2]
hrHigh <- univOut_sig.plot.univar[,3]
Hazard.ratio <- univOut_sig.plot.univar[,4]
pVal <- univOut_sig.plot.univar[,5]
pdf(file = "00_single_factorcox.pdf",height = 9,width = 9)
par(family="serif")
n <- nrow(univOut_sig.plot.univar)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
layout(matrix(c(1,1,1,2,2,2), 1, 6, byrow = TRUE))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=2
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.2+0.2,n:1,pVal,adj=1,cex=text.cex);text(1.2+0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
highlim <- max(as.numeric(hrHigh))+0.1
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(max(c(0,min(as.numeric(hrLow))))-0.1,highlim)
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio",cex.lab=2)
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1, lwd = 1.5,cex.axis=2,cex.lab=1.5)

dev.off()

#######################    多因素####
# library(survival)
# library(survminer)
# rt=read.table("uniSigExp.txt",header=T,sep="\t",check.names=F,row.names=1)
# multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
# multiCox=step(multiCox,direction = "both")
# multiCoxSum=summary(multiCox)
# 
# outTab=data.frame()
# outTab=cbind(
#   coef=multiCoxSum$coefficients[,"coef"],
#   HR=multiCoxSum$conf.int[,"exp(coef)"],
#   HR.95L=multiCoxSum$conf.int[,"lower .95"],
#   HR.95H=multiCoxSum$conf.int[,"upper .95"],
#   pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
# outTab=cbind(id=row.names(outTab),outTab)
# outTab=gsub("`","",outTab)
# write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)
# 
# riskScore=predict(multiCox,type="risk",newdata=rt)
# coxGene=rownames(multiCoxSum$coefficients)
# coxGene=gsub("`","",coxGene)
# outCol=c("futime","fustat",coxGene)
# write.table(cbind(id=rownames(cbind(rt[,outCol],riskScore)),cbind(rt[,outCol],riskScore)),
#             file="risk.train.txt",sep="\t",quote=F,row.names=F)
# 
# 
# p.value<-signif(multiCoxSum$coefficients[,5], digits=2)
# #wald.test<-signif(x$wald["test"], digits=2)
# Coefficient<-signif(multiCoxSum$coefficients[,1], digits=2);#coeficient beta
# HR <-signif(multiCoxSum$coefficients[,2], digits=2);#exp(beta)
# HR.confint.lower <- signif(multiCoxSum$conf.int[,"lower .95"], 2)
# HR.confint.upper <- signif(multiCoxSum$conf.int[,"upper .95"],2)
# z<-signif(multiCoxSum$coefficients[,4],2)
# HR.combine <- paste0(HR, " (",
#                      HR.confint.lower, "-", HR.confint.upper, ")")
# rescox.temp<-cbind(HR, HR.confint.lower,HR.confint.upper,HR.combine,p.value)
# names(rescox.temp)<-c("HR", "HR.confint.lower", "HR.confint.upper",'HR.combine',
#                       "p.value")
# rownames(rescox.temp) <- rownames(multiCoxSum$coefficients)
# 
# univOut_sig.plot <- rescox.temp
# gene <- rownames(univOut_sig.plot)
# hr <- univOut_sig.plot[,1]
# hrLow <- univOut_sig.plot[,2]
# hrHigh <- univOut_sig.plot[,3]
# Hazard.ratio <- univOut_sig.plot[,4]
# pVal <- univOut_sig.plot[,5]
# 
# n <- nrow(univOut_sig.plot)
# nRow <- n+1
# ylim <- c(1,nRow)
# layout(matrix(c(1,2),nc=2),width=c(3,2.5))
# xlim = c(0,3)
# par(mar=c(4,2.5,2,1))
# layout(matrix(c(1,1,1,2,2,2), 1, 6, byrow = TRUE))
# plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
# text.cex=1.3
# text(0,n:1,gene,adj=0,cex=text.cex)
# text(1.5+0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5+0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
# text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
# highlim <- max(as.numeric(hrHigh))+0.1
# par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
# xlim = c(max(c(0,min(as.numeric(hrLow))))-0.1,highlim)
# plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio",cex.lab=2)
# arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
# abline(v=1,col="black",lty=2,lwd=2)
# boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
# points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
# axis(1, lwd = 1.5,cex.axis=text.cex,cex.lab=1.5)
# 
# 
# ##########    ROC曲线####
# library(survivalROC)
# heatmap_train = read.table("risk.train.txt",header=T,sep="\t",
#                            check.names=F,row.names=1)
# #pdf(file="train.ROC.pdf",height = 8,width = 8)
# roc=survivalROC(Stime=heatmap_train$futime, status=heatmap_train$fustat,
#                 marker = heatmap_train$riskScore,predict.time =1*365, method="KM")
# plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='#FA8072',
#      xlab="False positive rate", ylab="True positive rate",
#      #main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
#      main="ROC curve",
#      lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
# abline(0,1)
# aucText=c()
# rocCol <- c('#FA8072','#63B8FF','#FFC1C1','#ADFF2F','#FFFF00')
# aucText=c(aucText,paste0("1 years"," (AUC=",sprintf("%.3f",roc$AUC),")"))
# j =0
# for (i in c(1,2,3,4,5)){
#   roc1=survivalROC(Stime=heatmap_train$futime, status=heatmap_train$fustat,
#                    marker = heatmap_train$riskScore,predict.time =i*365, method="KM")
#   j=j+1
#   aucText=c(aucText,paste0(i," years"," (AUC=",sprintf("%.3f",roc1$AUC),")"))
#   lines(roc1$FP, roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j+1],lwd = 3)
# }
# legend("bottomright", aucText,lwd=2,bty="n",col=rocCol,cex = 1)
# abline(0,1)
# #dev.off()
# 
# 
# ###################   KM曲线####
# library(survival)
# library("survminer")
# library(tidyr)
# risk.train = read.table("risk.train.txt",header=T,sep="\t",row.names=1)
# risk.train$futime=risk.train$futime/365
# # group <- as.vector(ifelse(risk.train$riskScore>bestthreshold.train,"high","low"))
# bestthreshold.surv <- surv_cutpoint(risk.train,time = "futime",
#                                     event = "fustat",
#                                     'riskScore',
#                                     minprop = 0.35,
#                                     progressbar = TRUE)
# bestthreshold.train <- bestthreshold.surv$cutpoint$cutpoint
# group <- ifelse(risk.train$riskScore > bestthreshold.train, 'high','low')
# 
# diff_train=survdiff(Surv(futime, fustat) ~group,data = risk.train)
# pValue=1-pchisq(diff_train$chisq,df=1)
# pValue=round(pValue,5)
# fit_train <- survfit(Surv(futime, fustat) ~ group, data = risk.train)
# summary(fit_train)
# risk.train.group <- cbind(risk.train,group)
# risk.train.group <- data.frame(id=row.names(risk.train.group),risk.train.group)
# write.table(risk.train.group,file="risk.train.group.txt",
#             sep="\t",quote=F,row.names=F)
# 
# ####        KM曲线
# #pdf(file="survival.coefficient.train.pdf",height = 9,width = 9)
# ggsurvplot(fit_train,
#            conf.int=TRUE,
#            pval=TRUE,
#            risk.table=TRUE,
#            xlab="time (years)",
#            #legend.labs=legend.labs,
#            legend.title="Risk",
#            #linetype = lty,
#            palette = c( "red", "blue"),
#            title="Kaplan-Meier Curve for Survival",
#            risk.table.height=.15) %>% print()
# #dev.off()
# 
# 
# ###############  风险曲线####
# loc_train = match(c('futime','fustat','riskScore'),colnames(heatmap_train))
# datalast_train = heatmap_train[,loc_train]
# datalast_train = data.frame(datalast_train)
# datalast_train$futime=datalast_train$futime/365
# #datalast_train$riskScore[which(datalast_train$riskScore >50)] <- 50
# 
# #pdf(file="生存曲线.pdf",height = 12,width = 12)
# par(mfrow=c(2,1))
# par(mar=c(2,4.5,1,1))
# col=c()
# col[sort(datalast_train$riskScore) <= bestthreshold.train]="blue"
# col[sort(datalast_train$riskScore) > bestthreshold.train]="red"
# plot(sort(datalast_train$riskScore),axes=F,xlab = NA,ylab = "Risk Score",col=col,
#      mgp=c(2.5,1,0),cex.lab=2)
# box(lwd=2)
# abline(v = length(which(sort(datalast_train$riskScore) <=  bestthreshold.train))+0.5,lty="dashed")
# abline(h =  bestthreshold.train,lty="dashed")
# text(50, bestthreshold.train+0.5,paste0('threshold = ',round( bestthreshold.train,4)),cex = 1.5,font=2)
# axis(2,seq(0,50,1),cex.axis=2)
# axis(1,seq(0,nrow(datalast_train),50),cex.axis=2)
# legend("topleft", c("High risk","Low risk"), pch=16:16, col=c("red","blue"),cex = 1.5)
# #图2
# datalastSORT = datalast_train[order(datalast_train[,"riskScore"]),]
# col=c()
# col[datalastSORT[,"fustat"]==1]= "red"
# col[datalastSORT[,"fustat"]==0]= "blue"
# par(mar=c(2,4.5,1,1))
# plot(datalastSORT[,"futime"],col=col,pch=16,axes=F,xlab = NA,mgp=c(2.5,1,0),cex.lab=1.8,
#      ylab = "Following up (years)")
# box(lwd=2)
# abline(v = length(which(sort(datalast_train$riskScore) <= bestthreshold.train))+0.5,lty="dashed")
# axis(2,seq(0,max(datalastSORT[,"futime"]),2),cex.axis=2)
# axis(1,seq(0,nrow(datalast_train),10),cex.axis=2)
# #dev.off()
# 



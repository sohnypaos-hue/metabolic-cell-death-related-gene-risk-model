rm(list = ls())
library(data.table)
setwd("/data/nas1/chenpeiru/44_GYZK-30212-7/04_lasso")
dat <- read.table("risk.train.txt",header=T,sep="\t",row.names=1,check.names=F)
head(dat)
data <- data.frame(
  New_Column_Name = rownames(dat),
  dat,
  row.names = NULL
)
head(data)
data <- data[,c(1:3,12)]
setwd("/data/nas1/chenpeiru/41_YQNC-10615-6/00_raw_data/01_TCGA")
phenotype_TCGA <- read.csv('phenotype_change_lable.csv',
                           header = T,sep = ',',row.names = 1)
head(phenotype_TCGA)
# phenotype_TCGA <- phenotype_TCGA[substr(rownames(phenotype_TCGA),13,23)=='BM_CD138pos',]
#phenotype_TCGA <- data.frame(id.small = rownames(phenotype_TCGA),
# futime = survival_TCGA$OS.time,
# fustat = survival_TCGA$OS)

phenotype_TCGA <- data.frame(
  New_Column_Name = rownames(phenotype_TCGA),
  phenotype_TCGA,
  row.names = NULL
)

#TCGA_phenotype <- read.csv("phenotype_change_lable.csv")
#head(TCGA_phenotype)
#colnames(TCGA_phenotype)[1] <- c("New_Column_Name")
head(phenotype_TCGA)
TCGA_phenotype_survival <- merge(phenotype_TCGA,data, id="New_Column_Name", all=FALSE)
head(TCGA_phenotype_survival)

rownames(TCGA_phenotype_survival) <- TCGA_phenotype_survival$New_Column_Name
TCGA_phenotype_survival <- TCGA_phenotype_survival[,-1]
head(TCGA_phenotype_survival)
TCGA_phenotype_survival$futime <- TCGA_phenotype_survival$futime/365
head(TCGA_phenotype_survival)

# TCGA_phenotype_survival <- na.omit(TCGA_phenotype_survival)

setwd("/data/nas1/chenpeiru/44_GYZK-30212-7/06_independent_prognosis")
save(TCGA_phenotype_survival,file = "TCGA_phenotype_survival.rda")#

####independent_prognosis###
rm(list = ls())
library(survival)
#setwd("/data/nas1/chenpeiru/41_YQNC-10615-6/08_independent_prognosis")
pFilter=0.05            #定义单因素显著性
################   单因素cox回归分析 ####
outTab=data.frame()
load("TCGA_phenotype_survival.rda")
dat <- read.table("/data/nas1/chenpeiru/44_GYZK-30212-7/04_lasso/risk.train.group.txt",
                  header=T,sep="\t",row.names=1,check.names=F)
head(dat)
data <- data.frame(
  New_Column_Name = rownames(dat),
  dat,
  row.names = NULL
)
head(data)
data <- data[,c(1,13)]
head(TCGA_phenotype_survival)
TCGA_phenotype_survival <- data.frame(
  New_Column_Name = rownames(TCGA_phenotype_survival),
  TCGA_phenotype_survival,
  row.names = NULL
)
head(TCGA_phenotype_survival)
head(data)
TCGA_phenotype_survival1 <- merge(TCGA_phenotype_survival,data, by.y="New_Column_Name",all=FALSE)
head(TCGA_phenotype_survival1)

# 
# #head(TCGA_phenotype_survival,20)
# TCGA_phenotype_survival1$Gender <- ifelse(TCGA_phenotype_survival1$Gender%in%"male",1,2)
# TCGA_phenotype_survival1$Stage <- factor(TCGA_phenotype_survival1$Stage,levels = c("I","II","III"),
#                                          labels = c(1,2,3))
# TCGA_phenotype_survival1$riskgroup <- factor(TCGA_phenotype_survival1$group,levels = c("high","low"),
#                                              labels = c(1,2))
# head(TCGA_phenotype_survival1)
# 
# TCGA_phenotype_survival1$Stage  <- as.numeric(TCGA_phenotype_survival1$Stage)
# TCGA_phenotype_survival1$riskgroup  <- as.numeric(TCGA_phenotype_survival1$riskgroup)
# rownames(TCGA_phenotype_survival1) <- TCGA_phenotype_survival1$New_Column_Name
# TCGA_phenotype_survival1 <- TCGA_phenotype_survival1[,c(-1,-7,-8)]
# head(TCGA_phenotype_survival1)
#class(TCGA_phenotype_survival$Stage)
head(TCGA_phenotype_survival1)
rownames(TCGA_phenotype_survival1) <- TCGA_phenotype_survival1$New_Column_Name
TCGA_phenotype_survival1 <- TCGA_phenotype_survival1[,-1]
head(TCGA_phenotype_survival1)
colnames(TCGA_phenotype_survival1)
univariate_rt <- TCGA_phenotype_survival1[,c(8,9,1,2,3,4,5,6,7,10)]
# head(univariate_rt)
# 
# colnames(univariate_rt)[10] <- "riskgroup"
# univariate_rt$riskgroup <- factor(univariate_rt$riskgroup,levels = c("low","high"),
#                                   labels = c(2,1))
head(univariate_rt)

outTab=data.frame()
for(i in colnames(univariate_rt[,3:ncol(univariate_rt)])){
  univar_cox_duli_ <- coxph(Surv(futime, fustat) ~ univariate_rt[,i], data = univariate_rt)
  #univar_cox_duli_=step(univar_cox_duli_,direction = "both")
  univar_cox_duli_ = summary(univar_cox_duli_)
  outTab=rbind(outTab,cbind(variable=i,
                            coef=univar_cox_duli_$coefficients[,"coef"],
                            HR=univar_cox_duli_$conf.int[,"exp(coef)"],
                            HR.95L=univar_cox_duli_$conf.int[,"lower .95"],
                            HR.95H=univar_cox_duli_$conf.int[,"upper .95"],
                            pvalue=univar_cox_duli_$coefficients[,"Pr(>|z|)"]))
}
outTab
# outTab1 <- outTab[-8,]
#setwd("/data/nas1/chenpeiru/41_YQNC-10615-6/08_independent_prognosis")
head(outTab)
# outTab1 <- outTab[-4,]#
write.table(outTab,file="univariateCox.Result.xls",sep="\t",row.names=F,quote=F)

###输出所有单因素的结果
outTab = outTab[is.na(outTab$pvalue)==FALSE,]
outTab=outTab[order(as.numeric(as.vector(outTab$pvalue))),]
head(outTab)
# outTab2 <- outTab[-4,]
write.table(outTab,file="uniCoxResult.txt",sep="\t",row.names=F,quote=F)
#输出单因素显著的结果
sigTab=outTab[as.numeric(as.vector(outTab$pvalue))<pFilter,]
head(sigTab)
# sigTab1 <- sigTab[-4,]
write.table(sigTab,file="uniCoxResult.Sig.txt",sep="\t",row.names=F,quote=F)
#输出单因素显著AS的PSI值，用于后续建模
sigclini=c("futime","fustat")
sigclini=c(sigclini,as.vector(sigTab[,1]))
uniSigExp=univariate_rt[,sigclini]
head(uniSigExp)
#PH 检验
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
head(uniSigExp)
#sigTab_filter <- sigTab[sigTab$id%in%ph_genes2,]
#head(sigTab_filter)
#write.table(sigTab_filter,file="uniCoxResult.Sig.txt",sep="\t",row.names=F,quote=F)
#head(sigTab)
#sigTab_filter <- sigTab[sigTab$id%in%ph_genes2,]
#head(sigTab_filter)
#write.table(sigTab_filter,file="2.input.txt",sep="\t",row.names=F,quote=F)

#uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
#write.table(uniSigExp,file="2.input.txt",sep="\t",row.names=F,quote=F)

###看一下单因素结果
# library(dplyr,warn.conflicts=F)
# library(kableExtra,warn.conflicts=F)
# sigTab %>% knitr::kable(format = "html",pad=0) %>%
#   kable_styling(bootstrap_options = "striped", full_width = F)

#####  单因素森林图
#setwd("/data/nas1/chenpeiru/41_YQNC-10615-6/08_independent_prognosis")
univariateCox_Result=read.table("uniCoxResult.Sig.txt",header=T,sep="\t",row.names=1,check.names=F)
head(univariateCox_Result)
pdf(file = "01_forestplot.pdf",height = 9,width = 12)
par(family="serif")
p.value<-signif(univariateCox_Result[,5], digits=3)
#wald.test<-signif(x$wald["test"], digits=2)
Coefficient<-signif(univariateCox_Result[,1], digits=3);#coeficient beta
HR <-signif(univariateCox_Result[,2], digits=4);#exp(beta)
HR.confint.lower <- signif(univariateCox_Result[,3], digits=3)
HR.confint.upper <- signif(univariateCox_Result[,4], digits=3)
HR.combine <- paste0(HR, " (",
                     HR.confint.lower, "-", HR.confint.upper, ")")
rescox.temp.1<-cbind(HR, HR.confint.lower,HR.confint.upper,HR.combine,p.value)
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
axis(1, lwd = 1.5,cex.axis=text.cex,cex.lab=1.5)
dev.off()
#只有"riskScore" "N_stage"   "Stage"     "T_stage"显著
############ 多因素cox回归分析 ####
multi_rt=read.table("uniSigExp.txt",header=T,sep="\t",row.names=1,check.names=F)
head(multi_rt)
cox <- coxph(Surv(futime, fustat) ~ ., data = multi_rt)
head(cox)
cox=step(cox,direction = "both")
cox
coxSum_multi=summary(cox)
coxSum_multi
outTab=data.frame()
outTab=cbind(
  coef=coxSum_multi$coefficients[,"coef"],
  HR=coxSum_multi$conf.int[,"exp(coef)"],
  HR.95L=coxSum_multi$conf.int[,"lower .95"],
  HR.95H=coxSum_multi$conf.int[,"upper .95"],
  pvalue=coxSum_multi$coefficients[,"Pr(>|z|)"])
outTab=as.data.frame(cbind(id=row.names(outTab),outTab))
outTab
#输出多因素显著的结果
outTab = outTab[as.numeric(as.vector(outTab$pvalue)) < pFilter,]
write.table(outTab,file="multiCox.Result.xls",sep="\t",row.names=F,quote=F)
head(outTab)
input3 <- multi_rt[,c("futime","fustat",outTab$id)]
head(input3)
write.table(input3,file="input3.txt",sep="\t",row.names=T,quote=F)

##看一下多因素结果
# library(dplyr,warn.conflicts=F)
# library(kableExtra,warn.conflicts=F)
# outTab %>% knitr::kable(format = "html",pad=0) %>%
#   kable_styling(bootstrap_options = "striped", full_width = F)


##多因素森林图
univariateCox_Result=read.table("multiCox.Result.xls",header=T,sep="\t",row.names=1,check.names=F)
head(univariateCox_Result)
pdf(file = "02_multi_forestplot.pdf",height = 8,width = 10)
par(family="serif")
p.value<-signif(univariateCox_Result[,5], digits=3)
#wald.test<-signif(x$wald["test"], digits=2)
Coefficient<-signif(univariateCox_Result[,1], digits=3);#coeficient beta
HR <-signif(univariateCox_Result[,2], digits=4);#exp(beta)
HR.confint.lower <- signif(univariateCox_Result[,3], digits=3)
HR.confint.upper <- signif(univariateCox_Result[,4], digits=3)
HR.combine <- paste0(HR, " (",
                     HR.confint.lower, "-", HR.confint.upper, ")")
rescox.temp.1<-cbind(HR, HR.confint.lower,HR.confint.upper,HR.combine,p.value)
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
axis(1, lwd = 1.5,cex.axis=text.cex,cex.lab=1.5)
dev.off()
# p.value<-signif(coxSum_multi$coefficients[,5], digits=3)
# #wald.test<-signif(x$wald["test"], digits=2)
# Coefficient<-signif(coxSum_multi$coefficients[,1], digits=3);#coeficient beta
# HR <-signif(coxSum_multi$coefficients[,2], digits=3);#exp(beta)
# HR.confint.lower <- signif(coxSum_multi$conf.int[,"lower .95"], 3)
# HR.confint.upper <- signif(coxSum_multi$conf.int[,"upper .95"],3)
# z<-signif(coxSum_multi$coefficients[,4],4)
# HR.combine <- paste0(HR, " (",
#                      HR.confint.lower, "-", HR.confint.upper, ")")
# rescox.temp<-cbind(HR, HR.confint.lower,HR.confint.upper,HR.combine,p.value)
# names(rescox.temp)<-c("HR", "HR.confint.lower", "HR.confint.upper",'HR.combine',
#                       "p.value")
# rownames(rescox.temp) <- rownames(coxSum_multi$coefficients)
# pdf(file = "02_multi_forestplot.pdf",height = 5,width = 9)
# par(family="serif")
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
# text.cex=2
# text(0,n:1,gene,adj=0,cex=text.cex)
# text(1.2+0.2,n:1,pVal,adj=1,cex=text.cex);text(1.2+0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
# text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
# highlim <- max(as.numeric(hrHigh))+0.1
# par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
# xlim = c(max(c(0,min(as.numeric(hrLow))))-1,highlim)
# plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio",cex.lab=2)
# arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
# abline(v=1,col="black",lty=2,lwd=2)
# boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
# points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
# axis(1, lwd = 1.5,cex.axis=text.cex,cex.lab=1.5)
# 
# dev.off()

########    列线图和校正曲线#### 

### 列线图 ####
# setwd("/data/nas1/chenpeiru/29_YQNJ-10132-5/07_independent_prognosis")
# pbc<-read.table("/data/nas1/chenpeiru/41_YQNC-10615-6/06_km_roc/risk.train.txt",
#                 header=TRUE,row.names=1)
# head(pbc)
# library(rms)
# dd<-datadist(pbc)
# # pbc$age <- pbc$age/365
# #pbc$futime <- pbc$futime*365
# options(datadist="dd")
# options(na.action="na.delete")
# summary(pbc$futime)
# colnames(pbc)
# coxpbc<-cph(formula = Surv(futime,fustat) ~  RAD9A + P4HB+PTGIS+G6PD+MYC,
#             data=pbc,x=T,y=T,surv = T,na.action=na.delete)
# print(coxpbc)
# surv<-Survival(coxpbc)
# #1,2,3
# surv1<-function(x) surv(1*365,x)
# surv2<-function(x) surv(2*365,x)
# surv3<-function(x) surv(3*365,x)

# x<-nomogram(coxpbc,fun = list(surv1,surv3,surv5),lp=T,
#             funlabel = c('1-year survival Probability','3-year survival Probability','5-year survival Probability'),
#             maxscale = 100,fun.at = c(0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1))
# pdf(file = "03_nomogram.pdf",height = 9,width = 12)
# plot(x, lplabel="Linear Predictor",
#      xfrac=.35,varname.label=TRUE, varname.label.sep="=", ia.space=.2,
#      tck=NA, tcl=-0.20, lmgp=0.3,
#      points.label='Points', total.points.label='Total Points',
#      total.sep.page=FALSE,
#      cap.labels=FALSE,cex.var = 1.6,cex.axis = 1.05,lwd=5,
#      label.every = 1,col.grid = gray(c(0.8, 0.95)))
# dev.off()

# nomo-new-----
rm(list = ls())
library(rms)
library(survival)
library(lattice)
library(Formula)
library(ggplot2)
library(Hmisc)
library(regplot)
# pbc <- rt
setwd("/data/nas1/chenpeiru/44_GYZK-30212-7/06_independent_prognosis")
pbc<-read.table("input3.txt",
                header=TRUE,row.names=1)
dd <- datadist(pbc)
head(pbc)

options(datadist="dd")
options(na.action="na.delete")
colnames(pbc)
# rt<-na.omit(rt)
# head(rt)
################列线图
# dd=datadist(rt)
# options(datadist="dd")
head(pbc)
pbccox <- coxph(Surv(futime,fustat) ~ riskScore+ N_stage,
                data = pbc) 

regplot(pbccox,plots = c('bean','boxes'),
        #对观测2的六个指标在列线图上进行计分展示
        observation=pbc[50,], #也可以不展示
        #预测3年和5年的死亡风险，此处单位是day
        failtime = c(1,2,3), 
        # prfail = TRUE, #cox回归中需要TRUE
        showP = T, #是否展示统计学差异
        droplines = T,#观测2示例计分是否画线
        # colors = mycol, #用前面自己定义的颜色
        rank="sd", #根据统计学差异的显著性进行变量的排序
        interval="confidence",
        subticks = T,
        dencol="#FF8C00", boxcol="yellow") #展示观测的可信区间

dev.copy2pdf(file = "01_Nomogram.pdf", width = 8,height = 6)


#####   校正曲线####
set.seed(38)
head(pbc)
pbc$futime <- pbc$futime*365
f1<-cph(formula = Surv(futime,fustat) ~  riskScore + N_stage,data=pbc,x=T,y=T,surv = T,
        na.action=na.delete, time.inc = 365)
#参数m=30表示每组30个样本进行重复计算
cal1<-calibrate(f1, cmethod="KM",method="boot",u=365,m=30,B=1000)

f2<-cph(formula = Surv(futime,fustat) ~ riskScore + N_stage,data=pbc,x=T,y=T,surv = T,
        na.action=na.delete, time.inc = 730)
cal2<-calibrate(f2, cmethod="KM",method="boot",u=730,m=30,B=1000)

f3<-cph(formula = Surv(futime,fustat) ~ riskScore + N_stage,data=pbc,x=T,y=T,surv = T,
        na.action=na.delete,time.inc = 1095)
cal3<-calibrate(f3, cmethod="KM", method="boot",u=1095,m=30,B=1000)



pdf(file = "04_calibration curve_123.pdf",height = 9,width = 12)
#lty = 0  不加误差线
plot(cal1,lwd = 2,lty = 1,errbar.col = c("#2166AC"),
     bty = "l", #只画左边和下边框
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#2166AC"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)
mtext("")

plot(cal2,lwd = 2,lty = 1,errbar.col = c("#B2182B"),
     xlim = c(0,1),ylim= c(0,1),col = c("#B2182B"),add = T)
lines(cal2[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#B2182B"), pch = 16)
mtext("")

plot(cal3,lwd = 2,lty = 1,errbar.col = c("yellow"),
     xlim = c(0,1),ylim= c(0,1),col = c("yellow"),add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("yellow"), pch = 16)
mtext("")

abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("topleft", #图例的位置
       legend = c("1-year","2-year","3-year"), #图例文字
       col =c("#2166AC","#B2182B","yellow"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.2,#图例字体大小
       bty = "n")#不显示图例边框

dev.off()


###预测生存率ROC####
############    ROC曲线####
setwd("/data/nas1/chenpeiru/44_GYZK-30212-7/06_independent_prognosis")

pbc<-read.table("input3.txt",
                header=TRUE,row.names=1)
head(pbc)
library(rms)
dd<-datadist(pbc)
# pbc$age <- pbc$age/365
head(pbc)
pbc$futime <- pbc$futime*365
options(datadist="dd")
options(na.action="na.delete")
ROC_data <- pbc
# ROC_data$futime <- ROC_data$futime/365
#模型的预测值
ROC_data$pred <- predict(pbccox,pbc)


pdf(file = "04_ROC123.pdf",height = 9,width = 9)
roc=survivalROC(Stime=ROC_data$futime, status=ROC_data$fustat,
                marker = ROC_data$pred,predict.time =1*365, method="KM")
head(roc)
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='#FA8072',
     xlab="False positive rate", ylab="True positive rate",
     #main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     main="ROC curve",
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
aucText=c()
rocCol <- c('#FA8072','#63B8FF','#FFC1C1','#ADFF2F','#FFFF00')
aucText=c(aucText,paste0("1 year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
j =0
for (i in c(2,3)){
  roc1=survivalROC(Stime=ROC_data$futime, status=ROC_data$fustat,
                   marker = ROC_data$pred,predict.time =i*365, method="KM")
  j=j+1
  aucText=c(aucText,paste0(i," year"," (AUC=",sprintf("%.3f",roc1$AUC),")"))
  lines(roc1$FP, roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j+1],lwd = 3)
}
legend("bottomright",aucText, lwd=2,bty="n",col=rocCol,cex = 1)
abline(0,1)
dev.off()


### 决策曲线#### 
library(survival)
library(ggplotify)
#library(magick)
library(prodlim)
#library(cmprsk)

library(rms)
library(regplot)
library(ggDCA)
#install.packages("ggDCA")
source("stdca.R")

pbc<-read.table("input3.txt",
                header=TRUE,row.names=1)
#pbc <- pbc[complete.cases(pbc),] #删掉缺失数据
head(pbc)
#pbc$futime <- pbc$futime/365
#DCA curve ####

riskScore  <- cph(Surv(futime, fustat) ~ riskScore,pbc)
N_stage <- cph(Surv(futime, fustat) ~ N_stage,pbc)

# SPC25 <- cph(Surv(futime, fustat) ~ SPC25 ,pbc)
# CXCL3 <- cph(Surv(futime, fustat) ~ CXCL3 ,pbc)
model <- cph(Surv(futime, fustat) ~ riskScore +N_stage,pbc)

dca_cph <- dca(riskScore, N_stage,model, 
               model.names = c("riskScore","N_stage","model"), times = c(1,2,3))
dca_p <- ggplot(dca_cph,linetype = F,lwd = 1)
dca_p
ggsave(filename = 'DCA.pdf',plot = dca_p,family="serif",height = 7,width = 12)


############################################################
## 16_risk_score_z_KM_ROC.R
## Risk score calculation with Z-score normalization,
## KM/ROC plotting, and Bootstrap evaluation
############################################################

rm(list=ls())
options(stringsAsFactors = FALSE)

############################################################
## 0. Packages
############################################################

library(dplyr)
library(survival)
library(survminer)
library(survivalROC)
library(readr)
library(tibble)
library(pec)

############################################################
## 1. Paths
############################################################

base_dir <- "/data/nas1/refinement9/yanxiuhang/GYZK-30122-11"
fig_dir <- file.path(base_dir, "16_risk_score_figures")
tab_dir <- file.path(base_dir, "16_risk_score_tables")
dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(tab_dir, recursive=TRUE, showWarnings=FALSE)

############################################################
## 2. Load RDA data
############################################################

# TCGA
tcga_file <- "/data/nas1/chenpeiru/44_GYZK-30212-7/00_raw_data/01_TCGA/dat_fpkm_tumor406_survival.rda"
load(tcga_file) # loads fpkm_tumor_p_t_survival
tcga_dat <- fpkm_tumor_p_t_survival
tcga_dat$futime <- as.numeric(tcga_dat$OS.time)  # 生存时间
tcga_dat$fustat <- as.numeric(tcga_dat$OS) 

# GEO
geo_file <- "/data/nas1/chenpeiru/44_GYZK-30212-7/00_raw_data/02_geo/expr_OS_geo.rda"
geo_dat <- get(load(geo_file))
geo_dat$futime <- geo_dat$futime * 30
geo_dat <- subset(geo_dat, futime > 30)

############################################################
## 3. Load Our MCD model coefficients
############################################################

# 1. Load LASSO gene IDs used in Our_MCD
lasso_geneids <- readRDS(file.path("/data/nas1/chenpeiru/44_GYZK-30212-7",
    "04_lasso/lasso_geneids.rds"))

# 2. Load LASSO coefficients (coef.min)
coef_min_raw <- readRDS(file.path("/data/nas1/chenpeiru/44_GYZK-30212-7",
    "04_lasso/coef.min.rds"))

# 3. Extract coefficients for selected genes
coef_min <- as.numeric(coef_min_raw[lasso_geneids, ])
names(coef_min) <- lasso_geneids

# 4. Create Our_MCD model list
all_models <- list(Our_MCD = coef_min)

# Optional: save for later use
# saveRDS(all_models, file=file.path("/data/nas1/refinement9/yanxiuhang/GYZK-30122-11", "all_models.rds"))

############################################################
## 4. Z-score normalization
############################################################

zscore_normalize <- function(dat, genes){
  mat <- dat[, genes, drop=FALSE]
  mat <- apply(mat, 2, function(x) as.numeric(scale(x, center=TRUE, scale=TRUE)))
  as.data.frame(mat)
}

############################################################
## 5. Risk score calculation
############################################################

calc_risk_zscore <- function(dat, coef_vec){
  genes_present <- intersect(names(coef_vec), colnames(dat))
  if(length(genes_present)==0) return(rep(NA, nrow(dat)))
  mat_z <- zscore_normalize(dat, genes_present)
  as.numeric(as.matrix(mat_z) %*% coef_vec[genes_present])
}

############################################################
## 6. Compute risk scores
############################################################

for(model in names(all_models)){
  coef_vec <- all_models[[model]]
  tcga_dat[[paste0("risk_",model)]] <- calc_risk_zscore(tcga_dat, coef_vec)
  geo_dat[[paste0("risk_",model)]] <- calc_risk_zscore(geo_dat, coef_vec)
}

############################################################
## 7. Determine cutoff (median)
############################################################

tcga_cutoffs <- sapply(names(all_models), function(m) median(tcga_dat[[paste0("risk_",m)]], na.rm=TRUE))
geo_cutoffs <- sapply(names(all_models), function(m) median(geo_dat[[paste0("risk_",m)]], na.rm=TRUE))

############################################################
## 8. Save risk scores
############################################################

write.csv(tcga_dat, file=file.path(tab_dir,"TCGA_risk_scores_z.csv"), row.names=TRUE)
write.csv(geo_dat, file=file.path(tab_dir,"GEO_risk_scores_z.csv"), row.names=TRUE)
write.csv(tcga_cutoffs, file=file.path(tab_dir,"TCGA_cutoffs_median.csv"))
write.csv(geo_cutoffs, file=file.path(tab_dir,"GEO_cutoffs_median.csv"))

############################################################
## 9. KM plot
############################################################

plot_km <- function(dat, score_col, cohort, model, cutoff, out_dir){
  dat$group <- ifelse(dat[[score_col]]>cutoff, "High", "Low")
  dat$group <- factor(dat$group, levels=c("Low","High"))
  
  fit <- survfit(Surv(futime/365,fustat) ~ group, data=dat)
  
  p <- ggsurvplot(fit, data=dat, pval=TRUE, conf.int=TRUE,
                  risk.table=TRUE, risk.table.height=0.25,
                  legend.title="Group", legend.labs=c("Low risk","High risk"),
                  palette=c("blue","red"),
                  xlab="Time (years)", ylab="Overall survival probability",
                  title=paste0(cohort," - ",model),
                  ggtheme=theme_classic(), tables.theme=theme_classic())
  
  pdf(file.path(out_dir,paste0("KM_",cohort,"_",model,".pdf")), width=6, height=6)
  print(p); dev.off()
  png(file.path(out_dir,paste0("KM_",cohort,"_",model,".png")), width=2400, height=2400, res=300)
  print(p); dev.off()
}

############################################################
## 10. ROC plot
############################################################

plot_roc <- function(dat, score_col, cohort, year, model, out_dir){
  tmp <- dat[, c("futime","fustat",score_col)]
  tmp <- tmp[complete.cases(tmp), ]
  if(nrow(tmp)<30) return(NULL)
  
  roc_obj <- survivalROC(Stime=tmp$futime, status=tmp$fustat, marker=tmp[[score_col]],
                         predict.time=year*365, method="KM")
  
  pdf(file.path(out_dir,paste0("ROC_",cohort,"_",model,"_",year,"y.pdf")), width=6, height=6)
  plot(roc_obj$FP, roc_obj$TP, type="l", lwd=2, col="red",
       xlab="False Positive Rate", ylab="True Positive Rate",
       main=paste0(cohort," ",model," ",year,"-year ROC"))
  abline(0,1,lty=2,col="grey"); dev.off()
  
  png(file.path(out_dir,paste0("ROC_",cohort,"_",model,"_",year,"y.png")), width=2200, height=2200, res=300)
  plot(roc_obj$FP, roc_obj$TP, type="l", lwd=2, col="red",
       xlab="False Positive Rate", ylab="True Positive Rate",
       main=paste0(cohort," ",model," ",year,"-year ROC"))
  abline(0,1,lty=2,col="grey"); dev.off()
}

############################################################
## 11. Generate KM and ROC
############################################################

years <- c(1,3,5)
for(model in names(all_models)){
  score_col <- paste0("risk_",model)
  
  # TCGA
  plot_km(tcga_dat, score_col, "TCGA", model, tcga_cutoffs[model], fig_dir)
  for(y in years) plot_roc(tcga_dat, score_col, "TCGA", y, model, fig_dir)
  
  # GEO
  plot_km(geo_dat, score_col, "GSE13507", model, geo_cutoffs[model], fig_dir)
  for(y in years) plot_roc(geo_dat, score_col, "GSE13507", y, model, fig_dir)
}

############################################################
## 12. Bootstrap AUC, C-index, IBS
############################################################

nboot <- 100  # 100次Bootstrap

# ---- AUC ----
calc_auc_boot <- function(dat, score_col, cohort, model, year, nboot=100){
  tmp <- dat[, c("futime","fustat",score_col)]; tmp <- tmp[complete.cases(tmp),]
  if(nrow(tmp)<30) return(NULL)
  aucs <- numeric(nboot); set.seed(123)
  for(i in 1:nboot){
    idx <- sample(1:nrow(tmp), replace=TRUE)
    boot_tmp <- tmp[idx,]
    roc_obj <- survivalROC(Stime=boot_tmp$futime, status=boot_tmp$fustat, marker=boot_tmp[,score_col],
                           predict.time=year*365, method="KM")
    aucs[i] <- roc_obj$AUC
  }
  data.frame(Cohort=cohort, Model=model, Year=year,
             AUC_mean=round(mean(aucs),4),
             AUC_lower=round(quantile(aucs,0.025),4),
             AUC_upper=round(quantile(aucs,0.975),4))
}

auc_res_boot <- data.frame()
for(model in names(all_models)){
  score_col <- paste0("risk_",model)
  for(y in years){
    x1 <- calc_auc_boot(tcga_dat, score_col, "TCGA", model, y, nboot)
    x2 <- calc_auc_boot(geo_dat, score_col, "GSE13507", model, y, nboot)
    auc_res_boot <- rbind(auc_res_boot, x1, x2)
  }
}
write.csv(auc_res_boot, file=file.path(tab_dir,"AUC_summary_boot.csv"), row.names=FALSE)

# ---- C-index ----
calc_cindex_boot <- function(dat, score_col, cohort, model, nboot=100){
  tmp <- dat[, c("futime","fustat",score_col)]; colnames(tmp) <- c("futime","fustat","score")
  tmp <- tmp[complete.cases(tmp),]
  if(nrow(tmp)<30 || length(unique(tmp$score))<2) return(NULL)
  cindexes <- numeric(nboot); set.seed(123)
  for(i in 1:nboot){
    idx <- sample(1:nrow(tmp), replace=TRUE)
    boot_tmp <- tmp[idx,]
    fit <- coxph(Surv(futime,fustat) ~ score, data=boot_tmp)
    cindexes[i] <- summary(fit)$concordance[1]
  }
  data.frame(Cohort=cohort, Model=model,
             Cindex_mean=round(mean(cindexes),4),
             Cindex_lower=round(quantile(cindexes,0.025),4),
             Cindex_upper=round(quantile(cindexes,0.975),4))
}

cindex_res_boot <- data.frame()
for(model in names(all_models)){
  score_col <- paste0("risk_",model)
  x1 <- calc_cindex_boot(tcga_dat, score_col, "TCGA", model, nboot)
  x2 <- calc_cindex_boot(geo_dat, score_col, "GSE13507", model, nboot)
  cindex_res_boot <- rbind(cindex_res_boot, x1, x2)
}
write.csv(cindex_res_boot, file=file.path(tab_dir,"Cindex_summary_boot.csv"), row.names=FALSE)

# ---- IBS ----
calc_ibs_boot <- function(dat, score_col, cohort, model, nboot=100){
  tmp <- dat[, c("futime","fustat",score_col)]; colnames(tmp) <- c("futime","fustat","score")
  tmp <- tmp[complete.cases(tmp),]
  if(nrow(tmp)<30 || length(unique(tmp$score))<2) return(NULL)
  ibs_vals <- numeric(nboot); set.seed(123)
  for(i in 1:nboot){
    idx <- sample(1:nrow(tmp), replace=TRUE)
    boot_tmp <- tmp[idx,]
    fit <- coxph(Surv(futime,fustat) ~ score, data=boot_tmp, x=TRUE, y=TRUE)
    pec_fit <- pec(fit, formula=Surv(futime,fustat) ~ 1, data=boot_tmp,
                   cens.model="cox", exact=FALSE, splitMethod="none")
    ibs_vals[i] <- crps(pec_fit)[2]
  }
  data.frame(Cohort=cohort, Model=model,
             IBS_mean=round(mean(ibs_vals),4),
             IBS_lower=round(quantile(ibs_vals,0.025),4),
             IBS_upper=round(quantile(ibs_vals,0.975),4))
}

ibs_res_boot <- data.frame()
for(model in names(all_models)){
  score_col <- paste0("risk_",model)
  x1 <- calc_ibs_boot(tcga_dat, score_col, "TCGA", model, nboot=100)
  x2 <- calc_ibs_boot(geo_dat, score_col, "GSE13507", model, nboot=100)
  ibs_res_boot <- rbind(ibs_res_boot, x1, x2)
}
write.csv(ibs_res_boot, file=file.path(tab_dir,"IBS_summary_boot.csv"), row.names=FALSE)

cat("\nAll analyses finished successfully.\n")

############################################################
## 15. Multi-time ROC plot for Our_MCD
############################################################

plot_multi_timeROC <- function(dat, score_col, cohort, years=c(1,3,5), out_prefix){
  
  # PDF 输出
  pdf(paste0(out_prefix,".pdf"), width=7, height=7)
  plot(0,0,type="n", xlim=c(0,1), ylim=c(0,1),
       xlab="False Positive Rate", ylab="True Positive Rate",
       main=paste0(cohort, " Our_MCD ROC"))
  abline(0,1,lty=2,col="grey")
  
  cols <- c("red", "blue", "green")
  legend_txt <- c()
  
  for(i in seq_along(years)){
    year <- years[i]
    tmp <- dat[, c("futime","fustat", score_col)]
    tmp <- tmp[complete.cases(tmp), ]
    if(nrow(tmp) < 30) next
    
    roc_obj <- survivalROC(
      Stime = tmp$futime,
      status = tmp$fustat,
      marker = tmp[, score_col],
      predict.time = year*365,
      method = "KM"
    )
    
    lines(roc_obj$FP, roc_obj$TP, col=cols[i], lwd=2)
    legend_txt <- c(legend_txt, paste0(year,"-year AUC=", round(roc_obj$AUC,3)))
  }
  
  legend("bottomright", legend=legend_txt, col=cols[1:length(legend_txt)], lwd=2, cex=0.8)
  dev.off()
  
  # PNG 输出
  png(paste0(out_prefix,".png"), width=2200, height=2200, res=300)
  plot(0,0,type="n", xlim=c(0,1), ylim=c(0,1),
       xlab="False Positive Rate", ylab="True Positive Rate",
       main=paste0(cohort, " Our_MCD ROC"))
  abline(0,1,lty=2,col="grey")
  
  for(i in seq_along(years)){
    year <- years[i]
    tmp <- dat[, c("futime","fustat", score_col)]
    tmp <- tmp[complete.cases(tmp), ]
    if(nrow(tmp) < 30) next
    
    roc_obj <- survivalROC(
      Stime = tmp$futime,
      status = tmp$fustat,
      marker = tmp[, score_col],
      predict.time = year*365,
      method = "KM"
    )
    
    lines(roc_obj$FP, roc_obj$TP, col=cols[i], lwd=2)
  }
  
  legend("bottomright", legend=legend_txt, col=cols[1:length(legend_txt)], lwd=2, cex=0.8)
  dev.off()
}

############################################################
## 16. Call multi-time ROC for TCGA and GEO
############################################################

# TCGA
plot_multi_timeROC(
  dat = tcga_dat,
  score_col = "risk_Our_MCD",
  cohort = "TCGA",
  years = c(1,3,5),
  out_prefix = file.path(fig_dir, "ROC_TCGA_Our_MCD")
)

# GEO
plot_multi_timeROC(
  dat = geo_dat,
  score_col = "risk_Our_MCD",
  cohort = "GSE13507",
  years = c(1,3,5),
  out_prefix = file.path(fig_dir, "ROC_GSE13507_Our_MCD")
)

############################################################
## 17. Use training set (TCGA) cutoff for GEO validation
############################################################

# 1. 训练集中位数 cutoff
train_cutoff <- median(tcga_dat$risk_Our_MCD, na.rm=TRUE)

# 2. TCGA 分组（训练集自己用自己的中位数）
tcga_dat$group <- ifelse(tcga_dat$risk_Our_MCD > train_cutoff, "High", "Low")
tcga_dat$group <- factor(tcga_dat$group, levels=c("Low","High"))

# 3. GEO 验证集使用训练集 cutoff
geo_dat$group <- ifelse(geo_dat$risk_Our_MCD > train_cutoff, "High", "Low")
geo_dat$group <- factor(geo_dat$group, levels=c("Low","High"))

############################################################
## 18. KM plot using training cutoff
############################################################

# TCGA
fit_tcga <- survfit(Surv(futime/365, fustat) ~ group, data=tcga_dat)
p_tcga <- ggsurvplot(fit_tcga, data=tcga_dat, pval=TRUE, conf.int=TRUE,
                     risk.table=TRUE, risk.table.height=0.25,
                     legend.title="Group", legend.labs=c("Low risk","High risk"),
                     palette=c("blue","red"),
                     xlab="Time (years)", ylab="Overall survival probability",
                     title="TCGA - Our_MCD (training cutoff)",
                     ggtheme=theme_classic(), tables.theme=theme_classic())

pdf(file.path(fig_dir,"KM_TCGA_Our_MCD_training_cutoff.pdf"), width=6, height=6)
print(p_tcga)
dev.off()

png(file.path(fig_dir,"KM_TCGA_Our_MCD_training_cutoff.png"), width=2400, height=2400, res=300)
print(p_tcga)
dev.off()

# GEO
fit_geo <- survfit(Surv(futime/365, fustat) ~ group, data=geo_dat)
p_geo <- ggsurvplot(fit_geo, data=geo_dat, pval=TRUE, conf.int=TRUE,
                    risk.table=TRUE, risk.table.height=0.25,
                    legend.title="Group", legend.labs=c("Low risk","High risk"),
                    palette=c("blue","red"),
                    xlab="Time (years)", ylab="Overall survival probability",
                    title="GSE13507 - Our_MCD (training cutoff)",
                    ggtheme=theme_classic(), tables.theme=theme_classic())

pdf(file.path(fig_dir,"KM_GSE13507_Our_MCD_training_cutoff.pdf"), width=6, height=6)
print(p_geo)
dev.off()

png(file.path(fig_dir,"KM_GSE13507_Our_MCD_training_cutoff.png"), width=2400, height=2400, res=300)
print(p_geo)
dev.off()

############################################################
## 19. Multi-time ROC using same training cutoff
############################################################

# 这里 ROC 本身使用连续 risk score，不受 cutoff 影响
# 但你 KM 分组可以用 training cutoff
plot_multi_timeROC(dat=tcga_dat, score_col="risk_Our_MCD",
                   cohort="TCGA", years=c(1,3,5),
                   out_prefix=file.path(fig_dir,"ROC_TCGA_Our_MCD_training_cutoff"))

plot_multi_timeROC(dat=geo_dat, score_col="risk_Our_MCD",
                   cohort="GSE13507", years=c(1,3,5),
                   out_prefix=file.path(fig_dir,"ROC_GSE13507_Our_MCD_training_cutoff"))

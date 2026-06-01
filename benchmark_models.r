############################################################
## 15_benchmark_models.R
## Benchmark prognostic models in BLCA
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

############################################################
## 0. Packages
############################################################

library(survival)
library(survminer)
library(survivalROC)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(Matrix)

############################################################
## 1. Paths
############################################################

old_base <- "/data/nas1/chenpeiru/44_GYZK-30212-7"

new_base <- "/data/nas1/refinement9/yanxiuhang/GYZK-30122-11"

outdir <- file.path(new_base,"15_benchmark_models")

dir.create(outdir,recursive = TRUE,showWarnings = FALSE)

setwd(outdir)

figdir <- file.path(outdir, "figures")
tabdir <- file.path(outdir, "tables")
rdsdir <- file.path(outdir, "rds")

dir.create(figdir, showWarnings = FALSE)
dir.create(tabdir, showWarnings = FALSE)
dir.create(rdsdir, showWarnings = FALSE)

############################################################
## 2. Published models
############################################################

published_models <- list(
  Sun_Ferroptosis = c(
    TFRC = 0.00144, G6PD = 0.00255, SLC38A1 = -0.00415, ZEB1 = 0.04419,
    SCD = 0.001382, SRC = -0.00382, PRDX6 = 0.001288),
  Hao_Ferro_Cupro = c( SCD = 0.18362,  DDR2 = 0.16029,  MT1A = 0.12655),
  Xiao_Cuproptosis = c( ZBTB41 = 0.692453, PRMT6 = 0.617085,  DDX10 = 0.543223,
    RPL17 = 0.3934, FANCF = -0.29078, MARS2 = -0.39922, HMGN4 = -0.42016,
    MRFAP1L1 = -0.43038, RBM34 = -0.60586, RSBN1L = -0.69634),
  Pan_Disulfidptosis = c( COL5A1 = 0.14, DIRAS3 = 0.28,  NKG7 = -0.23, POLR3G = 0.30)
)

############################################################
## 3. Our MCD model
############################################################

lasso_geneids <- readRDS( file.path(old_base, "04_lasso/lasso_geneids.rds"))

coef_min_raw <- readRDS( file.path( old_base,  "04_lasso/coef.min.rds"))

coef_min <- as.numeric( coef_min_raw[lasso_geneids, ])

names(coef_min) <- lasso_geneids

our_model <- list( Our_MCD = coef_min)

all_models <- c( our_model, published_models)

############################################################
## 4. Save coefficients
############################################################

coef_table <- do.call( rbind,
  lapply(names(all_models), function(m){
    data.frame(  Model = m,  Gene = names(all_models[[m]]),
      Coef = as.numeric(all_models[[m]])
    )
  })
)

write.table(coef_table, file = file.path( tabdir, "model_coefficients.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE)

############################################################
## 5. Helper functions
############################################################

standardize_surv_expr <- function(dat){
  
  dat <- as.data.frame( dat, check.names = FALSE)
  
  colnames(dat)[ colnames(dat) %in% c( "OS.time","OS_time")] <- "futime"
  
  colnames(dat)[ colnames(dat) %in% c("OS.status","OS_status")] <- "fustat"
  
  dat$futime <- as.numeric(as.character(dat$futime))
  
  dat$fustat <- as.numeric(as.character(dat$fustat))
  
  dat <- dat[!is.na(dat$futime) & !is.na(dat$fustat),]
  
  dat <- dat[ dat$futime > 0,]
  
  return(dat)
}

############################################################
## Risk score
############################################################

calc_risk <- function(dat, coef_vec){
  
  genes_present <- intersect(names(coef_vec),colnames(dat))
  
  if(length(genes_present) == 0){
    
    return(rep(NA, nrow(dat)))
    
  }
  
  expr_mat <- dat[, genes_present]
  
  expr_mat <- as.data.frame(expr_mat)
  
  expr_mat[] <- lapply(expr_mat, function(x){
      as.numeric(
        as.character(x)
      )
      
    }
  )
  
  expr_mat <- as.matrix(expr_mat)
  
  coef_use <- coef_vec[genes_present]
  
  risk <- as.numeric(expr_mat %*% coef_use)
  
  return(risk)
}

############################################################
## KM plot
############################################################

plot_km_one_model <- function(dat, score_col, cohort, model, out_prefix){
  
  tmp <- dat[, c("futime", "fustat", score_col)]
  
  colnames(tmp) <- c("futime", "fustat", "score")
  
  tmp <- tmp[complete.cases(tmp),]
  
  if(nrow(tmp) < 20){
    
    cat(model, cohort,": too few samples\n")
    
    return(NULL)
  }
  
  if(length(unique(tmp$fustat)) < 2){
    
    cat(model,
        cohort,": only one survival status\n")
    
    return(NULL)
  }
  
  if(length(unique(tmp$score)) < 2){
    
    cat(model,
        cohort,": score has only one value\n")
    
    return(NULL)
  }
  
  tmp$futime_year <- tmp$futime / 365
  
  cutoff <- median(tmp$score, na.rm = TRUE)
  
  tmp$group <- ifelse(tmp$score > cutoff, "High", "Low")
  
  tmp$group <- factor(tmp$group,levels = c("Low","High"))
  
  if(length(unique(tmp$group)) < 2){
    
    cat(model,
        cohort, ": only one risk group\n")
    
    return(NULL)
  }
  
  fit <- survfit( Surv(futime_year,fustat) ~ group, data = tmp)
  
  p <- ggsurvplot(fit, data = tmp, pval = TRUE, conf.int = FALSE,
    risk.table = FALSE, legend.title = model, legend.labs = c("Low risk","High risk"),
    xlab = "Time (years)",
    ylab = "Overall survival probability",
    title = paste0(
      cohort,
      " - ",
      model
    ),
    palette = c(
      "#3C8DBC",
      "#E74C3C"
    )
  )
  
  ##########################################################
  ## PDF
  ##########################################################
  
  pdf(paste0(out_prefix, ".pdf"),width = 7,height = 6)
  
  print(p)
  
  dev.off()
  
  ##########################################################
  ## PNG
  ##########################################################
  
  png(paste0(out_prefix, ".png"), width = 2200,height = 1800,res = 300)
  
  print(p)
  
  dev.off()
  
  cat("Finished:",cohort, model,"\n")
}

############################################################
## ROC plot
############################################################

plot_multiROC <- function(dat,cohort,predict_time,out_prefix){
  
  roc_list <- list()
  
  for(model in names(all_models)){
    
    score_col <- paste0("risk_",model)
    
    tmp <- dat[, c("futime","fustat",score_col)]
    
    tmp <- tmp[complete.cases(tmp),]
    
    if(nrow(tmp) < 30) next
    
    roc_obj <- survivalROC(Stime = tmp$futime, status = tmp$fustat,
      marker = tmp[, score_col], predict.time = predict_time, method = "KM")
    
    roc_list[[model]] <- roc_obj
  }
  
  ##########################################################
  ## PDF
  ##########################################################
  
  pdf(paste0(out_prefix, ".pdf"), width = 7,height = 7)
  
  plot(
    c(0,1),
    c(0,1),
    type = "n",
    xlab = "False Positive Rate",
    ylab = "True Positive Rate",
    main = paste0(
      cohort,
      " ",
      predict_time/365,
      "-year ROC"
    )
  )
  
  abline(0,1,lty=2,col="grey")
  
  cols <- c(
    "#E74C3C",
    "#3498DB",
    "#2ECC71",
    "#9B59B6",
    "#F39C12"
  )
  
  i <- 1
  
  legend_txt <- c()
  
  for(model in names(roc_list)){
    
    roc_obj <- roc_list[[model]]
    
    lines(
      roc_obj$FP,
      roc_obj$TP,
      col = cols[i],
      lwd = 2
    )
    
    legend_txt <- c(
      legend_txt,
      paste0(
        model,
        " AUC=",
        round(roc_obj$AUC,3)
      )
    )
    
    i <- i + 1
  }
  
  legend(
    "bottomright",
    legend = legend_txt,
    col = cols[1:length(legend_txt)],
    lwd = 2,
    cex = 0.8
  )
  
  dev.off()
  
  ##########################################################
  ## PNG
  ##########################################################
  
  png(
    paste0(out_prefix, ".png"),
    width = 2200,
    height = 2200,
    res = 300
  )
  
  plot(
    c(0,1),
    c(0,1),
    type = "n",
    xlab = "False Positive Rate",
    ylab = "True Positive Rate",
    main = paste0(
      cohort,
      " ",
      predict_time/365,
      "-year ROC"
    )
  )
  
  abline(0,1,lty=2,col="grey")
  
  i <- 1
  
  for(model in names(roc_list)){
    
    roc_obj <- roc_list[[model]]
    
    lines(
      roc_obj$FP,
      roc_obj$TP,
      col = cols[i],
      lwd = 2
    )
    
    i <- i + 1
  }
  
  legend(
    "bottomright",
    legend = legend_txt,
    col = cols[1:length(legend_txt)],
    lwd = 2,
    cex = 0.8
  )
  
  dev.off()
}

############################################################
## 6. Load GEO
############################################################

geo_file <- file.path(
  old_base,
  "00_raw_data/02_geo/expr_OS_geo.rda"
)

geo_dat <- get(
  load(geo_file)
)

geo_dat$futime <- geo_dat$futime * 30

geo_dat <- subset(
  geo_dat,
  futime > 30
)

############################################################
## 7. Load TCGA
############################################################

tcga_file <- file.path(
  old_base,
  "00_raw_data/01_TCGA/dat_fpkm_tumor406_survival.rda"
)

load(tcga_file)

tcga_dat <- fpkm_tumor_p_t_survival

colnames(tcga_dat)[colnames(tcga_dat) == "OS.time"] <- "futime"
colnames(tcga_dat)[colnames(tcga_dat) == "OS"] <- "fustat"

tcga_dat$futime <- as.numeric(as.character(tcga_dat$futime))
tcga_dat$fustat <- as.numeric(as.character(tcga_dat$fustat))

tcga_dat <- tcga_dat[
  !is.na(tcga_dat$futime) &
    !is.na(tcga_dat$fustat) &
    tcga_dat$futime > 0,
]

dim(tcga_dat)

setdiff(
  names(published_models$Xiao_Cuproptosis),
  colnames(tcga_dat)
)

setdiff(
  names(published_models$Pan_Disulfidptosis),
  colnames(tcga_dat)
)
############################################################
## 8. Calculate risk scores
############################################################

for(model in names(all_models)){
  
  coef_vec <- all_models[[model]]
  
  tcga_dat[[paste0(
    "risk_",
    model
  )]] <- calc_risk(
    tcga_dat,
    coef_vec
  )
  
  geo_dat[[paste0(
    "risk_",
    model
  )]] <- calc_risk(
    geo_dat,
    coef_vec
  )
}

############################################################
## 9. Save risk scores
############################################################

risk_cols <- paste0(
  "risk_",
  names(all_models)
)

write.table(
  tcga_dat[, c(
    "futime",
    "fustat",
    risk_cols
  )],
  file = file.path(
    tabdir,
    "benchmark_risk_TCGA.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE
)

write.table(
  geo_dat[, c(
    "futime",
    "fustat",
    risk_cols
  )],
  file = file.path(
    tabdir,
    "benchmark_risk_GSE13507.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE
)

############################################################
## 10. KM plots
## No loop version (stable)
############################################################

plot_km_one_model <- function(
    dat,
    score_col,
    cohort,
    model,
    out_prefix
){
  
  tmp <- dat[, c(
    "futime",
    "fustat",
    score_col
  )]
  
  colnames(tmp) <- c(
    "futime",
    "fustat",
    "score"
  )
  
  tmp <- tmp[
    complete.cases(tmp),
  ]
  
  if(nrow(tmp) < 20){
    cat(model, cohort, ": too few samples\n")
    return(NULL)
  }
  
  if(length(unique(tmp$fustat)) < 2){
    cat(model, cohort, ": only one survival status\n")
    return(NULL)
  }
  
  if(length(unique(tmp$score)) < 2){
    cat(model, cohort, ": score has only one value\n")
    return(NULL)
  }
  
  ##########################################################
  ## Time
  ##########################################################
  
  tmp$futime_year <- tmp$futime / 365
  
  ##########################################################
  ## Group
  ##########################################################
  
  cutoff <- median(
    tmp$score,
    na.rm = TRUE
  )
  
  tmp$group <- ifelse(
    tmp$score > cutoff,
    "highrisk",
    "lowrisk"
  )
  
  tmp$group <- factor(
    tmp$group,
    levels = c(
      "highrisk",
      "lowrisk"
    )
  )
  
  if(length(unique(tmp$group)) < 2){
    cat(model, cohort, ": only one risk group\n")
    return(NULL)
  }
  
  ##########################################################
  ## Survival fit
  ##########################################################
  
  fit <- survfit(
    Surv(
      futime_year,
      fustat
    ) ~ group,
    data = tmp
  )
  
  ##########################################################
  ## Plot
  ##########################################################
  
  p <- ggsurvplot(
    fit,
    data = tmp,
    
    conf.int = TRUE,
    
    pval = TRUE,
    
    risk.table = TRUE,
    
    risk.table.col = "strata",
    
    risk.table.height = 0.25,
    
    risk.table.y.text = TRUE,
    
    risk.table.y.text.col = TRUE,
    
    legend.title = "group",
    
    legend.labs = c(
      "highrisk",
      "lowrisk"
    ),
    
    palette = c(
      "red",
      "blue"
    ),
    
    xlab = "Time(years)",
    
    ylab = "Survival probability",
    
    title = paste0(
      "Kaplan-Meier Curve for ",
      model
    ),
    
    ggtheme = theme_classic(),
    
    tables.theme = theme_classic()
  )
  
  ##########################################################
  ## Save PDF
  ##########################################################
  
  pdf(
    paste0(out_prefix, ".pdf"),
    width = 6,
    height = 6
  )
  
  print(p)
  
  dev.off()
  
  ##########################################################
  ## Save PNG
  ##########################################################
  
  png(
    paste0(out_prefix, ".png"),
    width = 2400,
    height = 2400,
    res = 300
  )
  
  print(p)
  
  dev.off()
  
  cat(
    "Finished:",
    cohort,
    model,
    "\n"
  )
}

############################################################
## TCGA
############################################################

plot_km_one_model(
  tcga_dat,
  "risk_Our_MCD",
  "TCGA",
  "Our_MCD",
  file.path(
    figdir,
    "KM_TCGA_Our_MCD"
  )
)

plot_km_one_model(
  tcga_dat,
  "risk_Sun_Ferroptosis",
  "TCGA",
  "Sun_Ferroptosis",
  file.path(
    figdir,
    "KM_TCGA_Sun_Ferroptosis"
  )
)

plot_km_one_model(
  tcga_dat,
  "risk_Hao_Ferro_Cupro",
  "TCGA",
  "Hao_Ferro_Cupro",
  file.path(
    figdir,
    "KM_TCGA_Hao_Ferro_Cupro"
  )
)

plot_km_one_model(
  tcga_dat,
  "risk_Xiao_Cuproptosis",
  "TCGA",
  "Xiao_Cuproptosis",
  file.path(
    figdir,
    "KM_TCGA_Xiao_Cuproptosis"
  )
)

plot_km_one_model(
  tcga_dat,
  "risk_Pan_Disulfidptosis",
  "TCGA",
  "Pan_Disulfidptosis",
  file.path(
    figdir,
    "KM_TCGA_Pan_Disulfidptosis"
  )
)

############################################################
## GEO
############################################################

plot_km_one_model(
  geo_dat,
  "risk_Our_MCD",
  "GSE13507",
  "Our_MCD",
  file.path(
    figdir,
    "KM_GSE13507_Our_MCD"
  )
)

plot_km_one_model(
  geo_dat,
  "risk_Sun_Ferroptosis",
  "GSE13507",
  "Sun_Ferroptosis",
  file.path(
    figdir,
    "KM_GSE13507_Sun_Ferroptosis"
  )
)

plot_km_one_model(
  geo_dat,
  "risk_Hao_Ferro_Cupro",
  "GSE13507",
  "Hao_Ferro_Cupro",
  file.path(
    figdir,
    "KM_GSE13507_Hao_Ferro_Cupro"
  )
)

plot_km_one_model(
  geo_dat,
  "risk_Xiao_Cuproptosis",
  "GSE13507",
  "Xiao_Cuproptosis",
  file.path(
    figdir,
    "KM_GSE13507_Xiao_Cuproptosis"
  )
)

plot_km_one_model(
  geo_dat,
  "risk_Pan_Disulfidptosis",
  "GSE13507",
  "Pan_Disulfidptosis",
  file.path(
    figdir,
    "KM_GSE13507_Pan_Disulfidptosis"
  )
)

############################################################
## 11. ROC plots
############################################################

plot_multiROC(
  tcga_dat,
  "TCGA",
  365,
  file.path(
    figdir,
    "ROC_compare_TCGA_1year"
  )
)

plot_multiROC(
  tcga_dat,
  "TCGA",
  365*3,
  file.path(
    figdir,
    "ROC_compare_TCGA_3year"
  )
)

plot_multiROC(
  tcga_dat,
  "TCGA",
  365*5,
  file.path(
    figdir,
    "ROC_compare_TCGA_5year"
  )
)

plot_multiROC(
  geo_dat,
  "GSE13507",
  365,
  file.path(
    figdir,
    "ROC_compare_GSE13507_1year"
  )
)

plot_multiROC(
  geo_dat,
  "GSE13507",
  365*3,
  file.path(
    figdir,
    "ROC_compare_GSE13507_3year"
  )
)

plot_multiROC(
  geo_dat,
  "GSE13507",
  365*5,
  file.path(
    figdir,
    "ROC_compare_GSE13507_5year"
  )
)

cat("\nAll analyses finished successfully.\n")

############################################################
## 12. AUC summary table with Bootstrap
############################################################

nboot <- 100  # Bootstrap次数
auc_res_boot <- data.frame()

calc_auc_boot <- function(dat, score_col, cohort, model, year, nboot=1000){
  tmp <- dat[, c("futime","fustat",score_col)]
  tmp <- tmp[complete.cases(tmp), ]
  if(nrow(tmp) < 30) return(NULL)
  
  aucs <- numeric(nboot)
  set.seed(123)
  for(i in 1:nboot){
    idx <- sample(1:nrow(tmp), replace=TRUE)
    boot_tmp <- tmp[idx, ]
    roc_obj <- survivalROC(
      Stime=boot_tmp$futime,
      status=boot_tmp$fustat,
      marker=boot_tmp[,score_col],
      predict.time=year*365,
      method="KM"
    )
    aucs[i] <- roc_obj$AUC
  }
  
  data.frame(
    Cohort=cohort,
    Model=model,
    Year=year,
    AUC_mean=round(mean(aucs),4),
    AUC_lower=round(quantile(aucs,0.025),4),
    AUC_upper=round(quantile(aucs,0.975),4)
  )
}

for(model in names(all_models)){
  score_col <- paste0("risk_",model)
  for(year in c(1,3,5)){
    x1 <- calc_auc_boot(tcga_dat, score_col, "TCGA", model, year, nboot)
    x2 <- calc_auc_boot(geo_dat, score_col, "GSE13507", model, year, nboot)
    auc_res_boot <- rbind(auc_res_boot, x1, x2)
  }
}

write.csv(auc_res_boot, file=file.path(tabdir,"AUC_summary_boot.csv"), row.names=FALSE)

############################################################
## 13. C-index summary with Bootstrap
############################################################

calc_cindex_boot <- function(dat, score_col, cohort, model, nboot=1000){
  tmp <- dat[, c("futime","fustat",score_col)]
  colnames(tmp) <- c("futime","fustat","score")
  tmp <- tmp[complete.cases(tmp), ]
  if(nrow(tmp) < 30 || length(unique(tmp$score)) < 2) return(NULL)
  
  cindexes <- numeric(nboot)
  set.seed(123)
  for(i in 1:nboot){
    idx <- sample(1:nrow(tmp), replace=TRUE)
    boot_tmp <- tmp[idx, ]
    fit <- coxph(Surv(futime,fustat) ~ score, data=boot_tmp)
    cindexes[i] <- summary(fit)$concordance[1]
  }
  
  data.frame(
    Cohort=cohort,
    Model=model,
    Cindex_mean=round(mean(cindexes),4),
    Cindex_lower=round(quantile(cindexes,0.025),4),
    Cindex_upper=round(quantile(cindexes,0.975),4)
  )
}

cindex_res_boot <- data.frame()
for(model in names(all_models)){
  score_col <- paste0("risk_",model)
  x1 <- calc_cindex_boot(tcga_dat, score_col, "TCGA", model, nboot)
  x2 <- calc_cindex_boot(geo_dat, score_col, "GSE13507", model, nboot)
  cindex_res_boot <- rbind(cindex_res_boot, x1, x2)
}

write.csv(cindex_res_boot, file=file.path(tabdir,"Cindex_summary_boot.csv"), row.names=FALSE)

############################################################
## 14. IBS summary with Bootstrap
############################################################

library(pec)

calc_ibs_boot <- function(dat, score_col, cohort, model, nboot=100){
  tmp <- dat[, c("futime","fustat",score_col)]
  colnames(tmp) <- c("futime","fustat","score")
  tmp <- tmp[complete.cases(tmp), ]
  if(nrow(tmp) < 30 || length(unique(tmp$score)) < 2) return(NULL)
  
  ibs_vals <- numeric(nboot)
  set.seed(123)
  for(i in 1:nboot){
    idx <- sample(1:nrow(tmp), replace=TRUE)
    boot_tmp <- tmp[idx, ]
    fit <- coxph(Surv(futime,fustat) ~ score, data=boot_tmp, x=TRUE, y=TRUE)
    pec_fit <- pec(fit, formula=Surv(futime,fustat) ~ 1, data=boot_tmp,
                   cens.model="cox", exact=FALSE, splitMethod="none")
    ibs_vals[i] <- crps(pec_fit)[2]
  }
  
  data.frame(
    Cohort=cohort,
    Model=model,
    IBS_mean=round(mean(ibs_vals),4),
    IBS_lower=round(quantile(ibs_vals,0.025),4),
    IBS_upper=round(quantile(ibs_vals,0.975),4)
  )
}

ibs_res_boot <- data.frame()
for(model in names(all_models)){
  score_col <- paste0("risk_",model)
  x1 <- calc_ibs_boot(tcga_dat, score_col, "TCGA", model, nboot=100)  # IBS bootstrap次数可低于AUC/Cindex
  x2 <- calc_ibs_boot(geo_dat, score_col, "GSE13507", model, nboot=100)
  ibs_res_boot <- rbind(ibs_res_boot, x1, x2)
}

write.csv(ibs_res_boot, file=file.path(tabdir,"IBS_summary_boot.csv"), row.names=FALSE)

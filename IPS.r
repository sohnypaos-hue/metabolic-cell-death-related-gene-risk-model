# -----------------------------
# 0. 清理环境
# -----------------------------
rm(list = ls())

# -----------------------------
# 1. 载入必要包
# -----------------------------
library(readxl)
library(ggplot2)
library(dplyr)
library(tibble)
library(reshape2)
library(rstatix)

# -----------------------------
# 2. 读取表达矩阵和生存数据
# -----------------------------
load("./dat_fpkm_tumor406_survival.rda") 
# 假设对象名为 fpkm_tumor_p_t_survival，行是样本，列是基因
expr_matrix <- fpkm_tumor_p_t_survival

# -----------------------------
# 3. 读取 IPS 文件
# -----------------------------
ips <- read.table("/data/nas1/refinement9/yanxiuhang/GYZK-30122-11/TCIA-ClinicalData.tsv",
                  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# -----------------------------
# 4. 提取 Risk Score 所需的 8 个基因
# -----------------------------
genes8 <- c("WWTR1","SLC7A11","G6PD","ZEB1","EGR1","FLNA","FADS2","SLC39A14")
expr_8genes <- expr_matrix[, genes8, drop = FALSE]

# -----------------------------
# 5. 计算 Risk Score
# -----------------------------
# 假设 coef_vec 为 LASSO 模型的系数
coef_vec <- c(0.691452,0.08333,0.084101,0.083040,0.092249,0.015159,0.094294,0.001880)
risk_score <- as.numeric(as.matrix(expr_8genes) %*% coef_vec)

# -----------------------------
# 6. 构建 Risk Score 数据框
# -----------------------------
risk_df <- data.frame(sample = rownames(expr_8genes),
                      RiskScore = risk_score)

# -----------------------------
# 7. 处理 IPS 数据
# -----------------------------
ips_df <- ips %>%
  rename(sample = barcode) %>%   # 样本 ID 列统一
  filter(expData == 1)           # 只保留有表达数据的样本

# -----------------------------
# 8. 合并 Risk Score 与 IPS 数据（匹配 TCGA 患者 ID）
# -----------------------------

# 创建短 ID（前 12 位）以匹配患者
risk_df <- risk_df %>%
  mutate(sample_short = substr(sample, 1, 12))

ips_df <- ips_df %>%
  rename(sample = sample) %>%
  filter(expData == 1) %>%
  mutate(sample_short = substr(sample, 1, 12))

# 使用短 ID 进行合并
combined <- inner_join(risk_df, ips_df, by = "sample_short")

# 可选：保留原始 sample 名，并删除辅助列
combined <- combined %>%
  select(-sample_short)

# -----------------------------
# 9. 按 Risk Score 中位数划分高低风险组
# -----------------------------
median_score <- median(combined$RiskScore, na.rm = TRUE)
combined <- combined %>%
  mutate(RiskGroup = ifelse(RiskScore > median_score, "High", "Low"))

# -----------------------------
# 10. 绘制小提琴图比较 IPS（使用 BH 校正 p 值）
# -----------------------------
ips_cols <- c("ips_ctla4_neg_pd1_neg","ips_ctla4_neg_pd1_pos",
              "ips_ctla4_pos_pd1_neg","ips_ctla4_pos_pd1_pos")

# 检查 IPS 列是否存在
missing_cols <- ips_cols[!ips_cols %in% colnames(combined)]
if(length(missing_cols) > 0){
  stop(paste("数据中缺少以下 IPS 列:", paste(missing_cols, collapse = ", ")))
}

# 数据整理为长表格
plot_data <- combined[, c("RiskGroup", ips_cols)]
plot_data_melt <- melt(plot_data, id.vars = "RiskGroup")

# -----------------------------
# 11. Wilcoxon 检验 + BH 校正
# -----------------------------
wilcox_results <- plot_data_melt %>%
  group_by(variable) %>%
  wilcox_test(value ~ RiskGroup) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

print(wilcox_results)

# 准备显著性标记数据框用于绘图
signif_df <- wilcox_results %>%
  dplyr::select(variable, p.adj.signif) %>%
  mutate(ypos = max(plot_data_melt$value, na.rm = TRUE) * 1.05)

# -----------------------------
# 12. 绘图 + 校正显著性标记
# -----------------------------
ggplot(plot_data_melt, aes(x = variable, y = value, fill = RiskGroup)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA,
               position = position_dodge(0.9)) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2,
               position = position_dodge(0.9)) +
  geom_text(data = signif_df, aes(x = variable, y = ypos, label = p.adj.signif),
            inherit.aes = FALSE) +
  scale_fill_manual(values = c("Low" = "skyblue", "High" = "tomato")) +
  theme_bw(base_size = 14) +
  labs(x = "IPS Subtype", y = "IPS Score", fill = "Risk Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

# -----------------------------
# 13. 保存图形
# -----------------------------
ggsave("/data/nas1/refinement9/yanxiuhang/GYZK-30122-11/violin_ips_risk_signif.png", width = 8, height = 6, dpi = 300)
ggsave("/data/nas1/refinement9/yanxiuhang/GYZK-30122-11/violin_ips_risk_signif.pdf", width = 8, height = 6)

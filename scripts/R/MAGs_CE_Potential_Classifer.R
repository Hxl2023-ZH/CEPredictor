################################################################################
#                           CE-MAGs Predictor Version 02
################################################################################
# Xinliang Hu
# 20260424
# R version 4.5.2

### load R packages
library(tidyverse)
library(data.table)
library(reshape2)

library(pheatmap)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(vegan)

library(randomForest)
library(xgboost)
library(caret)
library(FactoMineR)
library(factoextra)
library(PRROC)
library(pROC)
library(ComplexHeatmap)
library(circlize)
library(uwot)

############################## 1.数据预处理模块 ################################
setwd('')

raw_df <- fread("MAGs_all_merged_annotation.txt")

#过滤低质量MAG
#raw_df_medium <- raw_df %>% filter(Quality != "Low") %>% select(-c('H2ase','BM','ACOCT','HadABC','CarC','HypCDEF'))

dim(raw_df_medium)
table(raw_df_medium$Quality)

raw_data <- raw_df_medium %>% select(-"Quality")
MAGs <- raw_data$MAGs
completeness <- raw_data$completeness

gene_table <- raw_data %>%
  select(-MAGs,-completeness)

dir.create("results",showWarnings=FALSE)
dir.create("results/tables",showWarnings=FALSE)
dir.create("results/figures",showWarnings=FALSE)

# 二进制转换
gene_binary <- gene_table

gene_binary[gene_binary >0] <- 1
gene_binary[is.na(gene_binary)] <- 0

gene_binary <- as.data.frame(gene_binary)

rownames(gene_binary) <- MAGs

print(table(unlist(gene_binary)))

gene_prev <- colSums(gene_binary)/nrow(gene_binary)

prev_table <- data.frame(
  gene=names(gene_prev),
  prevalence=gene_prev
)

write.table(prev_table, "results/tables/gene_prevalence.txt",sep="\t",row.names=FALSE)

################################# 2.评分与筛选模块 ###############################
### 模块评分
module_list <- list(
  
  condensation=c("ACAT"),
  reduction1=c("HAD"),
  dehydration=c("ECH"),
  reduction2=c("BCD","EtfAB"),
  product_release=c("ACT","CoAT"),
  lactate_pathway=c("LDH","LacP"),
  ethanol_pathway=c("ADH","ADA"),
  xylose_pathway=c("XI","XK","XylT")
  
)

module_score <- function(df,genes){
  
  valid <- intersect(genes,colnames(df))
  
  if(length(valid)==0){
    return(rep(0,nrow(df)))
  }
  
  mat <- as.matrix(df[,valid,drop=FALSE])
  
  score <- rowSums(mat)/length(genes)
  
  return(score)
  
}

module_matrix <- data.frame(row.names=MAGs)

for(i in names(module_list)){
  
  module_matrix[,i] <- module_score(gene_binary,module_list[[i]])
  
}

module_matrix$rBOX <- rowMeans(
  module_matrix[,c(
    "condensation",
    "reduction1",
    "dehydration",
    "reduction2"
  )]
)

core_rule <- (
  module_matrix$rBOX >= 0.625
)

core_genes <- c("ACAT","HAD","ECH","BCD","EtfAB")

core_count <- rowSums(gene_binary[,core_genes,drop=FALSE])

candidate_MAG <- core_count >=3

rule_table <- data.frame(
  MAG=MAGs,
  core_rule=core_rule,
  core_count=core_count
)

write.table(rule_table, "results/tables/rule_filter.txt", sep="\t",row.names=FALSE)


################################# 3.机器学习模型构建模块 #######################

feature_matrix <- as.matrix(gene_binary[,!colnames(gene_binary) %in% core_genes])

#feature_matrix <- as.matrix(gene_binary)
labels <- as.numeric(core_rule)

# 数据划分
set.seed(2026)

train_index <- createDataPartition(
  labels,
  p=0.8,
  list=FALSE
)

train_x <- feature_matrix[train_index,]
train_y <- labels[train_index]

test_x <- feature_matrix[-train_index,]
test_y <- labels[-train_index]

# 类别不平衡权重计算
pos <- sum(train_y==1)
neg <- sum(train_y==0)

scale_pos_weight <- neg/pos

print(paste("Positive:",pos))
print(paste("Negative:",neg))
print(paste("scale_pos_weight:",scale_pos_weight))
#"scale_pos_weight: 1.86614173228346"

# XGBoost
dtrain <- xgb.DMatrix(train_x,label=train_y)

dtest <- xgb.DMatrix(
  test_x,
  label=test_y
)

param <- list(
  objective="binary:logistic",
  eval_metric=c("auc","aucpr"),
  max_depth=5,
  eta=0.1,
  subsample=0.8,
  colsample_bytree=0.8,
  scale_pos_weight=scale_pos_weight
)

watchlist <- list(
  train=dtrain,
  eval=dtest
)

xgb_model <- xgb.train(
  
  params=param,
  data=dtrain,
  nrounds=500,
  watchlist=watchlist,
  early_stopping_rounds=30,
  verbose=1
  
)

# RandomForest
rf_model <- randomForest(
  
  x=train_x,
  y=as.factor(train_y),
  ntree=1000,
  classwt=c(
    "0"=1,
    "1"=scale_pos_weight
  ),
  importance=TRUE
  
)

# Ensemble predict （平均法）
prob_xgb <- predict(xgb_model,feature_matrix)
prob_rf <- predict(rf_model,feature_matrix,type="prob")[,2]

ml_prob <- (prob_xgb + prob_rf)/2

### 6.CEPS评分
CEPS <- 0.5*ml_prob +
  0.3*module_matrix$product_release +
  0.2*0.01*completeness

class <- rep("Low",length(CEPS))

class[CEPS > 0.60] <- "Medium"
class[CEPS > 0.80] <- "High"

result_table <- data.frame(
  MAG=MAGs,
  CEPS=CEPS,
  class=class,
  core_rule=core_rule,
  ML_prob=ml_prob 
)

table(result_table$class)

write.table(result_table,
            "results/tables/CEPS_prediction.txt",
            sep="\t",
            quote = FALSE,
            row.names=FALSE)


############################## 4.外部数据集验证模块 #############################
###### 4.1 构建参考特征矩阵
# 加载参考基因组的label信息与丰度表
reference_genomes <- read.delim("ref_genome_label.txt", header = TRUE, check.names = FALSE)
reference_gene_table_raw <- fread("ref_genome_matrix_counts.txt")

reference_gene_table <- reference_gene_table_raw %>% 
  select(-MAGs) %>% 
  select(-c('H2ase','ACOCT','HadABC','CarC','HypCDEF'))

ref_binary <- reference_gene_table

ref_binary[ref_binary > 0] <- 1
ref_binary[is.na(ref_binary)] <- 0

ref_binary <- as.data.frame(ref_binary)
rownames(ref_binary) <- reference_genomes$MAG

### 合并数据集
# check col names
cols_gene_binary <- names(gene_binary)
missing_cols <- setdiff(cols_gene_binary, names(ref_binary))
ref_binary[missing_cols] <- 0
ref_binary <- ref_binary[, cols_gene_binary]

# 检查列名是否一致
colnames(gene_binary) == colnames(ref_binary)

if (!identical(colnames(gene_binary), colnames(ref_binary))) {
  stop("Column names not match exactly (order or content)！！！")
}

combined_matrix <- rbind(
  gene_binary,
  ref_binary
)

combined_MAGs <- c(MAGs, reference_genomes$MAG)


###### 4.2 重新计算module score
combined_module <- data.frame(row.names = combined_MAGs)

for(i in names(module_list)){
  
  combined_module[,i] <- module_score(
    combined_matrix,
    module_list[[i]]
  )
  
}

combined_module$rBOX <- rowMeans(
  combined_module[,c(
    "condensation",
    "reduction1",
    "dehydration",
    "reduction2"
  )]
)

###### 4.3 基于已有模型进行预测

combined_feature <- as.matrix(combined_matrix[,!colnames(combined_matrix) %in% core_genes])

prob_xgb <- predict(xgb_model, combined_feature)
prob_rf  <- predict(rf_model, combined_feature, type="prob")[,2]

combined_ml <- (prob_xgb + prob_rf)/2

###### 4.4 计算CEPS
combined_CEPS <- 0.3*combined_module$product_release +
  0.5*combined_ml +
  0.2*0.01*c(completeness, reference_genomes$Completeness)


###提取参考基因组的结果
ref_result <- data.frame(
  MAG = reference_genomes$MAG,
  CEPS = combined_CEPS[(length(MAGs)+1):length(combined_CEPS)],
  ML_prob = combined_ml[(length(MAGs)+1):length(combined_ml)]
)

write.table(ref_result,"results/tables/Ref_genome_CEPS_ML_prob.txt",sep = "\t",quote = FALSE,row.names = FALSE)
write.table(combined_module,"results/tables/combined_module.txt",sep = "\t",quote = FALSE,row.names = TRUE)

###### 4.5 验证分析
### CEPS分布比对
ref_result$group <- "Reference"
result_table$group <- "MAG"

plot_df <- rbind(
  result_table[,c("MAG","CEPS","group")],
  ref_result[,c("MAG","CEPS","group")]
)

#添加参考基因组分组信息
tp3 <- merge(plot_df, 
             subset(reference_genomes, select = c("MAG", "Group")), 
             by = "MAG", all = TRUE)
tp3[is.na(tp3)] <- "Other"
tp3 <- rename(tp3, Type = "Group")

# extract annotation label
tp3_reference_points <- tp3[tp3$group=="Reference",]
reference_genomes$MAG == tp3_reference_points$MAG
tp3_reference_points$label_annotation <- reference_genomes$label_annotation

# box plot
tp3_reference_points$x_jitter <- 
  as.numeric(factor(tp3_reference_points$group)) + 1 +
  runif(nrow(tp3_reference_points),-0.15,0.15)

ceps_ds <- ggplot(tp3, aes(group, CEPS))+
  geom_boxplot(outlier.shape = NA, fill = NA, color = "gray50", linewidth = 0.6, width=0.7)+
  geom_jitter(aes(fill = Type),
              shape = 21,
              color = "black",  stroke = 0.5,
              width = 0.25, size = 1.5)+
  scale_fill_manual(
    name = "CE-Type",
    breaks = c("CE-C6","CE-C4", "non-CE", "Other"),
    values = c( "#CD5B45", "#efe6c1", '#C5B8C6', "#BBBBBB" ))+
    #values = c( "#CD5B45", "#efd7b6", '#00688B', "#BBBBBB" ))+
  #scale_shape_manual(values = c(0, 13, 4, 1)) + 
  
  geom_text_repel(
    data=tp3_reference_points,
    aes(x=x_jitter, y=CEPS,label=label_annotation),
    size=3,
    color="black",
    max.overlaps = 30
  ) +
  
  theme_bw(base_size = 16)+
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "",
    axis.title = element_text(face="bold")
  )+
  
  scale_x_discrete(labels = c("MAG" = "MAG\n(n=1365)",
                              "Reference" = "Reference\n(n=42) "
  ))+
  scale_y_continuous(breaks = c(0.00, 0.20, 0.40, 0.60, 0.80, 1.00))+

  geom_hline(yintercept = c(0.6, 0.8), 
             linetype = "dashed", color = "gray", linewidth = 0.8)+
  #xlab("")+
  ylab("Chain Elongation Potential Score")
ceps_ds

ggsave(filename = "results/figures/Figure01_CPES_distrubution_compare_boxplot2.pdf", 
       ceps_ds, width = 5, height = 6)


### 排名验证
rank_df <- data.frame(
  MAG = combined_MAGs,
  CEPS = combined_CEPS
)

rank_df <- rank_df[order(-rank_df$CEPS),]
head(rank_df,30)

###### ROC Curve
ref_binary_ml <- ref_binary[,!colnames(ref_binary) %in% core_genes]

prob_ref_xgb <- predict(xgb_model, as.matrix(ref_binary_ml))
prob_ref_rf  <- predict(rf_model, ref_binary_ml, type="prob")[,2]

ml_prob_ref <- (prob_ref_xgb + prob_ref_rf)/2

# label
labels_ref <- reference_genomes$label
roc_ref_ens <- roc(labels_ref, ml_prob_ref)
roc_ref_xgb <- roc(labels_ref, prob_ref_xgb)
roc_ref_rf <- roc(labels_ref, prob_ref_rf)

auc_ens_ref <- auc(roc_ref_ens)
auc_xgb_ref <- auc(roc_ref_xgb)
auc_rf_ref  <- auc(roc_ref_rf)

ref_roc_df <- rbind(
  
  data.frame(FPR = 1-roc_ref_xgb$specificities, TPR = roc_ref_xgb$sensitivities,model = "XGBoost"),
  data.frame(FPR = 1-roc_ref_rf$specificities,TPR = roc_ref_rf$sensitivities,model = "RandomForest"),
  data.frame(FPR = 1-roc_ref_ens$specificities,TPR = roc_ref_ens$sensitivities,model = "Ensemble")
  
)

## plot
p_roc_ref <- ggroc(
  list(
    XGBoost = roc_ref_xgb,
    RandomForest = roc_ref_rf,
    Ensemble = roc_ref_ens
  ),
  legacy.axes = TRUE,
  size = 1.3
) +
  
  scale_color_manual(
    values = c(
      "XGBoost"="#CD5B45",
      "RandomForest"="#00688B",
      "Ensemble"="#BBBBBB"
    )
  )+
  
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  )+
  
  annotate(
    "text",
    x = 0.65,
    y = 0.20,
    label = paste0(
      "XGB AUC=", round(auc_xgb_ref,3),
      "\nRF AUC=", round(auc_rf_ref,3),
      "\nENS AUC=", round(auc_ens_ref,3)
    ),
    size=5
  )+
  
  theme_bw(base_size = 16)+
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face="bold")
  )+
  
  labs(
    x="False Positive Rate",
    y="True Positive Rate",
    color="Model"
  )
p_roc_ref

ggsave("results/figures/Figure02_extern_validation_ROC_curve_3.pdf", p_roc_ref, width = 7, height = 5)

###### Precision Recall Curve
library(precrec)

mmdata_obj <- mmdata(
  
  scores = list(
    prob_ref_xgb,
    prob_ref_rf,
    ml_prob_ref
  ),
  
  labels = list(
    labels_ref,
    labels_ref,
    labels_ref
  ),
  
  modnames = c(
    "XGBoost",
    "RandomForest",
    "Ensemble"
  )
)

# 计算
eval_obj <- evalmod(mmdata_obj)

# 提取AUC
auc_df <- auc(eval_obj)

print(auc_df)

# 提取PR部分
auc_pr <- auc_df[auc_df$curvetypes=="PRC", ]

# 生成文字
pr_text <- paste0(
  "XGB AUPRC=",
  round(auc_pr$aucs[auc_pr$modnames=="XGBoost"],3),
  
  "\nRF AUPRC=",
  round(auc_pr$aucs[auc_pr$modnames=="RandomForest"],3),
  
  "\nENS AUPRC=",
  round(auc_pr$aucs[auc_pr$modnames=="Ensemble"],3)
)

# 绘图
p_pr_ref <- autoplot(eval_obj, curvetype="PRC")+
  scale_color_manual(
    values=c(
      "XGBoost"="#CD5B45",
      "RandomForest"="#00688B",
      "Ensemble"="#BBBBBB"
    )
  )+
  
  annotate(
    "text",
    x=0.2,
    y=0.2,
    label=pr_text,
    size=5
  )+
  
  theme_bw(base_size = 16)+
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face="bold")
  )+
  
  labs(
    x="Recall",
    y="Precision"
  )

p_pr_ref

ggsave("results/figures/Figure03_extern_validation_PRC_curve.pdf", p_pr_ref, width = 7, height = 5)


############################# 5.统计与绘图模块 #################################

#===============================6.2 PCA plot ===================================  
pca <- PCA(combined_matrix,graph=FALSE)
pca_coord <- as.data.frame(pca$ind$coord)
pca_coord$MAG <- rownames(pca_coord)
pca_coord <- merge(pca_coord, tp3, by = 'MAG', all = TRUE)

all_class <- rep("Low",length(pca_coord$CEPS))
all_class[tp3$CEPS > 0.6] <- "Medium"
all_class[tp3$CEPS > 0.8] <- "High"


pca_coord$class <- all_class

pca_df <- pca_coord %>% select(-c("Dim.3", "Dim.4", "Dim.5")) %>% 
  rename("Dim1"="Dim.1", "Dim2"="Dim.2")

reference_points <- pca_df[pca_df$group=="Reference",]
reference_genomes$MAG == reference_points$MAG
reference_points$label_annotation <- reference_genomes$label_annotation


mag_points <- pca_df[pca_df$group=="MAG",]

write.table(pca_df, "pca_df.txt", quote = FALSE, sep = '\t', row.names = TRUE)

p_pcoa <- ggplot(pca_df,aes(x = Dim1,y = Dim2,color = class))+
  
  geom_point(size = 1.5, alpha = 0.7)+
  
  scale_colour_manual(
    name = "CE-Potential",
    breaks = c("High","Medium", "Low"),
    values = c(
      "High"="#4575B4",
      "Medium"="#91BFDB",
      "Low"="grey75"#d9d9d9
    ))+
  
  theme_bw(base_size = 16)+
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face="bold")
  )+
  
  xlab("PCoA1")+
  ylab("PCoA2")+
  
  theme(legend.position="right")

p_pcoa


p_pcoa_add <- p_pcoa +
  
  geom_point(data = reference_points, alpha=1)+
  
  geom_label_repel(
    data = reference_points,
    aes(label = label_annotation),
    size = 2,
    color = "black",
    fill = "white",
    segment.color = "black",
    segment.size = 0.4,
    box.padding = 0.5,
    max.overlaps = 100
  )

p_pcoa_add

###=== 统计检验 ===
## 1、整体差异检验
pca_scores <- pca$ind$coord
dist_matrix <- dist(pca_scores)
adonis_res <- adonis2(dist_matrix ~ class, data = pca_df, permutations = 999)

# 准备注释文本
r2 <- round(adonis_res$R2[1], 3)
pval <- adonis_res$`Pr(>F)`[1]
p_label <- ifelse(pval < 0.001, "< 0.001", paste0("= ", round(pval, 3)))
anno_text <- paste0("PERMANOVA: R² = ", r2, ", p ", p_label)

###=== High/Medium vs Low ===
pca_df$class2 <- ifelse(pca_df$class %in% c("High", "Medium"), 
                        "High/Medium", "Low")
# 将 "High/Medium" 设为参考水平（顺序不影响 adonis，但便于查看）
pca_df$class2 <- factor(pca_df$class2, levels = c("High/Medium", "Low"))
pca_scores <- pca$ind$coord
dist_matrix <- dist(pca_scores)

adonis2_res <- adonis2(dist_matrix ~ class2, data = pca_df, permutations = 999)
adonis2_res              

pseudo_f <- round(adonis2_res$F[1], 3)
pval <- adonis2_res$`Pr(>F)`[1]
p_label <- ifelse(pval < 0.001, "< 0.001", paste0("= ", round(pval, 3)))
anno_text <- paste0("PERMANOVA (High/Medium vs Low): \nPseudo-F = ", pseudo_f, ", p ", p_label)

p_pcoa_add <- p_pcoa_add +
  annotate("text", x = Inf, y = Inf, label = anno_text,
           hjust = 1.1, vjust = 1.5, size = 4, color = "black")
p_pcoa_add

ggsave("results/figures/Figure05_PCoA_substrate.pdf",p_pcoa_add,width=7,height=5)

#=============================== 6.3.UMAP ======================================
umap_res <- umap(combined_matrix, n_neighbors = 15, min_dist = 0.1, metric = 'euclidean')

umap_df <- data.frame(
  UMAP1 = umap_res[,1],
  UMAP2 = umap_res[,2],
  class=all_class,
  #substrate_type=combined_substrate$substrate_type,
  dataset = pca_df$group
)

### 统计检验
# 创建合并分组
umap_df$class2 <- factor(ifelse(umap_df$class %in% c("High","Medium"), 
                                "High/Medium", "Low"),
                         levels = c("High/Medium", "Low"))
dist_matrix <- dist(combined_matrix)
adonis_res <- adonis2(dist_matrix ~ class2, data = umap_df, permutations = 999)

# 准备注释
ump_F <- round(adonis_res$F[1], 3)
pval <- adonis_res$`Pr(>F)`[1]
p_label <- ifelse(pval < 0.001, "< 0.001", paste0("= ", round(pval, 3)))
anno_text <- paste0("PERMANOVA: \nPseudo-F = ", ump_F, ", p ", p_label,
                    "\n(High & Medium vs. Low)")


reference_points_umap <- umap_df[umap_df$dataset=="Reference",]

p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = class))+
  
  geom_point(size=1.5,alpha=0.7)+
  
  scale_colour_manual(
    name = "CE-Potential",
    breaks = c("High","Medium", "Low"),
    values = c(
      "High"="#4575B4",
      "Medium"="#91BFDB",
      "Low"="grey75"
    ))+

  theme_bw(base_size = 16)+
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face="bold")
  )+
  
  xlab("UMAP1")+
  ylab("UMAP2")

p_umap

p_umap_add <- p_umap +
  
  geom_point(data = reference_points_umap, alpha=1)+
  
  geom_label_repel(
    data = reference_points_umap,
    aes(label = reference_points$label_annotation),
    size = 2,
    color = "black",
    fill = "white",
    segment.color = "black",
    box.padding = 0.5,
    max.overlaps = 100
  )

p_umap_add <- p_umap_add + annotate("text", x = -0.1, y = -4.0, label = anno_text,
                       hjust = 1.1, vjust = 1.5, size = 3.5)

p_umap_add

ggsave("results/figures/Figure06_UMAP_substrate.pdf", p_umap_add, width=7, height=5)

###### 6.4 模型一致性分析
agreement_df <- data.frame(
  RF = prob_rf,
  XGB = prob_xgb
)

correlation <- cor(agreement_df$RF, agreement_df$XGB)

diff <- abs(agreement_df$RF - agreement_df$XGB)
mean_diff <- mean(diff)

cat("Model correlation:", correlation, "\n")
#Model correlation: 0.9712345 
cat("Mean prediction difference:", mean_diff, "\n")
#Mean prediction difference: 0.07789901 

p_agreement <- ggplot(agreement_df, aes(x = RF, y = XGB)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE, color="#00688B") +
  theme_bw(base_size = 16)+
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face="bold")
  )+
  labs(
    #title = paste("Model agreement: Pearson R =", round(correlation,3)),
    x = "Random Forest probability",
    y = "XGBoost probability"
  )
p_agreement

p_agreement <- p_agreement + stat_cor(method = "pearson",
         label.x.npc = "left",
         label.y.npc = "top",
         aes(label = paste(after_stat(rr.label), 
                           after_stat(p.label), 
                           sep = "~`,`~")))

ggsave("results/figures/Figure07_model_agreement_RF_XGB.pdf", p_agreement, width=6, height=5)


###### 6.5 模型内部评估
### ROC curve
roc_xgb <- roc(test_y, predict(xgb_model,test_x))
roc_rf  <- roc(test_y, predict(rf_model,test_x,type="prob")[,2])

prob_ens_test <- (
  predict(xgb_model,test_x) +
    predict(rf_model,test_x,type="prob")[,2]
)/2

roc_ens <- roc(test_y, prob_ens_test)

auc_xgb <- pROC::auc(roc_xgb)
auc_rf  <- pROC::auc(roc_rf)
auc_ens <- pROC::auc(roc_ens)

roc_df <- rbind(
  
  data.frame(FPR = 1-roc_xgb$specificities, TPR = roc_xgb$sensitivities,model = "XGBoost"),
  data.frame(FPR = 1-roc_rf$specificities,TPR = roc_rf$sensitivities,model = "RandomForest"),
  data.frame(FPR = 1-roc_ens$specificities,TPR = roc_ens$sensitivities,model = "Ensemble")
  
)

p_roc <- ggplot(roc_df, aes(FPR,TPR,color=model))+
  
  geom_line(size=1.3)+
  
  geom_abline(
    slope=1,
    intercept=0,
    linetype="dashed"
  )+
  scale_color_manual(values = c("#BBBBBB", "#00688B", "#CD5B45")) +
  
  annotate(
    "text",
    x=0.65,
    y=0.20,
    label=paste0(
      "XGB AUC=",round(auc_xgb,3),
      "\nRF AUC=",round(auc_rf,3),
      "\nENS AUC=",round(auc_ens,3)
    )
  )+
  
  theme_bw(base_size = 16)+
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face="bold")
  )+
  
  xlab("False Positive Rate")+
  ylab("True Positive Rate")
p_roc

ggsave("results/figures/Figure08_ML_ROC.pdf", p_roc, width=7, height=5)


### Precision Recall Curve
# xgb
pr_xgb <- pr.curve(
  scores.class0=predict(xgb_model,test_x)[test_y==1],
  scores.class1=predict(xgb_model,test_x)[test_y==0],
  curve=TRUE
)

# rf
pr_rf <- pr.curve(
  scores.class0=predict(rf_model,test_x)[test_y==1],
  scores.class1=predict(rf_model,test_x)[test_y==0],
  curve=TRUE
)

# ensemble
prob_xgb_all <- predict(xgb_model, test_x)        # 返回概率向量
prob_rf_all  <- predict(rf_model, test_x, type = "prob")[, "1"]  # 取出正类概率

# 计算集成模型的平均预测概率
prob_ensemble <- (prob_xgb_all + prob_rf_all) / 2

pr_ens <- pr.curve(
  scores.class0 = prob_ensemble[test_y == 1],
  scores.class1 = prob_ensemble[test_y == 0],
  curve = TRUE
)

# 合并至数据框中
proc_df <- rbind(
  
  data.frame(Recall = pr_xgb$curve[,1], Precision = pr_xgb$curve[,2], model = "XGBoost"),
  data.frame(Recall = pr_rf$curve[,1], Precision = pr_rf$curve[,2], model = "RandomForest"),
  data.frame(Recall = pr_ens$curve[, 1], Precision = pr_ens$curve[, 2], model = "Ensemble")
  
)

p_pr <- ggplot(proc_df, aes(Recall,Precision,color=model))+
  
  geom_line(size=1)+
  
  scale_color_manual(values = c("gray", "#00688B", "#CD5B45")) +
  
  annotate(
    "text",
    x=0.35,
    y=0.50,
    label=paste0(
      "XGB AUC=",round(pr_xgb$auc.integral,3),
      "\nRF AUC=",round(pr_rf$auc.integral,3),
      "\nENS AUC=",round(pr_ens$auc.integral,3)
    )
  )+
  
  theme_bw(base_size = 16)+
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face="bold")
  )+
  
  xlab("Recall")+
  ylab("Precision")
p_pr

ggsave("results/figures/Figure09_ML_PRcurve.pdf", p_pr, width=7, height=5)


### Regression plot
df_reg <- data.frame(ML_probability = ml_prob, CEPS = CEPS)

p_reg <- ggplot(
  df_reg,
  aes(ML_probability,CEPS)
)+
  
  geom_point(alpha=0.5, size=2)+
  geom_smooth(method="lm", color="#00688B")+
  
  theme_bw(base_size = 16)+
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face="bold")
  )+
  
  xlab("ML predicted probability")+
  ylab("Chain Elongation Potential Score")

p_reg <- p_reg + stat_cor(method = "pearson",
           label.x.npc = "left",
           label.y.npc = "top",
           aes(label = paste(after_stat(rr.label), 
                             after_stat(p.label), 
                             sep = "~`,`~")))

p_reg 

ggsave("results/figures/Figure10_ML_regression.pdf",p_reg,width=5,height=5)

### Feature importance
rf_imp <- importance(rf_model)

imp_df <- data.frame(gene = rownames(rf_imp),importance = rf_imp[,1])

imp_df <- imp_df %>%
  arrange(desc(importance)) %>%
  slice(1:20)

p_imp <- ggplot(
  imp_df,
  aes(reorder(gene,importance),importance))+
  
  geom_bar(stat="identity",fill="#00688B")+
  coord_flip()+
  theme_bw(base_size=14)+
  scale_y_continuous(expand = c(0.008, 0.01))+
  
  xlab("")+
  ylab("Importance")
p_imp

ggsave("results/figures/Figure11_RF_importance.pdf",p_imp,width=5,height=5)

### XGBoost importance
xgb_imp <- xgb.importance(model=xgb_model)
xgb_imp <- xgb_imp[1:20,]

p_xgb <- ggplot(xgb_imp,
                aes(reorder(Feature,Gain),Gain))+
  geom_bar(stat="identity",fill="#00688B")+
  coord_flip()+
  scale_y_continuous(expand = c(0.0001, 0.001))+
  
  theme_bw(base_size = 16)+
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face="bold")
  )
p_xgb

ggsave("results/figures/Figure12_XGB_importance.pdf", p_xgb, width=5, height=5)

#################################### END #######################################

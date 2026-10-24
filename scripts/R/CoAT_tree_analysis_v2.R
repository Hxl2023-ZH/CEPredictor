################################################################################
#                             CoAT tree analysis
################################################################################

# phylogeny
library(ape)
library(phytools)
library(treeio)
library(ggtree)


# data
library(dplyr)
library(tidyr)
library(data.table)

# plotting
library(ggplot2)
library(pheatmap)
library(viridis)
library(ggrepel)
library(ggforce)

# clustering
library(cluster)
library(factoextra)
library(dynamicTreeCut)

set.seed(2026)
#================================== tree Analysis ==============================
setwd("~/tree_analysis/")

tree <- read.tree("data/merged_coat.treefile")
tree

meta <- fread("data/MAGs_coat_annotation_metadata.txt")
head(meta)


CE_refs <- meta %>%
  filter(CE_potential %in% c("CE-C6", "CE-C4")) %>%
  pull(SequenceID)

dist_mat <- cophenetic.phylo(tree)
write.csv(dist_mat, "results/CoAT_patristic_distance.csv")

refs <- meta %>% filter(CE_potential %in% c("CE-C6", "CE-C4")) %>% pull(SequenceID)
mags <- meta %>% filter(CE_potential %in% c("High", "Medium", "Low")) %>% pull(SequenceID)

nearest_df <- data.frame()

for(m in mags){
  
  d <- dist_mat[m, refs]
  nearest_ref <- names(which.min(d))
  nearest_distance <- min(d)
  nearest_df <- rbind(
    nearest_df,
    data.frame(
      MAG = m,
      nearest_ref = nearest_ref,
      phylo_distance = nearest_distance
    )
  )
}

# 归一化
dmin <- min(nearest_df$phylo_distance)
dmax <- max(nearest_df$phylo_distance)

nearest_df$EPS <- 1 -
  (nearest_df$phylo_distance - dmin) /
  (dmax - dmin)

write.csv(nearest_df, "results/CoAT_nearest_df.csv")

#========= CoAT进化邻近性与链延长潜力相关
nearest_df2 <- nearest_df %>% left_join(meta, by = c("MAG" = "SequenceID"))

table(is.na(nearest_df2$CE_potential))

aggregate(EPS ~ CE_potential, nearest_df2, mean)

# plot
library(ggpubr)

my_comparisons <- list( c("High", "Medium"), c("High", "Low"), c("Medium", "Low") )

p_eps <-
  
  ggplot(nearest_df2, aes(x = CE_potential, y = EPS))+
  geom_boxplot(outlier.shape = NA, fill = NA, color = "gray50", linewidth = 0.5)+
  
  geom_jitter(aes(colour = CE_potential), width = 0.15, size = 0.8)+
  
  scale_colour_manual(
    values = c(
      "High"="#4575B4",
      "Medium"="#91BFDB",
      "Low"="grey75"
    )
  )+
  
  theme_bw(base_size = 16)+
  
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )+
  
  labs(
    x = "",
    y = "Evolutionary Proximity Score (EPS)"
  )+
  
  scale_x_discrete(limits = c("High", "Medium", "Low"))+
  
  stat_compare_means(
    comparisons = my_comparisons,    
    method = "wilcox.test",               
    label = "p.signif",              
    hide.ns = FALSE, 
    bracket.size = 0.5,
    step.increase = 0.15 
  )
p_eps

ggsave(filename = "results/CoAT_Fig06_Evolutionary_Proximity_Score.pdf", 
       p_eps, width = 5, height = 5)

#================================
# 聚类分析
hc <- hclust(as.dist(dist_mat), method = "average")

# PCoA
pcoa <- cmdscale(as.dist(dist_mat), eig = TRUE, k = 2)

pcoa_df <- data.frame(
  Sample = rownames(pcoa$points),
  PC1 = pcoa$points[,1],
  PC2 = pcoa$points[,2]
)

pcoa_df <- left_join(pcoa_df, meta, by = c("Sample" = "SequenceID"))

###=== k值判断
fviz_nbclust(as.matrix(dist_mat), FUN = hcut, method = "silhouette")
fviz_nbclust(as.matrix(dist_mat), FUN = hcut, method = "gap_stat")

### Dynamic
clusters_dy <- cutreeDynamic(
  dendro = hc,
  distM = as.matrix(dist_mat),
  deepSplit = 3,
  pamRespectsDendro = TRUE,
  minClusterSize = 15
)


cluster_df_dy <- data.frame(SequenceID = names(clusters_dy), Clade = clusters_dy)
table(cluster_df_dy$Clade, meta$CE_potential)

### 功能clade识别
clusters <- cutree(hc, k = 8)
cluster_df <- data.frame(SequenceID = names(clusters), Clade = clusters)
table(cluster_df$Clade, meta$CE_potential)

# 统计检验
tbl <- table(cluster_df$Clade, meta$CE_potential)
chisq.test(tbl)

#富集分析
meta$Clade <- cluster_df$Clade
target <- "8"

tbl2 <- table(meta$Clade == target, meta$CE_potential == "High")
tbl2
fisher.test(tbl2)

df <- as.data.frame(tbl)

# plot
freq_bar <- ggplot(df,
       aes(x = Var1,
           y = Freq,
           fill = Var2)) +
  geom_bar(stat = "identity",
           position = "fill") +
  scale_fill_manual(
    values = c(
      "CE-C6" = "#CD5B45",
      "CE-C4" = "#efe6c1",
      "High" = "#2E86AB",
      "Medium" = "#90BFDB",
      "Low" = "grey70"
    )
  ) +
  
  scale_y_continuous(expand = c(0.01, 0.01),labels = scales::percent) +
  theme_bw(base_size = 14)
freq_bar

ggsave("results/CoAT_Fig01_clade_freq_bar.pdf", freq_bar, width = 5, height = 5)

## Basic phylogeny
phy <- read.tree("data/merged_coat.treefile")
class(phy)
meta <- fread("data/MAGs_coat_annotation_metadata.txt")
head(meta)

tree_df <- ggtree(phy)$data
tree_df <- left_join(tree_df, meta, by = c("label" = "SequenceID"))
tree_df$Ref <- FALSE
tree_df$Ref[tree_df$CE_potential %in% c("CE-C4","CE-C6")] <- TRUE

p <- ggtree(phy, layout = "rectangular", size = 0.35)
p$data <- tree_df


p1 <- p +
  
  geom_treescale(
    x = 0,
    y = -20,
    width = 0.5
  )+
  
  geom_tippoint(aes(color = CE_potential), size = 0.8,alpha = 0.7) +
  
  # reference strain
  geom_tippoint(
    data = subset(tree_df, Ref),
    aes(fill = CE_potential),
    shape = 21,
    color = "black",
    stroke = 0.5,
    size = 1.5
  )+
  scale_color_manual(
    values = c(
      "CE-C6" = "#CD5B45",
      "CE-C4" = "#efe6c1",
      "High" = "#2E86AB",
      "Medium" = "#90BFDB",
      "Low" = "grey70"
    )
  ) +
  scale_fill_manual(
    values = c(
      "CE-C6" = "#CD5B45",
      "CE-C4" = "#efe6c1",
      "High" = "#2E86AB",
      "Medium" = "#90BFDB",
      "Low" = "grey70"
    )
  ) 

p1

### add reference labels
ref_df <- subset(tree_df, CE_potential %in% c("CE-C4","CE-C6"))

p1 <- p1 + geom_tiplab(
    data = ref_df,
    aes(label = labels),
    size = 1.2,
    fontface = "bold.italic",
    color = "black",
    offset = 0.02
  )


p1

ggsave("results/CoAT_Fig02_phylogeny_tree3.pdf", p1, width = 9, height = 12)

#================================== PCoA Analysis ==============================
## PCoA
pcoa <- cmdscale(
  as.dist(dist_mat),
  eig = TRUE,
  k = 2
)

pcoa_df <- data.frame(
  SequenceID =
    rownames(pcoa$points),
  
  PC1 = pcoa$points[,1],
  PC2 = pcoa$points[,2]
)


pcoa_df <- left_join(
  pcoa_df,
  meta,
  by = 'SequenceID'
)


## Define major lineages
pcoa_df$Lineage <- "Lineage-IV"

pcoa_df$Lineage[
  pcoa_df$PC1 > 1.5
] <- "Lineage-I"

pcoa_df$Lineage[
  pcoa_df$PC2 > 1
] <- "Lineage-II"

pcoa_df$Lineage[
  pcoa_df$PC2 < -1.3
] <- "Lineage-III"


## Variance explained
eig_percent <- round(
  100 * pcoa$eig / sum(pcoa$eig),
  1
)

pc1_var <- eig_percent[1]
pc2_var <- eig_percent[2]

## Convex hull
hulls <- pcoa_df %>%
  group_by(Lineage) %>%
  slice(chull(PC1, PC2))

### 统计检验
library(vegan)
adonis2_res <- adonis2(dist_mat ~ Group, data = meta, permutations = 999)

# 准备注释文本
pseudo_f <- round(adonis2_res$F[1], 3)
pval <- adonis2_res$`Pr(>F)`[1]
p_label <- ifelse(pval < 0.001, "< 0.001", paste0("= ", round(pval, 3)))
anno_text <- paste0("PERMANOVA: \nPseudo-F = ", pseudo_f, ", p ", p_label)

## PCoA plot
ref_pcoa_df <- subset(pcoa_df, CE_potential %in% c("CE-C4","CE-C6"))

write.table(pcoa_df, "CoAT_pcoa_df.txt", row.names = FALSE, sep = '\t', quote = FALSE)
write.table(ref_pcoa_df, "ref_CoAT_pcoa_df.txt", row.names = FALSE, sep = '\t', quote = FALSE)

p_sci <- ggplot() +
  
  ## lineage hull
  geom_polygon(data = hulls, 
               aes(PC1, PC2, fill = Lineage, group = Lineage),
               alpha = 0.6, color = NA) +
  
  ## all points
  geom_point(data = pcoa_df,
    aes(PC1, PC2, color = CE_potential),
    size = 1.5,
    alpha = 0.85,
    stroke = 0.2,
    position = position_jitter(width=0.03, height=0.03)) +
  
  ## CE references highlighted
  geom_point(
    data = ref_pcoa_df,
    aes(PC1, PC2, fill = CE_potential),
    size = 2,
    shape = 21,
    #fill = "gold",
    #color = "black",
    stroke = 0.3
  ) +
  
  geom_label_repel(
    data = ref_pcoa_df,
    aes(
      x = PC1,
      y = PC2,
      label = labels
    ),
    size = 1.5,
    color = "black",
    fill = "white",
    segment.color = "black",
    segment.size = 0.15,
    box.padding = 0.3,
    max.overlaps = 100
  )+
  
  scale_color_manual(
    values = c(
      "CE-C6" = "#CD5B45",
      "CE-C4" = "#efe6c1",
      "High" = "#2E86AB",
      "Medium" = "#90BFDB",
      "Low" = "grey70"
    )
  ) +
  
  ## lineage labels
  scale_fill_manual(values = c(
    "Lineage-I" = "#F4A582",
    "Lineage-II" = "#92C5DE",
    "Lineage-III" = "#B2ABD2",
    "Lineage-IV" = "#D5D5D5",
    "CE-C6" = "#CD5B45",
    "CE-C4" = "#efe6c1",
    "High" = "#2E86AB",
    "Medium" = "#90BFDB",
    "Low" = "grey70"
  )) +
  
  theme_bw(base_size = 16) +
  
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 1.2),
    
    axis.title = element_text(face = "bold",size = 16),
    axis.text = element_text(color = "black"),
    
    plot.title = element_text(hjust = 0.5,face = "bold",size = 16),
    legend.title = element_text(face = "bold")
  ) +
  
  labs(
    title = "Phylogenetic ordination of CoAT proteins",
    x = paste0("PCoA1 (",pc1_var,"%)"),
    y = paste0("PCoA2 (",pc2_var,"%)")
  )

p_sci

p_sci <- p_sci +
  annotate("text", x = Inf, y = Inf, label = anno_text,
           hjust = 1.1, vjust = 1.5, size = 4, color = "black")
p_sci

ggsave("results/CoAT_Fig03_PCoA1.pdf", p_sci, width = 7, height = 5)


#================================== UMAP Analysis ==============================

#library(uwot)
library(dbscan)

set.seed(2026)


#------
library(umap)

config <- umap.defaults
config$input <- "dist"

umap_res <- umap(
  as.matrix(dist_mat),
  config=config
)

umap_df <- data.frame(
  SequenceID = rownames(dist_mat),
  UMAP1 = umap_res$layout[,1],
  UMAP2 = umap_res$layout[,2]
)

umap_df <- left_join(umap_df, meta, by = "SequenceID")
### 统计检验
#dist_matrix <- dist(dist_mat)
adonis_res <- adonis2(dist_mat ~ Group, data = umap_df, permutations = 999)
adonis_res

ump_F <- round(adonis_res$F[1], 3)
pval <- adonis_res$`Pr(>F)`[1]
p_label <- ifelse(pval < 0.001, "< 0.001", paste0("= ", round(pval, 3)))
anno_text <- paste0("PERMANOVA: \nPseudo-F = ", ump_F, ", p ", p_label)


###
cl <- dbscan(umap_df[,c("UMAP1","UMAP2")],eps = 0.25,minPts = 10)
umap_df$Cluster <- factor(cl$cluster)
table(umap_df$Cluster)

### 计算每个cluster中的CE富集
cluster_summary <- umap_df %>% 
  group_by(Cluster) %>%
  summarise(
    N = n(),
    CE_C4 = sum(CE_potential=="CE-C4"),
    CE_C6 = sum(CE_potential=="CE-C6"),
    High =sum(CE_potential=="High"),
    Medium =sum(CE_potential=="Medium"),
    Low =sum(CE_potential=="Low"))


cluster_summary$CE_ratio <- (cluster_summary$CE_C4 + 
                               cluster_summary$CE_C6 + 
                               cluster_summary$High) / cluster_summary$N

cluster_summary %>% arrange(desc(CE_ratio))

# 富集热图
mat <- cluster_summary %>% select(CE_C4, CE_C6, High, Medium, Low)
rownames(mat) <- cluster_summary$Cluster

pdf("results/CoAT_Fig04_CE_clusters_heatmap.pdf", width=5, height=7)
pheatmap(
  mat,
  scale="row",
  clustering_method="ward.D2",
  fontsize=9
)
dev.off()


# 绘制Convex Hull
ref_clusters <- umap_df %>% filter(Type == "isolate") %>% distinct(Cluster)

ref_clusters

cluster_centers <- umap_df %>%
  filter(Cluster %in% ref_clusters$Cluster) %>%
  group_by(Cluster) %>%
  summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), n = n())

cluster_centers$Label <- paste0("CL-", cluster_centers$Cluster)

# plot
p_umap <- ggplot() +
  
  #Cluster边界
  geom_mark_hull(
  data = subset(umap_df, Cluster %in% ref_clusters$Cluster),
  aes(UMAP1, UMAP2, group = Cluster),
  fill = NA,
  color = "grey50",
  linewidth = 0.5,
  expand = unit(1.5, "mm")
  ) +
  
  # 所有MAG背景
  geom_point(data =subset(umap_df,Type == "MAG"),
  aes(UMAP1, UMAP2, color=CE_potential),
  alpha = 0.35,
  size = 1.6) +
  
  # 参考菌高亮
  geom_point(
  data =subset(umap_df, Type == "isolate"),
  aes(UMAP1, UMAP2, fill = CE_potential),
  shape = 21,
  color = "black",
  stroke = 1.2,
  size = 2.0) +
  
  # 4. Cluster标签
  geom_label_repel(
  data = cluster_centers,
  aes(UMAP1, UMAP2, label = Label),
  size = 2,
  fontface = "bold",
  fill = "white",
  box.padding = 0.5,
  label.size = 0.2,
  max.overlaps = Inf) +
  
  # 配色
  scale_fill_manual(
    values = c(
      "CE-C6" = "#CD5B45",
      "CE-C4" = "#efe6c1",
      "High" = "#2E86AB",
      "Medium" = "#90BFDB",
      "Low" = "grey70"
    )
  ) +
  
  scale_color_manual(
    values = c(
      "CE-C6" = "#CD5B45",
      "CE-C4" = "#efe6c1",
      "High" = "#2E86AB",
      "Medium" = "#90BFDB",
      "Low" = "grey70"
    )
  ) +
  
  #主题
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 1.2),
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.title = element_text(face = "bold")
  ) +
  
  labs(
    title = "Phylogenetic UMAP of CoAT proteins",
    x = "UMAP1",
    y = "UMAP2"
  )


p_umap <- p_umap + annotate("text", x = 16, y = 16, label = anno_text,
                                    hjust = 1.1, vjust = 1.5, size = 3.5)


p_umap

ggsave("", p_umap, width = 7, height = 5)


write.table(umap_df, "", row.names = FALSE, sep = '\t', quote = FALSE)






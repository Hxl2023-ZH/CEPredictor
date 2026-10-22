################################################################################
#                      iTOL annotation generator
################################################################################

library(dplyr)
library(readr)
library(RColorBrewer)

setwd("~/Desktop/zhongke/Liu_metagenomic/Analysis/10_Tree/annotation_files/")

### 0. 输入数据
# metadata 格式要求：
# MAG | Phylum | Genus | Chain_elongation | Sample1 | Sample2 ...

meta <- read_tsv("MAGs_all_metadata.txt")
meta <- meta %>% filter(Quality != "Low")

# 样本列自动识别
#sample_cols <- grep("^A", colnames(meta), value = TRUE)
sample_cols1 <- c("ACAT","BCD","EtfAB","ECH","HAD")
sample_cols2 <- c("ACT","CoAT")

### 1. 预处理
meta$Genus[is.na(meta$Genus) | meta$Genus == ""] <- "Unclassified"

# 计算平均丰度 + log
#meta$Mean_abundance <- rowMeans(meta[, sample_cols])
#meta$Mean_abundance_log <- log10(meta$Mean_abundance + 1e-6)

### 2. 通用配色（你提供的调色板）
base_colors <- c(
  "#CCE8EF", "#4BB5CD", "#008FA4", "#00688B", "#003E5A", "#FDD0A2", "#FDAE69", "#FDAE4B", "#F46333", "#CD5B45"
)

#"#E6F598", "#ABDDA4", "#66C2A5", 


get_colors <- function(n) {
  if(n <= length(base_colors)){
    return(base_colors[1:n])
  } else {
    return(colorRampPalette(base_colors)(n))
  }
}

### 3. Phylum 注释
phylum_levels <- unique(meta$Phylum)
phylum_colors <- get_colors(length(phylum_levels))

phylum_map <- data.frame(
  Phylum = phylum_levels,
  Color = phylum_colors
)

meta <- left_join(meta, phylum_map, by = "Phylum")

write_itol_colorstrip <- function(df, label, file, legend_title){
  
  levels <- unique(df$Group)
  colors <- df %>% distinct(Group, Color)
  
  writeLines(c(
    "DATASET_COLORSTRIP",
    "SEPARATOR TAB",
    paste0("DATASET_LABEL\t", label),
    "COLOR\t#000000",
    "STRIP_WIDTH\t25",
    "MARGIN\t5",
    "BORDER_WIDTH\t1",
    
    paste0("LEGEND_TITLE\t", legend_title),
    paste0("LEGEND_SHAPES\t", paste(rep(1, nrow(colors)), collapse = "\t")),
    paste0("LEGEND_COLORS\t", paste(colors$Color, collapse = "\t")),
    paste0("LEGEND_LABELS\t", paste(colors$Group, collapse = "\t")),
    
    "DATA"
  ), file)
  
  write.table(
    df[, c("MAG", "Color")],
    file,
    append = TRUE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
}

phylum_df <- data.frame(
  MAG = meta$MAG,
  Group = meta$Phylum,
  Color = meta$Color
)

write_itol_colorstrip(phylum_df, "Phylum", "itol_phylum.txt", "Phylum")


### 4. Genus 注释（方案1：全部属）

genus_levels <- unique(meta$Family)
genus_colors <- get_colors(length(genus_levels))


genus_colors <- c(
  "#CD5B45", "#8E3B2C", "#E27F6D", "#B34E67", "#F2C0B6", '#D598C6', "#C5B8C6",
  "#00688B", "#00445C", "#4A6FA5", "#3D98B5", '#009DD0', '#00AFEB', "#9DCFE0", 
  "#008F7A", "#E6F598", "#ABDDA4", "#66C2A5", "#FDD0A2", "#FDAE69", "#FDAE4B","#D47B41",
  "#A89F91", "#90A4AE", "#B8C4B1"
                )

genus_colors <- c("#8E3B2C", "#E27F6D", '#D598C6', "#C5B8C6", 
                  "#FDD0A2", "#FDAE69", "#FDAE4B", "#D47B41",
                  "#00445C", "#3D98B5", '#009DD0', '#00AFEB', 
                  "#9DCFE0", "#B8C4B1")

tmp <- c("#EAEAEA")


genus_map <- data.frame(
  Family = genus_levels,
  Color = genus_colors
)

meta_full <- left_join(meta, genus_map, by = "Family")

genus_full_df <- data.frame(
  MAG = meta_full$MAG,
  Group = meta_full$Family,
  Color = meta_full$Color
)

write_itol_colorstrip(
  genus_full_df,
  "Family",
  "itol_family_all.txt",
  "Family"
)


### 5. Genus 注释（方案2：Top10 + Others）

genus_rank <- meta %>%
  group_by(Genus) %>%
  summarise(total = sum(Mean_abundance)) %>%
  arrange(desc(total))

top10 <- genus_rank$Genus[1:10]

meta$Genus2 <- ifelse(meta$Genus %in% top10, meta$Genus, "Others")

genus2_levels <- unique(meta$Genus2)
genus2_colors <- get_colors(length(genus2_levels))

genus2_map <- data.frame(
  Genus2 = genus2_levels,
  Color = genus2_colors
)

# 强制 Others 灰色
genus2_map$Color[genus2_map$Genus2 == "Others"] <- "gray"

meta_top10 <- left_join(meta, genus2_map, by = "Genus2")

genus_top10_df <- data.frame(
  MAG = meta_top10$MAG,
  Group = meta_top10$Genus2,
  Color = meta_top10$Color
)

write_itol_colorstrip(
  genus_top10_df,
  "Genus (Top10)",
  "itol_genus_top10.txt",
  "Genus (Top10)"
)


### 6. 碳链延长潜力（4分类）
ce_colors <-  c(
  "High"="#4575B4",
  "Medium"="#91BFDB",
  "Low"="#EAEAEA")

meta$CE_color <- ce_colors[meta$Chain_elongation]

ce_df <- data.frame(
  MAG = meta$MAG,
  Group = meta$Chain_elongation,
  Color = meta$CE_color
)

write_itol_colorstrip(
  ce_df,
  "Chain elongation",
  "itol_ce.txt",
  "Chain elongation potential"
)


### 7. 底物类型（4分类）

sb_colors <- c(
  "Ethanol-Type" = "#F46363",
  "Lactate-Type" = "#FDAE89",
  "Mixed" = "#FDAE4B",
  "Undefined" = "#FDD0A0",
  "Non-CE" = "#d9d9d9"
)

meta$SB_color <- sb_colors[meta$Subsrate]

ce_df <- data.frame(
  MAG = meta$MAG,
  Group = meta$Subsrate,
  Color = meta$SB_color
)

write_itol_colorstrip(
  ce_df,
  "Subsrate",
  "itol_sb.txt",
  "Subsrate Type"
)



### 8. 丰度（bar）

writeLines(c(
  "DATASET_SIMPLEBAR",
  "SEPARATOR TAB",
  "DATASET_LABEL\tMean abundance",
  "COLOR\t#ABDDA4",
  
  "LEGEND_TITLE\tMean abundance (log10)",
  "LEGEND_SHAPES\t1",
  "LEGEND_COLORS\t#ABDDA4",
  "LEGEND_LABELS\tRelative abundance",
  
  "DATA"
), "itol_abundance.txt")

write.table(
  meta[, c("MAG", "Mean_abundance")],
  "itol_abundance.txt",
  append = TRUE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

###### TREE_LABELS

out <- rbind(
  c("LABELS\n"),
  c("SEPARATOR TAB\n"),
  c("DATA"),
  cbind(meta$MAG, meta$Labe)
)

write.table(
  out,
  file = "itol_tree_labels.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


###############################
# 7. 丰度注释（Heatmap）
###############################

#sample_cols <- grep("^Sample", colnames(meta), value=TRUE)

sample_cols <- sample_cols2

n_samples <- length(sample_cols)

# 使用你提供的配色（连续型）
heat_colors <- c("#FFFFFF", "#C0E3FA")

file <- "itol_abundance_heatmap2.txt"

writeLines(c(
  "DATASET_HEATMAP",
  "SEPARATOR TAB",
  "DATASET_LABEL\tMAG abundance",
  "COLOR\t#000000",
  
  paste0("FIELD_LABELS\t", paste(sample_cols, collapse="\t")),
  
  # 渐变色
  paste0("COLOR_MIN\t", heat_colors[1]),
  paste0("COLOR_MAX\t", tail(heat_colors,1)),
  
  # 可选：中间颜色（更平滑）
  paste0("COLOR_MID\t", heat_colors[5]),
  
  # legend
  "LEGEND_TITLE\tRelative abundance",
  "LEGEND_SHAPES\t1",
  paste0("LEGEND_COLORS\t", heat_colors[5]),
  "LEGEND_LABELS\tLow → High",
  
  "DATA"
), file)

write.table(
  meta[, c("MAG", sample_cols)],
  file,
  append=TRUE,
  sep="\t",
  row.names=FALSE,
  col.names=FALSE,
  quote=FALSE
)

###############################
# DONE
###############################

cat("All iTOL annotation files generated successfully!\n")














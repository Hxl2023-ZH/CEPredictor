#################################################################################
#                              rex TFBS analysis
#################################################################################
# Xinliang Hu
# 20250428
# R 4.5.2

#================================= 研究方法 ====================================
#--- promoter --- TSS --- 5'UTR --- ATG --- CDS ---
# TSS ≈ CDS起始位点（start codon）

# 在细菌中，转录通常起始于起始密码子上游区域，且转录起始位点（TSS）在启动子区域内存在一定程度的变异。
#在缺乏实验验证TSS数据的情况下，本研究采用基因起始密码子作为参考位点定义启动子区域，这是一种在原核基因组分析中广泛采用的近似方法。

#在缺乏实验数据的情况下，使用基因起始位点作为TSS的近似替代是一种标准且被广泛接受的方法，
#适用于motif分布和调控趋势分析，但不适用于精确位点解析。

#由于宏基因组组装基因组（MAG）缺乏实验验证的转录起始位点（TSS）信息，本研究采用基因起始密码子位置作为TSS的近似替代，
#以此定义启动子区域。这种方法在原核生物基因组分析中被广泛采用。

# We transformed motif-level signals into biologically meaningful regulatory metrics 
# by integrating motif significance, positional information, and genome annotation, 
# enabling quantitative assessment of Rex-mediated regulation at the genome scale.

#=== 1. Rex结合位点识别与定位

#为系统评估转录因子Rex在宏基因组组装基因组（MAG）中的调控潜力，首先基于构建的Rex motif数据库，
#利用 FIMO（Find Individual Motif Occurrences） 对各MAG的启动子区域序列进行扫描。
#启动子区域定义为相对于基因起始位点（TSS）的上游300 bp至下游100 bp区间（-300 ~ +100 bp）。

#FIMO输出结果经过筛选（q-value < 0.01），并结合基因组注释文件（GFF）进行解析，
#以获取每个motif对应的基因ID及其在基因组中的精确位置。


#=== 2. Motif与基因及启动子区域的关联

#通过整合GFF注释与预先构建的启动子区域BED文件，将FIMO识别的motif位点映射至对应基因的启动子区域。
#对于每个motif，计算其相对于目标基因转录起始位点（TSS）的相对位置（rel_pos），用于后续空间分布分析。

#仅保留位于启动子窗口（-300 ~ +100 bp）内的motif用于后续分析，以减少非功能性背景信号的干扰。

#=== 3. Rex调控强度计算

#为量化每个motif的潜在调控能力，定义Rex调控强度（Rex regulatory strength）如下：

#motif显著性采用 −log10(p-value) 表示；
#结合位点距离TSS的衰减效应通过指数函数建模：position_weight = exp(-|rel_pos| / d)
#其中，d为距离衰减参数（通常设为100 bp）。
#最终，单个位点的调控强度定义为：rex_strength = -log10(p-value) × position_weight


#=== 4. 单个MAG层面的调控特征汇总

#在单个MAG层面，对所有符合条件的motif进行聚合分析，计算如下指标：

#总调控强度（rex_strength_sum）
#平均调控强度（rex_strength_mean）
#motif数量（n_motif）
#覆盖基因数量（n_gene）

#同时，根据功能注释（如ACAT、EtfAB、BCD等模块），进一步计算模块特异性的调控强度，以解析Rex在不同代谢模块中的调控差异。



#============== working directory
#project/
#├── fimo/
#│   ├── MAG01_fimo_filter.txt
#│   ├── MAG02_fimo_filter.txt
#│   └── ...
#├── gff/
#│   ├── MAG01.gff
#│   └── ...
#├── promoter/
#│   ├── MAG01_promoter.bed
#│   └── ...
#├── rbox/
#│   ├── MAG01_rbox_gene.txt
#│   └── ...
#└── output/

setwd('~/TFBS_analysis/local/')

library(dplyr)
library(stringr)
library(ggplot2)
library(GenomicRanges)

#============== core function
process_MAGs <- function(mag_id, base_dir='.') {
  
  message("Processing: ", mag_id)
  
  # ---------- 1. 读取 ----------
  fimo_file <- file.path(base_dir, "fimo", paste0(mag_id, "_fimo_filter.txt"))
  gff_file  <- file.path(base_dir, "gff", paste0(mag_id, ".gff"))
  bed_file  <- file.path(base_dir, "promoter", paste0(mag_id, "_promoter.bed"))
  rbox_file <- file.path(base_dir, "rbox", paste0(mag_id, "_rbox_gene.txt"))
  
  if(!file.exists(fimo_file)) return(NULL)
  
  fimo <- read.table(fimo_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  
  gff <- read.table(gff_file, sep="\t", comment.char="#", header=FALSE)
  colnames(gff) <- c("seqid","source","type","start","end",
                     "score","strand","phase","attr")
  
  promoter <- read.table(bed_file, header=FALSE)
  colnames(promoter) <- c("seqid","p_start","p_end","gene_id","score","strand")
  
  
  rbox <- read.table(rbox_file, header=FALSE)
  colnames(rbox) <- c("gene_id","module")
  
  
  genes <- gff %>%
    filter(type=="CDS") %>%
    mutate(
      gene_id = str_extract(attr, "ID=[^;]+") %>%
        str_replace("ID=", ""),
      product = str_extract(attr, "product=[^;]+") %>%
        str_replace("product=", "")
    ) %>%
    select(gene_id, seqid, start, end, strand, product)
  
  genes <- genes %>% rename_with(~ "strand_gene", .cols = "strand")
  genes <- genes %>% dplyr::rename(gene_start = start, gene_end = end)
  promoter <- promoter %>% rename_with(~ "strand_promoter", .cols = "strand")
  
  promoter <- promoter %>% distinct(gene_id, .keep_all = TRUE)
  
  genes <- genes %>% distinct(gene_id, .keep_all = TRUE)
  rbox <- rbox %>% distinct(gene_id, .keep_all = TRUE)
  
  dat <- fimo %>%
    left_join(genes, by="gene_id") %>%
    left_join(promoter, by="gene_id", suffix=c(".gene",".promoter")) %>%
    left_join(rbox, by="gene_id") %>%
    mutate(MAG = mag_id)
  
  # 坐标转换
  dat <- dat %>%
    mutate(
      
      TSS = gene_start,
      
      motif_genomic_start = ifelse(strand_gene == "+",
                                   p_start + start - 1,
                                   p_end - start + 1),
      
      rel_pos = ifelse(strand_gene == "+",
                       motif_genomic_start - TSS,
                       -(TSS - motif_genomic_start))
    )

  
  # 判断
  dat <- dat %>%
    mutate(
      in_promoter = rel_pos >= -300 & rel_pos <= 100
    )
  
  # 调控证据整合
  dat <- dat %>%
    mutate(
      strong_motif = q.value < 1e-5
    )
  
  # gene-level
  gene_summary <- dat %>%
    group_by(MAG, gene_id, product, module) %>%
    summarise(
      n_sites = n(),
      best_q = min(q.value),
      promoter_hit = any(in_promoter),
      strong_hit = any(strong_motif),
      .groups="drop"
    ) %>%
    mutate(
      regulated = promoter_hit & strong_hit
    )
  
  # MAG-level
  mag_summary <- gene_summary %>%
    summarise(
      MAG = mag_id,
      total_genes = n(),
      CE_genes = sum(!is.na(module)),
      CE_regulated = sum(!is.na(module) & regulated),
      regulated_genes = sum(regulated)
    ) %>%
    mutate(
      Rex_regulation = case_when(
        CE_regulated >= 2 ~ "Strong",
        CE_regulated == 1 ~ "Moderate",
        TRUE ~ "Weak"
      )
    )
  
  # save result
  write.table(dat,
              file.path(base_dir, "output", paste0(mag_id, "_motif.tsv")),
              sep="\t", quote=FALSE, row.names=FALSE)
  
  write.table(gene_summary,
              file.path(base_dir, "output", paste0(mag_id, "_gene.tsv")),
              sep="\t", quote=FALSE, row.names=FALSE)
  
  return(list(gene=gene_summary, mag=mag_summary, raw=dat))
  
}

#===================== Analysis
base_dir <- "project"

mag_list <- gsub("_fimo_filter.txt", "",
                 list.files(file.path(base_dir, "fimo"),
                            pattern="_fimo_filter.txt"))

results <- lapply(mag_list, process_MAGs, base_dir=base_dir)

# 合并结果
gene_all <- bind_rows(lapply(results, function(x) x$gene))
mag_all  <- bind_rows(lapply(results, function(x) x$mag))
raw_all  <- bind_rows(lapply(results, function(x) x$raw))

write.table(gene_all, "output/all_gene_regulation_local.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)

write.table(mag_all, "output/all_MAG_summary_local.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)

write.table(raw_all, "output/raw_all_local.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)


# 检查峰分布
ggplot(raw_all, aes(rel_pos)) + geom_histogram()


#==================== plot

theme_pub <- function(){
  theme_bw(base_size = 12) +
    theme(
      #panel.grid = element_blank(),
      axis.title = element_text(face="bold"),
      axis.text = element_text(color="black"),
      strip.background = element_rect(fill="grey90"),
      legend.position = "right"
    )
}

#raw_all <- raw_all %>%
#  mutate(module2 = ifelse(is.na(module),"Other",module))

cols_reg <- c("0"="#BDBDBD","1"="#D73027")  # 未调控 vs 调控
cols_module <- c("#BDBDBD", "#E6F598", "#66C2A5", "#CD5B45","#00688B")

p1 <- ggplot(raw_all,
             aes(x=rel_pos, fill=module)) +
  geom_histogram(binwidth=50, alpha=0.8, position="identity",color="black") +
  
  geom_vline(xintercept=0, linetype="dashed", size=0.6) +
  scale_fill_manual(values=cols_module) +
  coord_cartesian(xlim=c(-300,100)) +
  labs(x="Position relative to TSS (bp)",
       y="Rex TFBS Count",
       fill="Module") +
  theme_pub()
p1

ggsave("output/Fig1_motif_position_ref_raw.pdf", p1, width=6, height=4)

##

raw_plot <- raw_all %>%
  filter(rel_pos >= -300 & rel_pos <= 100)

raw_plot <- raw_plot %>%
  mutate(
    bin = cut(rel_pos,
              breaks = seq(-300, 100, by=50),
              include.lowest = TRUE)
  )

bin_stat <- raw_plot %>%
  group_by(bin, module) %>%
  summarise(
    count = n(),
    .groups="drop"
  ) %>%
  group_by(module) %>%
  mutate(
    proportion = count / sum(count)
  )

p1a <- ggplot(bin_stat,
              aes(x=bin, y=count, fill=module)) +
  geom_bar(stat="identity", position="stack", color="black", size=0.2) +
  
  geom_text(aes(label=count),
            position=position_stack(vjust=0.5),
            size=3) +
  
  scale_fill_manual(values=cols_module) +
  #coord_cartesian(xlim=c(-300,100)) +
  
  labs(x="Position relative to TSS (bp)",
       y="Rex TFBS count",
       fill="Module") +
  
  theme_bw() +
  geom_vline(xintercept=which(levels(bin_stat$bin) == "(-50,0]"),
             linetype="dashed")+
  annotate("text", x=6, y=max(bin_stat$count), label="TSS", size=4)

p1a

ggsave("output/Fig1_motif_bin_count_ref.pdf", p1a, width=7, height=4)


p1b <- ggplot(bin_stat,
              aes(x=bin, y=proportion, fill=module)) +
  geom_bar(stat="identity", position="fill", color="black") +
  
  scale_fill_manual(values=cols_module) +
  
  labs(x="Position relative to TSS (bp)",
       y="Proportion",
       fill="Module") +
  
  theme_bw()
p1b
ggsave("output/Fig1_motif_bin_prop.pdf", p1b, width=7, height=4)

#
gene_all <- gene_all %>%
  mutate(
    type = ifelse(is.na(module),"Non-CE","CE"),
    regulated = as.factor(regulated)
  )

p2 <- ggplot(gene_all,
             aes(x=type, fill=regulated)) +
  geom_bar(position="fill", color="black") +
  scale_fill_manual(values=cols_reg,
                    labels=c("No","Yes")) +
  labs(x="Gene category",
       y="Proportion",
       fill="Rex regulation") +
  theme_pub()
p2

ggsave("output/Fig2_CE_vs_nonCE.pdf", p2, width=4.5, height=4)


#
p3 <- ggplot(gene_all %>% filter(!is.na(module)),
             aes(x=module, fill=regulated)) +
  geom_bar(position="fill", color="black") +
  scale_fill_manual(values=cols_reg,
                    labels=c("No","Yes")) +
  labs(x="Functional module",
       y="Proportion",
       fill="Rex regulation") +
  theme_pub()
p3 

ggsave("output/Fig3_module_regulation.pdf", p3, width=5, height=4)


#
p4 <- ggplot(gene_all,
             aes(x=n_sites, fill=type)) +
  geom_histogram(binwidth=1, position="dodge", color="black") +
  scale_fill_manual(values=c("CE"="#E41A1C","Non-CE"="#377EB8")) +
  labs(x="Number of motifs per gene",
       y="Gene count",
       fill="Gene type") +
  theme_pub()
p4

ggsave("output/Fig4_motif_count.pdf", p4, width=5, height=4)


#
gene_all <- gene_all %>%
  left_join(raw_all %>%
              select(gene_id, start.y, end, MAG) %>%
              distinct(),
            by=c("gene_id","MAG")) %>%
  mutate(
    gene_len = abs(end - start.y),
    density = n_sites / (gene_len/1000)
  )

p5 <- ggplot(gene_all,
             aes(x=type, y=density, fill=type)) +
  geom_boxplot(outlier.shape=NA, alpha=0.8, width = 0.4) +
  geom_jitter(width=0.15, size=1, alpha=0.4) +
  scale_fill_manual(values=c("CE"="#E41A1C","Non-CE"="#377EB8")) +
  labs(x="Gene type",
       y="Motif density (per kb)") +
  theme_pub()
p5

ggsave("output/Fig5_density.pdf", p5, width=5, height=4)

###

meta_bin <- raw_plot %>%
  mutate(
    pos_bin = round(rel_pos / 5) * 5   # 5 bp resolution
  ) %>%
  group_by(MAG, module, pos_bin) %>%
  summarise(count = n(), .groups="drop")

meta_bin <- meta_bin %>%
  group_by(MAG, module) %>%
  mutate(
    density = count / sum(count)
  ) %>%
  ungroup()

meta_summary <- meta_bin %>%
  group_by(module, pos_bin) %>%
  summarise(
    mean_density = mean(density),
    sd_density   = sd(density),
    n_mag        = n(),
    se           = sd_density / sqrt(n_mag),
    ci_upper     = mean_density + 1.96 * se,
    ci_lower     = mean_density - 1.96 * se,
    .groups="drop"
  )

p_meta2 <- ggplot(meta_bin,
                  aes(x=pos_bin,
                      y=density,
                      color=module)) +
  
  geom_smooth(se=TRUE, method="loess", span=0.2) +
  
  geom_vline(xintercept=0, linetype="dashed") +
  
  labs(x="Position relative to TSS (bp)",
       y="Motif density",
       color="Module") +
  
  theme_bw()
p_meta2

ggsave("output/Fig_meta_loess.pdf", p_meta2, width=6, height=4)


##################################### END ######################################



################################################################################
#                            Rex regulation & CE-MAGs Analysis 
################################################################################
# Xinliang Hu
# 2026
# R 4.5.2

# Motif analyses were restricted to promoter regions (-300 to +100 bp relative to TSS) to focus on regulatory elements.

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(rstatix)
library(broom)

setwd('~/Desktop/zhongke/Liu_metagenomic/Analysis/09_Rex_TF/fimo/TFBS_analysis/')

dat <- read.table("MAG_rex-TFBS_merged_table.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

table(dat$module)
# ACAT   BCD   ECH EtfAB   HAD 
# 353   161     9   144   183 

dat <- dat %>% filter(module != "ECH") 


dat <- dat %>%
  mutate(
    Dataset = factor(Dataset, levels=c("reference","public","local")),
    CE_potential = factor(CE_potential, levels=c("High","Medium","Low")),
    module = factor(module, levels=c("ACAT","EtfAB","BCD","HAD")),
    logP = -log10(p.value),
    in_promoter = as.logical(in_promoter)
  )

dat_plot <- dat %>% filter(rel_pos >= -300, rel_pos <= 100)

#============================ Figure 01 =========================================

# build label dataset
reference_point <- dat_plot %>% filter(Dataset == "reference")
reference_genome <- read.delim("ref_genome_label.txt", header = TRUE, check.names = FALSE)
reference_genome <- reference_genome %>% filter(MAG %in% reference_point$MAG)

reference_data <- merge(reference_point, reference_genome, by = "MAG")

write.table(reference_data, file = "reference_rex_result.txt", quote = FALSE, sep = '\t', row.names = FALSE)
reference_data <-  read.delim("reference_rex_result.txt", header = TRUE, check.names = FALSE)

table(dat_plot$module)
# plot
#ACAT EtfAB   BCD   ECH   HAD 
# 353   143   161     9   181 

p1 <- ggplot(dat_plot,
             aes(x = rel_pos, y = logP, fill = module)) +
  
  geom_point(shape = 21, color = "gray60", alpha=0.6, size = 1.5) +
  
  # 高亮 reference
  geom_point(data = reference_data, shape = 21, color = "black", alpha=1, size = 1.5)+
  
  geom_label_repel(
    data = reference_data,
    aes(label = label_annotation),
    size = 1.5,
    color = "black",
    fill = "white",
    segment.color = "black",
    segment.size = 0.4,
    box.padding = 0.5,
    max.overlaps = 100
  )+
  
  coord_cartesian(xlim=c(-300,100)) +
  
  geom_vline(xintercept=0, linetype="dashed") +
  geom_vline(
    xintercept=c(-200, -50),
    linetype="dashed",
    color="#999999",
    linewidth=1
  )+
  
#  scale_colour_manual(
#    name = "CE-Potential",
#    breaks = c("High","Medium", "Low"),
#    values = c(
#      "High"="#4575B4",
#      "Medium"="#91BFDB",
#      "Low"="grey75"
#    ))+
  
  scale_fill_manual(name = "Module",
                     breaks = c("BCD", "ACAT", "HAD", "EtfAB"),
                     values = c("BCD" = "#b86161", "ACAT"="#efe6c1","HAD"="#ABDDA4", "EtfAB"="#aecae8"),
                     labels = c("BCD (n=161)", "ACAT (n=353)", "HAD (n=181)", "EtfAB (n=143)"))+

  scale_x_continuous(breaks = c(-300, -200, -100, -50, 0, 100),
                     labels = c("-300", "-200", "-100", "50", "0", "100")) +
  
  labs(
    x="Position relative to TSS (bp)",
    y=expression(-log[10](p.value)),
    fill="Module"
    #shape="CE potential"
  ) +
  
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face="bold")
  )

p1

ggsave("Rex_Fig01A_scatter_TSS.pdf", p1, width=9, height=5)

###
p1B <- ggplot(dat_plot,
             aes(x = rel_pos, y = logP, shape = CE_potential)) +
  
  geom_point(color = "black", alpha=0.4, size = 1.5) +
  
  # 高亮 reference
  geom_point(data = reference_data, color = "black", alpha=1, size = 1.5)+
  
  geom_label_repel(
    data = reference_data,
    aes(label = label_annotation),
    size = 1.5,
    color = "black",
    fill = "white",
    segment.color = "black",
    segment.size = 0.4,
    box.padding = 0.5,
    max.overlaps = 100
  )+
  
  coord_cartesian(xlim=c(-300,100)) +
  
  geom_vline(xintercept=0, linetype="dashed") +
  geom_vline(
    xintercept=c(-200, -50),
    linetype="dashed",
    color="#999999",
    linewidth=1
  )+
  
  scale_shape_manual(
    name = "CE-Potential",
      breaks = c("High","Medium", "Low"),
      values = c(0, 1, 2))+
  
  scale_x_continuous(breaks = c(-300, -200, -100, -50, 0, 100),
                     labels = c("-300", "-200", "-100", "50", "0", "100")) +
  
  labs(
    x="Position relative to TSS (bp)",
    y=expression(-log[10](p.value)),
    #fill="Module"
    shape="CE potential"
  ) +
  
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face="bold")
  )

p1B

ggsave("Rex_Fig01B_scatter_TSS.pdf", p1B, width=9, height=5)



#============================== Figure 02A =====================================
# “High CE genes show stronger enrichment near TSS”
# Rex motif 是否在不同 CE 潜力中表现出不同的空间分布模式？

# Rex motifs exhibit clear position-specific enrichment patterns across genomes 
#with different CE potentials, with the strongest and most focused enrichment 
#observed in high CE genomes, particularly for the ACAT module.

#High CE genomes show a strong and highly localized enrichment of Rex motifs 
#upstream of TSS, particularly for ACAT-associated genes, suggesting tight and 
#position-specific transcriptional regulation.

#In medium CE genomes, motif enrichment becomes more dispersed, with BCD-associated 
#motifs showing moderate enrichment, indicating a less specialized and potentially 
#transitional regulatory state.

#Low CE genomes exhibit a diffuse and flattened motif distribution, lacking clear 
#enrichment near TSS, suggesting weak or non-specific Rex-mediated regulation.

#These results reveal a progressive transition from highly localized and module-specific 
#Rex binding in high CE genomes to diffuse and non-specific distributions in low CE genomes, 
#suggesting that regulatory precision rather than mere motif presence is associated with 
#carbon chain elongation potential.

p2 <- ggplot(dat_plot, aes(x=rel_pos, color=module, fill=module)) +
  
  geom_density(alpha=0.2, size=1) +
  geom_vline(xintercept=0, linetype="dashed") +
  facet_wrap(~CE_potential, nrow=1) +
  
  scale_color_manual(values=c(
    "ACAT"="#efe6c1",
    "BCD"="#b86161",
    "ECH"="#6A3D9A",
    "EtfAB"="#aecae8",
    "HAD"="#ABDDA4"
  )) +
  
  scale_fill_manual(values=c(
    "ACAT"="#efe6c1",
    "BCD"="#b86161",
    "ECH"="#6A3D9A",
    "EtfAB"="#aecae8",
    "HAD"="#ABDDA4"
  )) +
  
  labs(x="Position relative to TSS (bp)",
       y="Density") +
  
  coord_cartesian(xlim=c(-300,100)) +
  theme_bw(base_size = 16)+
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face="bold")
  )

p2

ggsave("Rex_Fig02A_density_metagene.pdf", p2, width=8, height=5)

### 统计检验
kruskal.test(rel_pos ~ CE_potential, data=dat_plot)
# p-value = 0.5921
ks.test(
  dat$rel_pos[dat_plot$CE_potential=="High"],
  dat$rel_pos[dat_plot$CE_potential=="Low"]
)

pairwise.wilcox.test(dat_plot$rel_pos,
                     dat_plot$CE_potential,
                     p.adjust.method="BH")

#============================== Figure 02B =====================================

p_box <- ggplot(dat_plot,
                aes(x=CE_potential,
                    y=rel_pos,
                    fill=CE_potential)) +
  
  geom_boxplot(outlier.shape=NA, width=0.6) +
  geom_jitter(width=0.15, alpha=0.5, size=1) +
  
  stat_compare_means(method="kruskal.test") +
  
  stat_compare_means(
    comparisons=list(c("High","Medium"),
                     c("High","Low"),
                     c("Medium","Low")),
    method="wilcox.test",
    label="p.signif"
  ) +
  
  scale_fill_manual(
    name = "CE-Potential",
    breaks = c("High","Medium", "Low"),
    values = c(
      "High"="#4575B4",
      "Medium"="#91BFDB",
      "Low"="grey75"
    ))+
  
  scale_x_discrete(labels = c("High" = "High\n(n=89)",
                              "Medium" = "Medium\n(n=137)",
                              "Low" = "Low\n(n=135)"))+
  
  labs(x="CE potential",
       y="Motif position relative to TSS (bp)") +
  
  theme_bw(base_size = 14)+
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face="bold"),
    legend.position = ""
  )
p_box

ggsave("Rex_Fig2B_boxplot.pdf", p_box, width=5, height=5)


#============================== Figure 02C =====================================

dat_region <- dat_plot %>%
  mutate(region = case_when(
    rel_pos >= -50 & rel_pos <= 0 ~ "Core",
    rel_pos >= -200 & rel_pos < -50 ~ "Proximal",
    TRUE ~ "Distal"
  ))

region_stat <- dat_region %>%
  group_by(CE_potential, region) %>%
  summarise(n=n(), .groups="drop") %>%
  group_by(CE_potential) %>%
  mutate(prop = n/sum(n))

### 统计检验
# 1、不同 CE_potential 水平下，region 的分布是否相同？
#  卡方独立性检验
chisq.test(table(dat_region$region, dat_region$CE_potential))
# X-squared = 19.614, df = 4, p-value = 0.0005951
# 即不同 CE_potential 组的 region 分布有显著差异。

# 2、对每个region，比较不同 CE_potential 组的比例
# To assess whether the proportions of promoter regions differed among CE potential groups, 
# pairwise comparisons were performed using Fisher’s exact tests on 2 × 2 contingency tables 
# constructed from region counts and non-region counts within each group. 
# P values were adjusted for multiple testing using the Benjamini–Hochberg false discovery rate (FDR) method.

total_df <- region_stat %>%
  group_by(CE_potential) %>%
  summarise(total = sum(n))

df <- left_join(region_stat, total_df, by = "CE_potential")
#library(rstatix)

pairwise_results <- list()
regions <- unique(df$region)

for (r in regions) {
  
  sub <- df %>% filter(region == r)
  
  combs <- combn(levels(df$CE_potential), 2, simplify = FALSE)
  
  for (cc in combs) {
    
    g1 <- cc[1]
    g2 <- cc[2]
    
    d1 <- sub %>% filter(CE_potential == g1)
    d2 <- sub %>% filter(CE_potential == g2)
    
    # 2×2 contingency table
    mat <- matrix(
      c(d1$n,
        d1$total - d1$n,
        d2$n,
        d2$total - d2$n),
      nrow = 2,
      byrow = TRUE
    )
    
    p <- fisher.test(mat)$p.value
    
    pairwise_results[[length(pairwise_results)+1]] <-
      data.frame(
        region = r,
        group1 = g1,
        group2 = g2,
        p = p
      )
  }
}

pairwise_results <- bind_rows(pairwise_results)


##
pairwise_results <- pairwise_results %>%
  mutate(
    p.adj = p.adjust(p, method = "fdr"),
    p.signif = case_when(
      p.adj <= 0.001 ~ "***",
      p.adj <= 0.01  ~ "**",
      p.adj <= 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

pairwise_results

##
pairwise_results <- pairwise_results %>%
  group_by(region) %>%
  mutate(
    y.position = c(0.72, 0.80, 0.88)
  )

###
stat_df <- pairwise_results %>%
  mutate(
    
    # region对应x轴位置
    x_center = case_when(
      region == "Core" ~ 1,
      region == "Distal" ~ 2,
      region == "Proximal" ~ 3
    ),
    
    # group1横坐标
    xmin = case_when(
      group1 == "Low" ~ x_center - 0.3,
      group1 == "Medium" ~ x_center,
      group1 == "High" ~ x_center + 0.3
    ),
    
    # group2横坐标
    xmax = case_when(
      group2 == "Low" ~ x_center - 0.3,
      group2 == "Medium" ~ x_center,
      group2 == "High" ~ x_center + 0.3
    )
  )

stat_df <- stat_df %>%
  group_by(region) %>%
  mutate(
    y.position = c(0.72, 0.80, 0.88)
  )

#plot
p_region <- ggplot(region_stat, aes(x=region, y=prop, fill=CE_potential)) +
  
  geom_bar(stat="identity",
           position=position_dodge(0.9),
           color='black',
           width=0.8) +
  
  geom_text(aes(label=scales::percent(prop,
                                      accuracy=0.1)),
            position=position_dodge(0.9),
            vjust=-0.3,
            size=3.0) +
  
  scale_fill_manual(
    name = "CE-Potential",
    breaks = c("High","Medium", "Low"),
    values = c(
      "High"="#4575B4",
      "Medium"="#91BFDB",
      "Low"="grey75"
    ))+
  
  scale_x_discrete(labels = c(
    "Core" = "Core\n[-50,0]",
    "Distal" = "Distal\n(other)",
    "Proximal" = "Proximal\n[-200,-50)"
  )) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0,0.15))
  ) +
  
  labs(x="Promoter region",
       y="Proportion") +
  
  theme_bw(base_size = 16)+
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face="bold"),
    legend.position = ""
  )

p_region

p_region <- p_region +
  
  geom_segment(data=stat_df,
               aes(x=xmin,
                   xend=xmax,
                   y=y.position,
                   yend=y.position),
               inherit.aes = FALSE) +
  
  geom_segment(data=stat_df,
               aes(x=xmin,
                   xend=xmin,
                   y=y.position,
                   yend=y.position-0.015),
               inherit.aes = FALSE) +
  
  geom_segment(data=stat_df,
               aes(x=xmax,
                   xend=xmax,
                   y=y.position,
                   yend=y.position-0.015),
               inherit.aes = FALSE) +
  
  geom_text(data=stat_df,
            aes(x=(xmin+xmax)/2,
                y=y.position+0.015,
                label=p.signif),
            size=5,
            inherit.aes = FALSE)
p_region

ggsave("Rex_Fig2C_region.pdf", p_region, width=5, height=5)

#============================== Figure 03 ======================================
#构建调控强度指标
#To resolve inconsistencies, CE potential was harmonized at the MAG level
# using a maximum-ranking strategy (High > Medium > Low).

dat2 <- dat_plot %>%
  filter(in_promoter == TRUE) %>%
  mutate(
    score = -log10(p.value),
    pos_weight = exp(-abs(rel_pos)/100),# 距离TSS越近权重越高
    rex_strength = score * pos_weight
  )

dat2 <- dat2 %>%
  mutate(
    CE_rank = case_when(
      CE_potential == "High" ~ 3,
      CE_potential == "Medium" ~ 2,
      CE_potential == "Low" ~ 1
    )
  )

mag_ce <- dat2 %>%
  group_by(MAG) %>%
  summarise(
    CE_rank = max(CE_rank),
    .groups="drop"
  ) %>%
  mutate(
    CE_potential = case_when(
      CE_rank == 3 ~ "High",
      CE_rank == 2 ~ "Medium",
      TRUE ~ "Low"
    )
  )


dat2 <- dat2 %>%
  dplyr::select(-CE_potential) %>%
  left_join(mag_ce, by="MAG")

# 一致性检查
stopifnot(
  dat2 %>%
    group_by(MAG) %>%
    summarise(n=n_distinct(CE_potential)) %>%
    pull(n) %>%
    max() == 1
)

#
mag_stat <- dat2 %>%
  group_by(MAG) %>%
  summarise(
    rex_present = max(rex_present),
    CE_potential = dplyr::first(CE_potential),
    rex_strength_sum = sum(rex_strength),
    rex_strength_mean = mean(rex_strength),
    n_motif = n(),
    n_gene = n_distinct(sequence_name),
    .groups="drop"
  )

### ====== 机制模型 ======
### Logistic回归分析, Rex 是否影响 High CE？ 
#coef > 0 → Rex促进CE
#p < 0.05 → 显著
mag_stat <- mag_stat %>%
  mutate(
    CE_high = case_when(
      CE_potential == "High" ~ 1,
      CE_potential == "Medium" ~ 1,
      CE_potential == "Low" ~ 0
    )
  )

model1 <- glm(CE_high ~ rex_strength_sum + n_motif,
              data=mag_stat,
              family=binomial)

summary(model1)



# 控制module
# EtfAB-associated Rex regulation is significantly correlated with high CE potential?
module_stat <- dat2 %>%
  group_by(MAG, module) %>%
  summarise(
    strength = sum(rex_strength),
    .groups="drop"
  ) %>%
  pivot_wider(names_from=module, values_from=strength, values_fill=0)

mag_model <- mag_stat %>%
  left_join(module_stat, by="MAG")

model2 <- glm(CE_high ~ ACAT + EtfAB + BCD + HAD,
              data=mag_model,
              family=binomial)

summary(model2)


#=== 数据处理
# OR （Odds Ratio，优势比）本身是指数尺度：OR=e∧β，forest plot 通常用 log scale
# 某变量增加后，事件发生 odds 的变化倍数
# OR > 1（正关联），OR < 1（负关联）
# Odds ratios (ORs) and corresponding 95% confidence intervals (CIs) were estimated using logistic regression models.



coef_df0 <- tidy(model2) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    
    OR = exp(estimate),
    
    lower = exp(estimate - 1.96 * std.error),
    
    upper = exp(estimate + 1.96 * std.error),
    
    p_label = case_when(
      p.value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p.value)
    ),
    
    sig = p.value < 0.05,
    
    direction = ifelse(OR > 1,
                       "Positive",
                       "Negative")
  )

coef_df0 <- coef_df0 %>%
  arrange(OR) %>%
  mutate(
    term = factor(term,
                  levels = term)
  )

coef_df0 <- coef_df0 %>%
  mutate(
    OR_CI = sprintf(
      "%.2f (%.2f–%.2f)",
      OR, lower, upper
    )
  )

coef_df0 <- coef_df0 %>%
  mutate(
    fontface = ifelse(sig,
                      "bold",
                      "plain")
  )


### plot
forest_core <- ggplot(
  coef_df0,
  aes(x = term,
      y = OR)
) +
  
  geom_vline(
    xintercept = NULL
  ) +
  
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.7,
    color = "grey50"
  ) +
  
  geom_errorbar(
    aes(ymin = lower,
        ymax = upper,
        color = direction),
    width = 0.15,
    linewidth = 1.2
  ) +
  
  geom_point(
    aes(fill = direction),
    shape = 21,
    size = 4.5,
    stroke = 1.2,
    color = "black"
  ) +
  
  scale_y_log10() +
  
  scale_fill_manual(values = c(
    "Positive" = "#CD5B45",
    "Negative" = "#00688B"
  )) +
  
  scale_color_manual(values = c(
    "Positive" = "#CD5B45",
    "Negative" = "#00688B"
  )) +
  
  coord_flip() +
  
  labs(
    x = NULL,
    y = "Odds ratio (95% CI)"
  ) +
  
  theme_bw(base_size = 14) +
  
  theme(
    
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(
      face = "bold",
      color = "black"
    ),
    
    legend.position = "none",
    
    panel.border = element_rect(
      linewidth = 1
    )
  )

forest_core

forest_core <- forest_core +
  
  geom_text(
    aes(
      y = max(upper) * 1.1,
      label = OR_CI
    ),
    hjust = 0,
    size = 4
  )


forest_core <- forest_core +
  
  geom_text(
    aes(
      y = max(upper) * 3,
      label = paste0("P=", p_label),
      fontface = fontface
    ),
    hjust = 0,
    size = 4
  )

forest_core <- forest_core + expand_limits(y = max(coef_df0$upper) * 4)
forest_core

ggsave("Rex_Fig04_OR_forest2.pdf", forest_core, width=5, height=5)

#Coefficients:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)  0.169474   0.175566   0.965 0.334396    
#ACAT         0.064970   0.042240   1.538 0.124021    
#EtfAB       -0.053275   0.053809  -0.990 0.322135    
#BCD          0.306568   0.082415   3.720 0.000199 ***
#  HAD          0.003724   0.054490   0.068 0.945510  

# Rex regulatory strength associated with BCD was significantly positively 
# correlated with high CE potential (OR = 1.36, p = 0.035).
# Rex likely regulates CE through redox-linked electron transfer pathways.

# Module-specific regression analysis revealed that Rex regulation of EtfAB plays 
# a dominant role in determining carbon chain elongation potential, whereas ACAT shows a marginal contribution and BCD appears to be non-determinant.

# These results suggest that redox-linked regulation of electron transfer processes, 
# rather than core carbon condensation steps alone, is a key determinant of CE potential.


p3 <- ggplot(mag_stat,
             aes(x=rex_strength_sum,
                 y=as.numeric(CE_potential=="High"))) +
  
  geom_jitter(height=0.05, width=0, alpha=0.6) +
  
  geom_smooth(method="glm",
              method.args=list(family="binomial"),
              color="red",
              size=1) +
  
  labs(x="Rex regulatory strength (sum of weighted motif scores)",
       y="Probability of High CE") +
  
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face="bold")
  )
p3

ggsave("Fig03_model_curve.pdf", p3, width=5, height=4)

model1 <- glm(CE_high ~ rex_strength_sum,
              data=mag_stat,
              family=binomial)

summary(model1)


#============================== Figure 04 ======================================
#Figure 4 presents the odds ratios (OR) derived from logistic regression models assessing the contribution of different functional modules to CE potential. 
#The x-axis lists the modules, while the y-axis shows the estimated odds ratios with 95% confidence intervals.

#Modules with OR > 1 indicate a positive association with high CE potential, whereas OR < 1 suggests a negative association. 
#Statistically significant modules (p < 0.05) are marked with asterisks.

#Notably, the EtfAB-associated module exhibited the strongest positive association, highlighting its potential role in Rex-mediated regulation of carbon chain elongation.

library(MASS)

full_model <- glm(CE_high ~ ACAT + EtfAB + BCD + HAD,
                  family = binomial,
                  data = mag_model)


step_model <- stepAIC(full_model,direction = "both",trace = FALSE)

summary(step_model)
#Coefficients:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)  0.10814    0.14852   0.728 0.466556    
#ACAT         0.06985    0.04196   1.665 0.095993 .  
#BCD          0.31723    0.08225   3.857 0.000115 ***

exp(0.31723)
# 1.373318
# Rex regulatory strength associated with EtfAB was significantly positively correlated with CE potential (OR = 1.19, p = 0.0471).
# 逐步回归分析筛选出EtfAB和ACAT为与碳链延长潜力相关的关键模块。其中，EtfAB表现出显著正相关，而ACAT呈现边缘显著。
#结果表明，Rex调控主要通过与氧化还原相关的电子传递过程发挥作用，同时核心碳代谢过程可能具有辅助贡献。

coef_df <- tidy(step_model) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    OR = exp(estimate),
    lower = exp(estimate - 1.96*std.error),
    upper = exp(estimate + 1.96*std.error),
    sig = ifelse(p.value < 0.05, "*", "")
  )

p4 <- ggplot(coef_df,
             aes(x=term, y=OR)) +
  
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  geom_text(aes(label=sig, y=upper + 0.2), size=5) +
  
  geom_hline(yintercept=1, linetype="dashed") +
  
  coord_flip() +
  
  labs(x="Module",
       y="Odds ratio (logistic regression)") +
  
  theme_bw(base_size = 14)

p4

ggsave("Fig04_OR_forest.pdf", p4, width=3, height=5)

boot_fn <- function(data, index){
  fit <- glm(CE_high ~ ACAT + EtfAB,
             data=data[index,],
             family=binomial)
  coef(fit)
}

#============== 稳定性评估
### 置信区间
confint(step_model)
#Waiting for profiling to be done...
#2.5 %     97.5 %
#2.5 %    97.5 %
#  (Intercept) -0.18337193 0.3995110
# ACAT        -0.01026129 0.1546169
# BCD          0.16894906 0.4936321

### 模型预测能力
library(pROC)

prob <- predict(step_model, type="response")

roc_obj <- roc(mag_model$CE_high, prob)

plot(roc_obj)
auc(roc_obj)
# Area under the curve: 0.6576

#尽管EtfAB和ACAT在logistic回归模型中表现为正向系数，但其置信区间均跨越零，表明统计支持有限。模型预测能力较弱（AUC = 0.64），
#提示Rex调控可能对碳链延长潜力具有一定贡献，但并非决定性因素。
#该结果表明，Rex更可能作为一种情境依赖的调控因子，在基因表达调控中发挥精细调节作用，而非直接驱动碳链延长过程。



#============================== Figure 05 ======================================
#Figure 5 compares Rex regulatory strength across genomes with different CE potentials. 
#The x-axis represents CE potential categories, while the y-axis shows the aggregated Rex regulatory strength per genome.

#A significant difference in regulatory strength was observed among CE groups (Kruskal–Wallis test, p < 0.05). 
#Pairwise comparisons further revealed that genomes with high CE potential exhibited significantly higher Rex regulatory strength compared to medium and low groups.

#These results suggest that stronger Rex-mediated transcriptional regulation is associated with enhanced carbon chain elongation capacity, supporting a regulatory role of Rex in this metabolic process.


# 先做检验
stat_test <- mag_stat %>%
  kruskal_test(rex_strength_sum ~ CE_potential)

pairwise_test <- mag_stat %>%
  pairwise_wilcox_test(rex_strength_sum ~ CE_potential,
                       p.adjust.method="BH")

# 绘图
mag_stat_fc <- factor(mag_stat$CE_potential, levels = c("High", "Medium", "Low"))
p5 <- ggplot(mag_stat,
             aes(x=CE_potential, y=rex_strength_sum,)) +
  
  geom_boxplot(outlier.shape=NA, width=0.6, color = "gray50") +
  geom_jitter(aes(
    colour = CE_potential),
    width = 0.15, size = 1.5, alpha=0.6)+
  
  scale_colour_manual(
    name = "CE-Potential",
    breaks = c("High","Medium", "Low"),
    values = c(
      "High"="#4575B4",
      "Medium"="#91BFDB",
      "Low"="grey75"
    ))+
  
  scale_x_discrete(limits = c("High", "Medium", "Low"),
                   labels = c("High" = "High\n(n=89)",
                              "Medium" = "Medium\n(n=137)",
                              "Low" = "Low\n(n=135)"))+
  
  labs(x="CE-potential",
       y="Rex regulatory strength") +
  
  stat_compare_means(method="kruskal.test",
                     label.y=max(mag_stat$rex_strength_sum)*1.1) +
  
  stat_compare_means(comparisons=list(c("High","Medium"),
                                      c("High","Low"),
                                      c("Medium","Low")),
                     method="wilcox.test",
                     label="p.signif") +
  
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face="bold"),
    legend.position = ""
  )

p5

ggsave("Rex_Fig05_strength_vs_CE.pdf", p5, width=4, height=5)











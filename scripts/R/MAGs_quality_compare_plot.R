

library(ggplot2)
library(dplyr)
library(ggpubr)

setwd("/Users/hxl/Desktop/zhongke/Liu_metagenomic/Analysis/03_Bining/")

df <- read.delim("MAGs_score_public_in-house.txt", check.names = F, header = T)
dim(df)

df1 <- df %>% filter(Quality != "High")
dim(df1)

ggplot(df1, aes(x=score, fill=Group)) +
  geom_density(alpha=0.4) +
  theme_bw() +
  labs(x="MAG quality score",
       y="Density",
       title="Distribution of MAG quality score") +
  theme(
    text = element_text(size=14)
  )

ggsave("MAG_quality_score_compare.pdf", width = 6, height = 5)
ggsave("MAG_quality_score_compare.png", width = 6, height = 5, dpi = 300)

ggplot(df1, aes(x=Group, y=score, fill=Group)) +
  geom_violin(trim=FALSE, alpha=0.6) +
  geom_boxplot(width=0.15, fill="white") +
  theme_bw() +
  stat_compare_means(method="wilcox.test") +
  labs(y="MAG quality score")











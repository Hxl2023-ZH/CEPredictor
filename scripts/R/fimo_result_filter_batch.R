#################################################################################
#                              FIMO result filter
#################################################################################
# Xinliang Hu
# 20250427
# R 4.5.2

# FIMO hits were filtered using a q-value threshold of 0.001. 
# Overlapping hits within the same promoter region were merged, 
# retaining only the highest-scoring match to avoid redundancy. 
# Genes containing at least one high-confidence Rex-binding site were considered potentially regulated.

# 同一 sequence_name 出现多个非重叠 motif 是真实的多位点调控信号

library(dplyr)
library(GenomicRanges)

# 输入输出路径
input_dir <- "~/Desktop/zhongke/Liu_metagenomic/Analysis/09_Rex_TF/fimo/original_fimo_out/local_MAGs_rBOX_fimo_out/"
output_dir <- "~/Desktop/zhongke/Liu_metagenomic/Analysis/09_Rex_TF/fimo/fimo_filter/local"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

file_list <- list.files(input_dir, pattern = "\\.txt$", full.names = TRUE)

# ========== 核心处理函数 ==========
# 返回值：成功时返回结果数据框，跳过时返回 NULL，读取失败返回特殊标记（此处用 NA）
process_fimo <- function(file_path) {
  
  base_name <- tools::file_path_sans_ext(basename(file_path))
  cat("正在处理:", basename(file_path), "\n")
  
  # 预检文件大小
  if (file.size(file_path) == 0) {
    warning("文件 ", basename(file_path), " 大小为 0，已跳过。")
    return(NA)   # 用 NA 表示跳过（非成功、非 NULL）
  }
  
  # 安全读取
  fimo <- tryCatch(
    read.table(file_path, header = TRUE, sep = "\t", comment.char = "#"),
    error = function(e) {
      warning("读取失败: ", basename(file_path), " (", e$message, ")")
      return(NA)
    }
  )
  if (is.na(fimo)[1] && length(fimo) == 1) return(NA)  # 读取失败
  
  # 检查空数据
  if (nrow(fimo) == 0) {
    warning("文件 ", basename(file_path), " 无有效数据行，已跳过。")
    return(NA)
  }
  
  # ---------- 数据清洗 ----------
  fimo <- fimo %>% mutate(gene_id = sub("ID=", "", sequence_name))
  fimo_filt <- fimo %>% filter(q.value < 0.01)
  if (nrow(fimo_filt) == 0) {
    cat("  过滤后无显著 motif，跳过该文件。\n")
    return(NULL)   # 正常完成但无输出，不视为错误
  }
  
  fimo_filt <- fimo_filt %>%
    mutate(
      score   = as.numeric(score),
      q.value = as.numeric(q.value),
      p.value = as.numeric(p.value)
    )
  
  gr <- GRanges(
    seqnames = fimo_filt$gene_id,
    ranges   = IRanges(start = fimo_filt$start, end = fimo_filt$stop)
  )
  clusters <- reduce(gr)
  
  best_hits <- lapply(seq_along(clusters), function(i) {
    idx <- which(overlapsAny(gr, clusters[i]))
    hits <- fimo_filt[idx, , drop = FALSE]
    hits <- hits[order(hits$q.value, -hits$score), ]
    best <- hits[1, ]
    best$cluster_id <- i
    return(best)
  })
  
  motif_clean <- do.call(rbind, best_hits)
  
  # 写入输出文件
  out_file <- file.path(output_dir, paste0(base_name, "_clean.txt"))
  write.table(motif_clean, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  结果已保存至:", out_file, "\n\n")
  
  return(motif_clean)   # 成功，返回结果数据框
}

# ========== 批量处理 ==========
results <- lapply(file_list, process_fimo)

# ========== 收集跳过的文件 ==========
# NA 表示读取失败或无有效数据，NULL 表示过滤后无结果（通常不算错误）
skip_idx <- sapply(results, function(x) is.logical(x) && length(x) == 1 && is.na(x))
skip_files <- file_list[skip_idx]

if (length(skip_files) > 0) {
  writeLines(skip_files, file.path(output_dir, "skipped_files.txt"))
  cat("\n无法读取/无数据的文件记录已保存至:", file.path(output_dir, "skipped_files.txt"), "\n")
} else {
  cat("\n所有文件均正常读取，无跳过文件。\n")
}

# 如果你还需要汇总所有成功的结果，可以这样：
success_results <- results[!skip_idx & !sapply(results, is.null)]
# combined <- bind_rows(success_results)  # 若需要合并








#!/usr/bin/env Rscript
# Compare local cline angle metrics between focal and background loci
# Author: Ryosuke Ito

library(data.table)

set.seed(42)

### Input ###
fg_file <- "mvlmm-query_with_bgchm-0b.csv"
bg_file <- "bgchm_merged.csv"
outfile <- "fg_bg_angle_difference_summary.tsv"

n_boot <- 5000L
metric_cols <- c("v.med", "v.lower", "v.upper")

### Read data ###
fg_df <- fread(fg_file)[
  , .(chr, position = query_pos, v.med, v.lower, v.upper)
][
  , group := "FG"
]

bg_df <- fread(bg_file)[
  , .(chr, position, v.med, v.lower, v.upper)
][
  , group := "BG"
]

df <- rbindlist(list(fg_df, bg_df), use.names = TRUE, fill = TRUE)

### Prepare data ###
# Remove rows with missing values in any angle metric
df <- df[complete.cases(v.med, v.lower, v.upper)]

idx_fg <- which(df$group == "FG")
idx_bg <- which(df$group == "BG")

n_fg <- length(idx_fg)
n_bg <- length(idx_bg)

stopifnot(n_fg > 0)
stopifnot(n_bg >= n_fg)

### Define helper function ###
# Keep FG fixed and repeatedly sample BG to obtain a null distribution
get_diff_distribution <- function(vec_fg, vec_bg, n_fg, n_boot = 5000L) {
  diffs <- numeric(n_boot)

  for (b in seq_len(n_boot)) {
    bg_sub <- sample(vec_bg, n_fg, replace = FALSE)
    diffs[b] <- median(vec_fg) - median(bg_sub)
  }

  diffs
}

### Calculate difference distributions ###
dist_list <- lapply(metric_cols, function(col) {
  vec_fg <- df[[col]][idx_fg]
  vec_bg <- df[[col]][idx_bg]
  get_diff_distribution(vec_fg, vec_bg, n_fg = n_fg, n_boot = n_boot)
})

names(dist_list) <- metric_cols

### Summarize results ###
result_df <- rbindlist(lapply(metric_cols, function(col) {
  d <- dist_list[[col]]
  ci <- quantile(d, c(0.025, 0.975))

  data.table(
    metric = col,
    n = length(d),
    mean_diff = mean(d),
    median_diff = median(d),
    sd_diff = sd(d),
    se_diff = sd(d) / sqrt(length(d)),
    CI_lower = ci[[1]],
    CI_upper = ci[[2]],
    Pr_gt0 = mean(d > 0)
  )
}))

### Save output ###
fwrite(result_df, outfile, sep = "\t")
print(result_df)

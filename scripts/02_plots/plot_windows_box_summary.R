suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2)
})

per_sample <- read_tsv("results/windows_1Mb_summary_per_sample.tsv", show_col_types = FALSE)

dir.create("results/plots", showWarnings = FALSE, recursive = TRUE)

p1 <- ggplot(per_sample, aes(x=type, y=median_perc_mod)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, size = 2) +
  theme_minimal() +
  labs(title="Methylation (perc_mod) per sample (1Mb windows, cov≥5)",
       y="Median perc_mod per sample", x="Group")

ggsave("results/plots/BOX_median_per_sample.png", p1, width=7, height=5)

# También una tabla ordenada para comentar
per_sample %>%
  arrange(type, desc(median_perc_mod)) %>%
  write_tsv("results/windows_1Mb_summary_per_sample.sorted.tsv")

cat("OK: results/plots/BOX_median_per_sample.png\n")


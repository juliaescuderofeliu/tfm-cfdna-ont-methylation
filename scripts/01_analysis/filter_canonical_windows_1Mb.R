library(dplyr)
library(readr)

df <- read_tsv("results/windows_1Mb_diff_CRC_vs_CONTROL.tsv")

canon <- paste0("chr", c(1:22, "X", "Y"))

df_can <- df %>% filter(chr %in% canon)

write_tsv(df_can, "results/windows_1Mb_diff_CRC_vs_CONTROL_canonical.tsv")

cat("DONE: results/windows_1Mb_diff_CRC_vs_CONTROL_canonical.tsv\n")
cat("N windows:", nrow(df_can), "\n")

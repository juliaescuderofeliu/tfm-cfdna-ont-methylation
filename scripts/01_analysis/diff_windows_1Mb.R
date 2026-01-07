library(dplyr)
library(readr)

infile <- "results/windows_1Mb_all_samples_labeled.tsv"
outfile <- "results/windows_1Mb_diff_CRC_vs_CONTROL.tsv"

df <- read_tsv(infile, show_col_types = FALSE)

diff <- df %>%
  group_by(chr, win_start, win_end) %>%
  summarise(
    n_crc = sum(type == "CRC"),
    n_ctrl = sum(type == "CONTROL"),
    mean_crc = mean(perc_mod[type == "CRC"], na.rm = TRUE),
    mean_ctrl = mean(perc_mod[type == "CONTROL"], na.rm = TRUE),
    delta = mean_crc - mean_ctrl,
    p_value = ifelse(n_crc >= 2 & n_ctrl >= 2,
                     wilcox.test(perc_mod[type == "CRC"],
                                 perc_mod[type == "CONTROL"])$p.value,
                     NA_real_),
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
  arrange(p_adj)

write_tsv(diff, outfile)

cat("DONE:", outfile, "\n")


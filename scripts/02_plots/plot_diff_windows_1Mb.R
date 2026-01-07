library(dplyr)
library(readr)
library(ggplot2)

df <- read_tsv("results/windows_1Mb_diff_CRC_vs_CONTROL.tsv", show_col_types = FALSE)

p <- ggplot(df, aes(x = win_start/1e6, y = delta)) +
  geom_point(size = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  facet_wrap(~ chr, scales = "free_x", ncol = 6) +
  theme_bw() +
  labs(
    x = "Position (Mb)",
    y = "Delta methylation (CRC - CONTROL)",
    title = "Differential methylation by 1Mb window (CRC vs CONTROL)"
  )

ggsave("results/plots/DELTA_crc_minus_control_by_window.png",
       p, width = 14, height = 10)

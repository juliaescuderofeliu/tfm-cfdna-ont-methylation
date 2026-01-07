suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2)
})

w <- read_tsv("results/windows_1Mb_all_samples_labeled.tsv", show_col_types = FALSE)

# Promedio por grupo y ventana
g <- w %>%
  group_by(type, chr, win_start, win_end) %>%
  summarise(mean_perc_mod = mean(perc_mod, na.rm=TRUE), .groups="drop")

# Pivot a columnas CONTROL/CRC y delta
d <- g %>%
  tidyr::pivot_wider(names_from = type, values_from = mean_perc_mod) %>%
  mutate(delta_crc_minus_control = CRC - CONTROL) %>%
  arrange(chr, win_start)

dir.create("results/plots", showWarnings = FALSE, recursive = TRUE)

# Orden cromosómico básico (chr1..chr22, chrX, chrY) ignorando contigs raros si aparecen
keep <- paste0("chr", c(1:22,"X","Y"))
d2 <- d %>% filter(chr %in% keep) %>%
  mutate(chr = factor(chr, levels = keep),
         mid = (win_start + win_end)/2)

p <- ggplot(d2, aes(x=mid/1e6, y=delta_crc_minus_control)) +
  geom_point(size=0.6) +
  facet_wrap(~chr, scales="free_x", ncol=6) +
  theme_minimal() +
  labs(title="Delta methylation by 1Mb window (CRC - CONTROL), cov≥5",
       x="Position (Mb)", y="Δ perc_mod")

ggsave("results/plots/DELTA_crc_minus_control_by_window.png", p, width=14, height=8)

write_tsv(d2, "results/windows_1Mb_delta_crc_minus_control.tsv")
cat("OK: results/plots/DELTA_crc_minus_control_by_window.png\n")

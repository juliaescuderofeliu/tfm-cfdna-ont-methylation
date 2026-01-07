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

# Crear eje genómico continuo (solo cromosomas canónicos)
chr_offsets <- d2 %>%
  group_by(chr) %>%
  summarise(chr_max = max(mid), .groups = "drop") %>%
  mutate(offset = lag(cumsum(chr_max), default = 0))

d3 <- d2 %>%
  left_join(chr_offsets %>% select(chr, offset), by = "chr") %>%
  mutate(genomic_pos = mid + offset)

p <- ggplot(d3, aes(x = genomic_pos/1e6, y = delta_crc_minus_control)) +
  geom_point(size = 0.5, alpha = 0.35) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(
    title = "Genome-wide differential methylation across 1 Mb windows (CRC − CONTROL)",
    x = "Genomic position (hg38, canonical chromosomes)",
    y = "Delta methylation (CRC − CONTROL)"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

ggsave("results/plots/DELTA_crc_minus_control_by_window_canonical.png",
       p, width = 14, height = 5)

write_tsv(d2, "results/windows_1Mb_delta_crc_minus_control.tsv")
cat("OK: results/plots/DELTA_crc_minus_control_by_window.png\n")

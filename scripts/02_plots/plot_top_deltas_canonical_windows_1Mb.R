#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

infile <- "results/windows_1Mb_diff_CRC_vs_CONTROL_canonical.tsv"
outfile <- "results/plots/TOP_delta_windows1Mb_canonical.png"

if (!file.exists(infile)) {
  stop("No encuentro el archivo: ", infile)
}

df <- read_tsv(infile, show_col_types = FALSE)

# seguridad: columnas mínimas
needed <- c("chr","win_start","win_end","delta")
miss <- setdiff(needed, colnames(df))
if (length(miss) > 0) {
  stop("Faltan columnas en el TSV: ", paste(miss, collapse=", "))
}

df <- df %>%
  mutate(
    window = paste0(chr, ":", win_start, "-", win_end)
  ) %>%
  filter(!is.na(delta))

# top 15 hiper + top 15 hipo (según delta)
top_pos <- df %>% arrange(desc(delta)) %>% slice_head(n = 15)
top_neg <- df %>% arrange(delta) %>% slice_head(n = 15)

top <- bind_rows(top_pos, top_neg) %>%
  mutate(window = factor(window, levels = window[order(delta)]),
         direction = ifelse(delta > 0, "CRC > CONTROL", "CRC < CONTROL"))

p <- ggplot(top, aes(x = delta, y = window, fill = direction)) +
  geom_col() +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(
    title = "Top 1Mb windows by delta methylation (CRC - CONTROL) [canonical]",
    x = "Delta methylation (%)",
    y = "Window"
  ) +
  theme_bw()

ggsave(outfile, p, width = 8, height = 10, dpi = 200)

cat("OK. Plot guardado en:", outfile, "\n")


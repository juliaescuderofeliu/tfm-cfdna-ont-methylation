suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

infile  <- "results/windows_1Mb_all_samples_long.tsv"
metafile <- "metadata/samples.tsv"

stopifnot(file.exists(infile))
stopifnot(file.exists(metafile))

w <- read_tsv(infile, show_col_types = FALSE)
meta <- read_tsv(metafile, show_col_types = FALSE) %>%
  transmute(sample = barcode, type = type, notes = notes)

# Join
w2 <- w %>%
  left_join(meta, by = "sample")

# Chequeo: que no falte ninguna etiqueta
missing <- w2 %>% filter(is.na(type)) %>% distinct(sample)
if (nrow(missing) > 0) {
  write_tsv(missing, "results/ERROR_missing_labels.tsv")
  stop("Hay samples sin etiqueta. Mira: results/ERROR_missing_labels.tsv")
}

# Guardar archivo etiquetado
out_labeled <- "results/windows_1Mb_all_samples_labeled.tsv"
write_tsv(w2, out_labeled)

# Resumen por muestra
per_sample <- w2 %>%
  group_by(sample, type) %>%
  summarise(
    n_windows = n(),
    median_perc_mod = median(perc_mod, na.rm = TRUE),
    mean_perc_mod   = mean(perc_mod, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(type, sample)

write_tsv(per_sample, "results/windows_1Mb_summary_per_sample.tsv")

# Resumen por grupo
per_group <- per_sample %>%
  group_by(type) %>%
  summarise(
    n_samples = n(),
    median_of_medians = median(median_perc_mod),
    mean_of_means     = mean(mean_perc_mod),
    .groups = "drop"
  )

write_tsv(per_group, "results/windows_1Mb_summary_by_group.tsv")

cat("OK\n")
cat("Labeled:", out_labeled, "\n")
cat("Per-sample summary: results/windows_1Mb_summary_per_sample.tsv\n")
cat("Group summary:      results/windows_1Mb_summary_by_group.tsv\n")


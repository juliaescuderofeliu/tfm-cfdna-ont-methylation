#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
})

# ========= CONFIG =========
in_dir      <- "data/ch3"
out_dir     <- "results/windows_1Mb"
window_size <- 1e6     # 1 Mb
min_cov     <- 5
min_callprob <- 0.9
min_basequal <- 10
# ==========================

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ch3_files <- list.files(in_dir, pattern = "\\.ch3$", full.names = TRUE)
if (length(ch3_files) == 0) stop("No encuentro .ch3 en ", in_dir)

message("Encontrados ", length(ch3_files), " archivos .ch3")

for (f in ch3_files) {
  sample <- sub("-0\\.ch3$", "", basename(f))
  message("\n=== Procesando: ", sample, " ===")
  message("Archivo: ", f)

  ds <- open_dataset(f, format = "parquet")

  # Filtros mínimos (MUY importantes)
  ds_f <- ds %>%
    filter(ref_position >= 0) %>%                 # posiciones válidas
    filter(call_prob >= min_callprob) %>%
    filter(base_qual >= min_basequal) %>%
    filter(call_code %in% c("m", "h", "-"))       # por seguridad

  # Asignar ventana 1Mb y resumir
  win <- ds_f %>%
    mutate(
      chr = chrom,
      win_start = as.integer(floor(ref_position / window_size) * window_size),
      win_end   = win_start + as.integer(window_size)
    ) %>%
    group_by(chr, win_start, win_end) %>%
    summarise(
      cov   = n(),
      m     = sum(call_code == "m", na.rm = TRUE),
      h     = sum(call_code == "h", na.rm = TRUE),
      unmod = sum(call_code == "-", na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      mod = m + h,
      perc_m   = 100 * m / cov,
      perc_h   = 100 * h / cov,
      perc_mod = 100 * mod / cov
    ) %>%
    filter(cov >= min_cov) %>%
    arrange(chr, win_start)

  out_tsv <- file.path(out_dir, paste0(sample, ".windows1Mb.cov", min_cov, ".tsv"))
  message("Escribiendo: ", out_tsv)

  # Escritura (sin cargar todo lo bruto en memoria fuera de lo necesario)
  # collect() aquí es razonable porque ya estamos agregando por ventanas.
  win_df <- collect(win)
  write.table(win_df, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

  message("OK: ", nrow(win_df), " ventanas")
}

message("\nDONE.")

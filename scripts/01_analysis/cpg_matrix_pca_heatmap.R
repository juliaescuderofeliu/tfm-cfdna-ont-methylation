#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# ---------- CONFIG ----------
overlap_file <- "results/cpg_overlap_cov5_n12.tsv"   # o: results/cpg_overlap_cov5_ge10.tsv
ch3_dir      <- "data/ch3"
out_dir      <- "results/cpg_overlap"
plots_dir    <- "results/plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

controls <- c("barcode75", "barcode76")

# ---------- 1) Leer lista CpGs overlap ----------
cat("Leyendo overlap:", overlap_file, "\n")
ov <- read.delim(overlap_file, sep = "\t", stringsAsFactors = FALSE)
stopifnot(all(c("chr","pos") %in% colnames(ov)))
cat("CpGs en overlap:", nrow(ov), "\n")

ov_split <- split(ov$pos, ov$chr)

# ---------- 2) Detectar muestras disponibles ----------
ch3_files <- list.files(ch3_dir, pattern = "\\.ch3$", full.names = TRUE)
stopifnot(length(ch3_files) > 0)

samples <- sub("-0\\.ch3$", "", basename(ch3_files))
cat("Muestras detectadas:", length(samples), "\n")
cat("CRC:", length(setdiff(samples, controls)), " CONTROL:", length(intersect(samples, controls)), "\n")

# ---------- Helper: %mod por CpG ----------
calc_perc_mod_for_sample <- function(ch3_path, ov_split) {
  ds <- open_dataset(ch3_path, format="parquet")

  ds2 <- ds %>%
    filter(ref_position >= 0) %>%
    select(chrom, ref_position, call_code)

  res_list <- list()

  for (chr in names(ov_split)) {
    pos_vec <- ov_split[[chr]]

    sub <- ds2 %>%
      filter(chrom == chr, ref_position %in% pos_vec) %>%
      group_by(chrom, ref_position) %>%
      summarise(
        m = sum(call_code == "m", na.rm=TRUE),
        h = sum(call_code == "h", na.rm=TRUE),
        unmod = sum(call_code == "-", na.rm=TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        cov = m + h + unmod,
        perc_mod = ifelse(cov > 0, 100*(m+h)/cov, NA_real_)
      )

    res_list[[chr]] <- sub %>% collect()
  }

  bind_rows(res_list)
}

# ---------- 3) Construir long + wide ----------
all_tbl <- list()

for (i in seq_along(ch3_files)) {
  f <- ch3_files[[i]]
  s <- samples[[i]]
  cat("\nProcesando", s, "->", f, "\n")

  tmp <- calc_perc_mod_for_sample(f, ov_split) %>%
    transmute(chr = chrom, pos = ref_position, sample = s, perc_mod)

  out_one <- file.path(out_dir, paste0(s, ".cpg_overlap.perc_mod.tsv"))
  write.table(tmp, out_one, sep="\t", row.names=FALSE, quote=FALSE)
  cat("OK:", out_one, " filas:", nrow(tmp), "\n")

  all_tbl[[s]] <- tmp
}

long <- bind_rows(all_tbl)
out_long <- file.path(out_dir, "cpg_overlap_perc_mod_long.tsv")
write.table(long, out_long, sep="\t", row.names=FALSE, quote=FALSE)
cat("\nDONE long:", out_long, "\n")

wide <- long %>%
  mutate(id = paste(chr, pos, sep=":")) %>%
  select(id, sample, perc_mod) %>%
  pivot_wider(names_from = sample, values_from = perc_mod)

out_wide <- file.path(out_dir, "cpg_overlap_perc_mod_wide.tsv")
write.table(wide, out_wide, sep="\t", row.names=FALSE, quote=FALSE)
cat("DONE wide:", out_wide, "\n")

# ---------- 4) PCA (arreglado: quitar CpGs constantes) ----------
mat <- wide %>% select(-id)
mat <- as.data.frame(mat)
rownames(mat) <- wide$id

# quitar CpGs con NA (si los hubiese)
keep_na <- complete.cases(mat)
mat2 <- mat[keep_na, , drop=FALSE]
cat("CpGs sin NA:", nrow(mat2), "\n")

# quitar CpGs con varianza 0 en el conjunto de muestras (columna constante en PCA)
# OJO: PCA lo haremos sobre X = t(mat2), donde columnas = CpGs
sd_cols <- apply(mat2, 1, sd, na.rm=TRUE)  # sd por CpG (fila)
keep_var <- is.finite(sd_cols) & (sd_cols > 0)

cat("CpGs con varianza > 0:", sum(keep_var), "\n")
mat3 <- mat2[keep_var, , drop=FALSE]

X <- t(as.matrix(mat3))  # muestras x CpGs

pca <- prcomp(X, center=TRUE, scale.=TRUE)

scores <- as.data.frame(pca$x[,1:2, drop=FALSE])
scores$sample <- rownames(scores)
scores$type <- ifelse(scores$sample %in% controls, "CONTROL", "CRC")

p_pca <- ggplot(scores, aes(PC1, PC2, label=sample, shape=type)) +
  geom_point(size=3) +
  geom_text(vjust=-0.7, size=3) +
  theme_bw() +
  ggtitle("PCA on overlap CpGs (perc_mod, cov≥5)")

pca_png <- file.path(plots_dir, "PCA_cpg_overlap_cov5.png")
ggsave(pca_png, p_pca, width=8, height=6)
cat("OK PCA:", pca_png, "\n")

# ---------- 5) Heatmap top-var (300 CpGs) ----------
vars <- apply(mat3, 1, var, na.rm=TRUE)
topN <- 300
top_ids <- names(sort(vars, decreasing=TRUE))[1:min(topN, length(vars))]

hm <- mat3[top_ids, , drop=FALSE]

# z-score por CpG (maneja sd=0 ya filtrado, pero por seguridad)
hm_z <- t(scale(t(hm)))
hm_z[!is.finite(hm_z)] <- 0

hm_long <- as.data.frame(hm_z)
hm_long$cpg <- rownames(hm_long)
hm_long <- tidyr::pivot_longer(hm_long, -cpg, names_to="sample", values_to="z")

p_hm <- ggplot(hm_long, aes(sample, cpg, fill=z)) +
  geom_tile() +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  ggtitle(paste0("Heatmap top-var CpGs (n=", nrow(hm_z), ")"))

hm_png <- file.path(plots_dir, "HEATMAP_cpg_overlap_topvar300.png")
ggsave(hm_png, p_hm, width=8, height=10)
cat("OK Heatmap:", hm_png, "\n")

cat("\nALL DONE.\n")

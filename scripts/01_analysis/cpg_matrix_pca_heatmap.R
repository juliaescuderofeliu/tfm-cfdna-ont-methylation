#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

cat("== SCRIPT START ==\n")

# ---------- CONFIG ----------
overlap_file <- "results/cpg_overlap_cov5_n12.tsv"   # o: results/cpg_overlap_cov5_ge10.tsv
ch3_dir      <- "data/ch3"
out_dir      <- "results/cpg_overlap"
plots_dir    <- "results/plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

controls <- c("barcode75", "barcode76")

# ---------- 1) Leer lista CpGs overlap ----------
cat("Leyendo overlap: ", overlap_file, "\n")
ov <- read.delim(overlap_file, sep = "\t", stringsAsFactors = FALSE)
stopifnot(all(c("chr","pos") %in% colnames(ov)))

library(bit64)
ov$pos <- bit64::as.integer64(ov$pos)

cat("CpGs en overlap: ", nrow(ov), "\n")
ov_split <- split(ov$pos, ov$chr)


# ---------- 2) Detectar muestras disponibles ----------
cat("Buscando .ch3 en: ", ch3_dir, "\n")
ch3_files <- list.files(ch3_dir, pattern = "\\.ch3$", full.names = TRUE)
stopifnot(length(ch3_files) > 0)

samples <- sub("-0\\.ch3$", "", basename(ch3_files))
cat("Muestras detectadas: ", length(samples), "\n")
cat("CRC: ", length(setdiff(samples, controls)), " CONTROL: ", length(intersect(samples, controls)), "\n")

# ---------- Helper: %mod por CpG ----------
calc_perc_mod_for_sample <- function(ch3_path, ov_split) {
  cat("  open_dataset(): ", basename(ch3_path), "\n")
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
  cat("\nProcesando ", s, " -> ", basename(f), "\n")
  
  tmp <- calc_perc_mod_for_sample(f, ov_split) %>%
    transmute(chr = chrom, pos = ref_position, sample = s, perc_mod)
  
  out_one <- file.path(out_dir, paste0(s, ".cpg_overlap.perc_mod.tsv"))
  write.table(tmp, out_one, sep="\t", row.names=FALSE, quote=FALSE)
  cat("OK: ", out_one, " filas: ", nrow(tmp), "\n")
  
  all_tbl[[s]] <- tmp
}

long <- bind_rows(all_tbl)
out_long <- file.path(out_dir, "cpg_overlap_perc_mod_long.tsv")
write.table(long, out_long, sep="\t", row.names=FALSE, quote=FALSE)
cat("\nDONE long: ", out_long, "\n")

wide <- long %>%
  mutate(id = paste(chr, pos, sep=":")) %>%
  select(id, sample, perc_mod) %>%
  pivot_wider(names_from = sample, values_from = perc_mod)

out_wide <- file.path(out_dir, "cpg_overlap_perc_mod_wide.tsv")
write.table(wide, out_wide, sep="\t", row.names=FALSE, quote=FALSE)
cat("DONE wide: ", out_wide, "\n")

# ---------- 4) PCA ----------
mat <- wide %>% select(-id)
mat <- as.data.frame(mat)
rownames(mat) <- wide$id

keep_na <- complete.cases(mat)
mat2 <- mat[keep_na, , drop=FALSE]
cat("CpGs sin NA: ", nrow(mat2), "\n")

sd_cols <- apply(mat2, 1, sd, na.rm=TRUE)
keep_var <- is.finite(sd_cols) & (sd_cols > 0)
cat("CpGs con varianza > 0: ", sum(keep_var), "\n")
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
cat("OK PCA: ", pca_png, "\n")

# ---------- 5) Heatmap top-var (300 CpGs) ----------
vars <- apply(mat3, 1, var, na.rm=TRUE)
topN <- 300
top_ids <- names(sort(vars, decreasing=TRUE))[1:min(topN, length(vars))]

hm <- mat3[top_ids, , drop=FALSE]
hm_z <- t(scale(t(hm)))
hm_z[!is.finite(hm_z)] <- 0


# - escala tricolor
# - anotación CRC vs CONTROL
ann <- data.frame(
  type = ifelse(colnames(hm_z) %in% controls, "CONTROL", "CRC")
)
rownames(ann) <- colnames(hm_z)

# paleta tricolor (azul-blanco-rojo, muy típica en heatmaps)
pal <- colorRampPalette(c("#2c7bb6", "#ffffbf", "#d7191c"))(101)

hm_png <- file.path(plots_dir, "HEATMAP_cpg_overlap_topvar300.png")
png(hm_png, width=900, height=1200, res=150)
pheatmap(
  hm_z,
  color = pal,
  annotation_col = ann,
  show_rownames = FALSE,
  main = paste0("Top-variable overlap CpGs (n=", nrow(hm_z), "), z-score")
)
dev.off()

cat("OK Heatmap: ", hm_png, "\n")
cat("\nALL DONE.\n")


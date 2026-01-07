suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

infile <- "results/windows_1Mb_all_samples_labeled.tsv"
stopifnot(file.exists(infile))

w <- read_tsv(infile, show_col_types = FALSE)

# Matriz: filas = ventanas, columnas = muestras
# Key ventana = chr:win_start-win_end
mat <- w %>%
  mutate(window_id = paste(chr, win_start, win_end, sep=":")) %>%
  select(window_id, sample, perc_mod, type) %>%
  select(-type) %>%
  pivot_wider(names_from = sample, values_from = perc_mod)

# Quitar ventanas con NA en cualquier muestra (PCA simple)
m <- mat %>% drop_na()
X <- as.matrix(m %>% select(-window_id))
X <- scale(X)  # estandariza por muestra

p <- prcomp(t(X), center = TRUE, scale. = FALSE) # PCA de muestras (transpuesta)

scores <- as.data.frame(p$x[,1:2])
scores$sample <- rownames(scores)

meta <- w %>% distinct(sample, type)
scores <- scores %>% left_join(meta, by="sample")

dir.create("results/plots", showWarnings = FALSE, recursive = TRUE)

ggplot(scores, aes(PC1, PC2, label=sample, shape=type)) +
  geom_point(size=3) +
  geom_text(vjust=-0.8, size=3) +
  theme_minimal()

ggsave("results/plots/PCA_windows1Mb_cov5.png", width=7, height=5)
write_tsv(scores, "results/PCA_scores_PC1_PC2.tsv")

cat("OK PCA: results/plots/PCA_windows1Mb_cov5.png\n")

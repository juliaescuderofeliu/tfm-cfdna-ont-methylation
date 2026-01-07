library(dplyr)
library(readr)
library(purrr)

indir <- "results/windows_1Mb"
files <- list.files(indir, pattern = "\\.tsv$", full.names = TRUE)

message("Leyendo ", length(files), " archivos")

dfs <- lapply(files, function(f) {
  sample <- sub("\\.windows1Mb\\.cov5\\.tsv$", "", basename(f))
  read_tsv(f, show_col_types = FALSE) %>%
    select(chr, win_start, win_end, perc_mod) %>%
    mutate(sample = sample)
})

all <- bind_rows(dfs)

write_tsv(all, "results/windows_1Mb_all_samples_long.tsv")

message("DONE: results/windows_1Mb_all_samples_long.tsv")



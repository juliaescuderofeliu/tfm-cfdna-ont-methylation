library(arrow)
library(dplyr)
library(readr)
library(stringr)

files <- list.files("data/ch3", pattern="\\.ch3$", full.names=TRUE)

all <- lapply(files, function(f) {
  sample <- str_extract(basename(f), "barcode[0-9]+")
  ds <- open_dataset(f, format="parquet")
  
  ds %>%
    filter(ref_position >= 0, call_code %in% c("m","h")) %>%
    transmute(chr = chrom,
              pos = ref_position,
              sample = sample) %>%
    distinct() %>%
    collect()
})

df <- bind_rows(all)

overlap <- df %>%
  group_by(chr, pos) %>%
  summarise(n_samples = n_distinct(sample), .groups="drop") %>%
  arrange(desc(n_samples))

write_tsv(overlap, "results/cpg_overlap_cov5.tsv")

cat("DONE: results/cpg_overlap_cov5.tsv\n")


library(MethylSeqR)

setwd("~/crc-ont")

ch3_files <- list.files("data/ch3", pattern="\\.ch3$", full.names=TRUE)
ch3_70 <- ch3_files[grep("barcode70", ch3_files)]
stopifnot(length(ch3_70) == 1)

ch3_db <- make_pos_db(
  ch3_files = ch3_70,
  ch3_db    = "results/test70_chr22.duckdb",
  chrs      = "chr22",
  min_call_prob = 0.9,
  min_length    = 100,
  min_base_qual = 10
)

# crea ventanas 1kb (esto te da un output ya útil para reunión)
ch3_db <- summarize_windows(ch3_db, window_size=1000, step_size=1000)

dir.create("results/export_test70_chr22", showWarnings=FALSE, recursive=TRUE)
export_tables(ch3_db, out_dir="results/export_test70_chr22")

cat("OK\n")


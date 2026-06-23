suppressPackageStartupMessages({
  library(data.table)
})

folder_path <- getwd()

chip_file <- list.files(folder_path, pattern = "ChIP\\.xls$", full.names = TRUE)
hic_file  <- list.files(folder_path, pattern = "HiC\\.xls$",  full.names = TRUE)
stopifnot(length(chip_file) == 1, length(hic_file) == 1)

chip_path <- chip_file[1]
hic_path  <- hic_file[1]

message("ChIP: ", basename(chip_path))
message("HiC : ", basename(hic_path))

FDR_CUTOFF <- 0.05   # set to 1 to disable filtering

OVERLAP_MODE <- "any"

# Read tables
chip <- fread(chip_path, sep = "\t", header = TRUE, data.table = TRUE, fill = TRUE)
hic  <- fread(hic_path,  sep = "\t", header = TRUE, data.table = TRUE, fill = TRUE)


# Keep rows
chip <- chip[, .(chr, start, end, pileup)]
hic  <- hic[,  .(BIN1_CHR, BIN1_START, BIN1_END, BIN2_START, BIN2_END, FDR)]


chip[, `:=`(
  start  = as.integer(start),
  end    = as.integer(end),
  pileup = as.numeric(pileup)
)]

hic[, `:=`(
  BIN1_START = as.integer(BIN1_START),
  BIN1_END   = as.integer(BIN1_END),
  BIN2_START = as.integer(BIN2_START),
  BIN2_END   = as.integer(BIN2_END),
  FDR        = as.numeric(FDR)
)]

message("Rows: ChIP=", nrow(chip), "  HiC(before filter)=", nrow(hic))

# Optional Hi-C filter
if (!is.na(FDR_CUTOFF) && is.finite(FDR_CUTOFF) && FDR_CUTOFF < 1) {
  hic <- hic[FDR <= FDR_CUTOFF]
}
message("Rows: HiC(after  filter)=", nrow(hic))


setkey(chip, chr, start, end)

# ---- BIN1 anchor vs ChIP, project to BIN2
bin1 <- hic[, .(
  chr = BIN1_CHR,
  start = BIN1_START,
  end   = BIN1_END,
  BIN2_START, BIN2_END, FDR
)]
setkey(bin1, chr, start, end)

# Process Bin1
ov1 <- foverlaps(bin1, chip, type = "any", nomatch = 0L)

ov1 <- ov1[
  i.start <= end & i.end >= start
]

new_df1 <- ov1[, .(
  chr   = chr,
  start = BIN2_START,
  end   = BIN2_END,
  pileup = pileup,
  FDR   = FDR
)]

# process Bin2
# ---- BIN2 anchor vs ChIP, project to BIN1
bin2 <- hic[, .(
  chr = BIN1_CHR,
  start = BIN2_START,
  end   = BIN2_END,
  BIN1_START, BIN1_END, FDR
)]
setkey(bin2, chr, start, end)

# process Bin2
ov2 <- foverlaps(bin2, chip, type = "any", nomatch = 0L)

ov2 <- ov2[
  i.start <= end & i.end >= start
]

new_df2 <- ov2[, .(
  chr   = chr,
  start = BIN1_START,
  end   = BIN1_END,
  pileup = pileup,
  FDR   = FDR
)]

message("Projected rows: new_df1=", nrow(new_df1), "  new_df2=", nrow(new_df2))

# Merge
final_matrix <- rbindlist(list(new_df1, new_df2), use.names = TRUE)

message("Final unique rows: ", nrow(final_matrix))

# Rank 
rng01 <- function(x) {
  r <- range(x, na.rm = TRUE)
  if (isTRUE(all.equal(r[1], r[2]))) return(rep(0, length(x)))
  (x - r[1]) / (r[2] - r[1])
}

final_matrix[, pileup_norm := rng01(pileup)]
final_matrix[, fdr_norm    := rng01(FDR)]
final_matrix[, score       := pileup_norm - fdr_norm]
setorder(final_matrix, -score)

# Write output

out_file <- "final_ranked_output.xls"

write.table(
  final_matrix,
  file = out_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

message("Wrote ", out_file, " with rows: ", nrow(final_matrix))


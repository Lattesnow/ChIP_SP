# Remove chrX / chrY rows from ALL *HiC.xls files

hic_files <- list.files(
  path = getwd(),
  pattern = "HiC\\.xls$",
  full.names = TRUE
)

if (length(hic_files) == 0) {
  stop("No *HiC.xls files found in current directory.")
}

for (f in hic_files) {
  
  cat("Processing:", basename(f), "\n")
  
  hic <- read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
  
  # Detect chromosome columns
  hic_cols <- NULL
  if (all(c("BIN1_CHR", "BIN2_CHR") %in% colnames(hic))) {
    hic_cols <- c("BIN1_CHR", "BIN2_CHR")
  } else if (all(c("BIN1_CHROMOSOME", "BIN2_CHROMOSOME") %in% colnames(hic))) {
    hic_cols <- c("BIN1_CHROMOSOME", "BIN2_CHROMOSOME")
  } else {
    candidates <- c("BIN1_CHR", "BIN2_CHR",
                    "BIN1_CHROMOSOME", "BIN2_CHROMOSOME")
    hic_cols <- intersect(candidates, colnames(hic))
  }
  
  if (length(hic_cols) == 0) {
    warning("  No chromosome columns found — skipped")
    next
  }
  
  xy_regex <- "^(chr)?[XY]$"
  
  keep <- apply(hic[, hic_cols, drop = FALSE], 1, function(x) {
    x <- trimws(as.character(x))
    !any(grepl(xy_regex, x, ignore.case = TRUE) & !is.na(x))
  })
  
  n_removed <- sum(!keep)
  if (n_removed > 0) {
    cat("  Removed", n_removed, "rows containing chrX/chrY\n")
  } else {
    cat("  No chrX/chrY rows found\n")
  }
  
  hic_clean <- hic[keep, , drop = FALSE]
  
  # Overwrite file
  write.table(
    hic_clean,
    f,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

cat("Done cleaning all *HiC.xls files.\n")

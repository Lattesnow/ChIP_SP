suppressPackageStartupMessages({
  library(dplyr)
})

wd <- getwd()

merge_and_cleanup <- function(pattern, output_filename) {
  
  files <- list.files(path = wd, pattern = pattern, full.names = TRUE)
  
  if (length(files) == 0) {
    message("No files found for pattern: ", pattern)
    return(NULL)
  }
  
  message("Merging files:")
  message(basename(files), collapse = "\n")
  
  # Read all files
  dfs <- lapply(files, function(f) {
    read.delim(f, header = TRUE, stringsAsFactors = FALSE)
  })
  
  combined <- bind_rows(dfs)
  
  # Write merged file
  write.table(
    combined,
    file = output_filename,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  
  # delete original files
  if (file.exists(output_filename)) {
    file.remove(files)
    message("Merged file saved as: ", output_filename)
    message("Deleted original files.")
  } else {
    stop("Merged file was not created. Original files NOT deleted.")
  }
  
  invisible(combined)
}

# Merge HiC
merge_and_cleanup(
  pattern = "HiC\\.xls$",
  output_filename = "Combined_HiC.xls"
)

# Merge ChIP
merge_and_cleanup(
  pattern = "ChIP\\.xls$",
  output_filename = "Combined_ChIP.xls"
)

# Fix the chromosome chr missing in the chr column
# ---------- ChIP: chr column ----------
if (file.exists("Combined_ChIP.xls")) {
  chip <- read.delim("Combined_ChIP.xls", stringsAsFactors = FALSE)
  
  if ("chr" %in% colnames(chip)) {
    need_fix <- !is.na(chip$chr) & !grepl("^chr", chip$chr, ignore.case = TRUE)
    
    if (any(need_fix)) {
      chip$chr[need_fix] <- paste0("chr", chip$chr[need_fix])
      message("[ChIP] Added 'chr' prefix to ", sum(need_fix), " rows")
    } else {
      message("[ChIP] chr column already normalized")
    }
    
    write.table(chip, "Combined_ChIP.xls",
                sep = "\t", quote = FALSE, row.names = FALSE)
  } else {
    message("[ChIP] No 'chr' column found — skipped")
  }
}

# ---------- HiC: BIN1/BIN2 chr columns (CHR or CHROMOSOME) ----------
if (file.exists("Combined_HiC.xls")) {
  hic <- read.delim("Combined_HiC.xls", stringsAsFactors = FALSE)
  
  # detect which naming scheme is present
  hic_cols <- NULL
  if (all(c("BIN1_CHR", "BIN2_CHR") %in% colnames(hic))) {
    hic_cols <- c("BIN1_CHR", "BIN2_CHR")
  } else if (all(c("BIN1_CHROMOSOME", "BIN2_CHROMOSOME") %in% colnames(hic))) {
    hic_cols <- c("BIN1_CHROMOSOME", "BIN2_CHROMOSOME")
  } else {

    candidates <- c("BIN1_CHR", "BIN2_CHR", "BIN1_CHROMOSOME", "BIN2_CHROMOSOME")
    hic_cols <- intersect(candidates, colnames(hic))
    message("[HiC] Using available columns: ", paste(hic_cols, collapse = ", "))
  }
  
  for (col in hic_cols) {
    need_fix <- !is.na(hic[[col]]) & !grepl("^chr", hic[[col]], ignore.case = TRUE)
    
    if (any(need_fix)) {
      hic[[col]][need_fix] <- paste0("chr", hic[[col]][need_fix])
      message("[HiC] Added 'chr' prefix in ", col, " for ", sum(need_fix), " rows")
    } else {
      message("[HiC] ", col, " already normalized")
    }
  }
  
  write.table(hic, "Combined_HiC.xls",
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# Caliber the CHROMOSOME to CHR in the colnames of HiC file

if (file.exists("Combined_HiC.xls")) {
  
  hic <- read.delim("Combined_HiC.xls", stringsAsFactors = FALSE)
  
  old_names <- colnames(hic)
  new_names <- old_names
  
  # Force BIN1_CHR* -> BIN1_CHR
  new_names <- ifelse(grepl("^BIN1_CHR", new_names),
                      "BIN1_CHR",
                      new_names)
  
  # Force BIN2_CHR* -> BIN2_CHR
  new_names <- ifelse(grepl("^BIN2_CHR", new_names),
                      "BIN2_CHR",
                      new_names)
  
  # Report changes
  renamed <- old_names != new_names
  if (any(renamed)) {
    message("[HiC] Renamed columns:")
    for (i in which(renamed)) {
      message("  ", old_names[i], "  ->  ", new_names[i])
    }
  } else {
    message("[HiC] No BIN*_CHR column renaming needed")
  }
  
  colnames(hic) <- new_names
  
  write.table(
    hic,
    "Combined_HiC.xls",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

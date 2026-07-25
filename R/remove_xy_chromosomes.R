#' Remove chrX/chrY from Hi-C loop files
#'
#' Reads one or more Hi-C loop files, removes interactions involving
#' chromosome X or chromosome Y, and writes cleaned files in CSV format.
#' Supported input formats include CSV, TSV, TXT, TAB, BED, XLS, and XLSX.
#' Chromosome columns are detected automatically.
#'
#' @param files Character vector of Hi-C file paths.
#'
#' @return Invisibly returns a character vector containing the output file names.
#' @export
#'
#' @examples
#' hic_files <- list.files(
#'   pattern = "HiC\\.(csv|tsv|txt|xls|xlsx)$",
#'   full.names = TRUE
#' )
#'
#' removeXYChromosomes(hic_files)
removeXYChromosomes <- function(files) {

  if (length(files) == 0)
    stop("No Hi-C files supplied.")

  output_files <- character(length(files))

  for (i in seq_along(files)) {

    f <- files[i]

    message("Processing: ", basename(f))

    ext <- tolower(tools::file_ext(f))

    hic <- switch(
      ext,
      csv  = read.csv(f, stringsAsFactors = FALSE, check.names = FALSE),
      tsv  = read.delim(f, stringsAsFactors = FALSE, check.names = FALSE),
      txt  = read.delim(f, stringsAsFactors = FALSE, check.names = FALSE),
      tab  = read.delim(f, stringsAsFactors = FALSE, check.names = FALSE),
      bed  = read.delim(f,
                        stringsAsFactors = FALSE,
                        check.names = FALSE,
                        header = FALSE),
      xlsx = readxl::read_excel(f) |>
        as.data.frame(check.names = FALSE),
      xls = tryCatch(
        readxl::read_excel(f) |>
          as.data.frame(check.names = FALSE),
        error = function(e)
          read.delim(f,
                     stringsAsFactors = FALSE,
                     check.names = FALSE)
      ),
      stop("Unsupported file type: ", ext)
    )

    ## Detect chromosome columns
    if (all(c("BIN1_CHR", "BIN2_CHR") %in% colnames(hic))) {

      hic_cols <- c("BIN1_CHR", "BIN2_CHR")

    } else if (all(c("BIN1_CHROMOSOME", "BIN2_CHROMOSOME") %in% colnames(hic))) {

      hic_cols <- c("BIN1_CHROMOSOME", "BIN2_CHROMOSOME")

    } else {

      hic_cols <- intersect(
        c("BIN1_CHR",
          "BIN2_CHR",
          "BIN1_CHROMOSOME",
          "BIN2_CHROMOSOME"),
        colnames(hic)
      )
    }

    if (length(hic_cols) == 0) {
      warning("No chromosome columns detected in ", basename(f))
      next
    }

    xy_regex <- "^(chr)?[XY]$"

    keep <- apply(
      hic[, hic_cols, drop = FALSE],
      1,
      function(x) {
        x <- trimws(as.character(x))
        !any(grepl(xy_regex,
                   x,
                   ignore.case = TRUE) &
               !is.na(x))
      }
    )

    removed <- sum(!keep)

    message("Removed ", removed, " interactions on chrX/chrY")

    hic_clean <- hic[keep, , drop = FALSE]

    out_file <- file.path(
      dirname(f),
      paste0(
        tools::file_path_sans_ext(basename(f)),
        "_no_chrX_chrY.csv"
      )
    )

    write.csv(
      hic_clean,
      out_file,
      row.names = FALSE,
      quote = FALSE
    )

    message("Saved: ", basename(out_file))

    output_files[i] <- out_file
  }

  invisible(output_files)
}

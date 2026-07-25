#' Merge Hi-C loop outputs across replicates and resolutions
#'
#' @description
#' Reads and merges multiple Hi-C loop output files into a single data frame.
#' The function supports common delimited and Excel input formats, normalizes
#' Hi-C chromosome-column names, and optionally adds the `"chr"` prefix to
#' chromosome values.
#'
#' This function does not write output files or delete input files.
#'
#' @param hic_files Character vector containing paths to Hi-C loop files.
#'   Supported formats are CSV, TSV, TXT, TAB, XLS, and XLSX.
#' @param normalize_column_names Logical. If `TRUE`, chromosome columns named
#'   `"BIN1_CHROMOSOME"` and `"BIN2_CHROMOSOME"` are renamed to
#'   `"BIN1_CHR"` and `"BIN2_CHR"`. Default is `TRUE`.
#' @param add_chr_prefix Logical. If `TRUE`, chromosome values that do not
#'   begin with `"chr"` are given the `"chr"` prefix. Default is `TRUE`.
#'
#' @return A data.frame containing merged and normalized Hi-C loops.
#'
#' @examples
#' hic_files <- list.files(
#'   pattern = "HiC\\.(csv|tsv|txt)$",
#'   full.names = TRUE,
#'   ignore.case = TRUE
#' )
#'
#' hic_df <- mergeHiCLoops(hic_files)
#'
#' @export
mergeHiCLoops <- function(hic_files,
                          normalize_column_names = TRUE,
                          add_chr_prefix = TRUE) {

  if (length(hic_files) == 0) {
    stop("No Hi-C files provided.")
  }

  missing_files <- hic_files[!file.exists(hic_files)]

  if (length(missing_files) > 0) {
    stop(
      "The following Hi-C files do not exist: ",
      paste(missing_files, collapse = ", ")
    )
  }

  read_hic_file <- function(file) {

    extension <- tolower(tools::file_ext(file))

    hic <- switch(
      extension,

      "csv" = read.csv(
        file,
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),

      "tsv" = read.delim(
        file,
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),

      "txt" = read.delim(
        file,
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),

      "tab" = read.delim(
        file,
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),

      "xlsx" = as.data.frame(
        readxl::read_excel(file),
        check.names = FALSE
      ),

      "xls" = tryCatch(
        as.data.frame(
          readxl::read_excel(file),
          check.names = FALSE
        ),
        error = function(e) {
          read.delim(
            file,
            header = TRUE,
            stringsAsFactors = FALSE,
            check.names = FALSE
          )
        }
      ),

      stop(
        "Unsupported Hi-C file type: ",
        extension,
        ". Supported formats are csv, tsv, txt, tab, xls, and xlsx."
      )
    )

    hic
  }

  message("Merging ", length(hic_files), " Hi-C files:")

  for (file in hic_files) {
    message("  ", basename(file))
  }

  hic_list <- lapply(hic_files, read_hic_file)

  hic_df <- dplyr::bind_rows(hic_list)

  if (nrow(hic_df) == 0) {
    warning("The merged Hi-C data frame contains no rows.")
    return(hic_df)
  }

  # Normalize chromosome-column names
  if (normalize_column_names) {

    if ("BIN1_CHROMOSOME" %in% colnames(hic_df) &&
        !"BIN1_CHR" %in% colnames(hic_df)) {
      colnames(hic_df)[
        colnames(hic_df) == "BIN1_CHROMOSOME"
      ] <- "BIN1_CHR"
    }

    if ("BIN2_CHROMOSOME" %in% colnames(hic_df) &&
        !"BIN2_CHR" %in% colnames(hic_df)) {
      colnames(hic_df)[
        colnames(hic_df) == "BIN2_CHROMOSOME"
      ] <- "BIN2_CHR"
    }
  }

  chromosome_columns <- intersect(
    c(
      "BIN1_CHR",
      "BIN2_CHR",
      "BIN1_CHROMOSOME",
      "BIN2_CHROMOSOME"
    ),
    colnames(hic_df)
  )

  if (length(chromosome_columns) == 0) {
    warning(
      "No Hi-C chromosome columns were detected. ",
      "Expected BIN1_CHR/BIN2_CHR or ",
      "BIN1_CHROMOSOME/BIN2_CHROMOSOME."
    )
  }

  # Add chr prefix when missing
  if (add_chr_prefix && length(chromosome_columns) > 0) {

    for (column_name in chromosome_columns) {

      chromosome_values <- trimws(
        as.character(hic_df[[column_name]])
      )

      needs_prefix <- (
        !is.na(chromosome_values) &
        chromosome_values != "" &
        !grepl("^chr", chromosome_values, ignore.case = TRUE)
      )

      if (any(needs_prefix)) {
        chromosome_values[needs_prefix] <- paste0(
          "chr",
          chromosome_values[needs_prefix]
        )

        message(
          "Added 'chr' prefix in ",
          column_name,
          " for ",
          sum(needs_prefix),
          " rows."
        )
      }

      hic_df[[column_name]] <- chromosome_values
    }
  }

  rownames(hic_df) <- NULL

  message(
    "Merged Hi-C data: ",
    nrow(hic_df),
    " rows and ",
    ncol(hic_df),
    " columns."
  )

  hic_df
}

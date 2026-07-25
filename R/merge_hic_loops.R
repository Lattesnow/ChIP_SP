#' Merge Hi-C loop outputs across replicates and resolutions
#'
#' @description
#' Reads and merges multiple Hi-C loop output files into a single data frame.
#' Supported formats include CSV, TSV, TXT, TAB, XLS, and XLSX.
#'
#' Before merging, the function normalizes common Hi-C chromosome-column
#' names, optionally adds the \code{"chr"} prefix to chromosome values,
#' validates required loop-coordinate columns, and records the source file
#' for each row.
#'
#' The function does not write output files or delete input files.
#'
#' @param hic_files Character vector containing paths to Hi-C loop files.
#'   Supported formats are CSV, TSV, TXT, TAB, XLS, and XLSX.
#'
#' @param normalize_column_names Logical. If \code{TRUE},
#'   \code{BIN1_CHROMOSOME} and \code{BIN2_CHROMOSOME} are renamed to
#'   \code{BIN1_CHR} and \code{BIN2_CHR}, respectively. Default is
#'   \code{TRUE}.
#'
#' @param add_chr_prefix Logical. If \code{TRUE}, chromosome values that do
#'   not begin with \code{"chr"} are given the prefix. Default is
#'   \code{TRUE}.
#'
#' @param validate_columns Logical. If \code{TRUE}, the function checks that
#'   each Hi-C file contains the required chromosome and coordinate columns.
#'   Default is \code{TRUE}.
#'
#' @param add_source_file Logical. If \code{TRUE}, a \code{source_file}
#'   column is added to identify the original input file for each loop.
#'   Default is \code{TRUE}.
#'
#' @return A data.frame containing merged and normalized Hi-C loops.
#'
#' @examples
#' hic_files <- list.files(
#'   pattern = "HiC.*\\.(csv|tsv|txt|tab|xls|xlsx)$",
#'   full.names = TRUE,
#'   ignore.case = TRUE
#' )
#'
#' hic_df <- mergeHiCLoops(hic_files)
#'
#' write.csv(
#'   hic_df,
#'   "Combined_HiC.csv",
#'   row.names = FALSE,
#'   quote = FALSE
#' )
#'
#' @export
mergeHiCLoops <- function(hic_files,
                          normalize_column_names = TRUE,
                          add_chr_prefix = TRUE,
                          validate_columns = TRUE,
                          add_source_file = TRUE) {

  # ----------------------------------------------------------
  # 1. Validate input paths
  # ----------------------------------------------------------

  if (length(hic_files) == 0) {
    stop("No Hi-C files provided.")
  }

  if (!is.character(hic_files)) {
    stop("`hic_files` must be a character vector of file paths.")
  }

  hic_files <- unique(hic_files)

  missing_files <- hic_files[!file.exists(hic_files)]

  if (length(missing_files) > 0) {
    stop(
      "The following Hi-C files do not exist:\n",
      paste(missing_files, collapse = "\n")
    )
  }

  # Exclude a previously generated combined output
  hic_files <- hic_files[
    tolower(basename(hic_files)) != "combined_hic.csv"
  ]

  if (length(hic_files) == 0) {
    stop(
      "No valid Hi-C input files remain after excluding Combined_HiC.csv."
    )
  }

  # ----------------------------------------------------------
  # 2. Helper: read one Hi-C file
  # ----------------------------------------------------------

  read_hic_file <- function(file) {

    extension <- tolower(tools::file_ext(file))

    message("Reading: ", basename(file))

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
        {
          as.data.frame(
            readxl::read_excel(file),
            check.names = FALSE
          )
        },
        error = function(excel_error) {

          message(
            "  Excel reading failed; trying tab-delimited text: ",
            basename(file)
          )

          tryCatch(
            {
              read.delim(
                file,
                header = TRUE,
                stringsAsFactors = FALSE,
                check.names = FALSE
              )
            },
            error = function(text_error) {
              stop(
                "Could not read ",
                basename(file),
                " as either an Excel file or tab-delimited text.\n",
                "Excel error: ",
                conditionMessage(excel_error),
                "\nText error: ",
                conditionMessage(text_error)
              )
            }
          )
        }
      ),

      stop(
        "Unsupported Hi-C file type: ",
        extension,
        ". Supported formats are csv, tsv, txt, tab, xls, and xlsx."
      )
    )

    if (!is.data.frame(hic)) {
      hic <- as.data.frame(
        hic,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }

    hic
  }

  # ----------------------------------------------------------
  # 3. Helper: normalize column names
  # ----------------------------------------------------------

  normalize_hic_names <- function(hic) {

    old_names <- colnames(hic)
    new_names <- trimws(old_names)

    # Normalize common chromosome-column variants
    new_names[toupper(new_names) == "BIN1_CHROMOSOME"] <- "BIN1_CHR"
    new_names[toupper(new_names) == "BIN2_CHROMOSOME"] <- "BIN2_CHR"

    new_names[toupper(new_names) == "BIN1_CHR"] <- "BIN1_CHR"
    new_names[toupper(new_names) == "BIN2_CHR"] <- "BIN2_CHR"

    new_names[toupper(new_names) == "BIN1_START"] <- "BIN1_START"
    new_names[toupper(new_names) == "BIN1_END"] <- "BIN1_END"

    new_names[toupper(new_names) == "BIN2_START"] <- "BIN2_START"
    new_names[toupper(new_names) == "BIN2_END"] <- "BIN2_END"

    new_names[toupper(new_names) == "FDR"] <- "FDR"

    colnames(hic) <- new_names

    duplicated_names <- unique(
      colnames(hic)[duplicated(colnames(hic))]
    )

    if (length(duplicated_names) > 0) {
      stop(
        "Column-name normalization created duplicate columns: ",
        paste(duplicated_names, collapse = ", ")
      )
    }

    hic
  }

  # ----------------------------------------------------------
  # 4. Helper: normalize chromosome values
  # ----------------------------------------------------------

  normalize_chr_values <- function(x) {

    x <- trimws(as.character(x))

    missing_value <- is.na(x) | x == ""

    needs_prefix <- (
      !missing_value &
      !grepl("^chr", x, ignore.case = TRUE)
    )

    x[needs_prefix] <- paste0(
      "chr",
      x[needs_prefix]
    )

    # Standardize prefix capitalization
    x[!missing_value] <- sub(
      "^CHR",
      "chr",
      x[!missing_value],
      ignore.case = TRUE
    )

    x
  }

  # ----------------------------------------------------------
  # 5. Read and process each file before merging
  # ----------------------------------------------------------

  required_hic_columns <- c(
    "BIN1_CHR",
    "BIN1_START",
    "BIN1_END",
    "BIN2_CHR",
    "BIN2_START",
    "BIN2_END"
  )

  hic_list <- lapply(hic_files, function(file) {

    hic <- read_hic_file(file)

    if (normalize_column_names) {
      hic <- normalize_hic_names(hic)
    }

    if (validate_columns) {

      missing_columns <- setdiff(
        required_hic_columns,
        colnames(hic)
      )

      if (length(missing_columns) > 0) {
        stop(
          "Hi-C file ",
          basename(file),
          " is missing required columns: ",
          paste(missing_columns, collapse = ", ")
        )
      }
    }

    chromosome_columns <- intersect(
      c("BIN1_CHR", "BIN2_CHR"),
      colnames(hic)
    )

    if (add_chr_prefix && length(chromosome_columns) > 0) {

      for (column_name in chromosome_columns) {

        original_values <- as.character(
          hic[[column_name]]
        )

        hic[[column_name]] <- normalize_chr_values(
          hic[[column_name]]
        )

        n_changed <- sum(
          original_values != hic[[column_name]],
          na.rm = TRUE
        )

        if (n_changed > 0) {
          message(
            "  Added or standardized 'chr' prefix in ",
            column_name,
            " for ",
            n_changed,
            " rows."
          )
        }
      }
    }

    # Convert coordinates and FDR to numeric
    numeric_columns <- intersect(
      c(
        "BIN1_START",
        "BIN1_END",
        "BIN2_START",
        "BIN2_END",
        "FDR"
      ),
      colnames(hic)
    )

    for (column_name in numeric_columns) {

      original_values <- hic[[column_name]]

      converted_values <- suppressWarnings(
        as.numeric(original_values)
      )

      failed_conversion <- (
        !is.na(original_values) &
        trimws(as.character(original_values)) != "" &
        is.na(converted_values)
      )

      if (any(failed_conversion)) {
        warning(
          basename(file),
          ": ",
          sum(failed_conversion),
          " nonnumeric values detected in ",
          column_name,
          ". These values were converted to NA."
        )
      }

      hic[[column_name]] <- converted_values
    }

    # Remove rows with invalid required coordinates
    if (all(required_hic_columns %in% colnames(hic))) {

      valid_rows <- (
        !is.na(hic$BIN1_CHR) &
        hic$BIN1_CHR != "" &
        !is.na(hic$BIN1_START) &
        !is.na(hic$BIN1_END) &
        hic$BIN1_END > hic$BIN1_START &
        !is.na(hic$BIN2_CHR) &
        hic$BIN2_CHR != "" &
        !is.na(hic$BIN2_START) &
        !is.na(hic$BIN2_END) &
        hic$BIN2_END > hic$BIN2_START
      )

      n_invalid <- sum(!valid_rows)

      if (n_invalid > 0) {
        warning(
          basename(file),
          ": removed ",
          n_invalid,
          " rows with missing or invalid loop coordinates."
        )

        hic <- hic[valid_rows, , drop = FALSE]
      }
    }

    if (add_source_file) {
      hic$source_file <- basename(file)
    }

    rownames(hic) <- NULL

    message(
      "  Retained ",
      nrow(hic),
      " rows and ",
      ncol(hic),
      " columns."
    )

    hic
  })

  # ----------------------------------------------------------
  # 6. Merge files
  # ----------------------------------------------------------

  message(
    "Merging ",
    length(hic_list),
    " processed Hi-C files."
  )

  hic_df <- dplyr::bind_rows(hic_list)

  rownames(hic_df) <- NULL

  if (nrow(hic_df) == 0) {
    warning("The merged Hi-C data frame contains no rows.")
    return(hic_df)
  }

  # ----------------------------------------------------------
  # 7. Final validation
  # ----------------------------------------------------------

  if (validate_columns) {

    missing_final_columns <- setdiff(
      required_hic_columns,
      colnames(hic_df)
    )

    if (length(missing_final_columns) > 0) {
      stop(
        "The merged Hi-C table is missing required columns: ",
        paste(missing_final_columns, collapse = ", ")
      )
    }
  }

  message(
    "Merged Hi-C data: ",
    nrow(hic_df),
    " rows and ",
    ncol(hic_df),
    " columns."
  )

  hic_df
}

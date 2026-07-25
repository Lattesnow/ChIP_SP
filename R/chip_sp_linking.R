#' ChIP-SP Core Spatial Integration and Ranking
#'
#' @description
#' Integrates ChIP-seq peaks with Hi-C chromatin loops by identifying
#' overlaps between ChIP-seq peaks and either Hi-C loop anchor. When a
#' ChIP-seq peak overlaps one loop anchor, the partner anchor is reported
#' as the spatially linked ChIP-SP region.
#'
#' Spatially linked regions are ranked using normalized ChIP-seq pileup
#' and Hi-C loop false-discovery rate:
#'
#' \deqn{
#' score = normalized\ pileup - normalized\ FDR
#' }
#'
#' The function uses \code{data.table::foverlaps()} for efficient genomic
#' interval matching and avoids construction of a full Cartesian product.
#'
#' @param chip_file Character scalar. Path to a ChIP-seq peak file.
#'   Supported formats are BED, CSV, TSV, TXT, TAB, XLS, and XLSX.
#'   Delimited files must contain \code{chr}, \code{start}, and \code{end}
#'   columns. A \code{pileup} column is optional; when absent, all peaks
#'   are assigned a pileup value of 1.
#'
#' @param hic_df A data.frame or data.table containing Hi-C loops,
#'   typically returned by \code{mergeHiCLoops()}. Required columns are
#'   \code{BIN1_CHR}, \code{BIN1_START}, \code{BIN1_END},
#'   \code{BIN2_CHR}, \code{BIN2_START}, \code{BIN2_END}, and \code{FDR}.
#'
#' @param fdr_cutoff Numeric scalar specifying the maximum Hi-C loop FDR
#'   retained before overlap analysis. Default is \code{0.05}. Set to
#'   \code{1} or \code{NULL} to disable FDR filtering.
#'
#' @param overlap_mode Character scalar passed to
#'   \code{data.table::foverlaps(type = ...)}. Default is \code{"any"}.
#'
#' @param add_chr_prefix Logical. If \code{TRUE}, chromosome values lacking
#'   the \code{"chr"} prefix are normalized before overlap analysis.
#'   Default is \code{TRUE}.
#'
#' @return A data.frame containing spatially linked genomic regions ranked
#'   by ChIP-SP score. Output columns include:
#'
#' \itemize{
#'   \item \code{chr}, \code{start}, and \code{end}: partner-anchor coordinates;
#'   \item \code{pileup}: ChIP-seq peak signal;
#'   \item \code{FDR}: Hi-C loop confidence;
#'   \item \code{pileup_norm}: min-max normalized pileup;
#'   \item \code{fdr_norm}: min-max normalized FDR;
#'   \item \code{score}: final ChIP-SP ranking score;
#'   \item \code{source_anchor}: loop anchor overlapped by the ChIP-seq peak.
#' }
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
#' result <- chipSPLink(
#'   chip_file = "AR_ChIP.bed",
#'   hic_df = hic_df,
#'   fdr_cutoff = 0.05
#' )
#'
#' head(result)
#'
#' @export
chipSPLink <- function(chip_file,
                       hic_df,
                       fdr_cutoff = 0.05,
                       overlap_mode = "any",
                       add_chr_prefix = TRUE) {

  # ----------------------------------------------------------
  # 1. Validate inputs
  # ----------------------------------------------------------

  if (length(chip_file) != 1 || !is.character(chip_file)) {
    stop("`chip_file` must be a single file path.")
  }

  if (!file.exists(chip_file)) {
    stop("ChIP file not found: ", chip_file)
  }

  if (!is.data.frame(hic_df)) {
    stop(
      "`hic_df` must be a data.frame or data.table, ",
      "typically returned by mergeHiCLoops()."
    )
  }

  valid_overlap_modes <- c(
    "any",
    "within",
    "start",
    "end",
    "equal"
  )

  if (!overlap_mode %in% valid_overlap_modes) {
    stop(
      "`overlap_mode` must be one of: ",
      paste(valid_overlap_modes, collapse = ", ")
    )
  }

  # ----------------------------------------------------------
  # 2. Read ChIP-seq input
  # ----------------------------------------------------------

  read_chip_file <- function(file) {

    extension <- tolower(tools::file_ext(file))

    if (extension == "csv") {

      chip <- data.table::fread(
        file,
        sep = ",",
        header = TRUE,
        fill = TRUE,
        data.table = TRUE
      )

    } else if (extension %in% c("tsv", "txt", "tab")) {

      chip <- data.table::fread(
        file,
        sep = "\t",
        header = TRUE,
        fill = TRUE,
        data.table = TRUE
      )

    } else if (extension == "bed") {

      chip <- data.table::fread(
        file,
        sep = "\t",
        header = FALSE,
        fill = TRUE,
        data.table = TRUE
      )

      if (ncol(chip) < 3) {
        stop("BED file must contain at least three columns.")
      }

      data.table::setnames(
        chip,
        old = names(chip)[1:3],
        new = c("chr", "start", "end")
      )

      if (ncol(chip) >= 4 && !"pileup" %in% colnames(chip)) {
        data.table::setnames(
          chip,
          old = names(chip)[4],
          new = "pileup"
        )
      }

    } else if (extension %in% c("xls", "xlsx")) {

      chip <- data.table::as.data.table(
        readxl::read_excel(file)
      )

    } else {

      stop(
        "Unsupported ChIP file type: ",
        extension,
        ". Supported formats are bed, csv, tsv, txt, tab, xls, and xlsx."
      )
    }

    chip
  }

  chip <- read_chip_file(chip_file)
  hic <- data.table::as.data.table(hic_df)

  # ----------------------------------------------------------
  # 3. Normalize column names
  # ----------------------------------------------------------

  if ("BIN1_CHROMOSOME" %in% colnames(hic) &&
      !"BIN1_CHR" %in% colnames(hic)) {
    data.table::setnames(
      hic,
      "BIN1_CHROMOSOME",
      "BIN1_CHR"
    )
  }

  if ("BIN2_CHROMOSOME" %in% colnames(hic) &&
      !"BIN2_CHR" %in% colnames(hic)) {
    data.table::setnames(
      hic,
      "BIN2_CHROMOSOME",
      "BIN2_CHR"
    )
  }

  required_chip_columns <- c(
    "chr",
    "start",
    "end"
  )

  missing_chip_columns <- setdiff(
    required_chip_columns,
    colnames(chip)
  )

  if (length(missing_chip_columns) > 0) {
    stop(
      "Missing ChIP columns: ",
      paste(missing_chip_columns, collapse = ", ")
    )
  }

  if (!"pileup" %in% colnames(chip)) {
    chip[, pileup := 1]
    message(
      "No `pileup` column detected. ",
      "All ChIP peaks were assigned pileup = 1."
    )
  }

  required_hic_columns <- c(
    "BIN1_CHR",
    "BIN1_START",
    "BIN1_END",
    "BIN2_CHR",
    "BIN2_START",
    "BIN2_END",
    "FDR"
  )

  missing_hic_columns <- setdiff(
    required_hic_columns,
    colnames(hic)
  )

  if (length(missing_hic_columns) > 0) {
    stop(
      "Missing Hi-C columns: ",
      paste(missing_hic_columns, collapse = ", ")
    )
  }

  chip <- chip[, .(
    chr,
    start,
    end,
    pileup
  )]

  hic <- hic[, .(
    BIN1_CHR,
    BIN1_START,
    BIN1_END,
    BIN2_CHR,
    BIN2_START,
    BIN2_END,
    FDR
  )]

  # ----------------------------------------------------------
  # 4. Convert data types
  # ----------------------------------------------------------

  chip[, `:=`(
    chr = trimws(as.character(chr)),
    start = suppressWarnings(as.integer(start)),
    end = suppressWarnings(as.integer(end)),
    pileup = suppressWarnings(as.numeric(pileup))
  )]

  hic[, `:=`(
    BIN1_CHR = trimws(as.character(BIN1_CHR)),
    BIN1_START = suppressWarnings(as.integer(BIN1_START)),
    BIN1_END = suppressWarnings(as.integer(BIN1_END)),
    BIN2_CHR = trimws(as.character(BIN2_CHR)),
    BIN2_START = suppressWarnings(as.integer(BIN2_START)),
    BIN2_END = suppressWarnings(as.integer(BIN2_END)),
    FDR = suppressWarnings(as.numeric(FDR))
  )]

  chip <- chip[
    !is.na(chr) &
    chr != "" &
    !is.na(start) &
    !is.na(end) &
    end > start &
    !is.na(pileup)
  ]

  hic <- hic[
    !is.na(BIN1_CHR) &
    BIN1_CHR != "" &
    !is.na(BIN1_START) &
    !is.na(BIN1_END) &
    BIN1_END > BIN1_START &
    !is.na(BIN2_CHR) &
    BIN2_CHR != "" &
    !is.na(BIN2_START) &
    !is.na(BIN2_END) &
    BIN2_END > BIN2_START &
    !is.na(FDR)
  ]

  # ----------------------------------------------------------
  # 5. Normalize chromosome prefixes
  # ----------------------------------------------------------

  normalize_chr <- function(x) {

    x <- trimws(as.character(x))

    needs_prefix <- (
      !is.na(x) &
      x != "" &
      !grepl("^chr", x, ignore.case = TRUE)
    )

    x[needs_prefix] <- paste0(
      "chr",
      x[needs_prefix]
    )

    x
  }

  if (add_chr_prefix) {
    chip[, chr := normalize_chr(chr)]
    hic[, BIN1_CHR := normalize_chr(BIN1_CHR)]
    hic[, BIN2_CHR := normalize_chr(BIN2_CHR)]
  }

  message(
    "Rows: ChIP = ",
    nrow(chip),
    "; Hi-C before FDR filtering = ",
    nrow(hic)
  )

  # ----------------------------------------------------------
  # 6. Filter Hi-C loops by FDR
  # ----------------------------------------------------------

  if (!is.null(fdr_cutoff)) {

    if (length(fdr_cutoff) != 1 ||
        !is.numeric(fdr_cutoff) ||
        is.na(fdr_cutoff) ||
        !is.finite(fdr_cutoff) ||
        fdr_cutoff < 0) {
      stop(
        "`fdr_cutoff` must be NULL or one finite ",
        "non-negative numeric value."
      )
    }

    if (fdr_cutoff < 1) {
      hic <- hic[FDR <= fdr_cutoff]
    }
  }

  message(
    "Rows: Hi-C after FDR filtering = ",
    nrow(hic)
  )

  if (nrow(chip) == 0) {
    stop("No valid ChIP-seq peaks remain after input processing.")
  }

  if (nrow(hic) == 0) {
    warning("No Hi-C loops remain after input processing.")
    return(data.frame())
  }

  # ----------------------------------------------------------
  # 7. Overlap ChIP peaks with BIN1 and project to BIN2
  # ----------------------------------------------------------

  data.table::setkey(
    chip,
    chr,
    start,
    end
  )

  bin1 <- hic[, .(
    chr = BIN1_CHR,
    start = BIN1_START,
    end = BIN1_END,
    partner_chr = BIN2_CHR,
    partner_start = BIN2_START,
    partner_end = BIN2_END,
    FDR
  )]

  data.table::setkey(
    bin1,
    chr,
    start,
    end
  )

  overlap_bin1 <- data.table::foverlaps(
    bin1,
    chip,
    type = overlap_mode,
    nomatch = 0L
  )

  projected_from_bin1 <- overlap_bin1[, .(
    chr = partner_chr,
    start = partner_start,
    end = partner_end,
    pileup,
    FDR,
    source_anchor = "BIN1"
  )]

  # ----------------------------------------------------------
  # 8. Overlap ChIP peaks with BIN2 and project to BIN1
  # ----------------------------------------------------------

  bin2 <- hic[, .(
    chr = BIN2_CHR,
    start = BIN2_START,
    end = BIN2_END,
    partner_chr = BIN1_CHR,
    partner_start = BIN1_START,
    partner_end = BIN1_END,
    FDR
  )]

  data.table::setkey(
    bin2,
    chr,
    start,
    end
  )

  overlap_bin2 <- data.table::foverlaps(
    bin2,
    chip,
    type = overlap_mode,
    nomatch = 0L
  )

  projected_from_bin2 <- overlap_bin2[, .(
    chr = partner_chr,
    start = partner_start,
    end = partner_end,
    pileup,
    FDR,
    source_anchor = "BIN2"
  )]

  message(
    "Projected rows from BIN1 = ",
    nrow(projected_from_bin1),
    "; projected rows from BIN2 = ",
    nrow(projected_from_bin2)
  )

  # ----------------------------------------------------------
  # 9. Merge projected regions
  # ----------------------------------------------------------

  final_matrix <- data.table::rbindlist(
    list(
      projected_from_bin1,
      projected_from_bin2
    ),
    use.names = TRUE,
    fill = TRUE
  )

  if (nrow(final_matrix) == 0) {
    warning(
      "No ChIP-seq peaks overlapped retained Hi-C loop anchors."
    )
    return(as.data.frame(final_matrix))
  }

  # ----------------------------------------------------------
  # 10. Calculate ChIP-SP ranking score
  # ----------------------------------------------------------

  normalize_01 <- function(x) {

    x <- as.numeric(x)

    valid_values <- x[is.finite(x)]

    if (length(valid_values) == 0) {
      return(rep(NA_real_, length(x)))
    }

    value_range <- range(
      valid_values,
      na.rm = TRUE
    )

    if (isTRUE(all.equal(
      value_range[1],
      value_range[2]
    ))) {
      return(rep(0, length(x)))
    }

    (x - value_range[1]) /
      (value_range[2] - value_range[1])
  }

  final_matrix[, pileup_norm := normalize_01(pileup)]
  final_matrix[, fdr_norm := normalize_01(FDR)]
  final_matrix[, score := pileup_norm - fdr_norm]

  data.table::setorder(
    final_matrix,
    -score,
    FDR,
    -pileup
  )

  final_matrix[, rank := seq_len(.N)]

  message(
    "Final ChIP-SP output: ",
    nrow(final_matrix),
    " spatially linked rows."
  )

  as.data.frame(final_matrix)
}

#' Remove sex chromosomes from a genomic data.frame
#'
#' Filters rows corresponding to chromosome X and Y from a genomic
#' interval data.frame. This is useful as a preprocessing step for
#' ChIP-SP workflows when sex chromosomes should be excluded.
#'
#' @param df A data.frame containing genomic intervals.
#' @param chr_col Character scalar. Name of the chromosome column.
#'   Default is `"chr"`.
#' @param remove_chr Character vector of chromosome names to remove.
#'   Default is `c("chrX", "chrY", "X", "Y")`.
#' @param ignore_case Logical; whether chromosome matching should ignore case.
#'   Default is `FALSE`.
#'
#' @return A filtered data.frame with sex chromosomes removed.
#' @export
#'
#' @examples
#' x <- data.frame(
#'   chr = c("chr1", "chrX", "chr2", "chrY"),
#'   start = c(1, 10, 20, 30),
#'   end = c(5, 15, 25, 35)
#' )
#' removeXYChromosomes(x)
removeXYChromosomes <- function(df,
                                chr_col = "chr",
                                remove_chr = c("chrX", "chrY", "X", "Y"),
                                ignore_case = FALSE) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data.frame.")
  }
  
  if (!chr_col %in% colnames(df)) {
    stop("`chr_col` not found in `df`: ", chr_col)
  }
  
  chr_vals <- as.character(df[[chr_col]])
  
  if (ignore_case) {
    keep <- !(tolower(chr_vals) %in% tolower(remove_chr))
  } else {
    keep <- !(chr_vals %in% remove_chr)
  }
  
  out <- df[keep, , drop = FALSE]
  rownames(out) <- NULL
  out
}
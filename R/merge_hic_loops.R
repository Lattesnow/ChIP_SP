#' Merge Hi-C Loop Outputs Across Replicates and Resolutions
#'
#' @description
#' Merges multiple Hi-C loop output files into a single data frame.
#' This function performs no file deletion or writing by default.
#'
#' @param hic_files Character vector of Hi-C loop files (e.g. "*HiC.xls")
#'
#' @return A data.frame containing merged Hi-C loops
#'
#' @examples
#' hic_files <- list.files(pattern = "HiC\\.xls$", full.names = TRUE)
#' hic_df <- mergeHiCLoops(hic_files)
#'
#' @export
mergeHiCLoops <- function(hic_files) {
  
  if (length(hic_files) == 0) {
    stop("No Hi-C files provided.")
  }
  
  hic_list <- lapply(hic_files, function(f) {
    read.delim(f, stringsAsFactors = FALSE)
  })
  
  dplyr::bind_rows(hic_list)
}

#' ChIP-SP Core Spatial Integration and Ranking
#'
#' @description
#' Integrates ChIP-seq peaks with Hi-C chromatin loops by assigning
#' spatially linked regulatory regions and ranking interactions
#' based on ChIP pileup and Hi-C loop confidence.
#'
#' @param chip_file Character. ChIP-seq peak file (*ChIP.xls).
#' @param hic_df data.frame. Output from mergeHiCLoops().
#'
#' @return A data.frame of spatially linked regions ranked by ChIP-SP score.
#'
#' @examples
#' hic_df <- mergeHiCLoops(list.files(pattern = "HiC\\.xls$", full.names = TRUE))
#' res <- chipSPLink("TF_ChIP.xls", hic_df)
#'
#' @export
chipSPLink <- function(chip_file, hic_df) {
  
  if (!file.exists(chip_file)) {
    stop("ChIP file not found.")
  }
  
  if (!is.data.frame(hic_df)) {
    stop("hic_df must be a data.frame (output of mergeHiCLoops).")
  }
  
  file1 <- read.delim(chip_file, stringsAsFactors = FALSE)
  file2 <- hic_df
  
  # Cartesian product
  data1 <- file1[rep(seq_len(nrow(file1)), each = nrow(file2)), ]
  data2 <- file2[rep(seq_len(nrow(file2)), times = nrow(file1)), ]
  
  concatenated_df <- cbind(data1, data2)
  
  # BIN1 overlap
  filtered_BIN1 <- dplyr::filter(
    concatenated_df,
    chr == BIN1_CHR &
      ((BIN1_START < start & start < BIN1_END) |
         (BIN1_START < end & end < BIN1_END))
  )
  
  # BIN2 overlap
  filtered_BIN2 <- dplyr::filter(
    concatenated_df,
    chr == BIN2_CHR &
      ((BIN2_START < start & start < BIN2_END) |
         (BIN2_START < end & end < BIN2_END))
  )
  
  # Generate spatially linked regions
  new_df1 <- filtered_BIN1 %>%
    dplyr::select(chr, BIN2_START, BIN2_END, pileup, FDR) %>%
    dplyr::rename(start = BIN2_START, end = BIN2_END)
  
  new_df2 <- filtered_BIN2 %>%
    dplyr::select(chr, BIN1_START, BIN1_END, pileup, FDR) %>%
    dplyr::rename(start = BIN1_START, end = BIN1_END)
  
  final_matrix <- dplyr::bind_rows(new_df1, new_df2)
  
  # Ranking score (unchanged)
  final_matrix <- final_matrix %>%
    dplyr::mutate(
      pileup_norm = (pileup - min(pileup)) / (max(pileup) - min(pileup)),
      fdr_norm    = (FDR - min(FDR)) / (max(FDR) - min(FDR)),
      score       = pileup_norm - fdr_norm
    ) %>%
    dplyr::arrange(dplyr::desc(score))
  
  final_matrix
}

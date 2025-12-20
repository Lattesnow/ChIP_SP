# Step1 and Step2 for ChIP_SP
# Load packages
library(ChIPSP)

# Step 1: Merge Hi-C loop outputs
# Identify Hi-C loop files
hic_files <- list.files(
  pattern = "HiC\\.xls$",
  full.names = TRUE
)

# Merge Hi-C loops across replicates/resolutions
hic_df <- mergeHiCLoops(hic_files)

# Step 2: ChIP–Hi-C spatial integration and ranking
# Identify ChIP-seq peak file
chip_file <- list.files(
  pattern = "ChIP\\.xls$",
  full.names = TRUE
)

# Run ChIP-SP core algorithm
chipsp_results <- chipSPLink(
  chip_file = chip_file,
  hic_df    = hic_df
)

# Export results to BED format

exportChipSPBed(
  chipsp_results,
  file = "ChIPSP_final_ranked_output.bed"
)

############################################################
## Done
############################################################

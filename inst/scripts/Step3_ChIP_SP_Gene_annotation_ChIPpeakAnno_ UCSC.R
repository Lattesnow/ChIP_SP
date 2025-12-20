library(EnsDb.Hsapiens.v75)  
library(ChIPpeakAnno) 
library(org.Hs.eg.db) 
library(AnnotationDbi)  
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)
library(dplyr)

# Import the output from ChIP-seq peak calling
bed <- read.delim("r1881_hg19_test_peaks_ChIP.xls", sep = "\t")

if (!all(c("chr", "start", "end", "pileup", "fold_enrichment") %in% colnames(bed))) {
  stop("Required columns (chr, start, end, pileup, fold_enrichment) not found in the input file.")
}

# Convert peak data to GRanges object
ChIP_GR <- GRanges(seqnames = bed$chr,
              ranges = IRanges(start = bed$start, end = bed$end),
              pileup = bed$pileup,
              fold_enrichment = bed$fold_enrichment)  

# Create annotation data from UCSC hg19
annoData <- toGRanges(TxDb.Hsapiens.UCSC.hg19.knownGene, feature = "gene")

overlaps.annoA <- annotatePeakInBatch(ChIP_GR, 
                                      AnnotationData = annoData, 
                                      output = "both",
                                      bindingRegion = c(-5000, 5000))

overlaps.annoA <- addGeneIDs(overlaps.annoA,
                             "org.Hs.eg.db",
                             IDs2Add = "symbol",
                             feature_id_type = "entrez_id")

overlaps.annoA.df <- as.data.frame(overlaps.annoA)

overlaps.annoA.df <- merge(overlaps.annoA.df, 
                           as.data.frame(ChIP_GR)[, c("seqnames", "start", "end", "pileup", "fold_enrichment")], 
                           by = c("seqnames", "start", "end"), 
                           all.x = TRUE)

overlaps.annoA.df <- overlaps.annoA.df %>%
  filter(!is.na(symbol)) 

overlaps.annoA.df <- overlaps.annoA.df %>%
  arrange(desc(pileup.x)) %>%  
  mutate(rank = dense_rank(desc(pileup.x)))  

# Save
write.csv(overlaps.annoA.df, 
            file = "ChIP_anno_genes_5Kup5Kdown_UCSC_Both_Control.csv", 
            sep = "/t", 
            row.names = FALSE, 
            quote = FALSE)


# Import ChIPSP file

ChIP_SP_bed <- read.delim("final_ranked_output.xls", sep = "\t")

if (!all(c("chr", "start", "end", "pileup", "score") %in% colnames(ChIP_SP_bed))) {
  stop("Required columns (chr, start, end, pileup, score) not found in the input file.")
}

# Convert peak data to GRanges object
ChIP_SP_GR <- GRanges(seqnames = ChIP_SP_bed$chr,
                      ranges = IRanges(start = ChIP_SP_bed$start, end = ChIP_SP_bed$end),
                      pileup = ChIP_SP_bed$pileup,
                      score = ChIP_SP_bed$score)  

# Create annotation data from UCSC hg19
annoData <- toGRanges(TxDb.Hsapiens.UCSC.hg19.knownGene, feature = "gene")

overlaps.annoB <- annotatePeakInBatch(ChIP_SP_GR, 
                                      AnnotationData = annoData, 
                                      output = "both",
                                      bindingRegion = c(-5000, 5000))

overlaps.annoB <- addGeneIDs(overlaps.annoB,
                             "org.Hs.eg.db",
                             IDs2Add = "symbol",
                             feature_id_type = "entrez_id")

overlaps.annoB.df <- as.data.frame(overlaps.annoB)

overlaps.annoB.df <- merge(overlaps.annoB.df, 
                           as.data.frame(ChIP_SP_GR)[, c("seqnames", "start", "end", "pileup", "score")], 
                           by = c("seqnames", "start", "end"), 
                           all.x = TRUE)

overlaps.annoB.df <- overlaps.annoB.df %>%
  arrange(desc(pileup.x)) %>%  
  mutate(rank = dense_rank(desc(pileup.x)))  

overlaps.annoB.df <- overlaps.annoB.df %>%
  filter(!is.na(symbol)) 

# Save
write.csv(overlaps.annoB.df, 
          file = "ChIP_anno_genes_5Kup5Kdown_UCSC_Both_CHIPSP.csv", 
          sep = "/t", 
          row.names = FALSE, 
          quote = FALSE)



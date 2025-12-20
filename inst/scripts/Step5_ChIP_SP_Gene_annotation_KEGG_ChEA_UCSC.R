# plot KEGG from enrichr
# rename Adjusted.p.value to Adj

library(ggplot2)

kegg_results <- read.csv("KEGG_2021_Human_table_ChIP.txt", sep = "\t")  
chea_results <- read.csv("ChEA_2022_table_ChIP.txt", sep = "\t")  

kegg_results$Adj <- as.numeric(kegg_results$Adj)
chea_results$Adj <- as.numeric(chea_results$Adj)

kegg_results <- kegg_results[!is.na(kegg_results$Adj), ]
chea_results <- chea_results[!is.na(chea_results$Adj), ]

# Select up to 10 pathways (no sorting)
kegg_top <- head(kegg_results, 10)
chea_top <- head(chea_results, 10)

# Convert "Overlap" to numeric
kegg_top$GeneRatio <- as.numeric(sub("/.*", "", kegg_top$Overlap)) /
  as.numeric(sub(".*/", "", kegg_top$Overlap))

chea_top$GeneRatio <- as.numeric(sub("/.*", "", chea_top$Overlap)) /
  as.numeric(sub(".*/", "", chea_top$Overlap))

# function for plot
plot_enrichment <- function(df, title) {
  ggplot(df, aes(x = GeneRatio, y = reorder(Term, GeneRatio), 
                 size = Odds.Ratio, color = Adj)) +
    geom_point() +
    scale_color_gradientn(colors = c("red", "white", "lightblue", "blue")) +
    theme_minimal() +
    labs(title = title, x = "Gene Ratio", y = "Pathway", 
         color = "Adjusted P-value", size = "Odds Ratio") +
    theme(axis.text.y = element_text(size = 8))
}

kegg_plot <- plot_enrichment(kegg_top, "KEGG Pathway Enrichment")
chea_plot <- plot_enrichment(chea_top, "ChEA Transcription Factor Enrichment")

# Save
ggsave("KEGG_dotplot_ChIP.png", plot = kegg_plot, width = 8, height = 5, dpi = 1200)
ggsave("ChEA_dotplot_ChIP.png", plot = chea_plot, width = 8, height = 5, dpi = 1200)

# ChIP_SP
kegg_results <- read.csv("KEGG_2021_Human_table_ChIP_SP.txt", sep = "\t")  
chea_results <- read.csv("ChEA_2022_table_ChIP_SP.txt", sep = "\t")  

kegg_results$Adj <- as.numeric(kegg_results$Adj)
chea_results$Adj <- as.numeric(chea_results$Adj)

kegg_results <- kegg_results[!is.na(kegg_results$Adj), ]
chea_results <- chea_results[!is.na(chea_results$Adj), ]

# Select up to 10 pathways (no sorting)
kegg_top <- head(kegg_results, 10)
chea_top <- head(chea_results, 10)

# Convert "Overlap" to numeric
kegg_top$GeneRatio <- as.numeric(sub("/.*", "", kegg_top$Overlap)) /
  as.numeric(sub(".*/", "", kegg_top$Overlap))

chea_top$GeneRatio <- as.numeric(sub("/.*", "", chea_top$Overlap)) /
  as.numeric(sub(".*/", "", chea_top$Overlap))

# function for plot
plot_enrichment <- function(df, title) {
  ggplot(df, aes(x = GeneRatio, y = reorder(Term, GeneRatio), 
                 size = Odds.Ratio, color = Adj)) +
    geom_point() +
    scale_color_gradientn(colors = c("red", "white", "lightblue", "blue")) +
    theme_minimal() +
    labs(title = title, x = "Gene Ratio", y = "Pathway", 
         color = "Adjusted P-value", size = "Odds Ratio") +
    theme(axis.text.y = element_text(size = 8))
}

kegg_plot <- plot_enrichment(kegg_top, "KEGG Pathway Enrichment")
chea_plot <- plot_enrichment(chea_top, "ChEA Transcription Factor Enrichment")

# Save
ggsave("KEGG_dotplot_ChIP_SP.png", plot = kegg_plot, width = 8, height = 5, dpi = 1200)
ggsave("ChEA_dotplot_ChIP_SP.png", plot = chea_plot, width = 8, height = 5, dpi = 1200)

# ChIP_SP_Only
kegg_results <- read.csv("KEGG_2021_Human_table_Only.txt", sep = "\t")  
chea_results <- read.csv("ChEA_2022_table_Only.txt", sep = "\t")  

kegg_results$Adj <- as.numeric(kegg_results$Adj)
chea_results$Adj <- as.numeric(chea_results$Adj)

kegg_results <- kegg_results[!is.na(kegg_results$Adj), ]
chea_results <- chea_results[!is.na(chea_results$Adj), ]

# Select up to 10 pathways (no sorting)
kegg_top <- head(kegg_results, 10)
chea_top <- head(chea_results, 10)

# Convert "Overlap" to numeric
kegg_top$GeneRatio <- as.numeric(sub("/.*", "", kegg_top$Overlap)) /
  as.numeric(sub(".*/", "", kegg_top$Overlap))

chea_top$GeneRatio <- as.numeric(sub("/.*", "", chea_top$Overlap)) /
  as.numeric(sub(".*/", "", chea_top$Overlap))

# function for plot
plot_enrichment <- function(df, title) {
  ggplot(df, aes(x = GeneRatio, y = reorder(Term, GeneRatio), 
                 size = Odds.Ratio, color = Adj)) +
    geom_point() +
    scale_color_gradientn(colors = c("red", "white", "lightblue", "blue")) +
    theme_minimal() +
    labs(title = title, x = "Gene Ratio", y = "Pathway", 
         color = "Adjusted P-value", size = "Odds Ratio") +
    theme(axis.text.y = element_text(size = 8))
}

kegg_plot <- plot_enrichment(kegg_top, "KEGG Pathway Enrichment")
chea_plot <- plot_enrichment(chea_top, "ChEA Transcription Factor Enrichment")

# Save
ggsave("KEGG_dotplot_ChIP_SP_Only.png", plot = kegg_plot, width = 8, height = 5, dpi = 1200)
ggsave("ChEA_dotplot_ChIP_SP_Only.png", plot = chea_plot, width = 8, height = 5, dpi = 1200)
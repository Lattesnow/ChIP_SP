# ChIP_SP gene annotation Pathway analysis
library(clusterProfiler)
library(org.Hs.eg.db)  
library(enrichplot)    
library(ReactomePA)    
library(ggplot2)

# Import ChIP_SP genes
tab_ChIP <- read.csv("ChIP_anno_genes_upAdown_UCSC_Control.csv")

gene_list_ChIP <- unique(tab_ChIP$symbol)

# Convert gene symbols to Entrez IDs
gene_entrez_ids <- bitr(gene_list_ChIP, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Perform pathway analysis
go_bp <- enrichGO(gene = gene_entrez_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",  
                  pAdjustMethod = "BH", 
                  qvalueCutoff = 0.05, 
                  readable = TRUE)

go_cc <- enrichGO(gene = gene_entrez_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "CC",  
                  pAdjustMethod = "BH", 
                  qvalueCutoff = 0.05, 
                  readable = TRUE)

go_mf <- enrichGO(gene = gene_entrez_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "MF",  
                  pAdjustMethod = "BH", 
                  qvalueCutoff = 0.05, 
                  readable = TRUE)

kegg <- enrichKEGG(gene = gene_entrez_ids$ENTREZID,
                   organism = "hsa", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05)

reactome <- enrichPathway(gene = gene_entrez_ids$ENTREZID,
                          organism = "human", 
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05)

# Plot
go_bp_plot <- dotplot(go_bp, 
                      showCategory = 10, 
                      x = "GeneRatio", 
                      color = "p.adjust", 
                      size = "Count") +
  scale_color_gradientn(colors = c("red", "white","lightblue", "blue")) + 
  theme_minimal() + 
  ggtitle("GO Biological Process (BP)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

go_cc_plot <- dotplot(go_cc, 
                      showCategory = 10, 
                      x = "GeneRatio", 
                      color = "p.adjust", 
                      size = "Count") +
  scale_color_gradientn(colors = c("red", "white","lightblue", "blue")) + 
  theme_minimal() + 
  ggtitle("GO Cellular Component (CC)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

go_mf_plot <- dotplot(go_mf, 
                      showCategory = 10, 
                      x = "GeneRatio", 
                      color = "p.adjust", 
                      size = "Count") +
  scale_color_gradientn(colors = c("red", "white","lightblue", "blue")) + 
  theme_minimal() + 
  ggtitle("GO Molecular Function (MF)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

kegg_plot <- dotplot(kegg, 
                     showCategory = 10, 
                     x = "GeneRatio", 
                     color = "p.adjust", 
                     size = "Count") +
  scale_color_gradientn(colors = c("red", "white","lightblue", "blue")) +  
  theme_minimal() + 
  ggtitle("KEGG Pathway") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

reactome_plot <- dotplot(reactome, 
                         showCategory = 10, 
                         x = "GeneRatio", 
                         color = "p.adjust", 
                         size = "Count") +
  scale_color_gradientn(colors = c("red", "white","lightblue","blue")) +  
  theme_minimal() + 
  ggtitle("Reactome Pathway") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save
out_dir <- "GO_KEGG_Reactome_Plots_ChIP"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(out_dir, "go_bp_plot_ChIP.png"),
       plot = go_bp_plot, width = 8, height = 6, dpi = 1200)

ggsave(file.path(out_dir, "go_cc_plot_ChIP.png"),
       plot = go_cc_plot, width = 8, height = 6, dpi = 1200)

ggsave(file.path(out_dir, "go_mf_plot_ChIP.png"),
       plot = go_mf_plot, width = 8, height = 6, dpi = 1200)

ggsave(file.path(out_dir, "kegg_plot_ChIP.png"),
       plot = kegg_plot, width = 8, height = 6, dpi = 1200)

ggsave(file.path(out_dir, "reactome_plot_ChIP.png"),
       plot = reactome_plot, width = 8, height = 6, dpi = 1200)


#ChIP_SP
tab_CHIPSP <- read.csv("ChIP_anno_genes_upAdown_UCSC_CHIPSP.csv")

gene_list_ChIP_SP <- unique(tab_CHIPSP$symbol)

# Convert gene symbols to Entrez IDs
gene_entrez_ids <- bitr(gene_list_ChIP_SP, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Perform pathway analysis
go_bp <- enrichGO(gene = gene_entrez_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",  
                  pAdjustMethod = "BH", 
                  qvalueCutoff = 0.05, 
                  readable = TRUE)

go_cc <- enrichGO(gene = gene_entrez_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "CC",  
                  pAdjustMethod = "BH", 
                  qvalueCutoff = 0.05, 
                  readable = TRUE)

go_mf <- enrichGO(gene = gene_entrez_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "MF",  
                  pAdjustMethod = "BH", 
                  qvalueCutoff = 0.05, 
                  readable = TRUE)

head(kegg)

kegg <- enrichKEGG(gene = gene_entrez_ids$ENTREZID,
                   organism = "hsa")

reactome <- enrichPathway(gene = gene_entrez_ids$ENTREZID,
                          organism = "human", 
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05)

# Plot
go_bp_plot <- dotplot(go_bp, 
                      showCategory = 10, 
                      x = "GeneRatio", 
                      color = "p.adjust", 
                      size = "Count") +
  scale_color_gradientn(colors = c("red", "white","lightblue", "blue")) + 
  theme_minimal() + 
  ggtitle("GO Biological Process (BP)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

go_cc_plot <- dotplot(go_cc, 
                      showCategory = 10, 
                      x = "GeneRatio", 
                      color = "p.adjust", 
                      size = "Count") +
  scale_color_gradientn(colors = c("red", "white","lightblue", "blue")) + 
  theme_minimal() + 
  ggtitle("GO Cellular Component (CC)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

go_mf_plot <- dotplot(go_mf, 
                      showCategory = 10, 
                      x = "GeneRatio", 
                      color = "p.adjust", 
                      size = "Count") +
  scale_color_gradientn(colors = c("red", "white","lightblue", "blue")) + 
  theme_minimal() + 
  ggtitle("GO Molecular Function (MF)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

kegg_plot <- dotplot(kegg, 
                     showCategory = 10, 
                     x = "GeneRatio", 
                     color = "p.adjust", 
                     size = "Count") +
  scale_color_gradientn(colors = c("red", "white","lightblue", "blue")) +  
  theme_minimal() + 
  ggtitle("KEGG Pathway") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

reactome_plot <- dotplot(reactome, 
                         showCategory = 10, 
                         x = "GeneRatio", 
                         color = "p.adjust", 
                         size = "Count") +
  scale_color_gradientn(colors = c("red", "white","lightblue","blue")) +  
  theme_minimal() + 
  ggtitle("Reactome Pathway") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(head(go_bp))
print(head(go_cc))
print(head(go_mf))
print(head(kegg))
print(head(reactome))

# Save
out_dir <- "GO_KEGG_Reactome_Plots_ChIP_SP"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(out_dir, "go_bp_plot_ChIP_SP.png"),
       plot = go_bp_plot, width = 8, height = 6, dpi = 1200)

ggsave(file.path(out_dir, "go_cc_plot_ChIP_SP.png"),
       plot = go_cc_plot, width = 8, height = 6, dpi = 1200)

ggsave(file.path(out_dir, "go_mf_plot_ChIP_SP.png"),
       plot = go_mf_plot, width = 8, height = 6, dpi = 1200)

ggsave(file.path(out_dir, "kegg_plot_ChIP_SP.png"),
       plot = kegg_plot, width = 8, height = 6, dpi = 1200)

ggsave(file.path(out_dir, "reactome_plot_ChIP_SP.png"),
       plot = reactome_plot, width = 8, height = 6, dpi = 1200)

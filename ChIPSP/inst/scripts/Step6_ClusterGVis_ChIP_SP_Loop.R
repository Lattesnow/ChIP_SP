library(ClusterGVis)
library(readxl)
library(biomaRt)
library(dplyr)
# load data

expr<- read.delim("readcount_genename.xls")
exps <- expr %>% distinct(gene_name, .keep_all = TRUE)
rownames(exps) <- exps$gene_name
exps <- exps[,-1]

total_reads_per_sample <- colSums(exps)
exps <- sweep(exps, 2, total_reads_per_sample, "/") * 1e6  

write.csv(exps, file = "expression_normalized_RPM.csv")

ref <- read.delim("List of ChIP and ChIPSP identified genes.xls")
exp_df <- as.data.frame(exps)
exp_df$gene_id <- rownames(exp_df)

gene_sets <- c("ChIP_SP_Shared","ChIP", "ChIP_SP", "ChIP_SP_Unique" )

filtered_results <- list()

for (gene_set in gene_sets) {

# Convert rownames to a column to use dplyr functions


  exp_filtered_df <- exp_df %>%
    filter(gene_id %in% ref[[gene_set]])

rownames(exp_filtered_df) <- exp_filtered_df$gene_id
exp_filtered_df$gene_id <- NULL 

exp_filtered <- as.matrix(exp_filtered_df)
exps <- as.data.frame(exp_filtered)

exps$DaroR_1<- NULL
exps$DaroR_2<- NULL
# exps <- exps[rowSums(exps) > 200, ]

# check optimal cluster numbers
getClusters(exps)

# using mfuzz for clustering
cm <- clusterData(exps,
                  cluster.method = "mfuzz",
                  cluster.num = 6)
# using TCseq for clustering
ct <- clusterData(exps,
                  cluster.method = "TCseq",
                  cluster.num = 6)

# using kemans for clustering
ck <- clusterData(exps,
                  cluster.method = "kmeans",
                  cluster.num = 6)


# supply other aruguments passed by Heatmap function
visCluster(object = ck,
           plot.type = "heatmap",
           column_names_rot = 45)



markGenes = rownames(exps)[sample(1:nrow(exps),30,replace = F)]

pdf(paste0('addgene_CM_', gene_set, '.pdf'), height = 8, width = 6, onefile = FALSE)
visCluster(object = cm,
           plot.type = "heatmap",
           column_names_rot = 45,
           markGenes = markGenes,
           sample.group = c(rep("R1881", 2), rep("DMSO", 2), rep("Daro", 2)))
dev.off()


pdf(paste0('addgene_CK_', gene_set, '.pdf'), height = 8, width = 6, onefile = FALSE)
visCluster(object = ck,
           plot.type = "heatmap",
           column_names_rot = 45,
           markGenes = markGenes,
           sample.group = c(rep("R1881", 2), rep("DMSO", 2), rep("Daro", 2)))
dev.off()
}
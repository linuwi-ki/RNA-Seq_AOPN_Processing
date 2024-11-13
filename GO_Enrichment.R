library(pheatmap)
library(rafalib)
library(pvclust)
library(enrichR)
library(fgsea)
library(dplyr)
library(tibble)

# Read data for clustering and functional analysis, here we used filtered counts and metadata
data <- read.csv("output/outliers_removed_QC/gene_counts.csv",row.names = 1)
metadata <- read.csv("output/outliers_removed_QC/metadata_raw.csv",row.names = 1,stringsAsFactors=T)
degres <- read.csv("output/outliers_removed_QC/dge_results.csv")
degres2 <- degres[(degres$padj.PCB0_08_vs_Cd0 < 0.05 | degres$padj.PCB0_4_vs_Cd0 < 0.05 | 
                     degres$padj.PCB0_8_vs_Cd0 < 0.05 | degres$padj.Cd2_4_vs_Cd0 < 0.05 |
                     degres$padj.Cd12_vs_Cd0 < 0.05 | degres$padj.Cd24_vs_Cd0 < 0.05),]
write.table(degres2, "output/All_padj0.05_genes_tsv.txt",
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
degres2 <- degres[(abs(degres$log2FoldChange.PCB0_08_vs_Cd0) > 0.58 | abs(degres$log2FoldChange.PCB0_4_vs_Cd0) > 0.58 | 
                    abs(degres$log2FoldChange.PCB0_8_vs_Cd0) > 0.58 | abs(degres$log2FoldChange.Cd2_4_vs_Cd0) > 0.58 |
                    abs(degres$log2FoldChange.Cd12_vs_Cd0) > 0.58 | abs(degres$log2FoldChange.Cd24_vs_Cd0) > 0.58),]
write.table(degres2, "output/All_padj0.05+FC1.5_genes_tsv.txt",
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
rm(degres2)


# ***** ORA for GO terms *****

degres2 <- subset(degres, select = -gene_id)
rownames(degres2) <- NULL
degres2 <- column_to_rownames(degres2, "gene_symbol")
go_id_table <- data.frame(comp = c("PCB0_08_vs_Cd0", "PCB0_4_vs_Cd0", "PCB0_8_vs_Cd0", 
                                   "Cd2_4_vs_Cd0", "Cd12_vs_Cd0", "Cd24_vs_Cd0", 
                                   "Cd24_vs_PCB0_8"),
                            padj_id = c(6, 12, 18, 24, 30, 36, 42), 
                            fc_id = c(2, 8, 14, 20, 26, 32, 38), stringsAsFactors = FALSE)
for (i in 1:nrow(go_id_table)) {
  
  # Change values here to filter differently based on p-val and/or FC
  #selected_genes <- rownames(degres2[(degres2[,as.numeric(go_id_table[i,2])] < 0.05) &
  #                                       (abs(degres2[, go_id_table[i,3]]) > 0.58),])
  selected_genes <- rownames(degres2[(degres2[,as.numeric(go_id_table[i,2])] < 0.05),])
  
  # Create a heatmap to explore clustering
  cl <- pheatmap(data[rownames(data) %in% selected_genes,],scale="row",color=colorRampPalette(c("navy","white","firebrick"))(90),
                  border_color=NA,cluster_cols=F,cutree_rows=2)
  
  gene_clusters <- cutree(cl$tree_row,k=2)
  genes_cluster1 <- names(gene_clusters)[gene_clusters == 1]
  genes_cluster2 <- names(gene_clusters)[gene_clusters == 2]
  
  go_cluster <- enrichr(genes = genes_cluster2,databases = "GO_Biological_Process_2023")
  go_cluster <- go_cluster$GO_Biological_Process_2023
  go_cluster <- go_cluster[order(go_cluster$P.value),]
  
  png(file = paste("output/go_plots/", as.character(go_id_table[i,1]), "_GO_Enrich_noFC.png", sep=""),
      width = 1800, height = 1000, res = 150)
  {
    mypar(1,1,mar=c(4,30,2,2))
    barplot(-log10(go_cluster$P.value[15:1]),horiz=T,border=F,yaxs="i",
            names.arg=go_cluster$Term[15:1],las=1,cex.names=1, xlab = "-log10 P-value")
    abline(v=0,lwd=2)
  }
  dev.off()
  
  #Save all significant GO BPs as a table:
  go_cluster2 <- subset(go_cluster, P.value < 0.05)
  go_cluster2 <- subset(go_cluster2, select=(-c(Old.P.value, Old.Adjusted.P.value)))
  write.table(go_cluster2, paste("output/go_data/", as.character(go_id_table[i,1]), "_GO_tsv_noFC.txt", sep=""),
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  
}


# ***** GSEA for GO and Reactome*****

go_id_table <- data.frame(comp = c("PCB0_08_vs_Cd0", "PCB0_4_vs_Cd0", "PCB0_8_vs_Cd0", 
                                   "Cd2_4_vs_Cd0", "Cd12_vs_Cd0", "Cd24_vs_Cd0", 
                                   "Cd24_vs_PCB0_8"),
                          padj_id = c(6, 12, 18, 24, 30, 36, 42), 
                          fc_id = c(2, 8, 14, 20, 26, 32, 38), stringsAsFactors = FALSE)
# Pathways used for gene enrichment
enrich_pathwaysGO <- gmtPathways("MSigDB_files/c5.go.bp.v2023.1.Hs.symbols.gmt")
enrich_pathwaysGOMF <- gmtPathways("MSigDB_files/c5.go.mf.v2023.2.Hs.symbols.gmt")
enrich_pathwaysReact <- gmtPathways("MSigDB_files/c2.cp.reactome.v2023.1.Hs.symbols.gmt")

for (i in 1:nrow(go_id_table)) {
  #selected_genes <- rownames(degres2[(degres2[,as.numeric(go_id_table[i,2])] < 0.05) &
  #                                     (abs(degres2[, go_id_table[i,3]]) > 0.58),])
  selected_genes <- rownames(degres2[(degres2[,as.numeric(go_id_table[i,2])] < 0.05),])
  
  # GSEA - Create gene ranks based on expression
  gene_rank <- setNames(degres2[,as.numeric(go_id_table[i,3])],casefold(rownames(degres2),upper=T))
  
  # Load pathways to use for enrichment analysis
  # GO BP pathways
  fgseaRes <- fgsea(pathways=enrich_pathwaysGO, stats=gene_rank, minSize=15, maxSize=500, nproc=1)
  #fgseaRes <- fgsea(pathways=enrich_pathwaysGOMF, stats=gene_rank, minSize=15, maxSize=500, nproc=1)
  
  topPathwaysUp <- fgseaRes[NES > 0][head(order(pval), n=10), pathway]
  topPathwaysDown <- fgseaRes[NES < 0][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  
  # Nice summary table (shown as a plot)
  png(file = paste("output/go_plots/", as.character(go_id_table[i,1]), "_GSEA_GOBP_noFC.png", sep=""),
      type = "cairo", width = 2100, height = 1000, res = 150)
  gseaplot <- plotGseaTable(enrich_pathwaysGO[topPathways], gene_rank, fgseaRes, gseaParam = 0.5, c(7, 2.5, 0.8, 1.2, 1.2))
  print(gseaplot)
  dev.off()
  
  #Save all significant GSEA pathways as a table (GO only right now):
  fgseaRes <- subset(fgseaRes, pval < 0.05)
  fgseaRes <- fgseaRes[order(fgseaRes$pval),]
  fgseaRes$leadingEdge <- sapply(fgseaRes$leadingEdge, function(x) paste(x, collapse = ", "))
  write.table(fgseaRes, file = paste("output/go_data/", go_id_table[i,1], "_GSEA_GOBP_tsv_noFC.txt", sep=""),
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = TRUE)
  
  # # Reactome pathways
  fgseaRes <- fgsea(pathways=enrich_pathwaysReact, stats=gene_rank, minSize=15, maxSize=500, nproc=1)

  topPathwaysUp <- fgseaRes[NES > 0][head(order(pval), n=10), pathway]
  topPathwaysDown <- fgseaRes[NES < 0][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

  # # Nice summary table (shown as a plot)
  png(file = paste("output/go_plots/", as.character(go_id_table[i,1]), "_GSEA_React_noFC.png", sep=""),
      type = "cairo", width = 2100, height = 1000, res = 150)
  gseaplot <- plotGseaTable(enrich_pathwaysReact[topPathways], gene_rank, fgseaRes, gseaParam = 0.5, c(7, 2.5, 0.8, 1.2, 1.2))
  print(gseaplot)
  dev.off()
  
}

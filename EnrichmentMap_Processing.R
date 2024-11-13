library(dplyr)
library(tibble)

# Function to convert GO IDs to GO term descriptions

EM_formatting <- function(GOBP) {
  #Subset data and relocate rows to work in EnrichmentMap
  GOBP <- subset(GOBP, select = -c(log2err,ES,size))
  GOBP$description <- "NA"
  GOBP <- relocate(GOBP, "description", .after = pathway)
  colnames(GOBP) <- c("GO_pathway", "description", "pval", "FDR", "NES", "gene_list")
  
  return(GOBP)
}

# Read the enrichmentmap data from the tab-delimited file
enrichmentmap<-read.delim("output/go_data/No DEG FC filter/Cd2_4_vs_Cd0_GSEA_GOBP_tsv_noFC.txt",
                          header = TRUE, sep = "\t")
#Format data
results <- EM_formatting(enrichmentmap)
#Save table as tsv, this is then imported into Cytoscape.
write.table(results, file = "output/go_data/No DEG FC filter/EnrichmentMap files/GSEA_Cd2_4_for_EnrichmentMap_noFC.txt", sep = "\t", row.names = FALSE, quote = FALSE)
results <- results[results$FDR<0.05,]
write.table(results, file = "output/go_data/No DEG FC filter/EnrichmentMap files/GSEA_Cd2_4_for_EnrichmentMap_noFC_only_padj.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#Repeat for all substances/doses
enrichmentmap<-read.delim("output/go_data/No DEG FC filter/Cd12_vs_Cd0_GSEA_GOBP_tsv_noFC.txt",
                          header = TRUE, sep = "\t")
results <- EM_formatting(enrichmentmap)
write.table(results, file = "output/go_data/No DEG FC filter/EnrichmentMap files/GSEA_Cd12_for_EnrichmentMap_noFC.txt", sep = "\t", row.names = FALSE, quote = FALSE)
results <- results[results$FDR<0.05,]
write.table(results, file = "output/go_data/No DEG FC filter/EnrichmentMap files/GSEA_Cd12_for_EnrichmentMap_noFC_only_padj.txt", sep = "\t", row.names = FALSE, quote = FALSE)

enrichmentmap<-read.delim("output/go_data/No DEG FC filter/Cd24_vs_Cd0_GSEA_GOBP_tsv_noFC.txt",
                          header = TRUE, sep = "\t")
results <- EM_formatting(enrichmentmap)
write.table(results, file = "output/go_data/No DEG FC filter/EnrichmentMap files/GSEA_Cd24_for_EnrichmentMap_noFC.txt", sep = "\t", row.names = FALSE, quote = FALSE)
results <- results[results$FDR<0.05,]
write.table(results, file = "output/go_data/No DEG FC filter/EnrichmentMap files/GSEA_Cd24_for_EnrichmentMap_noFC_only_padj.txt", sep = "\t", row.names = FALSE, quote = FALSE)

enrichmentmap<-read.delim("output/go_data/No DEG FC filter/PCB0_08_vs_Cd0_GSEA_GOBP_tsv_noFC.txt",
                          header = TRUE, sep = "\t")
results <- EM_formatting(enrichmentmap)
write.table(results, file = "output/go_data/No DEG FC filter/EnrichmentMap files/GSEA_PCB0_08_for_EnrichmentMap_noFC.txt", sep = "\t", row.names = FALSE, quote = FALSE)
results <- results[results$FDR<0.05,]
write.table(results, file = "output/go_data/No DEG FC filter/EnrichmentMap files/GSEA_PCB0_08_for_EnrichmentMap_noFC_only_padj.txt", sep = "\t", row.names = FALSE, quote = FALSE)

enrichmentmap<-read.delim("output/go_data/No DEG FC filter/PCB0_4_vs_Cd0_GSEA_GOBP_tsv_noFC.txt",
                          header = TRUE, sep = "\t")
results <- EM_formatting(enrichmentmap)
write.table(results, file = "output/go_data/No DEG FC filter/EnrichmentMap files/GSEA_PCB0_4_for_EnrichmentMap_noFC.txt", sep = "\t", row.names = FALSE, quote = FALSE)
results <- results[results$FDR<0.05,]
write.table(results, file = "output/go_data/No DEG FC filter/EnrichmentMap files/GSEA_PCB0_4_for_EnrichmentMap_noFC_only_padj.txt", sep = "\t", row.names = FALSE, quote = FALSE)

enrichmentmap<-read.delim("output/go_data/No DEG FC filter/PCB0_8_vs_Cd0_GSEA_GOBP_tsv_noFC.txt",
                          header = TRUE, sep = "\t")
results <- EM_formatting(enrichmentmap)
write.table(results, file = "output/go_data/No DEG FC filter/EnrichmentMap files/GSEA_PCB0_8_for_EnrichmentMap_noFC.txt", sep = "\t", row.names = FALSE, quote = FALSE)
results <- results[results$FDR<0.05,]
write.table(results, file = "output/go_data/No DEG FC filter/EnrichmentMap files/GSEA_PCB0_8_for_EnrichmentMap_noFC_only_padj.txt", sep = "\t", row.names = FALSE, quote = FALSE)


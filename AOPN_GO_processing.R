library(dplyr)
library(tibble)
library(tidyr)
library(GO.db)

# Load GO terms from the AOPN
AOPN_GO <- read.delim("KE-list for automatic mapping.txt")
colnames(AOPN_GO) <- c("Description", "KE_ID", "Object_Term", "Process_Term", "Object_Ontology_ID",
                    "go_id")

# Load GSEA GO terms
Cd2_4_GSEA <- read.delim("output/go_data/NO DEG FC filter/Cd2_4_vs_Cd0_GSEA_GOBP_tsv_noFC.txt")
Cd12_GSEA <- read.delim("output/go_data/NO DEG FC filter/Cd12_vs_Cd0_GSEA_GOBP_tsv_noFC.txt")
Cd24_GSEA <- read.delim("output/go_data/NO DEG FC filter/Cd24_vs_Cd0_GSEA_GOBP_tsv_noFC.txt")
PCB0_08_GSEA <- read.delim("output/go_data/NO DEG FC filter/PCB0_08_vs_Cd0_GSEA_GOBP_tsv_noFC.txt")
PCB0_4_GSEA <- read.delim("output/go_data/NO DEG FC filter/PCB0_4_vs_Cd0_GSEA_GOBP_tsv_noFC.txt")
PCB0_8_GSEA <- read.delim("output/go_data/NO DEG FC filter/PCB0_8_vs_Cd0_GSEA_GOBP_tsv_noFC.txt")

# As GO terms from GSEA do not contain GO IDs, we need to
# convert the GOBP_ terms into GO ID to compare them to the AOPN
lt = as.list(GOTERM)
map = sapply(lt, function(x) Term(x))
map = map[sapply(lt, function(x) Ontology(x) == "BP")]
map = toupper(map)
map = gsub(" ", "_", map)
map2 = structure(names(map), names = map)

new_rn = map2[ gsub("^GOBP_", "", Cd2_4_GSEA$pathway) ]
l = !is.na(new_rn)
Cd2_4_GSEA = Cd2_4_GSEA[l, ]
rownames(Cd2_4_GSEA) = new_rn[l]
Cd2_4_GSEA <- rownames_to_column(Cd2_4_GSEA, var="go_id")

new_rn = map2[ gsub("^GOBP_", "", Cd12_GSEA$pathway) ]
l = !is.na(new_rn)
Cd12_GSEA = Cd12_GSEA[l, ]
rownames(Cd12_GSEA) = new_rn[l]
Cd12_GSEA <- rownames_to_column(Cd12_GSEA, var="go_id")

new_rn = map2[ gsub("^GOBP_", "", Cd24_GSEA$pathway) ]
l = !is.na(new_rn)
Cd24_GSEA = Cd24_GSEA[l, ]
rownames(Cd24_GSEA) = new_rn[l]
Cd24_GSEA <- rownames_to_column(Cd24_GSEA, var="go_id")

new_rn = map2[ gsub("^GOBP_", "", PCB0_08_GSEA$pathway) ]
l = !is.na(new_rn)
PCB0_08_GSEA = PCB0_08_GSEA[l, ]
rownames(PCB0_08_GSEA) = new_rn[l]
PCB0_08_GSEA <- rownames_to_column(PCB0_08_GSEA, var="go_id")

new_rn = map2[ gsub("^GOBP_", "", PCB0_4_GSEA$pathway) ]
l = !is.na(new_rn)
PCB0_4_GSEA = PCB0_4_GSEA[l, ]
rownames(PCB0_4_GSEA) = new_rn[l]
PCB0_4_GSEA <- rownames_to_column(PCB0_4_GSEA, var="go_id")

new_rn = map2[ gsub("^GOBP_", "", PCB0_8_GSEA$pathway) ]
l = !is.na(new_rn)
PCB0_8_GSEA = PCB0_8_GSEA[l, ]
rownames(PCB0_8_GSEA) = new_rn[l]
PCB0_8_GSEA <- rownames_to_column(PCB0_8_GSEA, var="go_id")

#Keep only significant genes based on padj < 0.05
Cd2_4_GSEA_sig <- Cd2_4_GSEA[Cd2_4_GSEA$padj < 0.05,]
Cd12_GSEA_sig <- Cd12_GSEA[Cd12_GSEA$padj < 0.05,]
Cd24_GSEA_sig <- Cd24_GSEA[Cd24_GSEA$padj < 0.05,]
PCB0_08_GSEA_sig <- PCB0_08_GSEA[PCB0_08_GSEA$padj < 0.05,]
PCB0_4_GSEA_sig <- PCB0_4_GSEA[PCB0_4_GSEA$padj < 0.05,]
PCB0_8_GSEA_sig <- PCB0_8_GSEA[PCB0_8_GSEA$padj < 0.05,]

# Get a list of all tables in the environment that match a pattern (padj significance)
all_tables <- mget(ls(pattern = "_GSEA_sig"))
result_df <- data.frame()
description_ke_id <- AOPN_GO[, c("go_id", "Description", "KE_ID")]

# Iterate over each additional table
for (i in seq_along(all_tables)) {
  # Attempt to merge the table with the existing result_df
  merged_df <- merge(all_tables[[i]], description_ke_id, by="go_id", all.x=TRUE)
  
  # Filter out rows where "Description" and "KE_ID" are NA
  merged_df <- merged_df[complete.cases(merged_df[, c("Description", "KE_ID")]), ]
  
  # Check if the merging was successful
  if (nrow(merged_df) > 0) {
    # Add an "exposure" column with the name of the original table
    merged_df$exposure <- names(all_tables)[i]
    
    # Update result_df with the merged data
    result_df <- rbind(result_df, merged_df[, c("go_id", "Description", "KE_ID", "pathway", "exposure")])
  } else {
    cat("Unable to merge with table", names(all_tables)[i], "\n")
  }
}

write.table(result_df, "output/Automatic_GO_AOPN_match_padj_sig_tsv.txt",
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


library(dplyr) # data wrangling
library(tibble) # data wrangling
library(tidyr) # data wrangling
library(pheatmap) # heatmap generation
library(ggplot2)
library(readr)

ipa_up_reg <- read.delim(file = "Cd+PCB_Upstream_reg_padj_and_zscore-filtered.txt", sep = "\t")
colnames(ipa_up_reg) <- c("Upstream.Regulator", "PCB 0.08", "PCB 0.4", "PCB 0.8", "Cd 2.4", "Cd 12", "Cd 24")
ipa_up_reg_all <- read.delim(file = "Cd+PCB_Upstream_reg_table_only_padj_sig.txt", sep = "\t")
ipa_up_reg_all <- subset(ipa_up_reg_all, select = c(Analysis, Upstream.Regulator, B.H.corrected.p.value))

ipa_up_reg_all <- ipa_up_reg_all %>%
  pivot_wider(names_from = Analysis, values_from = B.H.corrected.p.value)

# Merge the tables
merged_table <- merge(ipa_up_reg_all, ipa_up_reg, by = "Upstream.Regulator")
merged_table <- merged_table %>%
  mutate_all(~na_if(., "N/A"))

for (i in 2:7) {
  merged_table[, i] <- gsub("−", "-", merged_table[, i])
  merged_table[, i] <- gsub(",", ".", merged_table[, i])
}
# Function to convert scientific notation values to numeric decimal
convert_to_decimal <- function(x) {
  # Initialize vector to store converted values
  numeric_values <- numeric(length(x))
  
  # Convert each value
  for (i in seq_along(x)) {
    # Check if the value is in 'scientific notation' format
    if (grepl("×", x[i]) && grepl("-", x[i])) {
      # Split into coefficient and exponent parts
      parts <- unlist(strsplit(x[i], "×10^", fixed = TRUE))
      if (length(parts) == 2) {
        # Extract coefficient and exponent parts
        coefficient <- as.numeric(parts[1])
        exponent <- as.numeric(parts[2])
        numeric_values[i] <- coefficient * 10^exponent
      } else {
        numeric_values[i] <- NA  # Handle cases where the format is not as expected
      }
    } else {
      numeric_values[i] <- as.numeric(x[i])  # If not in scientific notation, convert directly
    }
  }
  return(numeric_values)
}

# Apply the conversion function to columns 2-7 of the merged_table
for (i in 2:7) {
  merged_table[, i] <- convert_to_decimal(merged_table[, i])
}

# convert the other columns to numeric to make sure all columns except the first one is numeric
merged_table[, 8:13] <- lapply(merged_table[, 8:13], as.numeric)

# Get the indices for columns ending in ".x" and ".y"
x_cols <- grep("\\.x$", names(merged_table))
y_cols <- grep("\\.y$", names(merged_table))

# Order the columns ending in ".x" and ".y" separately
ordered_x_cols <- names(merged_table)[x_cols][order(gsub("^.* ", "", names(merged_table)[x_cols]))]
ordered_y_cols <- names(merged_table)[y_cols][order(gsub("^.* ", "", names(merged_table)[y_cols]))]

# Combine and rename the columns
ordered_cols <- c("Upstream.Regulator", ordered_x_cols, ordered_y_cols)
new_names <- c("Upstream.Regulator", 
               gsub("\\.x$", " pval", ordered_x_cols), 
               gsub("\\.y$", " z-score", ordered_y_cols))

# Reorder and rename the columns in the merged_table
merged_table <- merged_table[, ordered_cols]
names(merged_table) <- new_names

merged_table <- merged_table %>%
  relocate("Cd 2.4 pval", .before = "Cd 12 pval")

merged_table <- merged_table %>%
  relocate("Cd 2.4 z-score", .before = "Cd 12 z-score")

# just to check which type of class each column is.
#sapply(merged_table, class)

# Select top 50 Upstream.Regulator values
top_50_regulators <- head(ipa_up_reg$Upstream.Regulator, 50)

# Filter merged_table using the top 50 Upstream.Regulator values
top_50_merged <- merged_table %>%
  filter(Upstream.Regulator %in% top_50_regulators)

heatmap_matrix <- top_50_merged
heatmap_matrix <- remove_rownames(heatmap_matrix)
heatmap_matrix <- column_to_rownames(heatmap_matrix, var = "Upstream.Regulator")
heatmap_matrix <- as.matrix(heatmap_matrix)

df <- as.data.frame(heatmap_matrix)
df <- df[order(df$`PCB 0.8 z-score`), ]
heatmap_matrix <- as.matrix(df)
remove(df)

#number_format = %.2f did not work in the pheatmap function, so we limit to 2 decimals here
heatmap_matrix[, 7:12] <- apply(heatmap_matrix[, 7:12], c(1, 2), function(x) as.numeric(sprintf("%.2f", x)))

# Replace NA values with 1 in the heatmap matrix
#heatmap_matrix[is.na(heatmap_matrix)] <- 1

# Create the heatmap with pheatmap
pheatmap(heatmap_matrix[,1:6], display_numbers = heatmap_matrix[,7:12],
         cluster_rows = FALSE, cluster_cols = FALSE,
         number_color = "black", na_col = "lightgrey",
         color = colorRampPalette(c("orangered2", "orange", "yellow1"))(50),
         labels_col = c("PCB 0.08", "PCB 0.4", "PCB 0.8", "Cd 2.4", "Cd 12", "Cd 24"),
         cellwidth = 60, cellheight = 25, fontsize = 22, fontsize_row = 18) 

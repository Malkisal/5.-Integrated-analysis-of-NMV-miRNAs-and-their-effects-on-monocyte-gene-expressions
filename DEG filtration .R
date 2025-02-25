# Load required libraries
library(readxl)
library(dplyr)

# Read files
mirna_targets <- read_excel("Final_miRNA_Atherosclerosis_Targets.xlsx")
nldl_genes <- read.csv("sig_diff_Pvalue_CTRL vs nLDL_NMVs3.csv")
oxldl_genes <- read.csv("sig_diff_Pvalue CTRL vs oxLDL_NMVs3.csv")

# Create a mapping of genes to their miRNAs
gene_mirna_map <- mirna_targets %>%
  select(miRNA = 1, `Target Gene`) %>%
  group_by(`Target Gene`) %>%
  summarise(miRNAs = paste(unique(miRNA), collapse = "; "))

# Filter and add miRNA information for nLDL dataset
nldl_filtered <- nldl_genes %>%
  filter(row %in% unique(mirna_targets$`Target Gene`)) %>%
  left_join(gene_mirna_map, by = c("row" = "Target Gene"))

# Filter and add miRNA information for oxLDL dataset
oxldl_filtered <- oxldl_genes %>%
  filter(row %in% unique(mirna_targets$`Target Gene`)) %>%
  left_join(gene_mirna_map, by = c("row" = "Target Gene"))

# Save results
write.csv(nldl_filtered, "filtered_nLDL_genes_with_mirna.csv", row.names = FALSE)
write.csv(oxldl_filtered, "filtered_oxLDL_genes_with_mirna.csv", row.names = FALSE)

# Print summary
cat("Matches in nLDL:", nrow(nldl_filtered), "\n")
cat("Matches in oxLDL:", nrow(oxldl_filtered), "\n")


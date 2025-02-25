# Load required libraries
library(VennDiagram)
library(readr)
library(dplyr)
library(scales)

# Read the filtered files
nldl_genes <- read_csv("filtered_nLDL_genes_with_mirna.csv")
oxldl_genes <- read_csv("filtered_oxLDL_genes_with_mirna.csv")

# Rename columns to match
names(nldl_genes)[names(nldl_genes) == "miRNAs in nLDL-NMVs"] <- "miRNAs"

# Extract gene lists
nldl_list <- nldl_genes$Gene
oxldl_list <- oxldl_genes$Gene

# Create Venn diagram with png instead of tiff
png("gene_overlap_venn.png", width = 600, height = 600 )
draw.pairwise.venn(
  area1 = length(nldl_list),
  area2 = length(oxldl_list),
  cross.area = length(intersect(nldl_list, oxldl_list)),
  category = c("nLDL", "oxLDL"),
  col = c("#440154FF", "#21908CFF"),
  fill = c(alpha("#440154FF", 0.3), alpha("#21908CFF", 0.3)),
  main = "Overlapping Genes between nLDL and oxLDL"
)
dev.off()

# Process gene lists
overlapped_genes <- intersect(nldl_list, oxldl_list)
unique_nldl <- setdiff(nldl_list, oxldl_list)
unique_oxldl <- setdiff(oxldl_list, nldl_list)

# Create separate dataframes for Excel export
overlapped_nldl <- nldl_genes[nldl_genes$Gene %in% overlapped_genes, ]
overlapped_oxldl <- oxldl_genes[oxldl_genes$Gene %in% overlapped_genes, ]
unique_nldl_data <- nldl_genes[nldl_genes$Gene %in% unique_nldl, ]
unique_oxldl_data <- oxldl_genes[oxldl_genes$Gene %in% unique_oxldl, ]

# Export to Excel
library(writexl)
write_xlsx(list(
  Overlapped_Genes_nLDL = overlapped_nldl,
  Overlapped_Genes_oxLDL = overlapped_oxldl,
  Unique_nLDL = unique_nldl_data,
  Unique_oxLDL = unique_oxldl_data
), path = "gene_lists_analysis.xlsx")

# Print summary
cat("Number of overlapping genes:", length(overlapped_genes), "\n")
cat("Number of unique nLDL genes:", length(unique_nldl), "\n")
cat("Number of unique oxLDL genes:", length(unique_oxldl), "\n")


# Set up custom colors
color1 <- "#4B0082"  # Deep indigo
color2 <- "#006400"  # Dark green

# Create high-resolution Venn diagram
png("gene_overlap_venn_publication.png", 
    width = 2000, 
    height = 2000, 
    res = 400)

draw.pairwise.venn(
  area1 = length(nldl_list),
  area2 = length(oxldl_list),
  cross.area = length(intersect(nldl_list, oxldl_list)),
  category = c("nLDL-NMV", "oxLDL-NMV"),
  lwd = 3,  # Border thickness
  lty = "solid",
  col = c(color1, color2),  # Border colors
  fill = c(alpha(color1, 0.4), alpha(color2, 0.4)),  # Fill colors
  cex = 2.5,  # Size of numbers
  cat.cex = 2,  # Size of category names
  cat.col = c(color1, color2),  # Category name colors
  cat.dist = 0.05,  # Distance of category names
  cat.pos = c(-10, 10),  # Position of category names
  euler.d = TRUE,
  scaled = TRUE,
  main = "Differential Gene Expression Overlap",
  main.cex = 3
)

dev.off()



# Load required libraries
library(readxl)
library(dplyr)
library(rtracklayer)
library(writexl)

# Function to standardize miRNA names
standardize_mirna_name <- function(name) {
  # Convert to lowercase for comparison
  name <- tolower(name)
  
  # Handle let-7 family
  if(grepl("^hsa-let", name)) {
    return(name)
  }
  
  # For all other miRNAs, keep only the base name without variant numbers
  base_name <- sub("(hsa-mir-[0-9]+[a-z]*).*$", "\\1", name)
  
  return(base_name)
}

# Function to extract miRNA lengths from GFF3 file
get_mirna_lengths <- function(gff_file) {
  print("Reading GFF3 file...")
  gff_data <- import(gff_file)
  
  # Get only primary transcript features
  precursor_features <- gff_data[gff_data$type == "miRNA_primary_transcript"]
  
  print(paste("Number of precursor features found:", length(precursor_features)))
  
  # Convert to data frame and calculate lengths
  mirna_lengths <- data.frame(
    miRNA = precursor_features$Name,
    Original_Name = precursor_features$Name,
    ID = precursor_features$ID,
    Length = width(precursor_features),
    stringsAsFactors = FALSE
  )
  
  # Add standardized names
  mirna_lengths$miRNA_std <- sapply(mirna_lengths$miRNA, standardize_mirna_name)
  
  # Handle duplicates by taking the longest variant for each miRNA
  dup_indices <- which(!duplicated(mirna_lengths$miRNA_std))
  if(length(dup_indices) < nrow(mirna_lengths)) {
    unique_names <- unique(mirna_lengths$miRNA_std)
    keep_indices <- sapply(unique_names, function(name) {
      rows <- which(mirna_lengths$miRNA_std == name)
      rows[which.max(mirna_lengths$Length[rows])]
    })
    mirna_lengths <- mirna_lengths[keep_indices, ]
  }
  
  print(paste("Number of unique precursor miRNAs:", nrow(mirna_lengths)))
  return(mirna_lengths)
}

# Function to calculate TPM using correct total reads
calculate_tpm <- function(counts, lengths, total_reads) {
  # Step 1: Calculate RPK (reads per kilobase)
  rpk <- counts / (lengths / 1000)
  
  # Print RPK range for debugging
  print("RPK range:")
  print(range(rpk, na.rm = TRUE))
  
  # Create the correctly ordered total reads vector to match sample order in counts
  ordered_total_reads <- c(
    total_reads[6],  # Donor 6
    total_reads[3],  # Donor 3
    total_reads[2],  # Donor 2
    total_reads[1],  # Donor 1
    total_reads[4],  # Donor 4
    total_reads[5]   # Donor 5
  )
  
  print("Using total reads per sample:")
  names(ordered_total_reads) <- colnames(counts)[1:6]
  print(ordered_total_reads)
  
  # Step 2: Calculate TPM for each sample (excluding Mean column)
  tpm <- matrix(0, nrow=nrow(rpk), ncol=ncol(rpk))
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  
  # Calculate TPM for each sample except Mean
  for(i in 1:6) {
    tpm[,i] <- rpk[,i] * 1e6 / ordered_total_reads[i]
  }
  
  # Calculate mean TPM for the last column
  tpm[,7] <- rowMeans(tpm[,1:6], na.rm=TRUE)
  
  return(tpm)
}

# Main script
main <- function() {
  # Read input files
  print("Reading input files...")
  mirna_counts <- read_excel("Monocyte_Dataset.xlsx")
  total_reads_df <- read_excel("Total reads.xlsx")
  
  # Extract total unique miRNA reads from row 5
  total_mirna_reads <- as.numeric(total_reads_df[5, -1])
  print("Total unique miRNA reads per sample:")
  print(total_mirna_reads)
  
  # Get miRNA lengths
  mirna_lengths <- get_mirna_lengths("hsa.gff3")
  
  # Prepare count matrix
  mirna_names <- mirna_counts[[1]]
  count_matrix <- as.matrix(mirna_counts[,-1])
  rownames(count_matrix) <- mirna_names
  
  # Standardize miRNA names
  std_names <- sapply(mirna_names, standardize_mirna_name)
  
  # Match lengths with standardized names
  print("Matching miRNA names with lengths...")
  matched_indices <- match(std_names, mirna_lengths$miRNA_std)
  print(paste("Number of matched miRNAs:", sum(!is.na(matched_indices))))
  
  if(sum(is.na(matched_indices)) > 0) {
    print(paste("Number of unmatched miRNAs:", sum(is.na(matched_indices))))
    print("Sample of unmatched miRNAs:")
    print(head(data.frame(
      Original = mirna_names[is.na(matched_indices)],
      Standardized = std_names[is.na(matched_indices)]
    )))
  }
  
  matched_lengths <- mirna_lengths$Length[matched_indices]
  
  # Remove rows with missing lengths
  valid_rows <- !is.na(matched_lengths)
  filtered_counts <- count_matrix[valid_rows, ]
  filtered_lengths <- matched_lengths[valid_rows]
  filtered_names <- mirna_names[valid_rows]
  
  print(paste("Processing", sum(valid_rows), "miRNAs with valid lengths"))
  
  # Calculate TPM with correct total reads
  print("Calculating TPM values...")
  tpm_matrix <- calculate_tpm(filtered_counts, filtered_lengths, total_mirna_reads)
  
  # Create output dataframe with original names
  tpm_df <- data.frame(
    miRNA = filtered_names,
    tpm_matrix,
    stringsAsFactors = FALSE
  )
  
  # Create matching information
  matching_info <- data.frame(
    Original_Name = mirna_names,
    Standardized_Name = std_names,
    GFF3_Name = mirna_lengths$Original_Name[matched_indices],
    Length = matched_lengths,
    Matched = !is.na(matched_indices),
    stringsAsFactors = FALSE
  )
  
  # Write output files
  print("Writing output files...")
  write_xlsx(tpm_df, "miRNA_TPM_normalized.xlsx")
  write_xlsx(matching_info, "miRNA_matching_info.xlsx")
  
  print("Process completed!")
  
  return(list(
    total_mirnas = length(mirna_names),
    matched_mirnas = sum(!is.na(matched_indices)),
    unmatched_mirnas = sum(is.na(matched_indices))
  ))
}

# Run the script
tryCatch({
  results <- main()
  print("\nSummary:")
  print(paste("Total miRNAs:", results$total_mirnas))
  print(paste("Successfully matched:", results$matched_mirnas))
  print(paste("Unmatched:", results$unmatched_mirnas))
}, error = function(e) {
  print(paste("Error occurred:", e$message))
})


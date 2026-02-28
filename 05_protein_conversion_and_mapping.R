#!/usr/bin/env Rscript
# scripts/05_protein_conversion_and_mapping_simple.R

cat("========================================\n")
cat("PROTEIN CONVERSION & MAPPING (Simple)\n")
cat("========================================\n\n")

# Load packages
if (!require("dplyr", quietly = TRUE)) install.packages("dplyr")
library(dplyr)

# Configuration
ranked_dir <- "C:/Users/DELL/Documents/my_project/ranked_lists"
output_dir <- "data/protein_mappings"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Find RNK files
rnk_files <- list.files(ranked_dir, pattern = "\\.rnk$", full.names = TRUE)
rnk_files <- rnk_files[!grepl("_genesymbols\\.rnk$", rnk_files)]

if (length(rnk_files) == 0) {
  cat("❌ No RNK files found\n")
  quit(status = 1)
}

cat("Found", length(rnk_files), "RNK files to process\n\n")

# Process each file
all_mappings <- list()

for (rnk_file in rnk_files) {
  cat("Processing:", basename(rnk_file), "\n")
  
  # Read protein IDs
  protein_data <- read.delim(rnk_file, header = FALSE,
                             col.names = c("protein_id", "score"),
                             stringsAsFactors = FALSE)
  
  # Ask user for mapping approach
  cat("  How do you want to map proteins to genes?\n")
  cat("  1. Use existing mapping file\n")
  cat("  2. Enter gene symbols manually (for few proteins)\n")
  cat("  3. Skip this file\n")
  
  choice <- readline("  Enter choice (1-3): ")
  
  if (choice == "1") {
    # Use mapping file
    mapping_file <- readline("  Enter path to mapping file (CSV with protein_id,gene_symbol): ")
    if (file.exists(mapping_file)) {
      mapping <- read.csv(mapping_file, stringsAsFactors = FALSE)
      merged <- merge(protein_data, mapping, by = "protein_id", all.x = TRUE)
    } else {
      cat("  Mapping file not found\n")
      next
    }
  } else if (choice == "2") {
    # Manual entry (for small datasets)
    merged <- protein_data
    merged$gene_symbol <- ""
    
    for (i in 1:min(10, nrow(merged))) {  # Limit to first 10
      cat("  Protein:", merged$protein_id[i], "\n")
      merged$gene_symbol[i] <- readline("  Enter gene symbol: ")
    }
  } else {
    next
  }
  
  # Save gene symbols file for GSEA
  if (sum(!is.na(merged$gene_symbol) & merged$gene_symbol != "") > 0) {
    gsea_data <- data.frame(
      gene_symbol = merged$gene_symbol,
      score = merged$score
    )
    gsea_data <- gsea_data[!is.na(gsea_data$gene_symbol) & gsea_data$gene_symbol != "", ]
    
    output_file <- sub("\\.rnk$", "_genesymbols.rnk", rnk_file)
    write.table(gsea_data, output_file, 
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    cat("  ✓ Saved gene symbols file:", basename(output_file), "\n")
    
    # Save mapping
    mapping_data <- merged %>%
      select(protein_id, gene_symbol, score) %>%
      filter(!is.na(gene_symbol) & gene_symbol != "")
    
    all_mappings[[basename(rnk_file)]] <- mapping_data
    cat("  ✓ Saved mapping for", nrow(mapping_data), "proteins\n")
  }
}

# Combine and save all mappings
if (length(all_mappings) > 0) {
  cat("\nCombining all mappings...\n")
  
  combined <- bind_rows(all_mappings)
  unique_mapping <- combined %>%
    select(protein_id, gene_symbol) %>%
    distinct() %>%
    arrange(protein_id, gene_symbol)
  
  # Save
  write.csv(unique_mapping, 
            file.path(output_dir, "protein_gene_mapping.csv"),
            row.names = FALSE)
  
  # Create lookup
  protein_to_gene <- setNames(unique_mapping$gene_symbol, unique_mapping$protein_id)
  saveRDS(protein_to_gene, 
          file.path(output_dir, "protein_to_gene_lookup.rds"))
  
  cat("✓ Saved comprehensive mapping:\n")
  cat("  - protein_gene_mapping.csv\n")
  cat("  - protein_to_gene_lookup.rds\n")
  cat("  Unique mappings:", nrow(unique_mapping), "\n")
}

cat("\n✅ Done! Now run GSEA with *_genesymbols.rnk files\n")
#!/usr/bin/env Rscript
# scripts/05_protein_gene_mapping_tracking.R
# Simple version - just extracts mapping from your conversion

cat("========================================\n")
cat("SIMPLE PROTEIN-GENE MAPPING EXTRACTION\n")
cat("========================================\n\n")

# ============================================================
# 1. FIND ALL ORIGINAL AND CONVERTED FILES
# ============================================================

ranked_dir <- "C:/Users/DELL/Documents/my_project/ranked_lists"

# Get all original .rnk files
original_files <- list.files(ranked_dir, pattern = "\\.rnk$", full.names = TRUE)
original_files <- original_files[!grepl("_genesymbols\\.rnk$", original_files)]

if (length(original_files) == 0) {
  cat("❌ No original .rnk files found\n")
  quit(status = 1)
}

cat("Found", length(original_files), "original RNK files\n")

# ============================================================
# 2. EXTRACT MAPPING FROM EACH FILE PAIR
# ============================================================

all_mappings <- list()

for (original_file in original_files) {
  # Find corresponding genesymbols file
  gene_file <- sub("\\.rnk$", "_genesymbols.rnk", original_file)
  
  if (!file.exists(gene_file)) {
    cat("Skipping", basename(original_file), "- no converted file\n")
    next
  }
  
  cat("Processing", basename(original_file), "...\n")
  
  # Read both files
  proteins <- read.delim(original_file, header = FALSE, 
                         col.names = c("protein_id", "score"),
                         stringsAsFactors = FALSE)
  
  genes <- read.delim(gene_file, header = FALSE,
                      col.names = c("gene_symbol", "score"),
                      stringsAsFactors = FALSE)
  
  # Check if they have the same number of rows
  if (nrow(proteins) == nrow(genes)) {
    # Perfect match - simple merge
    mapping <- data.frame(
      protein_id = proteins$protein_id,
      gene_symbol = genes$gene_symbol,
      score = proteins$score,
      contrast = gsub("\\.rnk$", "", basename(original_file)),
      stringsAsFactors = FALSE
    )
  } else {
    # Mismatch - merge by score
    mapping <- merge(proteins, genes, by = "score", 
                     suffixes = c("_protein", "_gene"))
    mapping$contrast <- gsub("\\.rnk$", "", basename(original_file))
  }
  
  # Clean gene symbols (take first if multiple)
  mapping$gene_clean <- sapply(strsplit(mapping$gene_symbol, ";|\\|"), function(x) x[1])
  mapping$gene_clean <- gsub("\\s+", "", mapping$gene_clean)
  
  all_mappings[[basename(original_file)]] <- mapping
  cat("  Extracted", nrow(mapping), "mappings\n")
}

# ============================================================
# 3. COMBINE AND SAVE
# ============================================================

if (length(all_mappings) == 0) {
  cat("❌ No mappings extracted\n")
  quit(status = 1)
}

# Combine all
combined <- do.call(rbind, all_mappings)

# Get unique mappings
unique_mapping <- combined %>%
  select(protein_id, gene_symbol = gene_clean) %>%
  filter(!is.na(gene_symbol) & gene_symbol != "") %>%
  distinct() %>%
  arrange(protein_id, gene_symbol)

cat("\nFinal mapping statistics:\n")
cat("Unique proteins:", length(unique(unique_mapping$protein_id)), "\n")
cat("Unique genes:", length(unique(unique_mapping$gene_symbol)), "\n")
cat("Total mappings:", nrow(unique_mapping), "\n")

# ============================================================
# 4. SAVE FILES
# ============================================================

output_dir <- "data/protein_mappings"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Save simple mapping
write.csv(unique_mapping, 
          file.path(output_dir, "protein_to_gene_mapping.csv"),
          row.names = FALSE)

# Create and save lookup
protein_to_gene <- setNames(unique_mapping$gene_symbol, unique_mapping$protein_id)
saveRDS(protein_to_gene, file.path(output_dir, "protein_to_gene_lookup.rds"))

cat("\n✅ Saved mapping files to", output_dir, "\n")
cat("1. protein_to_gene_mapping.csv\n")
cat("2. protein_to_gene_lookup.rds\n")
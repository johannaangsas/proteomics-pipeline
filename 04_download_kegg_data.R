#!/usr/bin/env Rscript
# scripts/04_download_kegg_data.R - FIXED VERSION
# Loads KEGG data from local GMT file for offline GSEA analysis

cat("========================================\n")
cat("KEGG DATA PREPARATION (Local GMT File)\n")
cat("========================================\n\n")

# ============================================================
# 1. CONFIGURATION
# ============================================================

KEGG_CONFIG <- list(
  organism = "hsa",                         # Homo sapiens
  output_dir = "data",                      # Where to save data
  local_gmt_file = "data/c2.cp.kegg_legacy.v2026.1.Hs.symbols.gmt"  # YOUR GMT FILE
)

# ============================================================
# 2. CREATE OUTPUT DIRECTORY
# ============================================================

cat("1. Setting up directories...\n")

if (!dir.exists(KEGG_CONFIG$output_dir)) {
  dir.create(KEGG_CONFIG$output_dir, recursive = TRUE, showWarnings = FALSE)
  cat("  Created directory:", KEGG_CONFIG$output_dir, "\n")
}

# ============================================================
# 3. DEFINE OUTPUT PATHS (CRITICAL FIX)
# ============================================================

cat("2. Defining output paths...\n")

output_path <- file.path(KEGG_CONFIG$output_dir, "kegg_human_data.rds")
csv_path <- file.path(KEGG_CONFIG$output_dir, "kegg_human_pathways.csv")
metadata_path <- file.path(KEGG_CONFIG$output_dir, "kegg_metadata.rds")
report_path <- file.path(KEGG_CONFIG$output_dir, "kegg_download_report.txt")

cat("  Main output file:", output_path, "\n")
cat("  CSV output file:", csv_path, "\n")

# ============================================================
# 4. LOAD KEGG DATA FROM LOCAL GMT FILE
# ============================================================

cat("3. Loading KEGG data from local GMT file...\n")

# Check if the GMT file exists
if (!file.exists(KEGG_CONFIG$local_gmt_file)) {
  cat("❌ ERROR: GMT file not found at:", KEGG_CONFIG$local_gmt_file, "\n")
  cat("  Please ensure the file exists in the correct location.\n")
  quit(status = 1)
}

cat("  Found GMT file:", basename(KEGG_CONFIG$local_gmt_file), "\n")

# Function to parse a GMT file into a list
parse_gmt <- function(file_path) {
  gmt_lines <- readLines(file_path)
  gene_sets <- list()
  
  for (line in gmt_lines) {
    fields <- strsplit(line, "\t")[[1]]
    pathway_name <- fields[1]
    gene_symbols <- fields[-(1:2)] # Skip pathway name and optional description
    gene_symbols <- gene_symbols[gene_symbols != ""] # Remove empty strings
    if (length(gene_symbols) > 0) {
      gene_sets[[pathway_name]] <- gene_symbols
    }
  }
  return(gene_sets)
}

# Parse the GMT file into a named list of gene sets
cat("  Parsing GMT file...\n")
kegg_gene_sets_list <- parse_gmt(KEGG_CONFIG$local_gmt_file)

if (length(kegg_gene_sets_list) == 0) {
  cat("❌ ERROR: No gene sets found in the GMT file.\n")
  quit(status = 1)
}

cat("  Successfully parsed", length(kegg_gene_sets_list), "pathways.\n")

# ============================================================
# 5. CONVERT TO CLUSTERPROFILER COMPATIBLE FORMAT
# ============================================================

cat("4. Converting to clusterProfiler compatible format...\n")

# Convert the list into a long-format data frame (KEGGPATHID2EXTID structure)
pathway_gene_pairs <- data.frame(
  from = rep(names(kegg_gene_sets_list), 
             times = lengths(kegg_gene_sets_list)),
  to = unlist(kegg_gene_sets_list, use.names = FALSE),
  stringsAsFactors = FALSE
)

# Remove any potential duplicates
pathway_gene_pairs <- unique(pathway_gene_pairs)

# Create the KEGGPATHID2NAME data frame (pathway ID to description)
pathway_names <- data.frame(
  from = names(kegg_gene_sets_list),
  to = names(kegg_gene_sets_list), # Using the same as ID for simplicity
  stringsAsFactors = FALSE
)

# Assemble the final kegg_data object in clusterProfiler's format
kegg_data <- list(
  KEGGPATHID2EXTID = pathway_gene_pairs,
  KEGGPATHID2NAME = pathway_names
)

cat("  Converted to", length(unique(kegg_data$KEGGPATHID2EXTID$from)), 
    "unique pathways\n")
cat("  Total gene-pathway associations:", nrow(kegg_data$KEGGPATHID2EXTID), "\n")

# ============================================================
# 6. SAVE DATA
# ============================================================

cat("5. Saving data for offline GSEA...\n")

# Save as RDS (efficient for R)
saveRDS(kegg_data, output_path)
cat("  ✓ Saved KEGG data to:", output_path, "\n")

# Also save as CSV for easy inspection
write.csv(kegg_data$KEGGPATHID2EXTID, csv_path, row.names = FALSE)
cat("  ✓ Saved pathway-gene mapping to:", csv_path, "\n")

# ============================================================
# 7. CREATE METADATA
# ============================================================

cat("6. Creating metadata...\n")

metadata <- list(
  creation_date = Sys.time(),
  data_source = paste("Local GMT file:", basename(KEGG_CONFIG$local_gmt_file)),
  organism = "Homo sapiens",
  pathways_count = length(unique(kegg_data$KEGGPATHID2EXTID$from)),
  genes_count = length(unique(kegg_data$KEGGPATHID2EXTID$to)),
  total_entries = nrow(kegg_data$KEGGPATHID2EXTID),
  gmt_file = KEGG_CONFIG$local_gmt_file
)

saveRDS(metadata, metadata_path)
cat("  ✓ Saved metadata to:", metadata_path, "\n")

# ============================================================
# 8. GENERATE SUMMARY REPORT
# ============================================================

cat("7. Generating summary report...\n")

report <- paste(
  "KEGG DATA PREPARATION REPORT (Local GMT File)",
  "==============================================",
  paste("Date:", Sys.time()),
  paste("Source file:", basename(KEGG_CONFIG$local_gmt_file)),
  paste("Organism:", metadata$organism),
  paste("Pathways loaded:", metadata$pathways_count),
  paste("Unique genes:", metadata$genes_count),
  paste("Total pathway-gene associations:", metadata$total_entries),
  "",
  "FILE LOCATIONS:",
  paste("1. Main KEGG data (RDS):", output_path),
  paste("2. Pathway-gene map (CSV):", csv_path),
  paste("3. Metadata (RDS):", metadata_path),
  paste("4. Original GMT file:", KEGG_CONFIG$local_gmt_file),
  "",
  "NEXT STEPS:",
  "1. Run protein ID conversion (05_protein_id_conversion.R)",
  "2. Run GSEA analysis (06_gsea_analysis.R) -- Will now work OFFLINE",
  "",
  "NOTE: This data is ready for offline GSEA analysis using clusterProfiler.",
  sep = "\n"
)

writeLines(report, report_path)
cat("  ✓ Saved report to:", report_path, "\n")

# ============================================================
# 9. FINAL SUMMARY
# ============================================================

cat("\n", paste0(rep("=", 60), collapse = ""), "\n", sep = "")
cat("✅ KEGG DATA PREPARATION COMPLETE\n")
cat(paste0(rep("=", 60), collapse = ""), "\n\n")

cat("Summary:\n")
cat("  Source: Local GMT file (", basename(KEGG_CONFIG$local_gmt_file), ")\n", sep = "")
cat("  Pathways:", metadata$pathways_count, "\n")
cat("  Unique genes:", metadata$genes_count, "\n")
cat("  Total associations:", metadata$total_entries, "\n\n")

cat("Data successfully saved in:", KEGG_CONFIG$output_dir, "\n")
cat("Your GSEA analysis pipeline is now fully offline.\n")
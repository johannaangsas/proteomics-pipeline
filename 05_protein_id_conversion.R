#!/usr/bin/env Rscript
# scripts/05_protein_id_conversion.R
# Converts UniProt IDs in .rnk files to gene symbols for GSEA
# Run after 03_de_analysis_for_gsea.R

cat("========================================\n")
cat("PROTEIN ID TO GENE SYMBOL CONVERSION\n")
cat("========================================\n\n")

# ============================================================
# 1. LOAD CONFIGURATION & PACKAGES
# ============================================================
# Add this as the FIRST line in 05_protein_id_conversion.R
ANALYSIS_CONFIG <- list(ranked_dir = "C:/Users/DELL/Documents/my_project/ranked_lists")


# Source your configuration (adjust path if needed)
config_file <- "C:/Users/DELL/Documents/my_project/config/analysis_parameters.R"
if (file.exists(config_file)) {
  source(config_file)
  cat("Loaded configuration from:", config_file, "\n")
} else {
  cat("Configuration file not found. Using defaults.\n")
  ANALYSIS_CONFIG <- list(
    ranked_dir = "ranked_lists",
    create_ranked_lists = TRUE
  )
}

# Load required packages
required_packages <- c("protti", "dplyr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# ============================================================
# 2. DEFINE CONVERSION FUNCTION
# ============================================================

convert_rnk_to_genesymbols <- function(rnk_file_path) {
  #' Convert UniProt IDs in .rnk files to gene symbols
  #'
  #' @param rnk_file_path Full path to the input .rnk file
  #' @return List with conversion statistics
  
  cat("  Processing:", basename(rnk_file_path), "\n")
  
  # 1. Read the .rnk file
  rnk_data <- tryCatch({
    read.delim(rnk_file_path, header = FALSE, 
               col.names = c("UniProt_ID", "Ranking_Score"))
  }, error = function(e) {
    cat("    ✗ ERROR reading file:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(rnk_data) || nrow(rnk_data) == 0) {
    cat("    ✗ File empty or not found\n")
    return(list(success = FALSE, n_total = 0, n_mapped = 0))
  }
  
  # 2. Get unique IDs (UniProt limits to 500 per query)
  protein_ids <- unique(rnk_data$UniProt_ID)
  cat("    Found", length(protein_ids), "unique protein IDs\n")
  
  # 3. Query UniProt in batches if needed
  cat("    Querying UniProt...\n")
  
  # Process in batches of 500 (UniProt API limit)
  batch_size <- 500
  n_batches <- ceiling(length(protein_ids) / batch_size)
  all_protein_info <- list()
  
  for (batch_idx in 1:n_batches) {
    start_idx <- (batch_idx - 1) * batch_size + 1
    end_idx <- min(batch_idx * batch_size, length(protein_ids))
    batch_ids <- protein_ids[start_idx:end_idx]
    
    cat("      Batch", batch_idx, "/", n_batches, 
        "(", length(batch_ids), "IDs)...\n")
    
    protein_info <- tryCatch({
      fetch_uniprot(
        uniprot_ids = batch_ids,
        columns = c("accession", "gene_names", "protein_name")
      )
    }, error = function(e) {
      cat("        Batch failed:", e$message, "\n")
      return(data.frame())
    })
    
    if (nrow(protein_info) > 0) {
      all_protein_info[[batch_idx]] <- protein_info
    }
    
    # Be nice to the UniProt server
    if (batch_idx < n_batches) Sys.sleep(0.5)
  }
  
  # Combine all batch results
  if (length(all_protein_info) == 0) {
    cat("    ✗ No data retrieved from UniProt\n")
    return(list(success = FALSE, 
                n_total = nrow(rnk_data), 
                n_mapped = 0))
  }
  
  protein_info_combined <- bind_rows(all_protein_info)
  
  # 4. Merge and create new .rnk file
  rnk_annotated <- rnk_data %>%
    left_join(protein_info_combined, by = c("UniProt_ID" = "accession")) %>%
    filter(!is.na(gene_names)) %>%
    select(Gene_Symbol = gene_names, Ranking_Score)
  
  # 5. Save new file
  new_path <- sub("\\.rnk$", "_genesymbols.rnk", rnk_file_path)
  write.table(rnk_annotated, new_path,
              sep = "\t", row.names = FALSE, 
              col.names = FALSE, quote = FALSE)
  
  # 6. Generate mapping statistics
  n_mapped <- nrow(rnk_annotated)
  n_total <- nrow(rnk_data)
  mapping_rate <- round(n_mapped / n_total * 100, 1)
  
  cat("    ✓ Successfully mapped", n_mapped, "/", n_total, 
      " (", mapping_rate, "%)\n", sep = "")
  cat("    ✓ Saved:", basename(new_path), "\n\n")
  
  return(list(
    success = TRUE,
    n_total = n_total,
    n_mapped = n_mapped,
    mapping_rate = mapping_rate,
    input_file = rnk_file_path,
    output_file = new_path
  ))
}

# ============================================================
# 3. MAIN EXECUTION: PROCESS ALL .RNK FILES
# ============================================================

cat("\n1. Locating .rnk files...\n")

# Find all .rnk files in the ranked_lists directory
if (!dir.exists(ANALYSIS_CONFIG$ranked_dir)) {
  cat("  ✗ Directory not found:", ANALYSIS_CONFIG$ranked_dir, "\n")
  cat("  Please run the DE analysis first to generate .rnk files.\n")
  quit(status = 1)
}

rnk_files <- list.files(ANALYSIS_CONFIG$ranked_dir, 
                        pattern = "\\.rnk$", 
                        full.names = TRUE)

# Filter out already converted files
rnk_files <- rnk_files[!grepl("_genesymbols\\.rnk$", rnk_files)]

if (length(rnk_files) == 0) {
  cat("  ✗ No .rnk files found to convert.\n")
  cat("  Files may have already been converted.\n")
} else {
  cat("  Found", length(rnk_files), ".rnk file(s) to convert:\n")
  for (f in rnk_files) cat("   -", basename(f), "\n")
  
  cat("\n2. Starting conversion...\n")
  
  # Process each file and collect statistics
  conversion_results <- list()
  successful_conversions <- 0
  
  for (i in seq_along(rnk_files)) {
    result <- convert_rnk_to_genesymbols(rnk_files[i])
    conversion_results[[i]] <- result
    
    if (result$success) {
      successful_conversions <- successful_conversions + 1
    }
  }
  
  # ============================================================
  # 4. GENERATE CONVERSION REPORT
  # ============================================================
  
  cat("\n3. Generating conversion report...\n")
  
  # Calculate overall statistics
  total_proteins <- sum(sapply(conversion_results, function(x) x$n_total))
  total_mapped <- sum(sapply(conversion_results, function(x) x$n_mapped))
  overall_rate <- round(total_mapped / total_proteins * 100, 1)
  
  # Create report
  report <- paste(
    "PROTEIN ID CONVERSION REPORT",
    "=============================",
    paste("Date:", Sys.time()),
    paste("Input directory:", ANALYSIS_CONFIG$ranked_dir),
    "",
    "SUMMARY:",
    paste("Files processed:", length(rnk_files)),
    paste("Successful conversions:", successful_conversions),
    paste("Total proteins processed:", total_proteins),
    paste("Total mapped to gene symbols:", total_mapped),
    paste("Overall mapping rate:", overall_rate, "%"),
    "",
    "DETAILED RESULTS:",
    sep = "\n"
  )
  
  # Add per-file details
  for (i in seq_along(conversion_results)) {
    r <- conversion_results[[i]]
    report <- paste0(report, 
                     "\n", basename(rnk_files[i]), ":",
                     "\n  Mapped: ", r$n_mapped, "/", r$n_total,
                     " (", r$mapping_rate, "%)",
                     "\n  Output: ", basename(r$output_file))
  }
  
  report <- paste0(report,
                   "\n\nNEXT STEPS:",
                   "\n1. Use *_genesymbols.rnk files for GSEA analysis",
                   "\n2. Original .rnk files are preserved for reference",
                   "\n3. Check unmapped proteins in UniProt if mapping rate is low")
  
  # Save report
  report_file <- file.path(ANALYSIS_CONFIG$ranked_dir, "conversion_report.txt")
  writeLines(report, report_file)
  cat("  ✓ Report saved:", report_file, "\n")
  
  # Print final summary
  cat("\n", paste0(rep("=", 50), collapse = ""), "\n", sep = "")
  cat("✅ CONVERSION COMPLETE\n")
  cat(paste0(rep("=", 50), collapse = ""), "\n\n")
  cat("Summary:\n")
  cat("  Files converted:", successful_conversions, "/", length(rnk_files), "\n")
  cat("  Proteins mapped:", total_mapped, "/", total_proteins, 
      " (", overall_rate, "%)\n", sep = "")
  cat("\nGSEA-ready files are available with '_genesymbols.rnk' suffix.\n")
}
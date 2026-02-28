#!/usr/bin/env Rscript
# scripts/check_mapping.R
# Minimal diagnostic for .rnk file mapping

cat("========================================\n")
cat("MINIMAL MAPPING DIAGNOSTIC\n")
cat("========================================\n\n")

# ============================================================
# 1. SETUP WITHOUT EXTERNAL PACKAGES
# ============================================================

# Use base R only - no package dependencies
dir_to_check <- "ranked_lists"
cat("Checking directory:", dir_to_check, "\n\n")

# Check if directory exists
if (!dir.exists(dir_to_check)) {
  cat("❌ ERROR: Directory not found\n")
  cat("Current working directory:", getwd(), "\n")
  cat("Files in current directory:\n")
  print(list.files())
  quit(status = 1)
}

# ============================================================
# 2. FIND AND CHECK FILES
# ============================================================

cat("Looking for .rnk files...\n")

# Get all files
all_files <- list.files(dir_to_check, full.names = TRUE)
cat("Total files in directory:", length(all_files), "\n")

# Find original .rnk files (not _genesymbols)
original_files <- all_files[grepl("\\.rnk$", all_files) & !grepl("_genesymbols", all_files)]
cat("Original .rnk files found:", length(original_files), "\n\n")

if (length(original_files) == 0) {
  cat("No original .rnk files found. Listing all files:\n")
  print(basename(all_files))
  quit(status = 0)
}

# ============================================================
# 3. SIMPLE DIAGNOSTIC FOR EACH FILE
# ============================================================

cat("FILE DIAGNOSTICS:\n")
cat(paste0(rep("-", 70), collapse = ""), "\n")

results <- data.frame()

for (orig_file in original_files) {
  file_name <- basename(orig_file)
  conv_file <- file.path(dir_to_check, 
                         sub("\\.rnk$", "_genesymbols.rnk", file_name))
  
  cat("File:", file_name, "\n")
  
  # Check if converted file exists
  if (!file.exists(conv_file)) {
    cat("  ❌ MISSING: No converted file found\n")
    next
  }
  
  # Try to read files
  tryCatch({
    # Read first few lines to check format
    orig_lines <- readLines(orig_file, n = 5)
    conv_lines <- readLines(conv_file, n = 5)
    
    # Count lines (simple row count)
    orig_count <- length(readLines(orig_file))
    conv_count <- length(readLines(conv_file))
    
    # Check first column format
    orig_first_col <- strsplit(orig_lines[1], "\t")[[1]][1]
    conv_first_col <- strsplit(conv_lines[1], "\t")[[1]][1]
    
    # Determine if likely converted
    is_uniprot_id <- grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]", orig_first_col)
    is_gene_symbol <- !grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]", conv_first_col)
    
    # Calculate mapping rate
    mapping_rate <- round(conv_count / orig_count * 100, 1)
    
    cat(sprintf("  Original: %4d rows, first ID: %s\n", 
                orig_count, substr(orig_first_col, 1, 15)))
    cat(sprintf("  Converted: %4d rows, first ID: %s\n", 
                conv_count, substr(conv_first_col, 1, 15)))
    cat(sprintf("  Mapping rate: %s%%\n", mapping_rate))
    
    # Check top 20 specifically
    if (orig_count > 20) {
      orig_top <- read.delim(orig_file, header = FALSE, nrows = 20)
      conv_top <- read.delim(conv_file, header = FALSE, nrows = 20)
      
      # Simple check: if first columns look different
      if (orig_top[1,1] != conv_top[1,1]) {
        cat("  ✓ Top proteins appear to be converted\n")
      } else {
        cat("  ⚠ WARNING: First IDs identical - may not be converted\n")
      }
    }
    
    # Store result
    results <- rbind(results, data.frame(
      File = file_name,
      Original_Rows = orig_count,
      Converted_Rows = conv_count,
      Mapping_Rate = mapping_rate,
      First_Original = substr(orig_first_col, 1, 15),
      First_Converted = substr(conv_first_col, 1, 15),
      stringsAsFactors = FALSE
    ))
    
  }, error = function(e) {
    cat("  ❌ ERROR reading files:", e$message, "\n")
  })
  
  cat(paste0(rep("-", 70), collapse = ""), "\n")
}

# ============================================================
# 4. SUMMARY
# ============================================================

cat("\nSUMMARY:\n")

if (nrow(results) > 0) {
  # Print table
  cat("\n")
  cat(sprintf("%-25s %-10s %-10s %-10s %-15s %-15s\n",
              "File", "Orig", "Conv", "Rate", "First Orig", "First Conv"))
  cat(paste0(rep("-", 85), collapse = ""), "\n")
  
  for (i in 1:nrow(results)) {
    cat(sprintf("%-25s %-10d %-10d %-10s %-15s %-15s\n",
                substr(results$File[i], 1, 24),
                results$Original_Rows[i],
                results$Converted_Rows[i],
                paste0(results$Mapping_Rate[i], "%"),
                results$First_Original[i],
                results$First_Converted[i]))
  }
  
  # Overall stats
  avg_rate <- mean(results$Mapping_Rate, na.rm = TRUE)
  cat("\n")
  cat(paste0(rep("=", 50), collapse = ""), "\n")
  cat("OVERALL: ", round(avg_rate, 1), "% average mapping rate\n", sep = "")
  cat(paste0(rep("=", 50), collapse = ""), "\n")
  
  # Check for low mapping rates
  low_rate <- results[results$Mapping_Rate < 70, ]
  if (nrow(low_rate) > 0) {
    cat("\n⚠ FILES WITH LOW MAPPING RATE (<70%):\n")
    for (i in 1:nrow(low_rate)) {
      cat("  ", low_rate$File[i], ": ", low_rate$Mapping_Rate[i], "%\n", sep = "")
    }
  }
} else {
  cat("No files successfully analyzed.\n")
}

cat("\n✅ Diagnostic complete\n")
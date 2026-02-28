#!/usr/bin/env Rscript
# scripts/08_analyze_leading_edge_PROTEIN.R
# Updated to work with the merged conversion+mapping script

cat("========================================\n")
cat("PROTEIN-LEVEL LEADING EDGE ANALYSIS\n")
cat("Using Saved Protein-Gene Mapping\n")
cat("========================================\n\n")

# ============================================================
# 1. LOAD REQUIRED PACKAGES
# ============================================================

cat("1. Loading packages...\n")

# Check for dplyr
if (!require("dplyr", quietly = TRUE)) {
  install.packages("dplyr", repos = "https://cloud.r-project.org")
  library(dplyr)
}

# ============================================================
# 2. LOAD SAVED PROTEIN-GENE MAPPING
# ============================================================

cat("2. Loading saved protein-gene mapping...\n")

# Try different possible locations/names
possible_mapping_files <- c(
  "data/protein_mappings/protein_mapping_lookup.rds",
  "data/protein_mappings/protein_gene_mapping.csv",
  "data/protein_mappings/protein_to_gene_lookup.rds",
  "data/protein_mappings/protein_to_gene_simple.rds"
)

mapping_file <- NULL
for (file in possible_mapping_files) {
  if (file.exists(file)) {
    mapping_file <- file
    break
  }
}

if (is.null(mapping_file)) {
  cat("❌ ERROR: No mapping file found.\n")
  cat("  Looking for files in: data/protein_mappings/\n")
  cat("  Available files:\n")
  all_files <- list.files("data/protein_mappings", full.names = TRUE)
  if (length(all_files) > 0) {
    for (f in all_files) cat("    -", f, "\n")
  } else {
    cat("    (directory empty or doesn't exist)\n")
  }
  cat("\n  Run scripts/05_protein_id_conversion_and_mapping.R first.\n")
  quit(status = 1)
}

cat("  ✓ Found mapping file:", mapping_file, "\n")

# Load based on file type
if (grepl("\\.rds$", mapping_file)) {
  # RDS file
  mapping_data <- readRDS(mapping_file)
  
  # Handle different RDS formats
  if (is.list(mapping_data) && "protein_to_gene" %in% names(mapping_data)) {
    # Complete lookup format
    protein_to_gene <- mapping_data$protein_to_gene
    gene_to_proteins <- mapping_data$gene_to_proteins
    cat("  ✓ Loaded comprehensive mapping object\n")
  } else if (is.character(mapping_data)) {
    # Simple named vector format
    protein_to_gene <- mapping_data
    # Create gene_to_proteins from it
    gene_to_proteins <- list()
    for (protein in names(protein_to_gene)) {
      gene <- protein_to_gene[protein]
      if (!gene %in% names(gene_to_proteins)) {
        gene_to_proteins[[gene]] <- character(0)
      }
      gene_to_proteins[[gene]] <- c(gene_to_proteins[[gene]], protein)
    }
    cat("  ✓ Loaded simple mapping vector\n")
  } else {
    cat("❌ ERROR: Unknown RDS format\n")
    quit(status = 1)
  }
} else if (grepl("\\.csv$", mapping_file)) {
  # CSV file
  mapping_df <- read.csv(mapping_file, stringsAsFactors = FALSE)
  
  # Check column names
  if ("protein_id" %in% colnames(mapping_df) && "gene_symbol" %in% colnames(mapping_df)) {
    protein_to_gene <- setNames(mapping_df$gene_symbol, mapping_df$protein_id)
    
    # Create gene_to_proteins
    gene_to_proteins <- list()
    for (i in 1:nrow(mapping_df)) {
      gene <- mapping_df$gene_symbol[i]
      protein <- mapping_df$protein_id[i]
      
      if (!gene %in% names(gene_to_proteins)) {
        gene_to_proteins[[gene]] <- character(0)
      }
      gene_to_proteins[[gene]] <- c(gene_to_proteins[[gene]], protein)
    }
    cat("  ✓ Loaded CSV mapping\n")
  } else {
    cat("❌ ERROR: CSV missing required columns (protein_id, gene_symbol)\n")
    quit(status = 1)
  }
} else {
  cat("❌ ERROR: Unknown file format\n")
  quit(status = 1)
}

cat("  ✓ Mapping covers", length(protein_to_gene), "proteins\n")
cat("  ✓ Mapping to", length(gene_to_proteins), "unique genes\n")

# Show examples
if (length(protein_to_gene) > 0) {
  cat("  Example mappings:\n")
  for (i in 1:min(3, length(protein_to_gene))) {
    prot <- names(protein_to_gene)[i]
    gene <- protein_to_gene[i]
    cat("    ", prot, " → ", gene, "\n", sep = "")
  }
}

# ============================================================
# 3. LOAD GSEA RESULTS
# ============================================================

cat("\n3. Loading GSEA results...\n")

# Try different possible GSEA result files
possible_gsea_files <- c(
  "gsea_results/tables/fgsea_all_results.csv",
  "gsea_results/fgsea_results.csv",
  "fgsea_results.csv"
)

gsea_file <- NULL
for (file in possible_gsea_files) {
  if (file.exists(file)) {
    gsea_file <- file
    break
  }
}

if (is.null(gsea_file)) {
  cat("❌ ERROR: No GSEA results file found.\n")
  cat("  Run GSEA analysis first.\n")
  quit(status = 1)
}

cat("  ✓ Loading from:", gsea_file, "\n")

gsea_results <- read.csv(gsea_file, stringsAsFactors = FALSE)

# Check if leading edge column exists
if (!"leadingEdge" %in% colnames(gsea_results)) {
  cat("❌ ERROR: GSEA results missing 'leadingEdge' column\n")
  quit(status = 1)
}

# Convert leading edge to list
gsea_results$leadingEdge_list <- strsplit(gsea_results$leadingEdge, ";")

cat("  ✓ Loaded", nrow(gsea_results), "pathway results\n")
cat("  ✓ Contrasts:", paste(unique(gsea_results$Contrast), collapse = ", "), "\n")

# ============================================================
# 4. MAP LEADING EDGE GENES TO PROTEINS
# ============================================================

cat("\n4. Mapping leading edge genes to proteins...\n")

# Create output directory
output_dir <- "gsea_results/protein_level"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Process each pathway
protein_results <- list()
pathway_count <- 0
total_proteins <- 0

for (i in 1:nrow(gsea_results)) {
  pathway <- gsea_results[i, ]
  leading_genes <- pathway$leadingEdge_list[[1]]
  
  # Skip if no leading edge genes
  if (length(leading_genes) == 0) next
  
  # Find proteins for these genes
  pathway_proteins <- character(0)
  protein_gene_pairs <- list()
  
  for (gene in leading_genes) {
    if (gene %in% names(gene_to_proteins)) {
      proteins <- gene_to_proteins[[gene]]
      pathway_proteins <- c(pathway_proteins, proteins)
      
      # Track which gene each protein came from
      for (protein in proteins) {
        protein_gene_pairs[[protein]] <- gene
      }
    }
  }
  
  # Skip if no proteins found
  if (length(pathway_proteins) == 0) next
  
  # Create result for this pathway
  pathway_result <- data.frame(
    pathway_id = pathway$pathway,
    pathway_name = pathway$pathway,
    contrast = pathway$Contrast,
    NES = pathway$NES,
    pval = pathway$pval,
    padj = ifelse("padj" %in% colnames(pathway), pathway$padj, NA),
    leading_edge_genes = length(leading_genes),
    leading_edge_proteins = length(pathway_proteins),
    protein_id = pathway_proteins,
    corresponding_gene = sapply(pathway_proteins, function(p) protein_gene_pairs[[p]]),
    stringsAsFactors = FALSE
  )
  
  protein_results[[i]] <- pathway_result
  pathway_count <- pathway_count + 1
  total_proteins <- total_proteins + length(pathway_proteins)
}

# Check if we found any mappings
if (pathway_count == 0) {
  cat("❌ ERROR: Could not map any leading edge genes to proteins\n")
  cat("  Check that your gene symbols match between GSEA and mapping.\n")
  cat("  Example GSEA gene:", gsea_results$leadingEdge_list[[1]][1], "\n")
  cat("  Example mapping genes:", paste(names(gene_to_proteins)[1:3], collapse = ", "), "\n")
  quit(status = 1)
}

# Combine all pathways
protein_level_results <- do.call(rbind, protein_results)

cat("  ✓ Mapped", pathway_count, "pathways to proteins\n")
cat("  ✓ Total protein-pathway associations:", total_proteins, "\n")
cat("  ✓ Unique proteins:", length(unique(protein_level_results$protein_id)), "\n")
cat("  ✓ Unique genes (from proteins):", length(unique(protein_level_results$corresponding_gene)), "\n")

# ============================================================
# 5. CREATE PROTEIN-CENTRIC SUMMARIES
# ============================================================

cat("\n5. Creating protein-centric summaries...\n")

# 5.1 Protein participation summary
protein_summary <- protein_level_results %>%
  group_by(protein_id, corresponding_gene) %>%
  summarise(
    n_pathways = n_distinct(pathway_id),
    pathways = paste(unique(pathway_id), collapse = ";"),
    n_contrasts = n_distinct(contrast),
    contrasts = paste(unique(contrast), collapse = ";"),
    avg_NES = mean(NES, na.rm = TRUE),
    min_pval = min(pval, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_pathways), desc(n_contrasts))

write.csv(protein_summary,
          file.path(output_dir, "protein_participation_summary.csv"),
          row.names = FALSE)

# 5.2 Per-contrast results
for (contrast_name in unique(protein_level_results$contrast)) {
  contrast_data <- protein_level_results[protein_level_results$contrast == contrast_name, ]
  
  contrast_summary <- contrast_data %>%
    group_by(protein_id, corresponding_gene) %>%
    summarise(
      n_pathways = n_distinct(pathway_id),
      pathways = paste(unique(pathway_id), collapse = ";"),
      avg_NES = mean(NES),
      min_pval = min(pval),
      .groups = "drop"
    ) %>%
    arrange(desc(n_pathways))
  
  safe_name <- gsub("[^A-Za-z0-9]", "_", contrast_name)
  write.csv(contrast_summary,
            file.path(output_dir, paste0("proteins_", safe_name, ".csv")),
            row.names = FALSE)
}

# 5.3 Save complete mapping
write.csv(protein_level_results,
          file.path(output_dir, "complete_protein_pathway_mapping.csv"),
          row.names = FALSE)

cat("  ✓ Saved all results to:", output_dir, "\n")

# ============================================================
# 6. GENERATE ANALYSIS REPORT
# ============================================================

cat("\n6. Generating analysis report...\n")

report_file <- file.path(output_dir, "protein_leading_edge_report.txt")

sink(report_file)
cat("PROTEIN-LEVEL LEADING EDGE ANALYSIS REPORT\n")
cat("==========================================\n")
cat("Date:", Sys.time(), "\n")
cat("GSEA results:", gsea_file, "\n")
cat("Protein mapping:", mapping_file, "\n\n")

cat("ANALYSIS OVERVIEW\n")
cat("----------------\n")
cat("Total GSEA pathways analyzed:", nrow(gsea_results), "\n")
cat("Pathways successfully mapped to proteins:", pathway_count, "\n")
cat("Total protein-pathway associations:", total_proteins, "\n")
cat("Unique proteins identified:", length(unique(protein_level_results$protein_id)), "\n")
cat("Mapping coverage:", round(pathway_count/nrow(gsea_results)*100, 1), "% of pathways\n\n")

cat("TOP PROTEINS (MOST PATHWAYS):\n")
cat("-----------------------------\n")
top_proteins <- head(protein_summary, 10)
for (i in 1:nrow(top_proteins)) {
  p <- top_proteins[i, ]
  cat(i, ". ", p$protein_id, " (", p$corresponding_gene, ")\n", sep = "")
  cat("    Pathways: ", p$n_pathways, "\n", sep = "")
  cat("    Contrasts: ", p$n_contrasts, "\n", sep = "")
}

cat("\nMAPPING STATISTICS:\n")
cat("------------------\n")
cat("Proteins appearing in 1 pathway:", sum(protein_summary$n_pathways == 1), "\n")
cat("Proteins appearing in 2-5 pathways:", sum(protein_summary$n_pathways >= 2 & protein_summary$n_pathways <= 5), "\n")
cat("Proteins appearing in >5 pathways:", sum(protein_summary$n_pathways > 5), "\n")

# Find multi-protein genes
multi_protein_genes <- protein_summary %>%
  group_by(corresponding_gene) %>%
  summarise(n_proteins = n()) %>%
  filter(n_proteins > 1) %>%
  arrange(desc(n_proteins))

if (nrow(multi_protein_genes) > 0) {
  cat("\nGENES WITH MULTIPLE PROTEIN ISOFORMS:\n")
  cat("-------------------------------------\n")
  for (i in 1:min(5, nrow(multi_protein_genes))) {
    gene <- multi_protein_genes[i, ]
    cat(gene$corresponding_gene, ": ", gene$n_proteins, " protein isoforms\n", sep = "")
  }
}

cat("\nOUTPUT FILES:\n")
cat("------------\n")
cat("1. protein_participation_summary.csv - All proteins across all contrasts\n")
cat("2. proteins_*.csv - Per-contrast protein results\n")
cat("3. complete_protein_pathway_mapping.csv - Complete mapping\n")
cat("4. This report\n")

cat("\nVALIDATION RECOMMENDATIONS:\n")
cat("--------------------------\n")
cat("1. Prioritize proteins appearing in multiple pathways\n")
cat("2. Focus on proteins with high |NES| values\n")
cat("3. Consider multiple isoforms for genes with several proteins\n")
cat("4. Use protein IDs (not gene names) for experimental design\n")

sink()

cat("  ✓ Report saved:", report_file, "\n")

# ============================================================
# 7. FINAL SUMMARY
# ============================================================

cat("\n", paste0(rep("=", 70), collapse = ""), "\n")
cat("✅ PROTEIN-LEVEL LEADING EDGE ANALYSIS COMPLETE\n")
cat(paste0(rep("=", 70), collapse = ""), "\n\n")

cat("Your analysis pipeline:\n")
cat("1. Proteins → Gene symbols (conversion script)\n")
cat("2. Gene symbols → GSEA pathway analysis\n")
cat("3. Leading edge genes → Original proteins (this script)\n\n")

cat("Results available in: gsea_results/protein_level/\n")
cat("\nKey files:\n")
cat("• protein_participation_summary.csv - Top proteins to validate\n")
cat("• proteins_[CONTRAST].csv - Protein results per comparison\n")
cat("• complete_protein_pathway_mapping.csv - Full traceability\n")
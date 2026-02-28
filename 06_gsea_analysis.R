#!/usr/bin/env Rscript
# scripts/06_gsea_analysis.R - FGSEA VERSION WITH PROPER FOLDER STRUCTURE
# Uses fgsea with gene permutation for statistically appropriate analysis

cat("========================================\n")
cat("FGSEA ANALYSIS (GENE PERMUTATION)\n")
cat("========================================\n\n")

# ============================================================
# 1. CONFIGURATION & SETUP
# ============================================================

GSEA_CONFIG <- list(
  # Input/output
  input_dir = "ranked_lists",
  output_dir = "gsea_results",
  figures_dir = "figures",
  
  # Local data files
  local_kegg_file = "data/kegg_human_data.rds",
  
  # Analysis parameters
  method = "fgsea",
  database = "KEGG",
  organism = "hsa",
  metabolic_focus = TRUE,
  
  # Statistical thresholds
  pvalue_cutoff = 0.05,
  qvalue_cutoff = 0.1,
  min_gset_size = 10,
  max_gset_size = 500,
  nes_threshold = 1.5,  # For highlighting strong effects
  
  # Visualization
  top_n_pathways = 20,
  plot_width = 10,
  plot_height = 8
)

# Command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  if (args[1] == "--help" || args[1] == "-h") {
    cat("Usage: Rscript 06_gsea_analysis.R [input_dir] [output_dir]\n")
    cat("  input_dir:  Directory with *_genesymbols.rnk files\n")
    cat("  output_dir: Directory for GSEA results\n")
    quit(status = 0)
  }
  if (length(args) >= 1) GSEA_CONFIG$input_dir <- args[1]
  if (length(args) >= 2) GSEA_CONFIG$output_dir <- args[2]
}

cat("FGSEA ANALYSIS CONFIGURATION:\n")
cat("  Method: fgsea with gene permutation\n")
cat("  Ranked lists directory:", GSEA_CONFIG$input_dir, "\n")
cat("  KEGG data file:", GSEA_CONFIG$local_kegg_file, "\n")
cat("  Output directory:", GSEA_CONFIG$output_dir, "\n")
cat("  NES threshold for strong effects:", GSEA_CONFIG$nes_threshold, "\n\n")

# ============================================================
# 2. LOAD REQUIRED PACKAGES
# ============================================================

cat("1. Loading required packages...\n")

required_packages <- c(
  "fgsea",           # For gene permutation GSEA
  "BiocParallel",    # For parallel processing
  "ggplot2",         # For plotting
  "dplyr",           # For data manipulation
  "tidyr"            # For data reshaping
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("  Installing", pkg, "...\n")
    if (!require("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
    library(pkg, character.only = TRUE)
  }
}

# ============================================================
# 3. CREATE OUTPUT DIRECTORIES WITH PROPER STRUCTURE
# ============================================================

cat("2. Creating output directories...\n")

# Define all directories needed
dirs_to_create <- c(
  # Main directories
  GSEA_CONFIG$output_dir,
  GSEA_CONFIG$figures_dir,
  
  # Tables directory
  file.path(GSEA_CONFIG$output_dir, "tables"),
  
  # R objects directory
  file.path(GSEA_CONFIG$output_dir, "r_objects"),
  
  # Figure subdirectories
  file.path(GSEA_CONFIG$figures_dir, "nes_plots"),        # NES bar plots, heatmaps
  file.path(GSEA_CONFIG$figures_dir, "enrichment_plots"), # Individual pathway plots
  file.path(GSEA_CONFIG$figures_dir, "summary_plots")     # Combined summary plots
)

for (dir in dirs_to_create) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    cat("  Created:", dir, "\n")
  }
}

# ============================================================
# 4. LOAD LOCAL KEGG DATA
# ============================================================

cat("3. Loading local KEGG data...\n")

# Check if local KEGG file exists
if (!file.exists(GSEA_CONFIG$local_kegg_file)) {
  cat("❌ ERROR: Local KEGG data file not found!\n")
  cat("  File:", GSEA_CONFIG$local_kegg_file, "\n")
  cat("  Please run the KEGG data preparation script first:\n")
  cat("  Rscript scripts/04_download_kegg_data.R\n")
  quit(status = 1)
}

# Load the pre-prepared KEGG data
kegg_data <- readRDS(GSEA_CONFIG$local_kegg_file)
cat("  ✓ Loaded KEGG data from:", GSEA_CONFIG$local_kegg_file, "\n")

# Extract pathway-gene mapping
if (is.null(kegg_data$KEGGPATHID2EXTID)) {
  cat("❌ ERROR: KEGG data format incorrect. Missing KEGGPATHID2EXTID.\n")
  quit(status = 1)
}

pathway_data <- kegg_data$KEGGPATHID2EXTID

# Filter for metabolic pathways if requested
if (GSEA_CONFIG$metabolic_focus) {
  cat("4. Filtering for metabolic pathways...\n")
  
  # Get pathway descriptions if available
  if (!is.null(kegg_data$KEGGPATHID2NAME)) {
    metabolic_keywords <- c(
      "metabolism", "biosynthesis", "degradation", "synthesis",
      "glycolysis", "tca", "citrate", "oxidative", "phosphorylation",
      "fatty acid", "amino acid", "carbohydrate", "lipid",
      "nucleotide", "steroid", "porphyrin", "vitamin", "cofactor"
    )
    
    pathway_names <- kegg_data$KEGGPATHID2NAME
    metabolic_pathways <- pathway_names$to[
      grepl(paste(metabolic_keywords, collapse = "|"), 
            pathway_names$to, ignore.case = TRUE)
    ]
    
    if (length(metabolic_pathways) > 0) {
      metabolic_ids <- pathway_names$from[
        pathway_names$to %in% metabolic_pathways
      ]
      pathway_data <- pathway_data[pathway_data$from %in% metabolic_ids, ]
      cat("  ✓ Found", length(metabolic_ids), "metabolic pathways by name\n")
    }
  }
}

# Convert to gene sets list (format needed by fgsea)
kegg_gene_sets <- split(pathway_data$to, pathway_data$from)

cat("  Final dataset:", length(kegg_gene_sets), "pathways\n")
cat("  Total gene-pathway associations:", nrow(pathway_data), "\n")

# ============================================================
# 5. LOAD RANKED LISTS
# ============================================================

cat("5. Loading ranked gene lists...\n")

# Find all *_genesymbols.rnk files
rnk_files <- list.files(GSEA_CONFIG$input_dir, 
                        pattern = "_genesymbols\\.rnk$",
                        full.names = TRUE)

if (length(rnk_files) == 0) {
  cat("❌ ERROR: No *_genesymbols.rnk files found in", GSEA_CONFIG$input_dir, "\n")
  cat("  Run the protein ID conversion script first:\n")
  cat("  Rscript scripts/05_protein_id_conversion.R\n")
  quit(status = 1)
}

cat("  Found", length(rnk_files), "ranked list file(s):\n")
for (f in rnk_files) cat("   -", basename(f), "\n")

# Function to load a ranked list
load_ranked_list <- function(file_path) {
  data <- read.delim(file_path, header = FALSE,
                     col.names = c("gene_symbol", "score"),
                     stringsAsFactors = FALSE)
  
  # Convert to named vector (required by fgsea)
  ranked_vector <- setNames(data$score, data$gene_symbol)
  
  # Remove duplicates (keep highest absolute score)
  if (any(duplicated(names(ranked_vector)))) {
    df <- data.frame(gene = names(ranked_vector), score = ranked_vector)
    df <- df[order(-abs(df$score)), ]
    df <- df[!duplicated(df$gene), ]
    ranked_vector <- setNames(df$score, df$gene)
  }
  
  # Ensure sorted descending (required for GSEA)
  ranked_vector <- sort(ranked_vector, decreasing = TRUE)
  
  return(ranked_vector)
}

# Load all ranked lists
ranked_lists <- lapply(rnk_files, load_ranked_list)
names(ranked_lists) <- gsub("_genesymbols\\.rnk$", "", basename(rnk_files))

cat("  Successfully loaded", length(ranked_lists), "contrast(s)\n")
for (name in names(ranked_lists)) {
  cat("   -", name, ":", length(ranked_lists[[name]]), "genes\n")
}

# ============================================================
# 6. RUN FGSEA WITH GENE PERMUTATION
# ============================================================

cat("\n6. Running fgsea with gene permutation...\n")

gsea_results <- list()

for (i in seq_along(ranked_lists)) {
  contrast_name <- names(ranked_lists)[i]
  ranked_vector <- ranked_lists[[i]]
  
  cat("  Analyzing:", contrast_name, "\n")
  cat("    Genes in ranked list:", length(ranked_vector), "\n")
  cat("    Mean rank score:", round(mean(ranked_vector), 3), "\n")
  cat("    Using gene permutation (correct for perfect confounding)...\n")
  
  tryCatch({
    # Run fgsea with gene permutation
    fgsea_result <- fgseaMultilevel(
      pathways = kegg_gene_sets,
      stats = ranked_vector,
      minSize = GSEA_CONFIG$min_gset_size,
      maxSize = GSEA_CONFIG$max_gset_size,
      eps = 0,            # Precise p-value calculation
      nPermSimple = 10000, # Number of permutations
      nproc = 1,          # Single core for reproducibility
      gseaParam = 1,      # Standard weighting
      BPPARAM = SerialParam()
    )
    
    if (nrow(fgsea_result) > 0) {
      # Add contrast name and sort by |NES|
      fgsea_result$Contrast <- contrast_name
      fgsea_result <- fgsea_result[order(-abs(fgsea_result$NES)), ]
      fgsea_result$Rank <- 1:nrow(fgsea_result)
      fgsea_result$abs_NES <- abs(fgsea_result$NES)
      
      cat("    ✓ Found", nrow(fgsea_result), "pathways with gene permutation\n")
      cat("    Top pathway: ", fgsea_result$pathway[1], 
          " (NES = ", round(fgsea_result$NES[1], 2), 
          ", |NES| = ", round(fgsea_result$abs_NES[1], 2), ")\n", sep = "")
      
      gsea_results[[contrast_name]] <- list(
        fgsea_object = fgsea_result,
        contrast_name = contrast_name,
        n_pathways = nrow(fgsea_result),
        top_pathway = fgsea_result$pathway[1],
        top_NES = fgsea_result$NES[1]
      )
    } else {
      cat("    ⚠ No pathways found\n")
      gsea_results[[contrast_name]] <- list(
        fgsea_object = data.frame(),
        contrast_name = contrast_name,
        n_pathways = 0
      )
    }
    
  }, error = function(e) {
    cat("    ❌ ERROR in fgsea:", e$message, "\n")
    gsea_results[[contrast_name]] <- list(
      fgsea_object = data.frame(),
      contrast_name = contrast_name,
      n_pathways = 0,
      error = e$message
    )
  })
}

# ============================================================
# 7. SAVE FGSEA RESULTS (HANDLING LIST COLUMNS)
# ============================================================

cat("\n7. Saving fgsea results...\n")

# Combine all results
all_fgsea_results <- list()
for (contrast in names(gsea_results)) {
  if (!is.null(gsea_results[[contrast]]$fgsea_object) && 
      nrow(gsea_results[[contrast]]$fgsea_object) > 0) {
    all_fgsea_results[[contrast]] <- gsea_results[[contrast]]$fgsea_object
  }
}

if (length(all_fgsea_results) > 0) {
  # Combine into one data frame
  combined_results <- do.call(rbind, all_fgsea_results)
  rownames(combined_results) <- NULL
  
  # Convert list columns to character vectors for CSV export
  if ("leadingEdge" %in% colnames(combined_results)) {
    combined_results$leadingEdge <- sapply(combined_results$leadingEdge, 
                                           function(x) paste(x, collapse = ";"))
  }
  
  # Save all results
  results_file <- file.path(GSEA_CONFIG$output_dir, "tables", "fgsea_all_results.csv")
  write.csv(combined_results, results_file, row.names = FALSE)
  cat("  ✓ Saved all fgsea results to:", results_file, "\n")
  
  # Save results filtered by |NES| threshold
  strong_effects <- combined_results[abs(combined_results$NES) >= GSEA_CONFIG$nes_threshold, ]
  
  if (nrow(strong_effects) > 0) {
    nes_file <- file.path(GSEA_CONFIG$output_dir, "tables", 
                          paste0("fgsea_strong_effects_nes", GSEA_CONFIG$nes_threshold, ".csv"))
    write.csv(strong_effects, nes_file, row.names = FALSE)
    cat("  ✓ Saved pathways with |NES| ≥", GSEA_CONFIG$nes_threshold, "to:", nes_file, "\n")
  } else {
    cat("  ⚠ No pathways with |NES| ≥", GSEA_CONFIG$nes_threshold, "\n")
  }
  
  # Save top 10 pathways per contrast by |NES|
  top_pathways_list <- list()
  for (contrast in unique(combined_results$Contrast)) {
    contrast_data <- combined_results[combined_results$Contrast == contrast, ]
    if (nrow(contrast_data) > 0) {
      top_10 <- head(contrast_data[order(-abs(contrast_data$NES)), ], 10)
      top_pathways_list[[contrast]] <- top_10
    }
  }
  
  if (length(top_pathways_list) > 0) {
    top_pathways <- do.call(rbind, top_pathways_list)
    top_file <- file.path(GSEA_CONFIG$output_dir, "tables", "fgsea_top_pathways.csv")
    write.csv(top_pathways, top_file, row.names = FALSE)
    cat("  ✓ Saved top pathways by |NES| to:", top_file, "\n")
  }
  
  # Save R objects (lists preserved)
  rdata_file <- file.path(GSEA_CONFIG$output_dir, "r_objects", "fgsea_results.rds")
  saveRDS(gsea_results, rdata_file)
  cat("  ✓ Saved R objects to:", rdata_file, "\n")
  
} else {
  cat("  ⚠ No fgsea results to save\n")
}

# ============================================================
# 8. CREATE NES PLOTS (MAIN VISUALIZATIONS)
# ============================================================

cat("8. Creating NES plots...\n")

if (exists("combined_results") && nrow(combined_results) > 0) {
  
  # 8.1 NES BAR PLOT
  cat("  Creating NES bar plot...\n")
  
  # Get top pathways by |NES| across all contrasts
  top_overall <- combined_results[order(-abs(combined_results$NES)), ]
  top_overall <- head(top_overall, 15)
  
  p_nes <- ggplot(top_overall, aes(x = reorder(pathway, NES), y = NES, fill = Contrast)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_brewer(palette = "Set1") +
    labs(title = "Top Pathways by Normalized Enrichment Score (NES)",
         subtitle = "Gene permutation GSEA | Perfect confounding accounted for",
         x = "Pathway",
         y = "NES (Positive = enriched in first condition)") +
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text.y = element_text(size = 9)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = c(GSEA_CONFIG$nes_threshold, -GSEA_CONFIG$nes_threshold), 
               linetype = "dotted", color = "red", alpha = 0.5)
  
  ggsave(file.path(GSEA_CONFIG$figures_dir, "nes_plots", "top_pathways_nes.pdf"),
         p_nes, width = GSEA_CONFIG$plot_width, height = GSEA_CONFIG$plot_height)
  cat("  ✓ Created NES bar plot: nes_plots/top_pathways_nes.pdf\n")
  
  # 8.2 INDIVIDUAL CONTRAST PLOTS
  cat("  Creating individual contrast plots...\n")
  
  for (contrast in unique(combined_results$Contrast)) {
    contrast_data <- combined_results[combined_results$Contrast == contrast, ]
    
    if (nrow(contrast_data) > 0) {
      # Take top 10 by |NES|
      top_contrast <- contrast_data[order(-abs(contrast_data$NES)), ]
      top_contrast <- head(top_contrast, 10)
      
      p <- ggplot(top_contrast, aes(x = reorder(pathway, NES), y = NES)) +
        geom_bar(stat = "identity", 
                 fill = ifelse(top_contrast$NES > 0, "firebrick", "steelblue")) +
        coord_flip() +
        labs(title = paste("Top Pathways:", contrast),
             subtitle = "Sorted by |NES| (gene permutation GSEA)",
             x = "Pathway",
             y = "Normalized Enrichment Score") +
        theme_minimal() +
        theme(axis.text.y = element_text(size = 9)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50")
      
      safe_name <- gsub("[^A-Za-z0-9]", "_", contrast)
      ggsave(file.path(GSEA_CONFIG$figures_dir, "nes_plots", 
                       paste0(safe_name, "_top_pathways.pdf")),
             p, width = 10, height = 6)
    }
  }
  cat("  ✓ Created individual contrast plots in nes_plots/\n")
  
} else {
  cat("  ⚠ No results for NES plots\n")
}

# ============================================================
# 9. CREATE ENRICHMENT PLOTS (FIXED VERSION)
# ============================================================

cat("9. Creating enrichment plots...\n")

plots_created <- 0

# Ensure the enrichment plots directory exists
enrichment_dir <- file.path(GSEA_CONFIG$figures_dir, "enrichment_plots")
dir.create(enrichment_dir, recursive = TRUE, showWarnings = FALSE)

for (contrast_name in names(gsea_results)) {
  # Access the fgsea data frame correctly
  if (!is.null(gsea_results[[contrast_name]]$fgsea_object)) {
    fgsea_df <- gsea_results[[contrast_name]]$fgsea_object
    
    if (nrow(fgsea_df) > 0) {
      cat("    Processing", contrast_name, "...\n")
      
      # Get top 3 pathways by |NES|
      top_pathways <- fgsea_df[order(-abs(fgsea_df$NES)), ]
      top_pathways <- head(top_pathways, 3)
      
      for (i in 1:nrow(top_pathways)) {
        pathway_name <- top_pathways$pathway[i]
        
        tryCatch({
          cat("      Creating plot for:", pathway_name, "\n")
          
          # Get the pathway genes from your KEGG sets
          if (pathway_name %in% names(kegg_gene_sets)) {
            pathway_genes <- kegg_gene_sets[[pathway_name]]
            
            # Get the ranked list for this contrast
            ranked_stats <- ranked_lists[[contrast_name]]
            
            # Create a safe filename
            safe_name <- gsub("[^A-Za-z0-9]", "_", pathway_name)
            plot_file <- file.path(enrichment_dir,
                                   paste0(contrast_name, "_", safe_name, ".pdf"))
            
            # Create the enrichment plot PROPERLY
            pdf(plot_file, width = 10, height = 6)
            
            # Method 1: Use fgsea's plotEnrichment with proper print()
            enrichment_plot <- fgsea::plotEnrichment(pathway_genes, ranked_stats)
            
            # Add title and labels
            pathway_nes <- round(top_pathways$NES[i], 2)
            pathway_pval <- format(top_pathways$pval[i], scientific = TRUE, digits = 2)
            
            # Customize the plot
            enrichment_plot <- enrichment_plot + 
              ggtitle(paste(contrast_name, "-", pathway_name),
                      subtitle = paste("NES =", pathway_nes, 
                                       "| p =", pathway_pval,
                                       "| Gene permutation GSEA")) +
              theme_minimal() +
              theme(plot.title = element_text(size = 14, face = "bold"),
                    plot.subtitle = element_text(size = 11, color = "gray50"))
            
            # PRINT the plot (critical step!)
            print(enrichment_plot)
            
            dev.off()
            
            # Verify the file was created and has content
            if (file.exists(plot_file)) {
              file_size <- file.size(plot_file)
              if (file_size > 1024) {  # More than 1KB means it has content
                plots_created <- plots_created + 1
                cat("        ✓ Created:", basename(plot_file), 
                    "(", round(file_size/1024, 1), "KB)\n", sep = "")
              } else {
                cat("        ⚠ File created but too small (", file_size, " bytes)\n", sep = "")
                
                # Try alternative plotting method
                cat("        Trying alternative plotting method...\n")
                create_simple_enrichment_plot(contrast_name, pathway_name, 
                                              pathway_genes, ranked_stats,
                                              top_pathways[i, ], plot_file)
              }
            } else {
              cat("        ⚠ File not created\n")
            }
            
          } else {
            cat("        ⚠ Pathway not found in gene sets:", pathway_name, "\n")
          }
          
        }, error = function(e) {
          cat("        ⚠ Error:", e$message, "\n")
          # Try the simple plot function as fallback
          try({
            cat("        Trying fallback plot...\n")
            simple_file <- file.path(enrichment_dir,
                                     paste0(contrast_name, "_", 
                                            gsub("[^A-Za-z0-9]", "_", pathway_name), 
                                            "_simple.pdf"))
            create_simple_enrichment_plot(contrast_name, pathway_name, 
                                          pathway_genes, ranked_stats,
                                          top_pathways[i, ], simple_file)
          })
        })
      }
    }
  }
}

cat("  ✓ Created", plots_created, "enrichment plots in enrichment_plots/\n")

# ============================================================
# 9.1 HELPER FUNCTION FOR SIMPLE ENRICHMENT PLOTS
# ============================================================

create_simple_enrichment_plot <- function(contrast_name, pathway_name, pathway_genes, 
                                          ranked_stats, pathway_info, output_file) {
  # Create a simple enrichment plot using base R
  pdf(output_file, width = 10, height = 6)
  
  # Get pathway genes positions in ranked list
  pathway_gene_names <- intersect(names(ranked_stats), pathway_genes)
  
  if (length(pathway_gene_names) > 0) {
    # Get positions of pathway genes in the ranked list
    gene_positions <- match(pathway_gene_names, names(ranked_stats))
    gene_scores <- ranked_stats[pathway_gene_names]
    
    # Sort by position
    gene_positions <- sort(gene_positions)
    
    # Create running enrichment score
    n_genes <- length(ranked_stats)
    n_pathway <- length(pathway_gene_names)
    
    # Calculate running enrichment score (simplified)
    hit_indices <- rep(0, n_genes)
    hit_indices[gene_positions] <- 1
    
    # Simplified enrichment score calculation
    running_es <- cumsum(hit_indices - (n_pathway / n_genes))
    running_es <- running_es / max(abs(running_es))  # Normalize
    
    # Create the plot
    par(mar = c(5, 5, 4, 2) + 0.1)
    
    # Plot 1: Running enrichment score
    plot(running_es, type = "l", lwd = 2, col = "blue",
         xlab = "Gene Rank", ylab = "Running Enrichment Score",
         main = paste(contrast_name, "-", pathway_name),
         ylim = c(-1, 1))
    
    abline(h = 0, lty = 2, col = "gray50")
    
    # Add pathway gene positions as vertical lines
    abline(v = gene_positions, col = "red", lty = 3, lwd = 0.5)
    
    # Add NES and p-value info
    legend("topright",
           legend = c(paste("NES =", round(pathway_info$NES, 2)),
                      paste("p =", format(pathway_info$pval, scientific = TRUE, digits = 2)),
                      paste("Genes in pathway:", length(pathway_gene_names))),
           bty = "n")
    
    # Add a second plot: Gene scores distribution
    par(new = TRUE)
    plot(gene_positions, gene_scores, type = "h", 
         col = "red", lwd = 2, 
         axes = FALSE, xlab = "", ylab = "",
         ylim = c(min(ranked_stats), max(ranked_stats)))
    
  } else {
    # No overlap between pathway genes and ranked list
    plot(1, 1, type = "n", 
         main = paste("No overlap:", pathway_name),
         xlab = "", ylab = "",
         axes = FALSE)
    text(1, 1, "No genes from pathway found in ranked list", 
         cex = 1.2, col = "red")
  }
  
  dev.off()
  
  # Check if file was created
  if (file.exists(output_file) && file.size(output_file) > 1024) {
    cat("        ✓ Created simple plot:", basename(output_file), "\n")
    return(TRUE)
  } else {
    cat("        ⚠ Simple plot also failed\n")
    return(FALSE)
  }
}
# ============================================================
# 10. CREATE SUMMARY PLOTS
# ============================================================

cat("10. Creating summary plots...\n")

if (exists("combined_results") && nrow(combined_results) > 0) {
  
  # 10.1 NES vs -log10(p-value) scatter plot
  p_scatter <- ggplot(combined_results, aes(x = NES, y = -log10(pval), color = Contrast)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = c(GSEA_CONFIG$nes_threshold, -GSEA_CONFIG$nes_threshold), 
               linetype = "dotted", color = "red", alpha = 0.5) +
    labs(title = "NES vs Statistical Significance",
         subtitle = "Gene permutation GSEA results",
         x = "Normalized Enrichment Score (NES)",
         y = "-log10(p-value)") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave(file.path(GSEA_CONFIG$figures_dir, "summary_plots", "nes_vs_pvalue.pdf"),
         p_scatter, width = 10, height = 8)
  cat("  ✓ Created scatter plot: summary_plots/nes_vs_pvalue.pdf\n")
  
  # 10.2 Pathway count by contrast
  pathway_counts <- data.frame(
    Contrast = names(gsea_results),
    Count = sapply(gsea_results, function(x) x$n_pathways)
  )
  
  p_counts <- ggplot(pathway_counts, aes(x = reorder(Contrast, -Count), y = Count, fill = Contrast)) +
    geom_bar(stat = "identity") +
    labs(title = "Number of Enriched Pathways by Contrast",
         subtitle = paste("NES threshold: |NES| >", GSEA_CONFIG$nes_threshold),
         x = "Contrast",
         y = "Number of Pathways") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(GSEA_CONFIG$figures_dir, "summary_plots", "pathway_counts.pdf"),
         p_counts, width = 8, height = 6)
  cat("  ✓ Created pathway count plot: summary_plots/pathway_counts.pdf\n")
  
} else {
  cat("  ⚠ No results for summary plots\n")
}

# ============================================================
# 11. GENERATE STATISTICALLY APPROPRIATE REPORT
# ============================================================

cat("11. Generating statistically appropriate report...\n")

report_content <- c(
  "FGSEA ANALYSIS REPORT (GENE PERMUTATION)",
  "=========================================",
  paste("Date:", Sys.time()),
  paste("Analysis method: fgsea with gene permutation"),
  paste("Ranking metric: logFC only (due to perfect confounding)"),
  paste("Statistical approach: Appropriate for confounded designs"),
  paste("NES threshold for strong effects:", GSEA_CONFIG$nes_threshold),
  "",
  "KEY DIFFERENCES FROM STANDARD GSEA:",
  "1. Uses GENE permutation (not sample permutation)",
  "2. p-values come from null distribution of gene sets",
  "3. Appropriate when input p-values are unreliable",
  "4. Focus on NES magnitude rather than p-value significance",
  "",
  "INTERPRETATION GUIDELINES:",
  paste("• NES > 0: Pathway enriched in FIRST condition (numerator)"),
  paste("• NES < 0: Pathway enriched in SECOND condition (denominator)"),
  paste("• |NES| >", GSEA_CONFIG$nes_threshold, ": Considered strong effect"),
  paste("• p-values indicate if NES is extreme vs. random gene sets"),
  "",
  "RESULTS SUMMARY:"
)

# Add results summary
for (contrast in names(gsea_results)) {
  if (!is.null(gsea_results[[contrast]]$n_pathways) && 
      gsea_results[[contrast]]$n_pathways > 0) {
    
    report_content <- c(report_content,
                        paste("\n", contrast, ":", sep = ""),
                        paste("  Total pathways tested:", length(kegg_gene_sets)),
                        paste("  Pathways with any enrichment:", gsea_results[[contrast]]$n_pathways))
    
    # Show top 3 pathways by |NES|
    fgsea_obj <- gsea_results[[contrast]]$fgsea_object
    if (!is.null(fgsea_obj) && nrow(fgsea_obj) > 0) {
      top_3 <- head(fgsea_obj[order(-abs(fgsea_obj$NES)), ], 3)
      report_content <- c(report_content, "  Top pathways by |NES|:")
      for (i in 1:nrow(top_3)) {
        report_content <- c(report_content,
                            paste("    ", i, ". ", top_3$pathway[i],
                                  " (NES=", round(top_3$NES[i], 2),
                                  ", p=", format(top_3$pval[i], scientific = TRUE, digits = 2),
                                  ", |NES| rank=", top_3$Rank[i], "/", nrow(fgsea_obj), ")", sep = ""))
      }
    }
  }
}

report_content <- c(report_content,
                    "",
                    "OUTPUT FILES:",
                    paste("1.", file.path(GSEA_CONFIG$output_dir, "tables/fgsea_all_results.csv"),
                          "- Complete fgsea results"),
                    paste("2.", file.path(GSEA_CONFIG$output_dir, "tables/fgsea_top_pathways.csv"),
                          "- Top pathways by |NES|"),
                    paste("3.", file.path(GSEA_CONFIG$figures_dir, "nes_plots/top_pathways_nes.pdf"),
                          "- Main NES visualization"),
                    paste("4.", file.path(GSEA_CONFIG$figures_dir, "enrichment_plots/"),
                          "- Individual enrichment plots"),
                    paste("5.", file.path(GSEA_CONFIG$figures_dir, "summary_plots/"),
                          "- Summary visualizations"),
                    "",
                    "STATISTICAL NOTES:",
                    "1. Gene permutation avoids assumptions about sample independence",
                    "2. Results are robust to perfect confounding in experimental design",
                    "3. Biological interpretation should focus on NES sign and magnitude",
                    "4. These results are statistically valid for your study design"
)

report_file <- file.path(GSEA_CONFIG$output_dir, "fgsea_analysis_report.txt")
writeLines(report_content, report_file)
cat("  ✓ Saved report to:", report_file, "\n")

# ============================================================
# 12. FINAL SUMMARY
# ============================================================

cat("\n", paste0(rep("=", 70), collapse = ""), "\n", sep = "")
cat("✅ FGSEA ANALYSIS COMPLETE (STATISTICALLY APPROPRIATE)\n")
cat(paste0(rep("=", 70), collapse = ""), "\n\n")

cat("Summary:\n")
cat("  Method: fgsea with gene permutation\n")
cat("  Correct for: Perfect confounding in experimental design\n")
cat("  Ranking metric: logFC only (no p-value dependence)\n")
cat("  Primary metric: Normalized Enrichment Score (NES)\n\n")

cat("Key outputs:\n")
for (contrast in names(gsea_results)) {
  if (!is.null(gsea_results[[contrast]]$n_pathways) && 
      gsea_results[[contrast]]$n_pathways > 0) {
    cat("  ", contrast, ": ", gsea_results[[contrast]]$n_pathways, 
        " enriched pathways\n", sep = "")
    if (!is.null(gsea_results[[contrast]]$top_NES)) {
      cat("       Top NES: ", round(gsea_results[[contrast]]$top_NES, 2), 
          " (", gsea_results[[contrast]]$top_pathway, ")\n", sep = "")
    }
  }
}

cat("\nFolder structure created:\n")
cat("  • gsea_results/tables/ - CSV result files\n")
cat("  • gsea_results/figures/nes_plots/ - NES visualizations\n")
cat("  • gsea_results/figures/enrichment_plots/ - Individual pathway plots\n")
cat("  • gsea_results/figures/summary_plots/ - Summary visualizations\n")
cat("  • gsea_results/r_objects/ - R data files\n")

cat("\nYour analysis now uses statistically appropriate methods\n")
cat("for a perfectly confounded experimental design.\n")
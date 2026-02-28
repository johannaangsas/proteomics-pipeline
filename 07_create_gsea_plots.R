#!/usr/bin/env Rscript
# scripts/07_create_gsea_plots.R
# Creates all plots from existing fgsea results

cat("========================================\n")
cat("CREATING GSEA PLOTS FROM EXISTING RESULTS\n")
cat("========================================\n\n")

# ============================================================
# 1. SETUP
# ============================================================

# Load required packages
library(fgsea)
library(ggplot2)
library(dplyr)

# Configuration
output_dir <- "gsea_results"
figures_dir <- "gsea_results/figures"

# Create all necessary directories
plot_dirs <- c(
  file.path(figures_dir, "nes_plots"),
  file.path(figures_dir, "enrichment_plots"),
  file.path(figures_dir, "summary_plots")
)

cat("1. Creating plot directories...\n")
for (dir in plot_dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    cat("  Created:", dir, "\n")
  }
}

# ============================================================
# 2. LOAD EXISTING DATA
# ============================================================

cat("2. Loading existing data...\n")

# Load fgsea results
results_file <- file.path(output_dir, "tables", "fgsea_all_results.csv")
if (!file.exists(results_file)) {
  cat("❌ ERROR: Results file not found:", results_file, "\n")
  quit(status = 1)
}

combined_results <- read.csv(results_file)
cat("  ✓ Loaded", nrow(combined_results), "results from", results_file, "\n")

# Load ranked lists (needed for enrichment plots)
cat("3. Loading ranked lists...\n")
ranked_dir <- "ranked_lists"
rnk_files <- list.files(ranked_dir, pattern = "_genesymbols\\.rnk$", full.names = TRUE)

ranked_lists <- list()
for (f in rnk_files) {
  contrast_name <- gsub("_genesymbols\\.rnk$", "", basename(f))
  data <- read.delim(f, header = FALSE, col.names = c("gene_symbol", "score"))
  ranked_vector <- setNames(data$score, data$gene_symbol)
  ranked_vector <- sort(ranked_vector, decreasing = TRUE)
  ranked_lists[[contrast_name]] <- ranked_vector
}
cat("  ✓ Loaded", length(ranked_lists), "ranked lists\n")

# Load KEGG gene sets
cat("4. Loading KEGG gene sets...\n")
kegg_data <- readRDS("data/kegg_human_data.rds")
pathway_data <- kegg_data$KEGGPATHID2EXTID
kegg_gene_sets <- split(pathway_data$to, pathway_data$from)
cat("  ✓ Loaded", length(kegg_gene_sets), "KEGG pathways\n")

# ============================================================
# 3. CREATE NES PLOTS
# ============================================================

cat("5. Creating NES plots...\n")

if (nrow(combined_results) > 0) {
  
  # Get top pathways by |NES|
  top_overall <- combined_results[order(-abs(combined_results$NES)), ]
  top_overall <- head(top_overall, 15)
  
  # Create NES bar plot
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
    geom_hline(yintercept = c(1.5, -1.5), linetype = "dotted", color = "red", alpha = 0.5)
  
  ggsave(file.path(figures_dir, "nes_plots", "top_pathways_nes.pdf"),
         p_nes, width = 10, height = 8)
  cat("  ✓ Created: nes_plots/top_pathways_nes.pdf\n")
  
  # Create individual contrast plots
  for (contrast in unique(combined_results$Contrast)) {
    contrast_data <- combined_results[combined_results$Contrast == contrast, ]
    
    if (nrow(contrast_data) > 0) {
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
      ggsave(file.path(figures_dir, "nes_plots", paste0(safe_name, "_top_pathways.pdf")),
             p, width = 10, height = 6)
    }
  }
  cat("  ✓ Created individual contrast plots\n")
  
} else {
  cat("  ⚠ No results for NES plots\n")
}

# ============================================================
# 4. CREATE ENRICHMENT PLOTS
# ============================================================

cat("6. Creating enrichment plots...\n")

plots_created <- 0

for (contrast in unique(combined_results$Contrast)) {
  contrast_data <- combined_results[combined_results$Contrast == contrast, ]
  
  if (nrow(contrast_data) > 0) {
    cat("    Processing", contrast, "...\n")
    
    # Get top 3 pathways by |NES|
    top_pathways <- contrast_data[order(-abs(contrast_data$NES)), ]
    top_pathways <- head(top_pathways, 3)
    
    for (i in 1:nrow(top_pathways)) {
      pathway_name <- top_pathways$pathway[i]
      
      tryCatch({
        cat("      Creating plot for:", pathway_name, "\n")
        
        if (pathway_name %in% names(kegg_gene_sets)) {
          pathway_genes <- kegg_gene_sets[[pathway_name]]
          
          # Create safe filename
          safe_name <- gsub("[^A-Za-z0-9]", "_", pathway_name)
          plot_file <- file.path(figures_dir, "enrichment_plots",
                                 paste0(contrast, "_", safe_name, ".pdf"))
          
          # Create PDF
          pdf(plot_file, width = 10, height = 6)
          
          # Create enrichment plot
          p <- plotEnrichment(pathway_genes, ranked_lists[[contrast]]) +
            ggtitle(paste(contrast, "-", pathway_name),
                    subtitle = paste("NES =", round(top_pathways$NES[i], 2),
                                     "| p =", format(top_pathways$pval[i], scientific = TRUE, digits = 2))) +
            theme_minimal()
          
          print(p)  # Critical: print the plot
          dev.off()
          
          # Check if file was created
          if (file.exists(plot_file) && file.size(plot_file) > 1024) {
            plots_created <- plots_created + 1
            cat("        ✓ Created:", basename(plot_file), "\n")
          }
        }
        
      }, error = function(e) {
        cat("        ⚠ Error:", e$message, "\n")
      })
    }
  }
}

cat("  ✓ Created", plots_created, "enrichment plots\n")

# ============================================================
# 5. FINAL SUMMARY
# ============================================================

cat("\n", paste0(rep("=", 60), collapse = ""), "\n", sep = "")
cat("✅ PLOT CREATION COMPLETE\n")
cat(paste0(rep("=", 60), collapse = ""), "\n\n")

cat("Files created in gsea_results/figures/:\n")
cat("  nes_plots/\n")
cat("    • top_pathways_nes.pdf\n")
cat("    • [contrast]_top_pathways.pdf\n")
cat("  enrichment_plots/\n")
cat("    • [contrast]_[pathway].pdf\n\n")

cat("Total plots created:", plots_created + length(unique(combined_results$Contrast)) + 1, "\n")
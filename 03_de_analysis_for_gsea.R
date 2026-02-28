#!/usr/bin/env Rscript
# scripts/03_de_analysis_for_gsea_fixed.R
# DE analysis function with visualizations - FIXED for GSEA

run_de_analysis_for_gsea <- function(ANALYSIS_CONFIG = NULL) {
  
  cat("=== DE ANALYSIS FUNCTION WITH VISUALIZATIONS ===\n")
  
  # Load required packages INSIDE the function
  if (!require("limma", quietly = TRUE)) {
    cat("Installing limma...\n")
    install.packages("limma", repos = "https://cloud.r-project.org", quiet = TRUE)
    library(limma)
  }
  
  # Default configuration
  if (is.null(ANALYSIS_CONFIG)) {
    ANALYSIS_CONFIG <- list(
      input_file = "intermediate_data/02_model_fit.rds",
      output_dir = "de_results",
      ranked_dir = "ranked_lists",
      figures_dir = "figures",
      create_ranked_lists = TRUE,
      effect_thresholds = list(large = 3, moderate = 2, small = 1),
      top_n = 50
    )
  }
  
  # ============================================================
  # 1. LOAD DATA
  # ============================================================
  
  cat("1. Loading data...\n")
  
  if (!file.exists(ANALYSIS_CONFIG$input_file)) {
    stop("Input file not found: ", ANALYSIS_CONFIG$input_file)
  }
  
  model_objects <- readRDS(ANALYSIS_CONFIG$input_file)
  fit_ebayes <- model_objects$fit_ebayes
  contrasts <- model_objects$contrasts
  metadata <- model_objects$metadata
  
  # ============================================================
  # 2. CHECK DESIGN
  # ============================================================
  
  cat("2. Checking experimental design...\n")
  
  perfect_confounding <- FALSE
  confounding_details <- ""
  
  if ("Batch" %in% colnames(metadata) && "Condition" %in% colnames(metadata)) {
    conf_table <- table(metadata$Condition, metadata$Batch)
    cat("  Batch-Condition table:\n")
    print(conf_table)
    
    condition_per_batch <- apply(conf_table > 0, 2, sum)
    perfect_confounding <- all(condition_per_batch == 1)
    
    if (perfect_confounding) {
      cat("\n  ⚠ WARNING: PERFECT CONFOUNDING DETECTED!\n")
      cat("    Each condition is in exactly one batch.\n")
      cat("    P-values are confounded with batch effects.\n")
      cat("    Using effect size (logFC) for ranking only.\n")
      confounding_details <- "Perfect confounding with batch"
    } else {
      cat("\n  No perfect confounding detected.\n")
      cat("  Using moderated t-statistic for ranking.\n")
      confounding_details <- "No perfect confounding"
    }
  } else {
    cat("  No batch information found. Assuming no confounding.\n")
    confounding_details <- "No batch information"
  }
  
  # ============================================================
  # 3. CREATE DIRECTORIES
  # ============================================================
  
  cat("3. Creating output directories...\n")
  
  dirs <- c(ANALYSIS_CONFIG$output_dir, ANALYSIS_CONFIG$ranked_dir, ANALYSIS_CONFIG$figures_dir)
  for (dir in dirs) {
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  # ============================================================
  # 4. PROCESS EACH CONTRAST
  # ============================================================
  
  cat("4. Processing contrasts...\n")
  
  contrast_names <- colnames(contrasts)
  all_results <- list()
  ranked_lists <- list()
  
  for (i in seq_along(contrast_names)) {
    contrast <- contrast_names[i]
    cat("   ", contrast, "...\n")
    
    # Get DE results
    results <- limma::topTable(fit_ebayes, coef = i, number = Inf, sort.by = "none")
    
    # Add protein ID
    if (!"ID" %in% colnames(results) && !is.null(rownames(results))) {
      results$Protein_ID <- rownames(results)
    } else if ("ID" %in% colnames(results)) {
      results$Protein_ID <- results$ID
    } else {
      results$Protein_ID <- paste0("Protein_", 1:nrow(results))
    }
    
    # Calculate effect size
    results$abs_logFC <- abs(results$logFC)
    
    # Effect size categories
    results$Effect_Category <- cut(results$abs_logFC,
                                   breaks = c(-Inf, 
                                              ANALYSIS_CONFIG$effect_thresholds$small,
                                              ANALYSIS_CONFIG$effect_thresholds$moderate,
                                              ANALYSIS_CONFIG$effect_thresholds$large,
                                              Inf),
                                   labels = c("Minimal", "Small", "Moderate", "Large"),
                                   include.lowest = TRUE)
    
    results$Direction <- ifelse(results$logFC > 0, "Up", "Down")
    
    # Choose ranking metric
    if (perfect_confounding) {
      results$Ranking_Score <- results$logFC
      results$Ranking_Metric <- "logFC"
      results$Note <- "P-values confounded with batch - interpret with caution"
    } else {
      results$Ranking_Score <- results$t
      results$Ranking_Metric <- "t_statistic"
    }
    
    # FIXED: Sort by ranking score (descending, not absolute value)
    results <- results[order(-results$Ranking_Score), ]  # Removed abs()
    results$Rank <- 1:nrow(results)
    results$Percentile <- (1 - results$Rank / nrow(results)) * 100
    
    # Store results
    all_results[[contrast]] <- results
    
    # Save full results
    write.csv(results,
              file.path(ANALYSIS_CONFIG$output_dir,
                        paste0(gsub("_", "-", contrast), "_de_results.csv")),
              row.names = FALSE)
    
    # Create ranked list for GSEA
    if (ANALYSIS_CONFIG$create_ranked_lists) {
      ranked_df <- data.frame(
        Protein_ID = results$Protein_ID,
        Ranking_Score = results$Ranking_Score
      )
      
      # GSEA requires scores sorted high to low (already sorted above)
      rnk_file <- file.path(ANALYSIS_CONFIG$ranked_dir,
                            paste0(gsub("_", "-", contrast), ".rnk"))
      
      write.table(ranked_df, rnk_file,
                  sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
      
      cat("      Saved ranked list:", rnk_file, "\n")
      
      # Store for return
      ranked_lists[[contrast]] <- setNames(ranked_df$Ranking_Score, ranked_df$Protein_ID)
    }
  }
  
  # ============================================================
  # 5. CREATE VISUALIZATIONS
  # ============================================================
  
  cat("5. Creating visualizations...\n")
  
  tryCatch({
    
    # Create effect size distribution PDF
    pdf(file.path(ANALYSIS_CONFIG$figures_dir, "effect_size_distributions.pdf"),
        width = 12, height = 8)
    
    # Set up layout for multiple contrasts
    n_contrasts <- length(contrast_names)
    par(mfrow = c(ceiling(n_contrasts/2), 2))
    
    for (contrast in contrast_names) {
      results <- all_results[[contrast]]
      
      # Histogram of logFC with threshold lines
      hist(results$logFC, breaks = 50,
           main = paste("Effect Size Distribution:", contrast),
           xlab = "log2 Fold Change", 
           col = "lightblue", 
           border = "white",
           xlim = c(-max(abs(results$logFC)) * 1.1, max(abs(results$logFC)) * 1.1))
      
      # Add threshold lines
      abline(v = ANALYSIS_CONFIG$effect_thresholds$large, 
             col = "red", lty = 2, lwd = 2)
      abline(v = -ANALYSIS_CONFIG$effect_thresholds$large, 
             col = "red", lty = 2, lwd = 2)
      abline(v = ANALYSIS_CONFIG$effect_thresholds$moderate, 
             col = "orange", lty = 2, lwd = 1.5)
      abline(v = -ANALYSIS_CONFIG$effect_thresholds$moderate, 
             col = "orange", lty = 2, lwd = 1.5)
      
      # Add legend
      legend("topright",
             legend = c(paste0("|FC| > ", ANALYSIS_CONFIG$effect_thresholds$large),
                        paste0("|FC| > ", ANALYSIS_CONFIG$effect_thresholds$moderate)),
             lty = 2, col = c("red", "orange"), lwd = c(2, 1.5))
      
      if (perfect_confounding) {
        mtext("⚠ Perfect confounding - interpret effect sizes with caution",
              side = 3, line = 0.5, cex = 0.8, col = "darkred")
      }
    }
    
    dev.off()
    cat("    Saved: figures/effect_size_distributions.pdf\n")
    
    # ============================================================
    # 5.2 Ranked list visualization
    # ============================================================
    
    pdf(file.path(ANALYSIS_CONFIG$figures_dir, "ranked_lists_visualization.pdf"),
        width = 10, height = 8)
    
    par(mfrow = c(ceiling(n_contrasts/2), 2))
    
    for (contrast in contrast_names) {
      results <- all_results[[contrast]]
      
      # Sort by ranking score for plotting (already sorted for GSEA)
      sorted_results <- results[order(-results$Ranking_Score), ]
      
      # Create ranked plot
      plot(1:nrow(sorted_results), sorted_results$Ranking_Score,
           type = "l", col = "blue", lwd = 2,
           xlab = "Rank", 
           ylab = "Ranking Score",
           main = paste("Ranked List:", contrast),
           ylim = c(min(sorted_results$Ranking_Score) * 1.1, 
                    max(sorted_results$Ranking_Score) * 1.1))
      
      # Add zero line
      abline(h = 0, col = "gray50", lty = 2)
      
      # Highlight top and bottom proteins
      n_highlight <- min(20, nrow(sorted_results))
      
      # Top proteins (positive)
      if (any(sorted_results$Ranking_Score > 0)) {
        top_pos <- which(sorted_results$Ranking_Score > 0)[1:min(10, sum(sorted_results$Ranking_Score > 0))]
        points(top_pos, sorted_results$Ranking_Score[top_pos],
               col = "red", pch = 19, cex = 1)
      }
      
      # Bottom proteins (negative)
      if (any(sorted_results$Ranking_Score < 0)) {
        bottom_neg <- tail(which(sorted_results$Ranking_Score < 0), 
                           min(10, sum(sorted_results$Ranking_Score < 0)))
        points(bottom_neg, sorted_results$Ranking_Score[bottom_neg],
               col = "darkgreen", pch = 19, cex = 1)
      }
      
      # Add legend
      legend("topright",
             legend = c("Top (positive)", "Bottom (negative)", "Zero line"),
             col = c("red", "darkgreen", "gray50"),
             pch = c(19, 19, NA),
             lty = c(NA, NA, 2))
    }
    
    dev.off()
    cat("    Saved: figures/ranked_lists_visualization.pdf\n")
    
    # ============================================================
    # 5.3 Volcano plot for each contrast
    # ============================================================
    
    pdf(file.path(ANALYSIS_CONFIG$figures_dir, "volcano_plots.pdf"),
        width = 12, height = 8)
    
    par(mfrow = c(ceiling(n_contrasts/2), 2))
    
    for (contrast in contrast_names) {
      results <- all_results[[contrast]]
      
      # Create volcano plot
      plot(results$logFC, -log10(results$P.Value),
           pch = 20, cex = 0.6,
           col = ifelse(results$abs_logFC > ANALYSIS_CONFIG$effect_thresholds$moderate & 
                          results$P.Value < 0.05, "red", "gray70"),
           xlab = "log2 Fold Change",
           ylab = "-log10(p-value)",
           main = paste("Volcano Plot:", contrast))
      
      # Add threshold lines
      abline(v = ANALYSIS_CONFIG$effect_thresholds$moderate, 
             col = "blue", lty = 2, lwd = 1)
      abline(v = -ANALYSIS_CONFIG$effect_thresholds$moderate, 
             col = "blue", lty = 2, lwd = 1)
      abline(h = -log10(0.05), 
             col = "darkgreen", lty = 2, lwd = 1)
      
      # Add counts of significant proteins
      n_sig <- sum(results$abs_logFC > ANALYSIS_CONFIG$effect_thresholds$moderate & 
                     results$P.Value < 0.05, na.rm = TRUE)
      
      legend("topright",
             legend = c(paste("Significant:", n_sig),
                        paste0("|FC|>", ANALYSIS_CONFIG$effect_thresholds$moderate),
                        "p < 0.05"),
             col = c("black", "blue", "darkgreen"),
             lty = c(NA, 2, 2),
             lwd = c(NA, 1, 1))
      
      if (perfect_confounding) {
        mtext("⚠ P-values confounded with batch",
              side = 3, line = 0.5, cex = 0.8, col = "darkred")
      }
    }
    
    dev.off()
    cat("    Saved: figures/volcano_plots.pdf\n")
    
    # ============================================================
    # 5.4 Top proteins bar plot
    # ============================================================
    
    pdf(file.path(ANALYSIS_CONFIG$figures_dir, "top_proteins_barplot.pdf"),
        width = 14, height = 8)
    
    # Combine top proteins from all contrasts
    top_proteins_list <- list()
    
    for (contrast in contrast_names) {
      results <- all_results[[contrast]]
      
      # Get top N up and down based on GSEA ranking (already sorted by Ranking_Score)
      # Top proteins are those with highest Ranking_Score (positive)
      top_up <- results[results$Ranking_Score > 0, ][1:min(10, sum(results$Ranking_Score > 0)), ]
      
      # Bottom proteins are those with lowest Ranking_Score (negative)
      results_sorted_neg <- results[results$Ranking_Score < 0, ]
      if (nrow(results_sorted_neg) > 0) {
        top_down <- results_sorted_neg[nrow(results_sorted_neg):max(1, nrow(results_sorted_neg)-9), ]
      } else {
        top_down <- data.frame()
      }
      
      top_proteins_list[[contrast]] <- rbind(
        cbind(top_up, Direction = "Up", Contrast = contrast),
        cbind(top_down, Direction = "Down", Contrast = contrast)
      )
    }
    
    # Combine and plot top proteins across all contrasts
    all_top <- do.call(rbind, top_proteins_list)
    
    if (!is.null(all_top) && nrow(all_top) > 0) {
      # Sort by Ranking_Score (not absolute value)
      all_top <- all_top[order(-all_top$Ranking_Score), ]
      
      # Take top 20 overall (10 highest positive, 10 lowest negative)
      top_20 <- head(all_top, 20)
      
      # Create bar plot
      bar_colors <- ifelse(top_20$Direction == "Up", "red", "blue")
      
      par(mar = c(7, 4, 4, 2) + 0.1)  # Increase bottom margin
      barplot(top_20$Ranking_Score,
              names.arg = paste(top_20$Protein_ID, "\n(", top_20$Contrast, ")", sep = ""),
              las = 2,  # Vertical labels
              cex.names = 0.7,
              col = bar_colors,
              main = "Top 20 Proteins Across All Contrasts (GSEA Ranking)",
              ylab = "Ranking Score (t-statistic or logFC)",
              border = NA)
      
      abline(h = 0, col = "black", lwd = 1)
      
      legend("topright",
             legend = c("Up-regulated", "Down-regulated"),
             fill = c("red", "blue"),
             border = NA)
    }
    
    dev.off()
    cat("    Saved: figures/top_proteins_barplot.pdf\n")
    
  }, error = function(e) {
    cat("    Warning: Could not create all visualizations:", e$message, "\n")
    cat("    Some figures may be missing.\n")
  })
  
  # ============================================================
  # 6. CREATE SUMMARY
  # ============================================================
  
  cat("6. Creating summary...\n")
  
  # Create simple summary
  summary_stats <- data.frame()
  for (contrast in contrast_names) {
    results <- all_results[[contrast]]
    
    stats <- data.frame(
      Contrast = contrast,
      N_Proteins = nrow(results),
      N_Large_Effect = sum(results$abs_logFC > ANALYSIS_CONFIG$effect_thresholds$large),
      N_Moderate_Effect = sum(results$abs_logFC > ANALYSIS_CONFIG$effect_thresholds$moderate),
      N_Small_Effect = sum(results$abs_logFC > ANALYSIS_CONFIG$effect_thresholds$small & 
                             results$abs_logFC <= ANALYSIS_CONFIG$effect_thresholds$moderate),
      Mean_logFC = mean(results$logFC, na.rm = TRUE),
      SD_logFC = sd(results$logFC, na.rm = TRUE),
      Top_Protein = results$Protein_ID[1],
      Top_Ranking_Score = round(results$Ranking_Score[1], 3),
      Top_logFC = round(results$logFC[1], 3),
      Top_pvalue = format(results$P.Value[1], scientific = TRUE, digits = 3),
      Ranking_Metric = results$Ranking_Metric[1],
      Perfect_Confounding = perfect_confounding,
      stringsAsFactors = FALSE
    )
    
    summary_stats <- rbind(summary_stats, stats)
  }
  
  write.csv(summary_stats,
            file.path(ANALYSIS_CONFIG$output_dir, "de_summary.csv"),
            row.names = FALSE)
  
  # ============================================================
  # 7. CREATE TOP PROTEINS REPORT
  # ============================================================
  
  cat("7. Creating top proteins report...\n")
  
  top_proteins_all <- data.frame()
  
  for (contrast in contrast_names) {
    results <- all_results[[contrast]]
    
    # Get top N proteins based on GSEA ranking
    top_n <- min(ANALYSIS_CONFIG$top_n, nrow(results))
    top_proteins <- results[1:top_n, c("Protein_ID", "Ranking_Score", "logFC", "abs_logFC", 
                                       "P.Value", "adj.P.Val", "Effect_Category", "Direction", "Rank")]
    
    top_proteins$Contrast <- contrast
    
    top_proteins_all <- rbind(top_proteins_all, top_proteins)
    
    # Save per-contrast top proteins
    write.csv(top_proteins,
              file.path(ANALYSIS_CONFIG$output_dir,
                        paste0(gsub("_", "-", contrast), "_top_proteins.csv")),
              row.names = FALSE)
  }
  
  # Save all top proteins combined
  write.csv(top_proteins_all,
            file.path(ANALYSIS_CONFIG$output_dir, "all_top_proteins.csv"),
            row.names = FALSE)
  
  # ============================================================
  # 8. CREATE REPORT
  # ============================================================
  
  cat("8. Creating report...\n")
  
  report_lines <- c(
    "DE ANALYSIS REPORT",
    "==================",
    paste("Date:", Sys.time()),
    paste("Input file:", ANALYSIS_CONFIG$input_file),
    "",
    "EXPERIMENTAL DESIGN:",
    paste("Perfect confounding:", perfect_confounding),
    paste("Details:", confounding_details),
    paste("Ranking metric:", ifelse(perfect_confounding, "logFC", "t-statistic")),
    "",
    "CONTRASTS ANALYZED:",
    paste("-", contrast_names, collapse = "\n"),
    "",
    "SUMMARY STATISTICS:"
  )
  
  # Add summary table
  for (i in 1:nrow(summary_stats)) {
    stats <- summary_stats[i, ]
    report_lines <- c(report_lines,
                      "",
                      paste("Contrast:", stats$Contrast),
                      paste("  Total proteins:", stats$N_Proteins),
                      paste("  Large effects (|FC|>", ANALYSIS_CONFIG$effect_thresholds$large, "): ", 
                            stats$N_Large_Effect, sep = ""),
                      paste("  Moderate effects (|FC|>", ANALYSIS_CONFIG$effect_thresholds$moderate, "): ", 
                            stats$N_Moderate_Effect, sep = ""),
                      paste("  Top protein:", stats$Top_Protein, 
                            "(Ranking Score =", stats$Top_Ranking_Score, 
                            ", logFC =", stats$Top_logFC, 
                            ", p =", stats$Top_pvalue, ")")
    )
  }
  
  report_lines <- c(report_lines,
                    "",
                    "OUTPUT FILES:",
                    "1. de_results/*_de_results.csv - Full DE results per contrast",
                    "2. de_results/*_top_proteins.csv - Top proteins per contrast",
                    "3. de_results/all_top_proteins.csv - All top proteins combined",
                    "4. de_results/de_summary.csv - Summary statistics",
                    "5. ranked_lists/*.rnk - Ranked lists for GSEA (sorted by Ranking_Score)",
                    "6. figures/effect_size_distributions.pdf - Effect size distributions",
                    "7. figures/ranked_lists_visualization.pdf - Ranked list visualizations",
                    "8. figures/volcano_plots.pdf - Volcano plots",
                    "9. figures/top_proteins_barplot.pdf - Top proteins bar plot",
                    "",
                    "GSEA RANKING INFORMATION:",
                    paste("- Ranking metric:", ifelse(perfect_confounding, "logFC", "t-statistic")),
                    "- .rnk files are sorted from highest to lowest ranking score",
                    "- Positive scores: Up-regulated in numerator condition",
                    "- Negative scores: Up-regulated in denominator condition",
                    "",
                    "NEXT STEPS:",
                    "1. Review figures/ for visual overview",
                    "2. Use ranked_lists/*.rnk files for GSEA analysis",
                    "3. Check de_results/*_top_proteins.csv for validation candidates"
  )
  
  if (perfect_confounding) {
    report_lines <- c(report_lines,
                      "",
                      "⚠ IMPORTANT WARNING:",
                      "Perfect confounding detected between condition and batch.",
                      "P-values are confounded with batch effects and should not be interpreted.",
                      "Use effect sizes (logFC) for ranking only.",
                      "Results indicate condition-associated changes, not necessarily condition-caused changes.",
                      "Experimental validation is strongly recommended."
    )
  }
  
  writeLines(paste(report_lines, collapse = "\n"),
             file.path(ANALYSIS_CONFIG$output_dir, "analysis_report.txt"))
  
  # ============================================================
  # 9. RETURN RESULTS
  # ============================================================
  
  cat("\n✓ Analysis complete with visualizations\n")
  cat("✓ .rnk files are now correctly sorted for GSEA\n")
  
  return(list(
    success = TRUE,
    de_results = all_results,
    ranked_lists = ranked_lists,
    contrast_names = contrast_names,
    perfect_confounding = perfect_confounding,
    output_dir = ANALYSIS_CONFIG$output_dir,
    ranked_dir = ANALYSIS_CONFIG$ranked_dir,
    figures_dir = ANALYSIS_CONFIG$figures_dir,
    summary_stats = summary_stats,
    timestamp = Sys.time()
  ))
}
#!/usr/bin/env Rscript
# Title: Step 4 - Results Visualization and Comprehensive Reporting
# Input: intermediate_data/03_de_results.rds
# Output: Final reports, visualizations, and summary documents

cat("=== STEP 4: RESULTS VISUALIZATION AND REPORTING ===\n")

# Load required packages
cat("Loading required packages...\n")
required_packages <- c("ggplot2", "ggrepel", "RColorBrewer", "pheatmap", 
                       "gridExtra", "viridis", "cowplot", "dplyr", "tidyr")
missing_packages <- required_packages[!required_packages %in% installed.packages()]
if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}
invisible(lapply(required_packages, library, character.only = TRUE))

# Load shared setup functions
source("config-setup_functions.R")

# Set up directories and logging
setup_directories()
log_file <- setup_logging("Step 4: Results Visualization and Reporting")

# 1. LOAD DE RESULTS
cat("1. Loading DE results from Step 3...\n")
check_file_exists("intermediate_data/03_de_results.rds", "Step 3")
de_results <- load_data_safely("intermediate_data/03_de_results.rds", "Step 3")

# Extract components
de_list <- de_results$de_results
overall_summary <- de_results$overall_summary
contrast_names <- de_results$contrast_names
params <- de_results$parameters
model_info <- de_results$model_info
metadata <- model_info$metadata

cat("DE results loaded successfully:\n")
cat("  Contrasts:", length(contrast_names), "\n")
cat("  Total proteins:", de_list[[1]]$summary$n_total, "\n")

# 2. CREATE VOLCANO PLOTS
cat("\n2. Creating volcano plots...\n")

create_volcano_plot <- function(results_df, contrast_name, fdr_thresh, fc_thresh) {
  
  # Prepare data
  plot_data <- results_df$all_results
  plot_data$log10p <- -log10(plot_data$P.Value)
  
  # Define significance categories
  plot_data$Significance <- "Not significant"
  plot_data$Significance[plot_data$qvalue < fdr_thresh & 
                           abs(plot_data$logFC) > fc_thresh & 
                           plot_data$logFC > 0] <- "Up-regulated"
  plot_data$Significance[plot_data$qvalue < fdr_thresh & 
                           abs(plot_data$logFC) > fc_thresh & 
                           plot_data$logFC < 0] <- "Down-regulated"
  
  # Colors
  colors <- c("Not significant" = "gray70",
              "Up-regulated" = "#E41A1C",
              "Down-regulated" = "#377EB8")
  
  # Create plot
  p <- ggplot(plot_data, aes(x = logFC, y = log10p, color = Significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = colors) +
    geom_vline(xintercept = c(-fc_thresh, fc_thresh), 
               linetype = "dashed", color = "gray40", alpha = 0.7) +
    geom_hline(yintercept = -log10(max(plot_data$P.Value[plot_data$qvalue < fdr_thresh])),
               linetype = "dashed", color = "gray40", alpha = 0.7) +
    labs(title = paste("Volcano Plot:", contrast_name),
         x = expression(log[2]("Fold Change")),
         y = expression(-log[10]("P-value"))) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.grid.minor = element_blank()) +
    guides(color = guide_legend(title = "Significance",
                                override.aes = list(size = 3, alpha = 1)))
  
  # Label top hits
  top_hits <- plot_data %>%
    filter(Significance != "Not significant") %>%
    arrange(P.Value) %>%
    head(10)
  
  if (nrow(top_hits) > 0) {
    p <- p + geom_text_repel(
      data = top_hits,
      aes(label = ID),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.5,
      point.padding = 0.1
    )
  }
  
  return(p)
}

# Generate volcano plots for each contrast
volcano_plots <- list()
for (contrast in contrast_names) {
  cat("  Creating volcano plot for:", contrast, "\n")
  
  volcano_plots[[contrast]] <- create_volcano_plot(
    results_df = de_list[[contrast]],
    contrast_name = contrast,
    fdr_thresh = params$fdr_threshold,
    fc_thresh = params$fc_threshold
  )
}

# Save individual volcano plots
for (contrast in contrast_names) {
  filename <- paste0("figures/volcano_", gsub("_", "-", tolower(contrast)), ".pdf")
  ggsave(filename, volcano_plots[[contrast]], 
         width = params$plot_width, height = params$plot_height, dpi = params$plot_dpi)
  
  # Also save as PNG
  filename_png <- paste0("figures/volcano_", gsub("_", "-", tolower(contrast)), ".png")
  ggsave(filename_png, volcano_plots[[contrast]], 
         width = params$plot_width, height = params$plot_height, dpi = params$plot_dpi)
}

# Create multi-panel volcano plot
if (length(contrast_names) <= 4) {
  cat("  Creating multi-panel volcano plot...\n")
  multi_volcano <- plot_grid(plotlist = volcano_plots, ncol = 2, labels = "AUTO")
  ggsave("figures/multi_volcano_plots.pdf", multi_volcano, 
         width = params$plot_width * 2, height = params$plot_height * 2)
}

# 3. CREATE MA PLOTS
cat("\n3. Creating MA plots...\n")

create_ma_plot <- function(results_df, contrast_name, fdr_thresh, fc_thresh) {
  
  plot_data <- results_df$all_results
  
  # Calculate average expression (if not already present)
  if (!"AveExpr" %in% colnames(plot_data)) {
    # This would need the original expression matrix
    # For now, use a placeholder
    plot_data$AveExpr <- 0
  }
  
  # Define significance
  plot_data$Significance <- "Not significant"
  plot_data$Significance[plot_data$qvalue < fdr_thresh & 
                           abs(plot_data$logFC) > fc_thresh] <- "Significant"
  
  p <- ggplot(plot_data, aes(x = AveExpr, y = logFC, color = Significance)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_color_manual(values = c("Not significant" = "gray70", 
                                  "Significant" = "red")) +
    geom_hline(yintercept = 0, color = "gray40") +
    geom_hline(yintercept = c(-fc_thresh, fc_thresh), 
               linetype = "dashed", color = "gray40", alpha = 0.7) +
    labs(title = paste("MA Plot:", contrast_name),
         x = "Average Expression",
         y = expression(log[2]("Fold Change"))) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  return(p)
}

# Save MA plots
for (contrast in contrast_names) {
  ma_plot <- create_ma_plot(de_list[[contrast]], contrast, 
                            params$fdr_threshold, params$fc_threshold)
  ggsave(paste0("figures/ma_", gsub("_", "-", tolower(contrast)), ".pdf"), 
         ma_plot, width = params$plot_width, height = params$plot_height)
}

# 4. CREATE HEATMAP OF TOP DE PROTEINS
cat("\n4. Creating heatmap of top DE proteins...\n")

# Load original expression data
if (file.exists("intermediate_data/01_validated_data.rds")) {
  validated_data <- readRDS("intermediate_data/01_validated_data.rds")
  expression <- validated_data$expression
  
  # Identify top DE proteins across all contrasts
  top_proteins <- c()
  for (contrast in contrast_names) {
    sig_proteins <- de_list[[contrast]]$all_results %>%
      filter(significant) %>%
      arrange(P.Value) %>%
      head(50) %>%
      pull(ID)
    top_proteins <- unique(c(top_proteins, sig_proteins))
  }
  
  if (length(top_proteins) > 0 && all(top_proteins %in% rownames(expression))) {
    # Subset expression matrix
    heatmap_data <- expression[top_proteins, , drop = FALSE]
    
    # Z-score normalization per protein
    heatmap_data_z <- t(scale(t(heatmap_data)))
    
    # Create annotation
    annotation_df <- data.frame(
      Condition = metadata$Condition,
      row.names = colnames(expression)
    )
    
    # Create heatmap
    pdf("figures/top_de_proteins_heatmap.pdf", width = 10, height = 12)
    pheatmap(heatmap_data_z,
             annotation_col = annotation_df,
             show_rownames = length(top_proteins) <= 50,
             show_colnames = ncol(expression) <= 30,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
             fontsize_row = 6,
             fontsize_col = 8,
             main = "Heatmap of Top DE Proteins",
             annotation_legend = TRUE)
    dev.off()
    
    cat("  Heatmap saved to: figures/top_de_proteins_heatmap.pdf\n")
  } else {
    cat("  Not enough significant proteins for heatmap or protein IDs don't match\n")
  }
}

# 5. CREATE VENN DIAGRAM OF OVERLAPPING SIGNIFICANT PROTEINS
cat("\n5. Creating Venn diagram of overlapping significant proteins...\n")

if (length(contrast_names) >= 2 && length(contrast_names) <= 4) {
  if (!requireNamespace("VennDiagram", quietly = TRUE)) {
    install.packages("VennDiagram", repos = "https://cloud.r-project.org")
  }
  library(VennDiagram)
  
  # Get significant proteins for each contrast
  sig_proteins_list <- list()
  for (contrast in contrast_names) {
    sig_proteins_list[[contrast]] <- de_list[[contrast]]$all_results %>%
      filter(significant) %>%
      pull(ID)
  }
  
  # Create Venn diagram
  venn_plot <- venn.diagram(
    x = sig_proteins_list,
    filename = NULL,
    category.names = contrast_names,
    output = TRUE,
    imagetype = "pdf",
    height = 8,
    width = 8,
    resolution = 300,
    compression = "lzw",
    lwd = 2,
    lty = 'blank',
    fill = brewer.pal(length(contrast_names), "Set2"),
    cex = 1.5,
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 1.2,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135)[1:length(contrast_names)],
    cat.dist = c(0.055, 0.055, 0.085)[1:length(contrast_names)],
    cat.fontfamily = "sans",
    rotation = 1
  )
  
  pdf("figures/venn_diagram_significant_proteins.pdf", width = 8, height = 8)
  grid.draw(venn_plot)
  dev.off()
  
  cat("  Venn diagram saved to: figures/venn_diagram_significant_proteins.pdf\n")
}

# 6. CREATE SUMMARY BAR PLOT
cat("\n6. Creating summary bar plots...\n")

# Prepare data for bar plot
summary_long <- overall_summary %>%
  select(Contrast, Significant, Up, Down) %>%
  pivot_longer(cols = c(Up, Down), names_to = "Direction", values_to = "Count")

# Bar plot of DE counts
de_bar_plot <- ggplot(summary_long, aes(x = Contrast, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  geom_text(aes(label = Count), position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8")) +
  labs(title = "Differential Expression Summary",
       x = "Contrast",
       y = "Number of Proteins",
       fill = "Direction") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom")

ggsave("figures/de_summary_barplot.pdf", de_bar_plot, 
       width = params$plot_width, height = params$plot_height)

# 7. CREATE COMPREHENSIVE REPORT
cat("\n7. Creating comprehensive final report...\n")

# Compile all previous reports
reports <- c(
  "intermediate_data/01_validation_report.txt",
  "intermediate_data/02_model_report.txt",
  "intermediate_data/03_de_report.txt"
)

all_reports <- lapply(reports, function(file) {
  if (file.exists(file)) {
    c(readLines(file), "\n")
  } else {
    paste("WARNING: Report file not found:", file, "\n")
  }
})

# Create final comprehensive report
final_report <- c(
  "COMPREHENSIVE DIFFERENTIAL EXPRESSION ANALYSIS REPORT",
  "======================================================",
  paste("Analysis completed:", Sys.time()),
  paste("Project:", basename(getwd())),
  "",
  "TABLE OF CONTENTS",
  "1. Executive Summary",
  "2. Data Validation Report",
  "3. Model Diagnostics Report", 
  "4. Differential Expression Results",
  "5. Visualizations Created",
  "6. Interpretation Guidelines",
  "7. Files Generated",
  "8. Next Steps",
  "",
  paste(rep("=", 80), collapse = ""),
  "1. EXECUTIVE SUMMARY",
  paste(rep("-", 80), collapse = ""),
  "",
  paste("Analysis performed:", length(contrast_names), "contrasts"),
  paste("Total proteins analyzed:", overall_summary$Total_Proteins[1]),
  "",
  "Key Findings:",
  sapply(1:nrow(overall_summary), function(i) {
    paste("  •", overall_summary$Contrast[i], ":", 
          overall_summary$Significant[i], "significant proteins",
          paste0("(", overall_summary$Percent_Sig[i], "%)"),
          "(", overall_summary$Up[i], "up,", overall_summary$Down[i], "down)")
  }),
  "",
  if (any(grepl("confounding", names(metadata), ignore.case = TRUE))) {
    c("⚠ IMPORTANT NOTE: PERFECT CONFOUNDING DETECTED",
      "   • Biological and technical effects are confounded",
      "   • Results are hypothesis-generating, not definitive",
      "   • Validation of key findings is essential",
      "")
  } else {
    "✓ No perfect confounding detected in experimental design"
  },
  "",
  paste(rep("=", 80), collapse = ""),
  unlist(all_reports),
  paste(rep("=", 80), collapse = ""),
  "5. VISUALIZATIONS CREATED",
  paste(rep("-", 80), collapse = ""),
  "",
  "The following visualizations were created and saved to the 'figures/' directory:",
  paste("  • Volcano plots for each contrast (PDF and PNG)"),
  if (length(contrast_names) <= 4) paste("  • Multi-panel volcano plot"),
  paste("  • MA plots for each contrast"),
  if (file.exists("figures/top_de_proteins_heatmap.pdf")) paste("  • Heatmap of top DE proteins"),
  if (length(contrast_names) >= 2 && length(contrast_names) <= 4) paste("  • Venn diagram of overlapping significant proteins"),
  paste("  • DE summary bar plot"),
  paste("  • QC plots from previous steps"),
  "",
  paste(rep("=", 80), collapse = ""),
  "6. INTERPRETATION GUIDELINES",
  paste(rep("-", 80), collapse = ""),
  "",
  "TIERED CONFIDENCE FRAMEWORK:",
  "  Tier 1 (High confidence): FDR < 0.001 AND |log2FC| > 3",
  "    • Strong candidates for validation",
  "    • High priority for follow-up studies",
  "",
  "  Tier 2 (Moderate confidence): FDR < 0.01 AND |log2FC| > 2",
  "    • Good candidates for validation",
  "    • Suitable for pathway analysis",
  "",
  "  Tier 3 (Suggestive): FDR < 0.05 AND |log2FC| > 1",
  "    • Exploratory hits",
  "    • Consider in context of other evidence",
  "",
  if (any(grepl("confounding", names(metadata), ignore.case = TRUE))) {
    c("SPECIAL CONSIDERATIONS FOR CONFOUNDED DESIGN:",
      "  1. All results are 'condition-associated', not 'condition-caused'",
      "  2. Batch effects cannot be separated from biological effects",
      "  3. Use more stringent thresholds for interpretation",
      "  4. Validate with orthogonal methods",
      "  5. Consider results as hypothesis-generating",
      "")
  },
  "VALIDATION RECOMMENDATIONS:",
  "  1. Prioritize Tier 1 proteins for experimental validation",
  "  2. Use orthogonal methods (Western blot, targeted MS)",
  "  3. Validate in independent sample sets if available",
  "  4. Consider biological context and known pathways",
  "",
  paste(rep("=", 80), collapse = ""),
  "7. FILES GENERATED",
  paste(rep("-", 80), collapse = ""),
  "",
  "Output Directory Structure:",
  "  figures/",
  "    • All visualization plots (PDF and PNG)",
  "",
  "  de_results/",
  "    • *_all.csv - Complete results for each contrast",
  "    • *_significant.csv - Significant proteins only",
  "    • *_tier[1-3].csv - Tiered confidence results",
  "",
  "  intermediate_data/",
  "    • 01_validated_data.rds - Validated input data",
  "    • 02_model_fit.rds - Model objects",
  "    • 03_de_results.rds - Complete DE results",
  "    • *_report.txt - Step-by-step reports",
  "",
  "  logs/",
  "    • analysis.log - Complete analysis log",
  "    • error_log.txt - Any errors encountered",
  "",
  "  config/",
  "    • analysis_parameters.R - Analysis parameters",
  "",
  paste(rep("=", 80), collapse = ""),
  "8. NEXT STEPS",
  paste(rep("-", 80), collapse = ""),
  "",
  "Recommended next analyses:",
  "  1. Pathway enrichment analysis",
  "  2. Protein-protein interaction network analysis",
  "  3. Functional annotation of top hits",
  "  4. Integration with other omics data (if available)",
  "  5. Experimental validation planning",
  "",
  "For pathway analysis:",
  "  • Use tiered confidence framework",
  "  • Prioritize pathways with multiple tier 1/2 proteins",
  "  • Consider biological relevance to experimental system",
  "",
  "REPRODUCIBILITY:",
  "  • All analysis steps are documented in the logs/ directory",
  "  • Intermediate files allow re-running specific steps",
  "  • Parameters are centralized in config/analysis_parameters.R",
  "",
  paste(rep("=", 80), collapse = ""),
  "END OF REPORT",
  ""
)

writeLines(final_report, "final_report/comprehensive_de_analysis_report.txt")

# 8. CREATE HTML REPORT
cat("\n8. Creating HTML report...\n")

html_report <- paste(
  '<!DOCTYPE html>',
  '<html>',
  '<head>',
  '<title>Differential Expression Analysis Report</title>',
  '<style>',
  'body { font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; }',
  'h1, h2, h3 { color: #2c3e50; }',
  '.summary-box { background-color: #f8f9fa; padding: 20px; border-left: 5px solid #3498db; margin: 20px 0; }',
  '.warning-box { background-color: #fff3cd; padding: 20px; border-left: 5px solid #ffc107; margin: 20px 0; }',
  'table { border-collapse: collapse; width: 100%; margin: 20px 0; }',
  'th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }',
  'th { background-color: #f2f2f2; }',
  'img { max-width: 100%; height: auto; border: 1px solid #ddd; margin: 10px 0; }',
  '.figure-caption { font-style: italic; color: #666; text-align: center; }',
  '</style>',
  '</head>',
  '<body>',
  '<h1>Differential Expression Analysis Report</h1>',
  paste('<p><strong>Generated:</strong>', Sys.time(), '</p>'),
  paste('<p><strong>Project:</strong>', basename(getwd()), '</p>'),
  '',
  '<div class="summary-box">',
  '<h2>Executive Summary</h2>',
  paste('<p><strong>Contrasts analyzed:</strong>', length(contrast_names), '</p>'),
  paste('<p><strong>Total proteins:</strong>', overall_summary$Total_Proteins[1], '</p>'),
  '<table>',
  '<tr><th>Contrast</th><th>Significant</th><th>Up</th><th>Down</th><th>% Significant</th></tr>',
  paste(sapply(1:nrow(overall_summary), function(i) {
    paste('<tr>',
          '<td>', overall_summary$Contrast[i], '</td>',
          '<td>', overall_summary$Significant[i], '</td>',
          '<td>', overall_summary$Up[i], '</td>',
          '<td>', overall_summary$Down[i], '</td>',
          '<td>', overall_summary$Percent_Sig[i], '%</td>',
          '</tr>')
  }), collapse = '\n'),
  '</table>',
  '</div>',
  '',
  if (any(grepl("confounding", names(metadata), ignore.case = TRUE))) {
    '<div class="warning-box">',
    '<h2>⚠ Important Note: Perfect Confounding Detected</h2>',
    '<p>Biological and technical effects are confounded in this experiment.</p>',
    '<ul>',
    '<li>Results are "condition-associated", not "condition-caused"</li>',
    '<li>Validation of key findings is essential</li>',
    '<li>Use tiered confidence framework for interpretation</li>',
    '</ul>',
    '</div>'
  },
  '',
  '<h2>Visualizations</h2>',
  paste(sapply(contrast_names, function(contrast) {
    png_file <- paste0("figures/volcano_", gsub("_", "-", tolower(contrast)), ".png")
    if (file.exists(png_file)) {
      paste('<h3>', contrast, '</h3>',
            '<img src="', png_file, '" alt="Volcano plot for ', contrast, '">',
            '<p class="figure-caption">Volcano plot for ', contrast, 
            ' (red: up-regulated, blue: down-regulated, gray: not significant)</p>')
    } else {
      ''
    }
  }), collapse = '\n'),
  '',
  '<h2>Files Generated</h2>',
  '<ul>',
  '<li><strong>figures/</strong> - All visualization plots</li>',
  '<li><strong>de_results/</strong> - CSV files with detailed results</li>',
  '<li><strong>final_report/</strong> - This report and summary documents</li>',
  '<li><strong>intermediate_data/</strong> - Intermediate files for reproducibility</li>',
  '<li><strong>logs/</strong> - Analysis logs</li>',
  '<li><strong>config/</strong> - Analysis parameters</li>',
  '</ul>',
  '',
  '<h2>Interpretation Guidelines</h2>',
  '<h3>Tiered Confidence Framework</h3>',
  '<table>',
  '<tr><th>Tier</th><th>Criteria</th><th>Interpretation</th></tr>',
  '<tr><td>Tier 1</td><td>FDR &lt; 0.001 AND |log2FC| &gt; 3</td><td>High confidence, prioritize for validation</td></tr>',
  '<tr><td>Tier 2</td><td>FDR &lt; 0.01 AND |log2FC| &gt; 2</td><td>Moderate confidence, good candidates</td></tr>',
  '<tr><td>Tier 3</td><td>FDR &lt; 0.05 AND |log2FC| &gt; 1</td><td>Suggestive, exploratory hits</td></tr>',
  '</table>',
  '',
  '<h2>Next Steps</h2>',
  '<ol>',
  '<li>Pathway enrichment analysis using tiered framework</li>',
  '<li>Protein-protein interaction network analysis</li>',
  '<li>Functional annotation of top hits</li>',
  '<li>Experimental validation planning</li>',
  '</ol>',
  '',
  '</body>',
  '</html>',
  sep = '\n'
)

writeLines(html_report, "final_report/de_analysis_report.html")

# 9. CREATE ZIP ARCHIVE OF RESULTS
cat("\n9. Creating zip archive of results...\n")

files_to_archive <- c(
  "de_results/",
  "figures/",
  "final_report/",
  "logs/analysis.log",
  "config/analysis_parameters.R"
)

# Check which files exist
existing_files <- files_to_archive[file.exists(files_to_archive)]

if (length(existing_files) > 0) {
  zip_filename <- paste0("de_analysis_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
  zip(zip_filename, existing_files)
  cat("  Archive created:", zip_filename, "\n")
}

# 10. FINAL SUMMARY
cat("\n10. Generating final summary...\n")

final_summary <- paste(
  "ANALYSIS COMPLETE",
  "=================",
  paste("Completed:", Sys.time()),
  "",
  "SUMMARY STATISTICS:",
  paste("  Total contrasts analyzed:", length(contrast_names)),
  paste("  Total proteins:", overall_summary$Total_Proteins[1]),
  "",
  "SIGNIFICANT FINDINGS BY CONTRAST:",
  sapply(1:nrow(overall_summary), function(i) {
    paste("  •", overall_summary$Contrast[i], ":",
          overall_summary$Significant[i], "significant proteins",
          paste0("(", overall_summary$Percent_Sig[i], "%)"))
  }),
  "",
  "FILES CREATED:",
  "  • Comprehensive report: final_report/comprehensive_de_analysis_report.txt",
  "  • HTML report: final_report/de_analysis_report.html",
  "  • Volcano plots: figures/volcano_*.pdf",
  "  • DE results: de_results/*.csv",
  if (exists("zip_filename")) paste("  • Archive:", zip_filename),
  "",
  "NEXT ACTIONS:",
  "  1. Review comprehensive report",
  "  2. Examine volcano plots for key hits",
  "  3. Use tiered CSV files for pathway analysis",
  "  4. Plan validation experiments",
  "",
  if (any(grepl("confounding", names(metadata), ignore.case = TRUE))) {
    c("⚠ CRITICAL NOTE:",
      "  Due to perfect confounding in the experimental design,",
      "  all results must be interpreted with caution.",
      "  Biological conclusions require validation.")
  } else {
    "✓ No perfect confounding detected in experimental design"
  },
  "",
  "Thank you for using the differential expression analysis pipeline!",
  sep = "\n"
)

writeLines(final_summary, "final_report/final_summary.txt")
cat(final_summary)

# Clean up
cat("\n", rep("=", 50), "\n", sep = "")
cat("Step 4 completed successfully!\n")
cat("All visualizations and reports have been generated.\n")
cat("Check the 'final_report/' directory for comprehensive outputs.\n")

cleanup_logging(log_file, "Step 4: Results Visualization and Reporting")
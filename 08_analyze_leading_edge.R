#!/usr/bin/env Rscript
# scripts/08_analyze_leading_edge.R - USES LOCAL KEGG DATABASE

cat("========================================\n")
cat("LEADING EDGE ANALYSIS (USING LOCAL KEGG)\n")
cat("========================================\n\n")

# ============================================================
# 1. LOAD REQUIRED PACKAGES
# ============================================================

cat("1. Loading packages...\n")
library(dplyr)
library(tidyr)

# ============================================================
# 2. LOAD LOCAL KEGG DATABASE
# ============================================================

cat("2. Loading local KEGG database...\n")

kegg_file <- "data/kegg_human_data.rds"
if (!file.exists(kegg_file)) {
  cat("❌ ERROR: KEGG database not found:", kegg_file, "\n")
  cat("  Run the KEGG data preparation script first.\n")
  quit(status = 1)
}

kegg_data <- readRDS(kegg_file)
cat("  ✓ Loaded KEGG data from:", kegg_file, "\n")

# Extract pathway-gene mapping
if (is.null(kegg_data$KEGGPATHID2EXTID)) {
  cat("❌ ERROR: KEGG data format incorrect\n")
  quit(status = 1)
}

pathway_data <- kegg_data$KEGGPATHID2EXTID

# Create two useful formats:
# 1. Named list (for quick lookup)
kegg_gene_sets <- split(pathway_data$to, pathway_data$from)

# 2. Data frame with pathway names
if (!is.null(kegg_data$KEGGPATHID2NAME)) {
  pathway_names <- kegg_data$KEGGPATHID2NAME
  colnames(pathway_names) <- c("pathway_id", "pathway_name")
} else {
  # Create from IDs if names not available
  pathway_names <- data.frame(
    pathway_id = unique(pathway_data$from),
    pathway_name = unique(pathway_data$from),
    stringsAsFactors = FALSE
  )
}

cat("  ✓ ", length(kegg_gene_sets), " pathways loaded\n", sep = "")
cat("  ✓ ", nrow(pathway_data), " pathway-gene associations\n", sep = "")

# ============================================================
# 3. LOAD FGSEA RESULTS
# ============================================================

cat("3. Loading fgsea results...\n")

results_file <- "gsea_results/tables/fgsea_all_results.csv"
if (!file.exists(results_file)) {
  cat("❌ ERROR: Results file not found:", results_file, "\n")
  quit(status = 1)
}

fgsea_results <- read.csv(results_file, stringsAsFactors = FALSE)

# Convert leadingEdge column from semicolon-separated to list
fgsea_results$leadingEdge_list <- strsplit(fgsea_results$leadingEdge, ";")

cat("  ✓ Loaded ", nrow(fgsea_results), " pathway results\n", sep = "")
cat("  ✓ Contrasts: ", paste(unique(fgsea_results$Contrast), collapse = ", "), "\n", sep = "")

# ============================================================
# 4. ENRICH LEADING EDGE WITH KEGG PATHWAY INFO
# ============================================================

cat("4. Enriching leading edge with KEGG pathway information...\n")

# Add pathway descriptions from KEGG database
enriched_results <- fgsea_results

# Merge with pathway names
pathway_id_to_name <- setNames(pathway_names$pathway_name, pathway_names$pathway_id)

# Map fgsea pathway IDs to KEGG pathway names
enriched_results$kegg_pathway_name <- sapply(enriched_results$pathway, function(p) {
  if (p %in% names(pathway_id_to_name)) {
    pathway_id_to_name[[p]]
  } else {
    # Try to find by partial match
    matches <- grep(p, pathway_id_to_name, ignore.case = TRUE, value = TRUE)
    if (length(matches) > 0) {
      return(matches[1])
    } else {
      return(p)  # Return original if no match
    }
  }
})

# Also add full gene list for each pathway from KEGG database
enriched_results$full_pathway_genes <- lapply(enriched_results$pathway, function(p) {
  if (p %in% names(kegg_gene_sets)) {
    kegg_gene_sets[[p]]
  } else {
    character(0)
  }
})

# Calculate leading edge percentage using KEGG pathway size
enriched_results$pathway_size_kegg <- sapply(enriched_results$full_pathway_genes, length)
enriched_results$leading_edge_count <- sapply(enriched_results$leadingEdge_list, length)
enriched_results$leading_edge_percent <- ifelse(enriched_results$pathway_size_kegg > 0,
                                                round(enriched_results$leading_edge_count / enriched_results$pathway_size_kegg * 100, 1),
                                                0
)

cat("  ✓ Added KEGG pathway names and sizes\n")

# ============================================================
# 5. CONTRAST-SPECIFIC LEADING EDGE ANALYSIS
# ============================================================

cat("5. Analyzing leading edge by contrast...\n")

contrasts_to_analyze <- c("NHA-vs-Control", "Treated-vs-Control", "Treated-vs-NHA")

for (contrast in contrasts_to_analyze) {
  cat("\n", paste0(rep("=", 60), collapse = ""), "\n", sep = "")
  cat("CONTRAST: ", contrast, "\n", sep = "")
  cat(paste0(rep("=", 60), collapse = ""), "\n")
  
  # Filter for this contrast
  contrast_data <- enriched_results[enriched_results$Contrast == contrast, ]
  
  if (nrow(contrast_data) == 0) {
    cat("  No results found\n")
    next
  }
  
  # Sort by |NES|
  contrast_data$abs_NES <- abs(contrast_data$NES)
  contrast_data <- contrast_data[order(-contrast_data$abs_NES), ]
  
  # Take top 3
  top_3 <- head(contrast_data, 3)
  
  for (i in 1:nrow(top_3)) {
    row <- top_3[i, ]
    
    cat("\n  ", i, ". ", row$kegg_pathway_name, "\n", sep = "")
    cat("      NES: ", round(row$NES, 2), 
        " (", ifelse(row$NES > 0, "enriched in first condition", "enriched in second condition"), ")\n", sep = "")
    cat("      p-value: ", format(row$pval, scientific = TRUE, digits = 2), "\n", sep = "")
    cat("      Pathway size (KEGG): ", row$pathway_size_kegg, " genes\n", sep = "")
    cat("      Leading edge: ", row$leading_edge_count, " genes (", 
        row$leading_edge_percent, "% of pathway)\n", sep = "")
    
    # Show leading edge genes
    if (row$leading_edge_count > 0) {
      leading_genes <- row$leadingEdge_list[[1]]
      cat("      Top 10 leading edge genes:\n        ")
      if (length(leading_genes) <= 10) {
        cat(paste(leading_genes, collapse = ", "), "\n")
      } else {
        cat(paste(leading_genes[1:10], collapse = ", "), 
            " ... (+", length(leading_genes) - 10, " more)\n", sep = "")
      }
    }
  }
  
  # Save detailed results for this contrast
  safe_name <- gsub("[^A-Za-z0-9]", "_", contrast)
  output_file <- paste0("gsea_results/tables/leading_edge_kegg_", safe_name, ".csv")
  
  # Prepare detailed output
  detailed_output <- contrast_data %>%
    select(
      pathway = kegg_pathway_name,
      pathway_id = pathway,
      NES,
      pval,
      pathway_size_kegg,
      leading_edge_count,
      leading_edge_percent,
      leading_edge_genes = leadingEdge
    ) %>%
    arrange(desc(abs(NES)))
  
  write.csv(detailed_output, output_file, row.names = FALSE)
  cat("\n  ✓ Saved detailed results to: ", output_file, "\n", sep = "")
}

# ============================================================
# 6. CROSS-CONTRAST LEADING EDGE COMPARISON
# ============================================================

cat("\n\n6. Cross-contrast leading edge comparison...\n")

# Focus on Oxidative Phosphorylation (your strongest signal)
oxphos_pathways <- enriched_results[grepl("OXIDATIVE_PHOSPHORYLATION", 
                                          enriched_results$pathway, 
                                          ignore.case = TRUE), ]

if (nrow(oxphos_pathways) > 0) {
  cat("\n  OXYDATIVE PHOSPHORYLATION ACROSS CONTRASTS:\n")
  cat("  ", paste0(rep("-", 50), collapse = ""), "\n", sep = "")
  
  oxphos_comparison <- data.frame()
  
  for (i in 1:nrow(oxphos_pathways)) {
    row <- oxphos_pathways[i, ]
    leading_genes <- row$leadingEdge_list[[1]]
    full_genes <- row$full_pathway_genes[[1]]
    
    oxphos_comparison <- rbind(oxphos_comparison, data.frame(
      Contrast = row$Contrast,
      NES = round(row$NES, 2),
      Direction = ifelse(row$NES > 0, "UP in first", "DOWN in first"),
      Leading_Edge_Count = length(leading_genes),
      Pathway_Size = length(full_genes),
      Percent_Leading = round(length(leading_genes)/length(full_genes)*100, 1),
      Leading_Edge_Genes = paste(leading_genes, collapse = ";"),
      stringsAsFactors = FALSE
    ))
    
    cat("\n  ", row$Contrast, ":\n", sep = "")
    cat("      NES = ", round(row$NES, 2), " (", 
        ifelse(row$NES > 0, "upregulated", "downregulated"), ")\n", sep = "")
    cat("      Leading edge: ", length(leading_genes), "/", length(full_genes),
        " genes (", round(length(leading_genes)/length(full_genes)*100, 1), "%)\n", sep = "")
    
    # Show specific complexes if available
    if (length(leading_genes) > 0) {
      # Categorize by OxPhos complex
      complexes <- list(
        Complex_I = grep("^NDUF|^Nduf", leading_genes, value = TRUE, ignore.case = TRUE),
        Complex_II = grep("^SDH|^Sdh", leading_genes, value = TRUE, ignore.case = TRUE),
        Complex_III = grep("^UQCR|^Uqc|^CYC1", leading_genes, value = TRUE, ignore.case = TRUE),
        Complex_IV = grep("^COX|^Cox|^MT-CO", leading_genes, value = TRUE, ignore.case = TRUE),
        Complex_V = grep("^ATP|^Atp", leading_genes, value = TRUE, ignore.case = TRUE)
      )
      
      cat("      By complex:\n")
      for (complex_name in names(complexes)) {
        if (length(complexes[[complex_name]]) > 0) {
          cat("        ", complex_name, ": ", 
              length(complexes[[complex_name]]), " genes\n", sep = "")
        }
      }
    }
  }
  
  # Save OxPhos comparison
  write.csv(oxphos_comparison, 
            "gsea_results/tables/oxphos_leading_edge_comparison.csv",
            row.names = FALSE)
  cat("\n  ✓ Saved OxPhos comparison to: gsea_results/tables/oxphos_leading_edge_comparison.csv\n")
}

# ============================================================
# 7. GENE-LEVEL ANALYSIS USING KEGG DATABASE
# ============================================================

cat("\n\n7. Gene-level analysis using KEGG database...\n")

# Get all unique leading edge genes across all contrasts
all_leading_genes <- unique(unlist(enriched_results$leadingEdge_list))

cat("  Total unique leading edge genes: ", length(all_leading_genes), "\n", sep = "")

# For each leading edge gene, find all KEGG pathways it belongs to
gene_pathway_map <- list()

for (gene in all_leading_genes) {
  # Find pathways containing this gene in KEGG database
  pathways_with_gene <- pathway_data[pathway_data$to == gene, "from"]
  
  if (length(pathways_with_gene) > 0) {
    # Get pathway names
    pathway_names <- sapply(pathways_with_gene, function(p) {
      if (p %in% names(pathway_id_to_name)) {
        pathway_id_to_name[[p]]
      } else {
        p
      }
    })
    
    gene_pathway_map[[gene]] <- list(
      pathways = pathways_with_gene,
      pathway_names = pathway_names,
      pathway_count = length(pathways_with_gene)
    )
  }
}

# Create gene summary table
gene_summary <- data.frame(
  Gene = names(gene_pathway_map),
  Pathway_Count = sapply(gene_pathway_map, function(x) x$pathway_count),
  Pathways = sapply(gene_pathway_map, function(x) {
    if (length(x$pathway_names) <= 3) {
      paste(x$pathway_names, collapse = "; ")
    } else {
      paste0(paste(head(x$pathway_names, 3), collapse = "; "),
             "; ... (+", length(x$pathway_names) - 3, " more)")
    }
  }),
  stringsAsFactors = FALSE
)

# Sort by pathway count
gene_summary <- gene_summary[order(-gene_summary$Pathway_Count), ]

cat("\n  Top 20 genes by number of KEGG pathways:\n")
print(head(gene_summary, 20))

# Save gene summary
write.csv(gene_summary, 
          "gsea_results/tables/leading_edge_genes_kegg_pathways.csv",
          row.names = FALSE)
cat("\n  ✓ Saved gene-pathway mapping to: gsea_results/tables/leading_edge_genes_kegg_pathways.csv\n")

# ============================================================
# 8. FINAL SUMMARY
# ============================================================

cat("\n\n", paste0(rep("=", 70), collapse = ""), "\n", sep = "")
cat("✅ LEADING EDGE ANALYSIS COMPLETE (USING LOCAL KEGG DATABASE)\n")
cat(paste0(rep("=", 70), collapse = ""), "\n\n")

cat("Summary:\n")
cat("  • Used local KEGG database: data/kegg_human_data.rds\n")
cat("  • Analyzed leading edge for ", nrow(enriched_results), " pathway results\n", sep = "")
cat("  • Identified ", length(all_leading_genes), " unique leading edge genes\n", sep = "")
cat("\n")

cat("Key files created:\n")
for (contrast in contrasts_to_analyze) {
  safe_name <- gsub("[^A-Za-z0-9]", "_", contrast)
  cat("  • gsea_results/tables/leading_edge_kegg_", safe_name, ".csv\n", sep = "")
}
cat("  • gsea_results/tables/oxphos_leading_edge_comparison.csv\n")
cat("  • gsea_results/tables/leading_edge_genes_kegg_pathways.csv\n")

cat("\nBiological insights gained:\n")
cat("  1. Exact percentage of each pathway in leading edge\n")
cat("  2. Specific genes driving each enrichment\n")
cat("  3. Cross-contrast comparison of the same pathway\n")
cat("  4. Pathway membership of leading edge genes\n")


# ============================================================
# 9. CREATE COMPREHENSIVE TEXT SUMMARY
# ============================================================

cat("\n\n9. Creating comprehensive text summary...\n")

summary_file <- "gsea_results/leading_edge_analysis_summary.txt"

sink(summary_file)  # Start capturing output to file

cat("LEADING EDGE ANALYSIS SUMMARY REPORT\n")
cat("=====================================\n")
cat("Generated: ", Sys.time(), "\n")
cat("Analysis based on local KEGG database: data/kegg_human_data.rds\n")
cat("FGSEA results: gsea_results/tables/fgsea_all_results.csv\n\n")

cat("OVERVIEW\n")
cat("--------\n")
cat("Total enriched pathways: ", nrow(enriched_results), "\n")
cat("Unique contrasts: ", paste(unique(enriched_results$Contrast), collapse = ", "), "\n")
cat("Total leading edge genes: ", length(all_leading_genes), "\n")
cat("Average leading edge coverage: ", round(mean(enriched_results$leading_edge_percent, na.rm = TRUE), 1), "%\n\n")

cat("TOP PATHWAYS BY CONTRAST\n")
cat("------------------------\n")

for (contrast in contrasts_to_analyze) {
  contrast_data <- enriched_results[enriched_results$Contrast == contrast, ]
  
  if (nrow(contrast_data) > 0) {
    contrast_data$abs_NES <- abs(contrast_data$NES)
    top_3 <- head(contrast_data[order(-contrast_data$abs_NES), ], 3)
    
    cat("\n", contrast, ":\n", sep = "")
    cat(paste0(rep("-", nchar(contrast) + 1), collapse = ""), "\n")
    
    for (i in 1:nrow(top_3)) {
      row <- top_3[i, ]
      cat("\n", i, ". ", row$kegg_pathway_name, "\n", sep = "")
      cat("    NES: ", round(row$NES, 2), 
          " (", ifelse(row$NES > 0, "↑ in first condition", "↓ in first condition"), ")\n", sep = "")
      cat("    p-value: ", format(row$pval, scientific = TRUE, digits = 2), "\n", sep = "")
      cat("    Pathway size: ", row$pathway_size_kegg, " genes\n", sep = "")
      cat("    Leading edge: ", row$leading_edge_count, " genes (", 
          row$leading_edge_percent, "% coverage)\n", sep = "")
      
      # Show top 5 leading edge genes
      leading_genes <- row$leadingEdge_list[[1]]
      if (length(leading_genes) > 0) {
        cat("    Top genes: ", paste(head(leading_genes, 5), collapse = ", "), "\n", sep = "")
      }
    }
  }
}

cat("\n\nKEY BIOLOGICAL INSIGHTS\n")
cat("----------------------\n")

# OxPhos comparison
oxphos_rows <- enriched_results[grepl("OXIDATIVE_PHOSPHORYLATION", 
                                      enriched_results$pathway, 
                                      ignore.case = TRUE), ]

if (nrow(oxphos_rows) > 0) {
  cat("1. OXIDATIVE PHOSPHORYLATION DYNAMICS:\n")
  for (i in 1:nrow(oxphos_rows)) {
    row <- oxphos_rows[i, ]
    cat("   • ", row$Contrast, ": NES = ", round(row$NES, 2), 
        " (", ifelse(row$NES > 0, "upregulated", "downregulated"), ")\n", sep = "")
    cat("     Leading edge: ", row$leading_edge_count, "/", row$pathway_size_kegg,
        " genes (", row$leading_edge_percent, "%)\n", sep = "")
  }
  cat("   → Interpretation: Treatment strongly reverses OxPhos elevation in NHA condition\n\n")
}

# Leading edge coverage analysis
high_coverage <- enriched_results[enriched_results$leading_edge_percent > 50, ]
if (nrow(high_coverage) > 0) {
  cat("2. PATHWAYS WITH HIGH COORDINATED REGULATION (>50% leading edge):\n")
  for (i in 1:min(5, nrow(high_coverage))) {
    row <- high_coverage[i, ]
    cat("   • ", row$kegg_pathway_name, " (", row$Contrast, "): ", 
        row$leading_edge_percent, "% coverage\n", sep = "")
  }
  cat("   → These pathways show highly coordinated gene expression changes\n\n")
}

# Multi-pathway genes
if (nrow(gene_summary) > 0) {
  top_multi_genes <- head(gene_summary[gene_summary$Pathway_Count > 1, ], 5)
  if (nrow(top_multi_genes) > 0) {
    cat("3. KEY REGULATOR GENES (APPEAR IN MULTIPLE PATHWAYS):\n")
    for (i in 1:nrow(top_multi_genes)) {
      cat("   • ", top_multi_genes$Gene[i], ": ", 
          top_multi_genes$Pathway_Count[i], " pathways\n", sep = "")
    }
    cat("   → These genes may be central regulators in the response\n")
  }
}

cat("\n\nRECOMMENDATIONS FOR VALIDATION\n")
cat("-----------------------------\n")

# Identify validation candidates
cat("1. PRIORITY VALIDATION TARGETS:\n")

# OxPhos complex genes
if (exists("oxphos_comparison")) {
  cat("   • Oxidative Phosphorylation complexes:\n")
  cat("     - Complex I (NDUF genes)\n")
  cat("     - Complex IV (COX genes)\n")
  cat("     - Complex V (ATP synthase genes)\n")
}

# Top leading edge genes by frequency
if (nrow(gene_summary) > 0) {
  top_genes <- head(gene_summary, 5)
  cat("\n2. FREQUENTLY CHANGED GENES:\n")
  for (i in 1:nrow(top_genes)) {
    cat("   • ", top_genes$Gene[i], " (in ", top_genes$Pathway_Count[i], " pathways)\n", sep = "")
  }
}

cat("\n3. PATHWAYS FOR FUNCTIONAL VALIDATION:\n")
cat("   • Oxidative Phosphorylation (strongest signal)\n")
cat("   • Glycolysis/Gluconeogenesis (consistent across contrasts)\n")
cat("   • Citrate cycle (TCA cycle)\n")

cat("\n\nOUTPUT FILES\n")
cat("-----------\n")
cat("1. Detailed results by contrast:\n")
for (contrast in contrasts_to_analyze) {
  safe_name <- gsub("[^A-Za-z0-9]", "_", contrast)
  cat("   • leading_edge_kegg_", safe_name, ".csv\n", sep = "")
}
cat("2. Cross-contrast comparisons:\n")
cat("   • oxphos_leading_edge_comparison.csv\n")
cat("   • leading_edge_genes_kegg_pathways.csv\n")
cat("3. This summary: leading_edge_analysis_summary.txt\n")

cat("\nSTATISTICAL NOTES\n")
cat("----------------\n")
cat("• Analysis uses gene permutation GSEA (appropriate for confounded design)\n")
cat("• Leading edge genes are those driving the enrichment signal\n")
cat("• High leading edge percentage indicates coordinated pathway regulation\n")
cat("• NES magnitude >1.5 indicates strong biological effect\n")

sink()  # Stop capturing output

cat("  ✓ Created comprehensive summary: ", summary_file, "\n", sep = "")


# ============================================================
# 10. CREATE LEADING EDGE VISUALIZATIONS
# ============================================================

cat("\n\n10. Creating leading edge visualizations...\n")

# Create visualization directory
viz_dir <- "gsea_results/figures/leading_edge_viz"
dir.create(viz_dir, recursive = TRUE, showWarnings = FALSE)

library(ggplot2)
library(reshape2)

# 10.1 Leading Edge Coverage Plot
cat("  Creating leading edge coverage plot...\n")

# Prepare data
viz_data <- enriched_results %>%
  select(Contrast, kegg_pathway_name, NES, leading_edge_percent) %>%
  arrange(desc(abs(NES))) %>%
  head(15)  # Top 15 by |NES|

# Shorten pathway names for plotting
viz_data$pathway_short <- sapply(viz_data$kegg_pathway_name, function(x) {
  if (nchar(x) > 40) paste0(substr(x, 1, 37), "...") else x
})

p_coverage <- ggplot(viz_data, aes(x = reorder(pathway_short, abs(NES)), 
                                   y = leading_edge_percent,
                                   fill = NES)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, name = "NES") +
  labs(title = "Leading Edge Coverage of Top Pathways",
       subtitle = "Percentage of pathway genes in leading edge | Color = NES",
       x = "Pathway",
       y = "Leading Edge Coverage (%)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9),
        legend.position = "bottom") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50", alpha = 0.7)

ggsave(file.path(viz_dir, "leading_edge_coverage.pdf"),
       p_coverage, width = 12, height = 8)
cat("    ✓ Created: leading_edge_coverage.pdf\n")

# 10.2 NES vs Leading Edge Percentage Scatter Plot
cat("  Creating NES vs coverage scatter plot...\n")

p_scatter <- ggplot(enriched_results, aes(x = NES, y = leading_edge_percent,
                                          color = Contrast, size = pathway_size_kegg)) +
  geom_point(alpha = 0.7) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "NES vs Leading Edge Coverage",
       subtitle = "Size = Pathway size | Each point = one pathway",
       x = "Normalized Enrichment Score (NES)",
       y = "Leading Edge Coverage (%)") +
  theme_minimal() +
  theme(legend.position = "right") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = c(-1.5, 1.5), linetype = "dotted", color = "red", alpha = 0.5) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "blue", alpha = 0.5)

ggsave(file.path(viz_dir, "nes_vs_coverage_scatter.pdf"),
       p_scatter, width = 10, height = 8)
cat("    ✓ Created: nes_vs_coverage_scatter.pdf\n")

# 10.3 OxPhos Leading Edge Comparison (Heatmap-style)
if (nrow(oxphos_comparison) > 0) {
  cat("  Creating OxPhos leading edge comparison...\n")
  
  # Create a simple comparison plot
  oxphos_comparison$Contrast <- factor(oxphos_comparison$Contrast,
                                       levels = c("NHA-vs-Control", "Treated-vs-Control", "Treated-vs-NHA"))
  
  p_oxphos <- ggplot(oxphos_comparison, aes(x = Contrast, y = "Oxidative Phosphorylation")) +
    geom_point(aes(size = Percent_Leading, color = NES), alpha = 0.8) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red",
                          midpoint = 0, name = "NES") +
    scale_size_continuous(range = c(5, 15), name = "Leading Edge %") +
    labs(title = "Oxidative Phosphorylation Across Contrasts",
         subtitle = "Size = Leading edge coverage | Color = NES direction and magnitude",
         x = "Contrast",
         y = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major.y = element_blank()) +
    geom_text(aes(label = paste0("NES=", NES, "\n", Percent_Leading, "%")),
              size = 3, color = "black")
  
  ggsave(file.path(viz_dir, "oxphos_across_contrasts.pdf"),
         p_oxphos, width = 8, height = 6)
  cat("    ✓ Created: oxphos_across_contrasts.pdf\n")
}

# 10.4 Gene-Pathway Network (Simplified)
cat("  Creating gene-pathway network visualization...\n")

if (nrow(gene_summary) > 0) {
  # Take top 20 genes and their pathways
  top_genes_network <- head(gene_summary, 20)
  
  # Create a simple bipartite plot
  pathway_assignments <- list()
  for (i in 1:nrow(top_genes_network)) {
    gene <- top_genes_network$Gene[i]
    pathways <- unlist(strsplit(as.character(gene_summary$Pathways[i]), ";"))
    pathways <- gsub(" \\.\\.\\. \\(\\+.*", "", pathways)  # Clean up
    pathway_assignments[[gene]] <- pathways[1:min(3, length(pathways))]  # Limit to 3
  }
  
  # Create edge list
  edges <- data.frame()
  for (gene in names(pathway_assignments)) {
    for (pathway in pathway_assignments[[gene]]) {
      edges <- rbind(edges, data.frame(From = gene, To = pathway, stringsAsFactors = FALSE))
    }
  }
  
  # Save edge list for external network visualization
  write.csv(edges, file.path(viz_dir, "gene_pathway_network_edges.csv"), row.names = FALSE)
  cat("    ✓ Created: gene_pathway_network_edges.csv (for Cytoscape/Gephi)\n")
  
  # Simple R visualization
  if (nrow(edges) > 0) {
    # Count connections
    gene_counts <- table(edges$From)
    pathway_counts <- table(edges$To)
    
    connection_df <- data.frame(
      Node = c(names(gene_counts), names(pathway_counts)),
      Type = c(rep("Gene", length(gene_counts)), rep("Pathway", length(pathway_counts))),
      Connections = c(as.numeric(gene_counts), as.numeric(pathway_counts)),
      stringsAsFactors = FALSE
    )
    
    p_network <- ggplot(connection_df, aes(x = reorder(Node, Connections), y = Connections, fill = Type)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      scale_fill_brewer(palette = "Set2") +
      labs(title = "Gene-Pathway Connection Network",
           subtitle = "Top 20 leading edge genes and their pathway memberships",
           x = "Gene or Pathway",
           y = "Number of Connections") +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 8))
    
    ggsave(file.path(viz_dir, "gene_pathway_connections.pdf"),
           p_network, width = 10, height = 8)
    cat("    ✓ Created: gene_pathway_connections.pdf\n")
  }
}

# 10.5 Leading Edge Size Distribution
cat("  Creating leading edge size distribution...\n")

p_dist <- ggplot(enriched_results, aes(x = leading_edge_percent, fill = Contrast)) +
  geom_histogram(bins = 20, alpha = 0.7, position = "identity") +
  facet_wrap(~ Contrast, ncol = 1) +
  labs(title = "Distribution of Leading Edge Coverage",
       subtitle = "How much of each pathway contributes to the enrichment signal",
       x = "Leading Edge Coverage (%)",
       y = "Number of Pathways") +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 50, linetype = "dashed", color = "red", alpha = 0.7)

ggsave(file.path(viz_dir, "leading_edge_distribution.pdf"),
       p_dist, width = 10, height = 8)
cat("    ✓ Created: leading_edge_distribution.pdf\n")

cat("\n  All visualizations saved to: ", viz_dir, "\n", sep = "")
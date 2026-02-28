#!/usr/bin/env Rscript
# Title: Step 2 - Model Design and Linear Model Fitting
# Input: intermediate_data/01_validated_data.rds
# Output: Fitted model objects with design matrices

cat("=== STEP 2: MODEL SETUP AND FITTING ===\n")

# Load required packages
cat("Loading required packages...\n")
required_packages <- c("limma", "matrixStats")
missing_packages <- required_packages[!required_packages %in% installed.packages()]
if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}
invisible(lapply(required_packages, library, character.only = TRUE))

# Load shared setup functions
source("config-setup_functions.R")

# Set up directories and logging
setup_directories()
log_file <- setup_logging("Step 2: Model Setup and Fitting")

# 1. LOAD VALIDATED DATA
cat("1. Loading validated data from Step 1...\n")
check_file_exists("intermediate_data/01_validated_data.rds", "Step 1")
validated_data <- load_data_safely("intermediate_data/01_validated_data.rds", "Step 1")

expression <- validated_data$expression
metadata <- validated_data$metadata
sanity_checks <- validated_data$sanity_checks

# Verify data is still valid
if (!sanity_checks$has_condition) {
  stop("ERROR: Metadata does not contain 'Condition' column!")
}

cat("Data loaded successfully:\n")
cat("  Proteins:", nrow(expression), "\n")
cat("  Samples:", ncol(expression), "\n")
cat("  Conditions:", paste(levels(metadata$Condition), collapse = ", "), "\n")

# 2. CHECK SAMPLE SIZES
cat("\n2. Checking sample sizes per condition...\n")

condition_counts <- table(metadata$Condition)
min_samples <- min(condition_counts)

if (min_samples < 3) {
  warning(paste("WARNING: Some conditions have fewer than 3 samples!",
                "This may reduce statistical power."))
}

cat("Samples per condition:\n")
print(condition_counts)

# 3. CREATE DESIGN MATRIX
cat("\n3. Creating design matrix...\n")

# Check for perfect confounding
if ("Batch" %in% colnames(metadata)) {
  confounding_table <- table(metadata$Condition, metadata$Batch)
  confounding_score <- sum(apply(confounding_table > 0, 2, sum) == 1) / ncol(confounding_table)
  
  cat("\nCONFOUNDING STATUS:\n")
  if (confounding_score == 1) {
    cat("  PERFECT CONFOUNDING DETECTED!\n")
    cat("  Each condition is completely confounded with batch.\n")
    cat("  Cannot separate biological from technical effects.\n")
    cat("  Model: expression ~ Condition (batch effects confounded)\n")
    
    # Note in metadata
    metadata$confounding_note <- "Perfect confounding: Condition = Batch"
  } else if (confounding_score > 0.5) {
    cat("  HIGH CONFOUNDING DETECTED!\n")
    cat("  Most conditions are confounded with batch.\n")
  } else {
    cat("  Minimal or no confounding detected.\n")
  }
}

# Create design matrix
cat("\nCreating design matrix with conditions...\n")
design <- model.matrix(~ 0 + Condition, data = metadata)
colnames(design) <- gsub("Condition", "", colnames(design))

# Check design matrix properties
cat("\nDesign matrix properties:\n")
cat("  Dimensions:", dim(design), "\n")
cat("  Rank:", qr(design)$rank, "\n")
cat("  Condition matrix:\n")
print(round(crossprod(design), 2))

# Sanity check: ensure columns sum to sample sizes
col_sums <- colSums(design)
cat("\nColumn sums (should equal group sizes):\n")
print(col_sums)

# 4. CHECK FOR COLLINEARITY
cat("\n4. Checking for design matrix issues...\n")

# Calculate condition correlations
if (ncol(design) > 1) {
  condition_cor <- cor(design)
  diag(condition_cor) <- NA
  
  max_cor <- max(abs(condition_cor), na.rm = TRUE)
  cat("  Maximum condition correlation:", round(max_cor, 3), "\n")
  
  if (max_cor > 0.8) {
    warning("WARNING: High correlation between conditions in design matrix!")
  }
}

# 5. FIT LINEAR MODELS
cat("\n5. Fitting linear models...\n")

# Check for extreme values that might affect model fitting
expression_matrix <- as.matrix(expression)
expression_stats <- list(
  range = range(expression_matrix, na.rm = TRUE),
  mean = mean(expression_matrix, na.rm = TRUE),
  sd = sd(expression_matrix, na.rm = TRUE)
)

cat("Expression data statistics:\n")
cat("  Range:", round(expression_stats$range[1], 2), "to", 
    round(expression_stats$range[2], 2), "\n")
cat("  Mean:", round(expression_stats$mean, 2), "\n")
cat("  SD:", round(expression_stats$sd, 2), "\n")

# Check for variance heterogeneity
protein_vars <- rowVars(expression_matrix, na.rm = TRUE)
cat("Protein variance statistics:\n")
cat("  Min:", round(min(protein_vars, na.rm = TRUE), 3), "\n")
cat("  Median:", round(median(protein_vars, na.rm = TRUE), 3), "\n")
cat("  Max:", round(max(protein_vars, na.rm = TRUE), 3), "\n")

# Flag proteins with very low variance
low_var_proteins <- sum(protein_vars < 0.001, na.rm = TRUE)
if (low_var_proteins > 0) {
  cat("  WARNING:", low_var_proteins, "proteins have variance < 0.001\n")
}

# Fit linear model
cat("\nFitting linear models with limma...\n")
fit <- tryCatch({
  lmFit(expression, design)
}, error = function(e) {
  stop("ERROR in lmFit: ", e$message)
})

cat("  Model fitting successful\n")
cat("  Coefficients:", ncol(fit$coefficients), "\n")

# 6. CREATE CONTRASTS
cat("\n6. Setting up contrasts...\n")

# Load parameters
source("config/analysis_parameters.R")
contrast_defs <- PARAMS$contrasts

# Create contrast matrix
contrast_names <- names(contrast_defs)
contrasts <- matrix(0, nrow = ncol(design), ncol = length(contrast_names))
colnames(contrasts) <- contrast_names
rownames(contrasts) <- colnames(design)

for (i in seq_along(contrast_names)) {
  contrast_name <- contrast_names[i]
  groups <- contrast_defs[[contrast_name]]
  
  if (!all(groups %in% colnames(design))) {
    stop(paste("ERROR: Contrast", contrast_name, 
               "contains undefined group(s):", 
               paste(setdiff(groups, colnames(design)), collapse = ", ")))
  }
  
  contrasts[groups[1], i] <- 1
  contrasts[groups[2], i] <- -1
}

cat("Contrast matrix:\n")
print(contrasts)

# Fit contrasts
cat("\nFitting contrasts...\n")
fit_contrasts <- tryCatch({
  contrasts.fit(fit, contrasts)
}, error = function(e) {
  stop("ERROR in contrasts.fit: ", e$message)
})

# Apply empirical Bayes moderation
cat("Applying empirical Bayes moderation...\n")
fit_ebayes <- tryCatch({
  eBayes(fit_contrasts)
}, error = function(e) {
  stop("ERROR in eBayes: ", e$message)
})

cat("  Empirical Bayes moderation successful\n")

# 7. MODEL DIAGNOSTICS
cat("\n7. Running model diagnostics...\n")

# Check residuals
cat("  Calculating model residuals...\n")
residuals <- residuals.MArrayLM(fit_ebayes, expression)
residual_stats <- list(
  mean = mean(residuals, na.rm = TRUE),
  sd = sd(residuals, na.rm = TRUE),
  qq = quantile(residuals, probs = c(0.01, 0.05, 0.5, 0.95, 0.99), na.rm = TRUE)
)

cat("  Residual statistics:\n")
cat("    Mean:", round(residual_stats$mean, 4), "(should be ~0)\n")
cat("    SD:", round(residual_stats$sd, 4), "\n")
cat("    Quantiles:", round(residual_stats$qq, 4), "\n")

# Check variance inflation (if applicable)
if (ncol(design) > 2) {
  vif_diag <- diag(solve(cor(design)))
  cat("  Variance Inflation Factors (VIF):\n")
  print(round(vif_diag, 2))
  
  if (any(vif_diag > 10)) {
    warning("WARNING: High VIF (>10) detected - potential multicollinearity!")
  }
}

# 8. SAVE MODEL OBJECTS
cat("\n8. Saving model objects...\n")

model_objects <- list(
  design = design,
  contrasts = contrasts,
  fit = fit,
  fit_contrasts = fit_contrasts,
  fit_ebayes = fit_ebayes,
  condition_counts = condition_counts,
  expression_stats = expression_stats,
  protein_vars = protein_vars,
  residuals = residuals,
  metadata = metadata,
  parameters = PARAMS,
  timestamp = Sys.time()
)

saveRDS(model_objects, "intermediate_data/02_model_fit.rds")

# 9. GENERATE MODEL DIAGNOSTIC REPORT
cat("\n9. Generating model diagnostic report...\n")

model_report <- paste(
  "MODEL SETUP AND DIAGNOSTICS REPORT",
  "===================================",
  paste("Date:", Sys.time()),
  "",
  "EXPERIMENTAL DESIGN:",
  paste("  Conditions:", paste(colnames(design), collapse = ", ")),
  paste("  Samples per condition:", paste(paste(names(condition_counts), condition_counts, sep = ": "), collapse = ", ")),
  "",
  "DESIGN MATRIX:",
  paste("  Dimensions:", paste(dim(design), collapse = " × ")),
  paste("  Rank:", qr(design)$rank),
  paste("  Column sums:", paste(round(col_sums, 2), collapse = ", ")),
  "",
  "CONFOUNDING STATUS:",
  if (exists("confounding_score")) {
    c(paste("  Confounding score:", round(confounding_score, 3)),
      if (confounding_score == 1) {
        "  STATUS: PERFECT CONFOUNDING - Biological and technical effects cannot be separated"
      } else if (confounding_score > 0.5) {
        "  STATUS: HIGH CONFOUNDING - Results may be strongly influenced by batch effects"
      } else {
        "  STATUS: Minimal confounding"
      }
    )
  } else {
    "  Batch information not available for confounding assessment"
  },
  "",
  "EXPRESSION DATA:",
  paste("  Range:", round(expression_stats$range[1], 2), "to", round(expression_stats$range[2], 2)),
  paste("  Mean ± SD:", round(expression_stats$mean, 2), "±", round(expression_stats$sd, 2)),
  paste("  Protein variance: min=", round(min(protein_vars, na.rm = TRUE), 4), 
        " median=", round(median(protein_vars, na.rm = TRUE), 4),
        " max=", round(max(protein_vars, na.rm = TRUE), 4)),
  "",
  "CONTRASTS DEFINED:",
  paste("  ", paste(names(contrast_defs), collapse = ", ")),
  "",
  "MODEL FITTING:",
  "  ✓ Linear models fitted successfully",
  "  ✓ Contrasts defined and fitted",
  "  ✓ Empirical Bayes moderation applied",
  "",
  "RESIDUAL DIAGNOSTICS:",
  paste("  Residual mean:", round(residual_stats$mean, 4), "(ideal: 0)"),
  paste("  Residual SD:", round(residual_stats$sd, 4)),
  paste("  Residual quantiles (1%, 5%, 50%, 95%, 99%):", 
        paste(round(residual_stats$qq, 4), collapse = ", ")),
  "",
  "MODEL ASSUMPTIONS:",
  if (exists("vif_diag") && any(vif_diag > 10)) {
    paste("  ⚠  Multicollinearity detected (VIF > 10):", 
          paste(names(vif_diag)[vif_diag > 10], collapse = ", "))
  } else {
    "  ✓ No severe multicollinearity detected"
  },
  if (low_var_proteins > 0) {
    paste("  ⚠ ", low_var_proteins, "proteins have very low variance (< 0.001)")
  } else {
    "  ✓ All proteins have reasonable variance"
  },
  "",
  "OUTPUT FILES:",
  "  intermediate_data/02_model_fit.rds - All model objects",
  "",
  "NEXT STEP:",
  "  Run Step 3 for differential expression analysis",
  sep = "\n"
)

writeLines(model_report, "intermediate_data/02_model_report.txt")

# 10. CREATE DESIGN VISUALIZATION
cat("\n10. Creating design visualization...\n")

if (PARAMS$create_plots) {
  if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)
  
  pdf("figures/02_design_matrix.pdf", width = 8, height = 6)
  par(mar = c(6, 6, 4, 2))
  
  # Heatmap of design matrix
  image(1:ncol(design), 1:nrow(design), t(design[nrow(design):1, ]),
        col = c("white", "steelblue"),
        xlab = "", ylab = "Samples", main = "Design Matrix",
        xaxt = "n", yaxt = "n")
  axis(1, at = 1:ncol(design), labels = colnames(design), las = 2, cex.axis = 0.8)
  axis(2, at = seq(1, nrow(design), length.out = min(10, nrow(design))), 
       labels = seq(nrow(design), 1, length.out = min(10, nrow(design))), las = 1)
  box()
  
  # Add condition sample sizes
  mtext(paste("Total samples:", nrow(design)), side = 3, line = 0.5, cex = 0.8)
  
  dev.off()
  
  cat("  Design matrix plot saved to: figures/02_design_matrix.pdf\n")
}

# Clean up
cat("\nStep 2 completed successfully!\n")
cat("Model objects saved to: intermediate_data/02_model_fit.rds\n")
cat("Diagnostic report saved to: intermediate_data/02_model_report.txt\n")

cleanup_logging(log_file, "Step 2: Model Setup and Fitting")
#!/usr/bin/env Rscript
# Title: Step 1 - Data Loading and Validation
# Input: step2b_data_imputed.rds
# Output: Validated data structure for downstream analysis

cat("=== STEP 1: DATA VALIDATION ===\n")

# Load required packages
cat("Loading required packages...\n")
required_packages <- c("limma", "qvalue", "data.table", "matrixStats")
missing_packages <- required_packages[!required_packages %in% installed.packages()]
if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}
invisible(lapply(required_packages, library, character.only = TRUE))

# Load shared setup functions
source("config-setup_functions.R")

# Set up directories and logging
setup_directories()
log_file <- setup_logging("Step 1: Data Validation")

# 1. LOAD DATA
cat("1. Loading input data...\n")
if (!file.exists("step2b_data_imputed.rds")) {
  stop("ERROR: Input file 'step2b_data_imputed.rds' not found!")
}

data <- tryCatch({
  readRDS("step2b_data_imputed.rds")
}, error = function(e) {
  stop("ERROR: Failed to load RDS file: ", e$message)
})

# 2. SANITY CHECKS
cat("\n2. Performing data sanity checks...\n")

# Check structure
check_results <- list()

# Check 1: Required components
check_results$has_expression <- "expression" %in% names(data)
check_results$has_metadata <- "metadata" %in% names(data)

if (!all(unlist(check_results[c("has_expression", "has_metadata")]))) {
  stop("ERROR: Data object missing required components!")
}

expression <- data$expression
metadata <- data$metadata

# Check 2: Dimensions match
check_results$n_samples_match <- ncol(expression) == nrow(metadata)
check_results$sample_names_match <- all(colnames(expression) == rownames(metadata))

# Check 3: Missing values
check_results$missing_expression <- sum(is.na(expression))
check_results$missing_metadata <- sum(is.na(metadata))

# Check 4: Metadata structure
check_results$has_condition <- "Condition" %in% colnames(metadata)
if (check_results$has_condition) {
  check_results$condition_levels <- levels(metadata$Condition)
  check_results$condition_counts <- table(metadata$Condition)
}

# Print check results
cat("\nSanity check results:\n")
cat("  Samples in expression:", ncol(expression), "\n")
cat("  Samples in metadata:", nrow(metadata), "\n")
cat("  Proteins:", nrow(expression), "\n")
cat("  Missing values in expression:", check_results$missing_expression, "\n")
cat("  Sample names match:", check_results$sample_names_match, "\n")

if (check_results$has_condition) {
  cat("\nCondition distribution:\n")
  print(check_results$condition_counts)
}

# 3. DATA QUALITY METRICS
cat("\n3. Calculating data quality metrics...\n")

quality_metrics <- list(
  expression_stats = list(
    mean_intensity = mean(as.matrix(expression), na.rm = TRUE),
    median_intensity = median(as.matrix(expression), na.rm = TRUE),
    cv_proteins = apply(expression, 1, function(x) sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)),
    cv_samples = apply(expression, 2, function(x) sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE))
  ),
  metadata_quality = list(
    complete_cases = sum(complete.cases(metadata)),
    unique_conditions = length(unique(metadata$Condition))
  )
)

cat("  Mean protein intensity:", round(quality_metrics$expression_stats$mean_intensity, 2), "\n")
cat("  Median protein intensity:", round(quality_metrics$expression_stats$median_intensity, 2), "\n")
cat("  Average protein CV:", round(mean(quality_metrics$expression_stats$cv_proteins, na.rm = TRUE), 3), "\n")

# 4. DETECT CONFOUNDING
cat("\n4. Checking for confounding...\n")

if (check_results$has_condition && "Batch" %in% colnames(metadata)) {
  confounding_table <- table(metadata$Condition, metadata$Batch)
  cat("\nCondition × Batch table:\n")
  print(confounding_table)
  
  # Calculate confounding score
  condition_per_batch <- apply(confounding_table > 0, 2, sum)
  batch_per_condition <- apply(confounding_table > 0, 1, sum)
  
  confounding_score <- sum(condition_per_batch == 1) / ncol(confounding_table)
  
  cat("\nConfounding analysis:\n")
  cat("  Conditions per batch:", paste(condition_per_batch, collapse = ", "), "\n")
  cat("  Batches per condition:", paste(batch_per_condition, collapse = ", "), "\n")
  cat("  Confounding score:", round(confounding_score, 3), "\n")
  
  if (confounding_score == 1) {
    cat("  WARNING: PERFECT CONFOUNDING DETECTED!\n")
    cat("  Each condition is completely confounded with batch.\n")
  }
}

# 5. SAVE VALIDATED DATA
cat("\n5. Saving validated data...\n")

validated_data <- list(
  expression = expression,
  metadata = metadata,
  sanity_checks = check_results,
  quality_metrics = quality_metrics,
  timestamp = Sys.time()
)

saveRDS(validated_data, "intermediate_data/01_validated_data.rds")

# 6. GENERATE VALIDATION REPORT
cat("\n6. Generating validation report...\n")

validation_report <- paste(
  "DATA VALIDATION REPORT",
  "=====================",
  paste("Date:", Sys.time()),
  paste("Input file: step2b_data_imputed.rds"),
  "",
  "DATA STRUCTURE:",
  paste("  Expression matrix: ", nrow(expression), "proteins ×", ncol(expression), "samples"),
  paste("  Metadata: ", nrow(metadata), "samples with", ncol(metadata), "variables"),
  paste("  Sample names match: ", check_results$sample_names_match),
  "",
  "QUALITY METRICS:",
  paste("  Missing values in expression: ", check_results$missing_expression),
  paste("  Mean protein intensity: ", round(quality_metrics$expression_stats$mean_intensity, 2)),
  paste("  Median protein intensity: ", round(quality_metrics$expression_stats$median_intensity, 2)),
  paste("  Average protein CV: ", round(mean(quality_metrics$expression_stats$cv_proteins, na.rm = TRUE), 3)),
  "",
  "CONFOUNDING ASSESSMENT:",
  if (exists("confounding_table")) {
    c("  Condition × Batch table:",
      capture.output(print(confounding_table)),
      paste("  Confounding score: ", round(confounding_score, 3)),
      if (confounding_score == 1) "  STATUS: PERFECT CONFOUNDING DETECTED" else "  STATUS: Partial or no confounding"
    )
  } else {
    "  Batch information not available for confounding check"
  },
  "",
  "VALIDATION STATUS:",
  if (all(unlist(check_results[1:2]))) "  PASS: Data structure is valid" else "  FAIL: Data structure issues",
  if (check_results$n_samples_match) "  PASS: Sample counts match" else "  FAIL: Sample count mismatch",
  if (check_results$sample_names_match) "  PASS: Sample names match" else "  FAIL: Sample name mismatch",
  if (check_results$missing_expression == 0) "  PASS: No missing values" else paste("  WARNING:", check_results$missing_expression, "missing values"),
  "",
  sep = "\n"
)

writeLines(validation_report, "intermediate_data/01_validation_report.txt")

# Clean up
cat("\nValidation complete!\n")
cat("Output saved to: intermediate_data/01_validated_data.rds\n")
cat("Report saved to: intermediate_data/01_validation_report.txt\n")

cleanup_logging(log_file, "Step 1: Data Validation")
# Run RIbench benchmark evaluation for RINet
#
# Bypasses evaluateBiomarkerTestSets() to avoid spawning a subprocess per
# test set (which reloads the model each time). Instead, loads the model
# once, reads all test set CSVs, batch-predicts, and writes .Rdata files
# in the format evaluateAlgorithmResults() expects.
#
# Prerequisites:
#   install.packages("RIbench")
#   # rinet must be installed (or loaded via devtools::load_all())
#
# Usage:
#   Rscript benchmark/run_benchmark_rinet.R [subset] [--outliers] [--no-log-normal] [--workdir PATH] [--outdir PATH]
#
#   subset (optional): controls which test sets to run
#     - "all" (default)    — run all 5,760 test sets
#     - "test"             — quick test with 3 test sets per biomarker
#     - integer (e.g. 10)  — run N test sets per biomarker
#     - biomarker name     — e.g. "Ca", "AST", "Hb"
#     - distribution type  — e.g. "normal", "skewed", "heavilySkewed", "shifted"
#
#   --outliers (optional): apply Tukey 1.5x IQR outlier removal on log-scale
#     Can appear anywhere in the argument list.
#
#   --no-log-normal (optional): disable log-scaling for analytes with Distribution == "normal"
#     Can appear anywhere in the argument list.
#
#   --workdir PATH (required): path to working directory containing Data/ subfolder
#
#   --outdir PATH (optional): path for Results/ and Evaluation/ output
#     Defaults to the same as --workdir.

library(RIbench)
library(rinet)
library(data.table)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

# Separate flag arguments from positional arguments
all_flags <- c("--outliers", "--tukey", "--iqr", "--no-log-normal", "--workdir", "--outdir")
use_outlier_removal <- any(tolower(args) %in% c("--outliers", "--tukey", "--iqr"))
no_log_normal <- any(args == "--no-log-normal")
workdir_pos <- match("--workdir", args)
outdir_pos <- match("--outdir", args)
positional_args <- args[!args %in% all_flags & !seq_along(args) %in% c(workdir_pos + 1, outdir_pos + 1)]

# Parse subset argument (first positional)
subset_arg <- if (length(positional_args) >= 1) positional_args[1] else "all"
if (subset_arg == "test") {
  subset_val <- 3L
  subset_arg <- "3 (test)"
} else {
  subset_val <- suppressWarnings(as.integer(subset_arg))
  if (is.na(subset_val)) {
    subset_val <- subset_arg  # character subset (biomarker or distribution)
  }
}

if (is.na(workdir_pos)) {
  stop("--workdir is required. Usage: Rscript benchmark/run_benchmark_rinet.R --workdir PATH [--outdir PATH] [subset] [--outliers] [--no-log-normal]")
}
working_dir <- normalizePath(args[workdir_pos + 1], mustWork = TRUE)
output_dir <- if (!is.na(outdir_pos)) normalizePath(args[outdir_pos + 1], mustWork = FALSE) else working_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

algo_name   <- "RINet"
percentiles <- c(0.025, 0.975)

cat("=== RIbench Evaluation for RINet ===\n")
cat("Working directory:", working_dir, "\n")
cat("Output directory:", output_dir, "\n")
cat("Subset:", subset_arg, "\n")
cat("Outlier removal:", if (use_outlier_removal) "ON (Tukey 1.5x IQR on log-scale)" else "OFF", "\n")
cat("Log-scale for normal:", if (no_log_normal) "OFF" else "ON", "\n\n")

# ---------------------------------------------------------------------------
# Step 1: Determine test set definitions
# ---------------------------------------------------------------------------
testsets <- loadTestsetDefinition()

if (is.numeric(subset_val)) {
  testsets <- RIbench:::defineSubset(testsets, N = subset_val)
  testsets <- testsets[testsets$Subset == 1, ]
} else if (is.character(subset_val) && subset_val != "all") {
  if (subset_val %in% c("normal", "skewed", "heavilySkewed", "shifted")) {
    testsets <- testsets[testsets$Distribution == subset_val, ]
  } else if (subset_val %in% unique(testsets$Analyte)) {
    testsets <- testsets[testsets$Analyte == subset_val, ]
  }
}

# ---------------------------------------------------------------------------
# Step 2: Check benchmark data
# ---------------------------------------------------------------------------
data_dir <- file.path(working_dir, "Data")
n_existing <- length(list.files(data_dir, pattern = "\\.csv$", recursive = TRUE))
n_expected <- nrow(testsets)

if (n_existing < n_expected) {
  stop(sprintf(
    "Expected %d test set CSV files in %s, found %d.\nGenerate them first with RIbench::generateBiomarkerTestSets(workingDir = '%s').",
    n_expected, data_dir, n_existing, working_dir
  ))
}
cat(sprintf("Step 1: Found %d test set CSV files.\n\n", n_existing))

# ---------------------------------------------------------------------------
# Step 3: Run RINet
# ---------------------------------------------------------------------------
cat("Step 2: Running RINet on benchmark test sets (in-process, single model load)...\n")

# Create results directory structure
for (analyte in unique(testsets$Analyte)) {
  dir.create(file.path(output_dir, "Results", algo_name, analyte),
             recursive = TRUE, showWarnings = FALSE)
}

# Trigger model load with a dummy prediction
dummy <- tryCatch(
  predict_rinet_1d(exp(rnorm(100, 2, 0.5)), log_scale = TRUE, percentiles = percentiles),
  error = function(e) NULL
)
cat("Model loaded. Processing", nrow(testsets), "test sets...\n")

t_start <- proc.time()

n_success <- 0
n_fail    <- 0

for (i in seq_len(nrow(testsets))) {
  row <- testsets[i, ]
  csv_file <- file.path(working_dir, "Data", row$Analyte,
                        paste0(row$Index, "_", row$Analyte, "_seed_", row$startSeed, ".csv"))

  rdata_file <- file.path(output_dir, "Results", algo_name, row$Analyte,
                          paste0(row$Index, "_", row$Analyte, "_seed_", row$startSeed, "_", algo_name, ".Rdata"))

  # Skip if already computed
  if (file.exists(rdata_file)) {
    n_success <- n_success + 1
    next
  }

  # Determine whether to use log-scale for this test set
  use_log_scale <- !(no_log_normal && row$Distribution == "normal")

  # Read data; remove non-positive values only when using log-scale
  data_vec <- fread(csv_file)$V1
  if (use_log_scale) data_vec <- data_vec[data_vec > 0]

  # Optional outlier removal using Tukey's fences (1.5x IQR) on log-scale
  if (use_outlier_removal && use_log_scale) {
    log_vec <- log(data_vec)
    q1 <- quantile(log_vec, 0.25)
    q3 <- quantile(log_vec, 0.75)
    iqr <- q3 - q1
    keep <- log_vec >= (q1 - 1.5 * iqr) & log_vec <= (q3 + 1.5 * iqr)
    data_vec <- data_vec[keep]
  }

  # CRP: use 95th percentile (one-sided upper RI); others: 2.5th and 97.5th
  pctiles <- if (row$Analyte == "CRP") c(0.05, 0.95) else percentiles

  # Predict
  t_iter <- system.time({
    result <- tryCatch({
      res <- predict_rinet_1d(data_vec, log_scale = use_log_scale, percentiles = pctiles)
      ri <- res[[1]]$reference_interval
      data.frame(Percentile = pctiles, PointEst = as.numeric(ri))
    }, error = function(e) {
      data.frame(Percentile = pctiles, PointEst = NA_real_)
    })
  })

  # RIbench expects CRP to have only the upper RI (1 row).
  # getRIsAllwithoutModel() does c(0, model$PointEst), so PointEst must be scalar.
  fit <- result
  if (row$Analyte == "CRP") {
    fit <- data.frame(Percentile = 0.95, PointEst = result$PointEst[2])
    fit$Runtime <- t_iter[3]
  } else {
    fit$Runtime <- c(t_iter[3], NA)
  }

  # Save
  save(fit, file = rdata_file)

  if (any(is.na(result$PointEst))) {
    n_fail <- n_fail + 1
  } else {
    n_success <- n_success + 1
  }

  # Progress
  if (i %% 10 == 0 || i == nrow(testsets)) {
    cat(sprintf("\r  Progress: %d/%d (%.0f%%)", i, nrow(testsets), 100 * i / nrow(testsets)))
  }
}

t_elapsed <- (proc.time() - t_start)[3]
cat(sprintf("\nDone in %.1f seconds. Success: %d, Failed: %d\n\n", t_elapsed, n_success, n_fail))

# ---------------------------------------------------------------------------
# Step 4: Evaluate results and compute benchmark score
# ---------------------------------------------------------------------------
cat("Step 3: Evaluating results and generating plots...\n")
col_rinet <- rgb(20, 130, 250, maxColorValue = 255)

benchmarkScore <- evaluateAlgorithmResults(
  workingDir = output_dir,
  algoNames  = algo_name,
  subset     = subset_val,
  cols       = col_rinet
)

cat("\n=== Benchmark Score ===\n")
print(benchmarkScore)
cat("\nPlots saved to:", file.path(output_dir, "Evaluation"), "\n")

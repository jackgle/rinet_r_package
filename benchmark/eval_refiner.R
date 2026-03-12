#!/usr/bin/env Rscript
# Evaluate refineR predictions using RIbench.
#
# Usage:
#   Rscript benchmark/eval_refiner.R <predictions_dir> [--workdir PATH] [--outdir PATH] [subset]
#
# Example:
#   Rscript benchmark/eval_refiner.R ../rinet_v1/evaluation/ribench/refineR_predictions \
#     --workdir ../rinet_v1/data/RIbench --outdir benchmark/results_refiner
#
# Arguments:
#   predictions_dir: directory containing refineR CSV prediction files
#   --workdir PATH:  directory containing Data/ subfolder (default: benchmark/results)
#   --outdir PATH:   directory for Results/ and Evaluation/ output (default: same as --workdir)
#   subset:          "all" (default), "test", integer, biomarker name, or distribution type

suppressPackageStartupMessages({
  library(RIbench)
})

args <- commandArgs(trailingOnly = TRUE)

# Parse flags
workdir_pos <- match("--workdir", args)
outdir_pos  <- match("--outdir", args)
flag_indices <- c(
  which(args %in% c("--workdir", "--outdir")),
  if (!is.na(workdir_pos)) workdir_pos + 1,
  if (!is.na(outdir_pos)) outdir_pos + 1
)
positional_args <- if (length(flag_indices) > 0) args[-flag_indices] else args

if (length(positional_args) < 1) {
  stop("Usage: Rscript eval_refiner.R <predictions_dir> [--workdir PATH] [--outdir PATH] [subset]")
}

predictions_dir <- normalizePath(positional_args[1], mustWork = TRUE)
subset_arg <- if (length(positional_args) >= 2) positional_args[2] else "all"

if (subset_arg == "test") {
  subset_val <- 3L
  subset_arg <- "3 (test)"
} else {
  subset_val <- suppressWarnings(as.integer(subset_arg))
  if (is.na(subset_val)) {
    subset_val <- subset_arg
  }
}

working_dir <- if (!is.na(workdir_pos)) normalizePath(args[workdir_pos + 1], mustWork = FALSE) else file.path(getwd(), "benchmark", "results")
output_dir  <- if (!is.na(outdir_pos))  normalizePath(args[outdir_pos + 1], mustWork = FALSE) else working_dir

algo_name <- "refineR"

cat("=== RIbench Evaluation for refineR ===\n")
cat("Predictions dir:", predictions_dir, "\n")
cat("Working directory:", working_dir, "\n")
cat("Output directory:", output_dir, "\n")
cat("Subset:", subset_arg, "\n\n")

# ---------------------------------------------------------------------------
# Step 1: Convert refineR CSVs to .Rdata
# ---------------------------------------------------------------------------
cat("Step 1: Converting refineR prediction CSVs to .Rdata...\n")

csv_files <- list.files(predictions_dir, pattern = "\\.csv$", full.names = TRUE)
cat(sprintf("  Found %d CSV files.\n", length(csv_files)))

n_converted <- 0
for (csv_file in csv_files) {
  base  <- tools::file_path_sans_ext(basename(csv_file))
  parts <- strsplit(base, "_")[[1]]
  index   <- parts[1]
  analyte <- parts[2]
  seed    <- parts[4]

  df <- read.csv(csv_file, stringsAsFactors = FALSE)

  if (analyte == "CRP") {
    fit <- data.frame(
      Percentile = 0.975,
      PointEst   = df$PointEst[df$Percentile == 0.975]
    )
  } else {
    fit <- data.frame(
      Percentile = df$Percentile,
      PointEst   = df$PointEst
    )
  }
  fit$Runtime <- NA

  out_subdir <- file.path(output_dir, "Results", algo_name, analyte)
  dir.create(out_subdir, recursive = TRUE, showWarnings = FALSE)
  rdata_file <- file.path(out_subdir,
    paste0(index, "_", analyte, "_seed_", seed, "_", algo_name, ".Rdata"))
  save(fit, file = rdata_file)
  n_converted <- n_converted + 1
}
cat(sprintf("  Converted %d files.\n\n", n_converted))

# ---------------------------------------------------------------------------
# Step 2: Generate benchmark test sets if needed
# ---------------------------------------------------------------------------
cat("Step 2: Checking benchmark test sets...\n")
testsets <- loadTestsetDefinition()

if (is.numeric(subset_val)) {
  testsets <- testsets[testsets$NTCperBiomarker == subset_val, ]
} else if (subset_val != "all") {
  if (subset_val %in% testsets$Analyte) {
    testsets <- testsets[testsets$Analyte == subset_val, ]
  } else if (subset_val %in% testsets$Distribution) {
    testsets <- testsets[testsets$Distribution == subset_val, ]
  }
}

data_dir <- file.path(working_dir, "Data")
n_existing <- length(list.files(data_dir, pattern = "\\.csv$", recursive = TRUE))
n_expected <- nrow(testsets)

if (n_existing >= n_expected) {
  cat(sprintf("  Already exist (%d files), skipping generation.\n\n", n_existing))
} else {
  cat(sprintf("  Generating %d test sets...\n", n_expected))
  generateBiomarkerTestSets(
    workingDir = working_dir,
    subset     = subset_val
  )
  cat("  Done.\n\n")
}

# ---------------------------------------------------------------------------
# Step 3: Evaluate
# ---------------------------------------------------------------------------
cat("Step 3: Evaluating results and generating plots...\n")
col_refiner <- rgb(180, 60, 60, maxColorValue = 255)

benchmarkScore <- evaluateAlgorithmResults(
  workingDir = output_dir,
  algoNames  = algo_name,
  subset     = subset_val,
  cols       = col_refiner
)

cat("\n=== Benchmark Score ===\n")
print(benchmarkScore)
cat("\nPlots saved to:", file.path(output_dir, "Evaluation"), "\n")

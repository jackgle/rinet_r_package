# Run RIbench benchmark evaluation for RINet
#
# This script generates the RIbench test sets, runs RINet on each one,
# and evaluates the results to produce benchmark scores and plots.
#
# Prerequisites:
#   install.packages("RIbench")
#   # rinet must be installed (or loaded via devtools::load_all())
#
# Usage:
#   Rscript benchmark/run_benchmark.R [subset]
#
#   subset (optional): controls which test sets to run
#     - integer (e.g. 3)  — run N test sets per biomarker (good for quick test)
#     - "all"              — run all 5,760 test sets (takes hours)
#     - biomarker name     — e.g. "Ca", "AST", "Hb"
#     - distribution type  — e.g. "normal", "skewed", "heavilySkewed", "shifted"
#   Default: 3 (quick test run)

library(RIbench)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

# Parse subset argument
subset_arg <- if (length(args) >= 1) args[1] else "3"
subset_val <- suppressWarnings(as.integer(subset_arg))
if (is.na(subset_val)) {
  subset_val <- subset_arg  # character subset (biomarker or distribution)
}

working_dir <- file.path(getwd(), "benchmark", "results")
dir.create(working_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== RIbench Evaluation for RINet ===\n")
cat("Working directory:", working_dir, "\n")
cat("Subset:", subset_arg, "\n\n")

# ---------------------------------------------------------------------------
# Step 1: Generate benchmark test sets
# ---------------------------------------------------------------------------
cat("Step 1: Generating benchmark test sets...\n")
generateBiomarkerTestSets(workingDir = working_dir, subset = subset_val)
cat("Done.\n\n")

# ---------------------------------------------------------------------------
# Step 2: Evaluate test sets using RINet
# ---------------------------------------------------------------------------
cat("Step 2: Running RINet on benchmark test sets...\n")
wrapper_path <- file.path(getwd(), "benchmark", "rinet_wrapper.R")

progress <- evaluateBiomarkerTestSets(
  workingDir         = working_dir,
  algoName           = "RINet",
  algoFunction       = "estimateRIs",
  libs               = c("rinet"),
  sourceFiles        = list(wrapper_path),
  requirePercentiles = TRUE,
  subset             = subset_val,
  timeLimit          = 14400
)

# Report any failures
if (!is.null(progress) && nrow(progress) > 0) {
  cat("\nWarning: Some test sets failed:\n")
  print(progress)
} else {
  cat("All test sets completed successfully.\n")
}
cat("\n")

# ---------------------------------------------------------------------------
# Step 3: Evaluate results and compute benchmark score
# ---------------------------------------------------------------------------
cat("Step 3: Evaluating results and generating plots...\n")
col_rinet <- rgb(20, 130, 250, maxColorValue = 255)

benchmarkScore <- evaluateAlgorithmResults(
  workingDir = working_dir,
  algoNames  = "RINet",
  subset     = subset_val,
  cols       = col_rinet
)

cat("\n=== Benchmark Score ===\n")
print(benchmarkScore)
cat("\nPlots saved to:", file.path(working_dir, "Evaluation"), "\n")

#!/usr/bin/env Rscript
# Compute accuracy of RINet and refineR on RIbench at a given z-score deviation threshold.
#
# Usage:
#   Rscript benchmark/accuracy.R [--rinet-dir PATH] [--refiner-dir PATH] [--workdir PATH] [--threshold T]
#
# Defaults:
#   --rinet-dir   benchmark/results_rinet
#   --refiner-dir benchmark/results_refiner
#   --workdir     ../RIbench
#   --threshold   0.1

suppressPackageStartupMessages(library(RIbench))

args <- commandArgs(trailingOnly = TRUE)

parse_flag <- function(flag, default) {
  pos <- match(flag, args)
  if (!is.na(pos) && pos < length(args)) args[pos + 1] else default
}

rinet_dir   <- parse_flag("--rinet-dir",   "benchmark/results_rinet")
refiner_dir <- parse_flag("--refiner-dir", "benchmark/results_refiner")
workdir     <- parse_flag("--workdir",     "../RIbench")
threshold   <- as.numeric(parse_flag("--threshold", "0.1"))

cat(sprintf("Threshold: %.2f (abs z-score deviation)\n\n", threshold))

# Load results via RIbench
algos <- list(
  RINet   = rinet_dir,
  refineR = refiner_dir
)

results <- list()
for (name in names(algos)) {
  res <- readResultsAndComputeErrors(
    workingDir = algos[[name]],
    algoName   = name,
    subset     = "all"
  )
  results[[name]] <- res$perfDf
}

# Compute accuracy per distribution type and overall
dist_types <- c("normal", "skewed", "heavilySkewed", "shifted")

rows <- list()
for (name in names(results)) {
  df <- results[[name]]

  for (dt in dist_types) {
    sub <- df[df$Distribution == dt, ]
    valid <- sub$zzDevAbs_Ov[!is.na(sub$zzDevAbs_Ov)]
    n <- length(valid)
    correct <- sum(valid <= threshold)
    rows <- c(rows, list(data.frame(
      Method = name, Distribution = dt,
      Correct = correct, Total = n,
      Accuracy = if (n > 0) correct / n else NA
    )))
  }

  # Overall
  valid <- df$zzDevAbs_Ov[!is.na(df$zzDevAbs_Ov)]
  n <- length(valid)
  correct <- sum(valid <= threshold)
  rows <- c(rows, list(data.frame(
    Method = name, Distribution = "Overall",
    Correct = correct, Total = n,
    Accuracy = if (n > 0) correct / n else NA
  )))
}

out <- do.call(rbind, rows)
out$Accuracy <- sprintf("%.1f%%", out$Accuracy * 100)

cat("=== Accuracy by Distribution (|z-dev| <= threshold) ===\n")
print(out, row.names = FALSE)

# Breakdown by analyte
cat("\n=== Accuracy by Analyte ===\n")
analyte_rows <- list()
for (name in names(results)) {
  df <- results[[name]]
  for (a in unique(df$Analyte.x)) {
    sub <- df[df$Analyte.x == a, ]
    valid <- sub$zzDevAbs_Ov[!is.na(sub$zzDevAbs_Ov)]
    n <- length(valid)
    correct <- sum(valid <= threshold)
    med_z <- median(valid)
    analyte_rows <- c(analyte_rows, list(data.frame(
      Method = name, Analyte = a,
      Correct = correct, Total = n,
      Accuracy = if (n > 0) sprintf("%.1f%%", 100 * correct / n) else "NA",
      Median_zDev = sprintf("%.3f", med_z)
    )))
  }
}
analyte_out <- do.call(rbind, analyte_rows)
print(analyte_out, row.names = FALSE)

# Breakdown by analyte x distribution (shifted only, to diagnose 0%)
cat("\n=== Shifted distribution by Analyte ===\n")
shifted_rows <- list()
for (name in names(results)) {
  df <- results[[name]]
  sub <- df[df$Distribution == "shifted", ]
  for (a in unique(sub$Analyte.x)) {
    s <- sub[sub$Analyte.x == a, ]
    valid <- s$zzDevAbs_Ov[!is.na(s$zzDevAbs_Ov)]
    n <- length(valid)
    correct <- sum(valid <= threshold)
    shifted_rows <- c(shifted_rows, list(data.frame(
      Method = name, Analyte = a,
      Correct = correct, Total = n,
      Accuracy = if (n > 0) sprintf("%.1f%%", 100 * correct / n) else "NA",
      Median_zDev = sprintf("%.3f", median(valid)),
      Mean_zDev = sprintf("%.3f", mean(valid))
    )))
  }
}
shifted_out <- do.call(rbind, shifted_rows)
print(shifted_out, row.names = FALSE)

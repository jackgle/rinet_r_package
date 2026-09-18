#!/usr/bin/env Rscript
# Plot a grid of RIbench test sets with RINet and refineR predictions.
#
# Usage:
#   Rscript benchmark/plot_comparison.R [--rinet PATH] [--refiner PATH] [--workdir PATH] [--outfile PATH] [--n N] [--seed S] [--analyte NAME]
#
# Defaults:
#   --rinet   benchmark/results_rinet/Results/RINet
#   --refiner benchmark/results_refiner/Results/refineR
#   --workdir ../RIbench
#   --outfile benchmark/comparison_plot.png
#   --n       16
#   --seed    42
#   --log     (flag) plot histograms in log-space

suppressPackageStartupMessages(library(RIbench))

args <- commandArgs(trailingOnly = TRUE)

parse_flag <- function(flag, default) {
  pos <- match(flag, args)
  if (!is.na(pos) && pos < length(args)) args[pos + 1] else default
}

rinet_dir   <- parse_flag("--rinet",   "benchmark/results_rinet/Results/RINet")
refiner_dir <- parse_flag("--refiner", "benchmark/results_refiner/Results/refineR")
workdir     <- parse_flag("--workdir", "../RIbench")
outfile     <- parse_flag("--outfile", "benchmark/comparison_plot.png")
n_plots     <- as.integer(parse_flag("--n", "16"))
seed        <- as.integer(parse_flag("--seed", "42"))
analyte_filter <- parse_flag("--analyte", NULL)
use_log        <- "--log" %in% args

spec_file <- file.path(workdir, "SpecificationTestSets.csv")
data_dir  <- file.path(workdir, "Data")

stopifnot(file.exists(spec_file), dir.exists(data_dir))

spec <- read.csv(spec_file, stringsAsFactors = FALSE)
spec$key <- paste0(spec$Index, "_", spec$Analyte, "_seed_", spec$startSeed)

# Find test sets that have both RINet and refineR results
rinet_files  <- list.files(rinet_dir, pattern = "\\.Rdata$", recursive = TRUE)
refiner_files <- list.files(refiner_dir, pattern = "\\.Rdata$", recursive = TRUE)

strip_algo <- function(f) sub("_(RINet|refineR)\\.Rdata$", "", basename(f))
rinet_keys   <- strip_algo(rinet_files)
refiner_keys <- strip_algo(refiner_files)
common_keys  <- intersect(rinet_keys, refiner_keys)
common_keys  <- common_keys[common_keys %in% spec$key]

# Optional analyte filter
if (!is.null(analyte_filter)) {
  get_analyte <- function(k) strsplit(k, "_")[[1]][2]
  common_keys <- common_keys[sapply(common_keys, get_analyte) == analyte_filter]
}

cat(sprintf("Found %d RINet, %d refineR results, %d in common%s.\n",
            length(rinet_keys), length(refiner_keys), length(common_keys),
            if (!is.null(analyte_filter)) paste0(" (", analyte_filter, ")") else ""))

set.seed(seed)
n_plots <- min(n_plots, length(common_keys))
sampled <- sample(common_keys, n_plots)

# Load an .Rdata result file and return the PointEst values
load_result <- function(dir, key, algo) {
  parts <- strsplit(key, "_")[[1]]
  analyte <- parts[2]
  f <- file.path(dir, analyte, paste0(key, "_", algo, ".Rdata"))
  if (!file.exists(f)) return(NULL)
  env <- new.env()
  load(f, envir = env)
  fit <- env$fit
  if (any(is.na(fit$PointEst))) return(NULL)
  fit$PointEst
}

# Z-score deviation using the Box-Cox transform (same as RIbench::computePerfMeas)
z_score_dev <- function(pred, row) {
  if (is.null(pred)) return(c(NA, NA))
  lambda <- row$nonp_lambda
  mu     <- row$nonp_mu
  sigma  <- row$nonp_sigma
  shift  <- row$nonp_shift

  bc <- function(x) {
    x <- pmax(x - shift, 1e-20)
    if (lambda == 0) log(x) else (x^lambda - 1) / lambda
  }

  gt_z_lower <- (bc(row$GT_LRL) - mu) / sigma
  gt_z_upper <- (bc(row$GT_URL) - mu) / sigma

  if (row$Analyte == "CRP") {
    est_z_upper <- (bc(pred[2]) - mu) / sigma
    return(c(lower = NA, upper = abs(gt_z_upper - est_z_upper)))
  }

  est_z_lower <- (bc(pred[1]) - mu) / sigma
  est_z_upper <- (bc(pred[2]) - mu) / sigma
  c(lower = abs(gt_z_lower - est_z_lower),
    upper = abs(gt_z_upper - est_z_upper))
}

# Grid layout
nc <- 4
nr <- ceiling(n_plots / nc)

png(outfile, width = 3.2 * nc, height = 2.8 * nr, units = "in", res = 150)
par(mfrow = c(nr, nc), mar = c(2.5, 1.5, 2, 0.5), mgp = c(1.2, 0.4, 0))

for (key in sampled) {
  row <- spec[spec$key == key, ]
  analyte <- row$Analyte
  gt_lower <- row$GT_LRL
  gt_upper <- row$GT_URL

  csv_file <- file.path(data_dir, analyte, paste0(key, ".csv"))
  vals <- scan(csv_file, quiet = TRUE)

  rinet_pred   <- load_result(rinet_dir, key, "RINet")
  refiner_pred <- load_result(refiner_dir, key, "refineR")

  # For CRP, only upper is predicted; prepend 0 for plotting
  if (analyte == "CRP") {
    if (length(rinet_pred) == 1)   rinet_pred   <- c(0, rinet_pred)
    if (length(refiner_pred) == 1) refiner_pred <- c(0, refiner_pred)
  }

  # Apply log-transform for plotting if requested
  plot_vals <- vals
  plot_gt_lower  <- gt_lower
  plot_gt_upper  <- gt_upper
  plot_rinet     <- rinet_pred
  plot_refiner   <- refiner_pred

  if (use_log) {
    plot_vals    <- log(plot_vals[plot_vals > 0])
    plot_gt_lower <- log(max(gt_lower, 1e-20))
    plot_gt_upper <- log(max(gt_upper, 1e-20))
    if (!is.null(plot_rinet))   plot_rinet   <- log(pmax(plot_rinet, 1e-20))
    if (!is.null(plot_refiner)) plot_refiner <- log(pmax(plot_refiner, 1e-20))
  }

  # Trim tails for plotting only
  plot_vals <- plot_vals[plot_vals >= quantile(plot_vals, 0.005) &
                         plot_vals <= quantile(plot_vals, 0.995)]

  # Histogram
  breaks <- pretty(range(plot_vals), n = 80)
  x_label <- paste0(analyte, " (", row$Distribution, ")")
  if (use_log) x_label <- paste0("log(", analyte, ") (", row$Distribution, ")")
  hist(plot_vals, breaks = breaks, col = "grey85", border = "grey70",
       main = "", xlab = x_label,
       ylab = "", yaxt = "n", freq = FALSE, cex.axis = 0.7, cex.lab = 0.8)

  # Ground truth
  abline(v = plot_gt_lower, col = "black", lwd = 2, lty = 2)
  abline(v = plot_gt_upper, col = "black", lwd = 2, lty = 2)

  # RINet
  if (!is.null(plot_rinet)) {
    abline(v = plot_rinet[1], col = "dodgerblue", lwd = 2)
    abline(v = plot_rinet[2], col = "dodgerblue", lwd = 2)
  }

  # refineR
  if (!is.null(plot_refiner)) {
    abline(v = plot_refiner[1], col = "firebrick", lwd = 2)
    abline(v = plot_refiner[2], col = "firebrick", lwd = 2)
  }

  # Title: absolute z-score deviation (mean of lower + upper, or upper-only for CRP)
  ri_z  <- z_score_dev(rinet_pred, row)
  rf_z  <- z_score_dev(refiner_pred, row)
  ri_zz <- mean(ri_z, na.rm = TRUE)
  rf_zz <- mean(rf_z, na.rm = TRUE)
  title(sprintf("R:%.2f | rR:%.2f  (|z-dev|)", ri_zz, rf_zz), cex.main = 0.75)
}

# Legend on a blank panel if there's space
if (n_plots %% nc != 0) {
  plot.new()
  legend("center", legend = c("Ground truth", "RINet", "refineR"),
         col = c("black", "dodgerblue", "firebrick"),
         lwd = 2, lty = c(2, 1, 1), cex = 1.1, bty = "n")
}

dev.off()
cat("Saved to", outfile, "\n")

# Diagnostic plots for RINet predictions on RIbench test sets
#
# For each distribution type, picks a representative test set and plots:
#   - Histogram of the data
#   - Ground truth reference interval (green)
#   - RINet predicted reference interval (blue)
#   - Predicted reference distribution overlay
#
# Usage:
#   Rscript benchmark/plot_diagnostics.R [subset]
#   Default subset: "test" (3 per biomarker)

library(RIbench)
library(rinet)
library(data.table)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
subset_arg <- if (length(args) >= 1) args[1] else "test"
if (subset_arg == "test") {
  subset_val <- 3L
} else {
  subset_val <- suppressWarnings(as.integer(subset_arg))
  if (is.na(subset_val)) subset_val <- subset_arg
}

working_dir <- file.path(getwd(), "benchmark", "results")
percentiles <- c(0.025, 0.975)

# Load test set definitions
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

# Warm up model
dummy <- tryCatch(
  predict_rinet_1d(exp(rnorm(100, 2, 0.5)), log_scale = TRUE, percentiles = percentiles),
  error = function(e) NULL
)

# ---------------------------------------------------------------------------
# Pick representative test sets: one per (distribution type x biomarker)
# Also vary pathological fraction to see different difficulty levels
# ---------------------------------------------------------------------------
dist_types <- c("normal", "skewed", "heavilySkewed", "shifted")
# For each distribution type, pick one low-pathol and one high-pathol example
examples <- do.call(rbind, lapply(dist_types, function(dt) {
  sub <- testsets[testsets$Distribution == dt, ]
  if (nrow(sub) == 0) return(NULL)
  # Pick from different biomarkers if possible
  analytes <- unique(sub$Analyte)
  picks <- NULL
  for (a in analytes) {
    asub <- sub[sub$Analyte == a, ]
    # Pick lowest and highest pathological fraction available
    asub <- asub[order(asub$fractionPathol), ]
    low <- head(asub, 1)
    high <- tail(asub, 1)
    picks <- rbind(picks, low, high)
  }
  picks
}))
# Remove duplicates
examples <- examples[!duplicated(examples$Index), ]

cat("Plotting", nrow(examples), "diagnostic cases...\n")

# ---------------------------------------------------------------------------
# Create diagnostic plots
# ---------------------------------------------------------------------------
out_dir <- file.path(working_dir, "Diagnostics")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pdf(file.path(out_dir, "diagnostic_plots.pdf"), width = 12, height = 8)

for (i in seq_len(nrow(examples))) {
  row <- examples[i, ]

  csv_file <- file.path(working_dir, "Data", row$Analyte,
                        paste0(row$Index, "_", row$Analyte, "_seed_", row$startSeed, ".csv"))
  if (!file.exists(csv_file)) {
    cat("  Skipping", basename(csv_file), "(not found)\n")
    next
  }

  data_vec <- fread(csv_file)$V1
  data_vec <- data_vec[is.finite(data_vec)]
  data_pos <- data_vec[data_vec > 0]

  if (length(data_pos) < 10) {
    cat("  Skipping", row$Analyte, row$Index, "(too few positive values)\n")
    next
  }

  # Get RINet prediction
  pred <- tryCatch({
    res <- predict_rinet_1d(data_pos, log_scale = TRUE, percentiles = percentiles)
    res[[1]]
  }, error = function(e) NULL)

  # Ground truth RI and Box-Cox params from test set definition
  gt_lower <- row$GT_LRL
  gt_upper <- row$GT_URL
  gt_lambda <- row$nonp_lambda
  gt_mu     <- row$nonp_mu
  gt_sigma  <- row$nonp_sigma
  gt_shift  <- row$nonp_shift

  # Box-Cox transform (matching RIbench's implementation)
  boxcox <- function(x, lambda) {
    if (abs(lambda) < 1e-20) log(x) else (x^lambda - 1) / lambda
  }

  # Compute z-scores for GT and predicted RIs
  zz_gt_lower  <- (boxcox(gt_lower - gt_shift, gt_lambda) - gt_mu) / gt_sigma
  zz_gt_upper  <- (boxcox(gt_upper - gt_shift, gt_lambda) - gt_mu) / gt_sigma

  # --- Plot layout: histogram + info panel ---
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

  # Left panel: Histogram with RI overlay
  title_str <- sprintf("%s [%s] | N=%s | pathol=%.0f%%",
                       row$Analyte, row$Distribution, format(row$N, big.mark = ","),
                       row$fractionPathol * 100)

  # Determine x-range from finite values only
  all_vals <- data_vec
  if (!is.null(pred) && !is.null(pred$reference_interval)) {
    all_vals <- c(all_vals, as.numeric(pred$reference_interval))
  }
  all_vals <- c(all_vals, gt_lower, gt_upper)
  all_vals <- all_vals[is.finite(all_vals)]
  xlim <- range(all_vals)

  hist(data_vec, breaks = 100, col = "grey85", border = "grey60",
       main = title_str, xlab = row$Analyte, freq = FALSE, xlim = xlim)

  # Ground truth RI
  abline(v = gt_lower, col = "forestgreen", lwd = 2.5, lty = 2)
  abline(v = gt_upper, col = "forestgreen", lwd = 2.5, lty = 2)

  if (!is.null(pred)) {
    ri <- pred$reference_interval

    # Predicted RI
    abline(v = ri["lower"], col = "dodgerblue", lwd = 2.5)
    abline(v = ri["upper"], col = "dodgerblue", lwd = 2.5)

    # Compute z-scores for predicted RIs
    zz_pred_lower <- (boxcox(ri["lower"] - gt_shift, gt_lambda) - gt_mu) / gt_sigma
    zz_pred_upper <- (boxcox(ri["upper"] - gt_shift, gt_lambda) - gt_mu) / gt_sigma
    zz_dev_lower  <- abs(zz_gt_lower - zz_pred_lower)
    zz_dev_upper  <- abs(zz_gt_upper - zz_pred_upper)
    zz_dev_ov     <- (zz_dev_lower + zz_dev_upper) / 2

    # Overlay predicted lognormal density
    pred_mean <- pred$mean   # log-scale
    pred_std  <- pred$std    # log-scale
    x_seq <- seq(max(0.001, xlim[1]), xlim[2], length.out = 500)
    y_dens <- dlnorm(x_seq, meanlog = pred_mean, sdlog = pred_std)
    # Scale by reference fraction
    y_dens <- y_dens * pred$reference_fraction
    lines(x_seq, y_dens, col = "dodgerblue", lwd = 2)

    legend("topright",
           legend = c(
             sprintf("GT RI: [%.2f, %.2f]", gt_lower, gt_upper),
             sprintf("Pred RI: [%.2f, %.2f]", ri["lower"], ri["upper"]),
             sprintf("Pred ref frac: %.1f%%", pred$reference_fraction * 100),
             sprintf("|zDev| LRL=%.3f  URL=%.3f  Ov=%.3f", zz_dev_lower, zz_dev_upper, zz_dev_ov),
             "Pred lognormal density"
           ),
           col = c("forestgreen", "dodgerblue", NA, NA, "dodgerblue"),
           lty = c(2, 1, NA, NA, 1), lwd = c(2.5, 2.5, NA, NA, 2),
           cex = 0.7, bg = "white")
  } else {
    legend("topright", legend = c(
      sprintf("GT RI: [%.2f, %.2f]", gt_lower, gt_upper),
      "Prediction FAILED"
    ), col = c("forestgreen", "red"), lty = c(2, NA), lwd = 2.5,
    cex = 0.75, bg = "white")
  }

  # Right panel: Log-scale histogram
  if (any(data_pos > 0)) {
    hist(log(data_pos), breaks = 100, col = "grey85", border = "grey60",
         main = paste(title_str, "(log scale)"),
         xlab = paste0("log(", row$Analyte, ")"), freq = FALSE)

    # Ground truth in log-scale
    abline(v = log(gt_lower), col = "forestgreen", lwd = 2.5, lty = 2)
    abline(v = log(gt_upper), col = "forestgreen", lwd = 2.5, lty = 2)

    if (!is.null(pred)) {
      # Predicted RI in log-scale
      abline(v = log(ri["lower"]), col = "dodgerblue", lwd = 2.5)
      abline(v = log(ri["upper"]), col = "dodgerblue", lwd = 2.5)

      # Overlay predicted normal density on log-scale
      x_log <- seq(min(log(data_pos)), max(log(data_pos)), length.out = 500)
      y_norm <- dnorm(x_log, mean = pred_mean, sd = pred_std)
      y_norm <- y_norm * pred$reference_fraction
      lines(x_log, y_norm, col = "dodgerblue", lwd = 2)
    }
  }
}

dev.off()
cat("Saved to:", file.path(out_dir, "diagnostic_plots.pdf"), "\n")

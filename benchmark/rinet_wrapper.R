# RINet wrapper function for RIbench benchmark evaluation
#
# Adapts rinet's predict_rinet_1d() to the interface expected by
# RIbench's evaluateBiomarkerTestSets() with requirePercentiles = TRUE.
#
# RINet assumes lognormal distributions (log_scale = TRUE) for all biomarkers.

library(rinet)

estimateRIs <- function(Data = NULL, percentiles = c(0.025, 0.975), ...) {

  # Remove non-positive values (rounding artifacts) for log-scale
  Data <- Data[Data > 0]

  result <- tryCatch(
    predict_rinet_1d(
      data        = Data,
      log_scale   = TRUE,
      percentiles = percentiles
    ),
    error = function(e) NULL
  )

  # If prediction failed, return NAs
  if (is.null(result) || length(result) == 0) {
    return(data.frame(Percentile = percentiles, PointEst = NA_real_))
  }

  ri <- result[[1]]$reference_interval
  point_estimates <- as.numeric(ri)

  data.frame(Percentile = percentiles, PointEst = point_estimates)
}

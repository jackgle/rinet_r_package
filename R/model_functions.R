#' @importFrom stats sd qnorm qchisq quantile
NULL

# Internal environment to cache loaded models
.models <- new.env(parent = emptyenv())

# Internal helper to silence TensorFlow/keras chatter (stdout/stderr/messages)
.silence_tf <- function(expr) {
  # R-level silencing
  tf_log <- file(tempfile(), open = "w")
  on.exit({
    sink(NULL, type = "message")
    sink(NULL, type = "output")
    close(tf_log)
  }, add = TRUE)
  sink(tf_log, type = "message")
  sink(tf_log, type = "output")
  
  # C-level silencing via Python if available
  if (reticulate::py_available()) {
    tryCatch({
      reticulate::py_run_string("import os; _orig_stderr = os.dup(2); _devnull = os.open(os.devnull, os.O_WRONLY); os.dup2(_devnull, 2)")
      on.exit({
        reticulate::py_run_string("os.dup2(_orig_stderr, 2); os.close(_orig_stderr); os.close(_devnull)")
      }, add = TRUE)
    }, error = function(e) {})
  }

  withCallingHandlers(expr,
    warning = function(w) invokeRestart("muffleWarning"),
    message = function(m) invokeRestart("muffleMessage")
  )
}

#' Internal function to load model and scaler
#' @keywords internal
.load_model <- function(ndim) {
  model_key <- paste0("model_", ndim, "d")
  
  if (is.null(.models[[model_key]])) {
    # Ensure TensorFlow is initialized
    if (!reticulate::py_module_available("tensorflow")) {
      stop("TensorFlow is not available. Please install TensorFlow:\n",
           "  library(keras)\n",
           "  install_keras()")
    }
    helpers <- tryCatch(get(".py_silence", envir = parent.env(environment()), inherits = TRUE), error = function(e) NULL)
    
    model_file <- paste0("rinet_", ndim, "d.keras")
    model_path <- system.file("models", model_file, package = "rinet")
    
    if (!file.exists(model_path)) {
      stop("Model file not found: ", model_path, 
           "\nMake sure ", model_file, " is in inst/models/")
    }
    
    message("Loading ", ndim, "D model (this may take a moment on first use)...")
    keras <- reticulate::import("keras")
    .models[[model_key]] <- keras$models$load_model(model_path)
    message("Model loaded.")
  }
  
  return(.models[[model_key]])
}

#' Internal function to load scaler
#' @keywords internal
.load_scaler <- function(ndim) {
  scaler_key <- paste0("scaler_", ndim, "d")
  
  if (is.null(.models[[scaler_key]])) {
    scaler_file <- paste0("scaler_", ndim, "d.pkl")
    scaler_path <- system.file("models", scaler_file, package = "rinet")
    
    if (!file.exists(scaler_path)) {
      stop("Scaler file not found: ", scaler_path,
           "\nMake sure ", scaler_file, " is in inst/models/")
    }
    
    .silence_tf({
      reticulate::py_run_string("import pickle")
      .models[[scaler_key]] <- reticulate::py_run_string(
        sprintf("scaler = pickle.load(open('%s', 'rb'))", scaler_path)
      )$scaler
    })
  }
  
  return(.models[[scaler_key]])
}

#' Extract histogram features from standardized data
#' @keywords internal
.extract_features <- function(data, ndim, feature_grid_range = c(-4, 4), 
                              feature_grid_nbins = 100) {
  data <- as.matrix(data)
  
  if (ndim == 1) {
    # Check for standardization (using population SD to match Python)
    n <- length(data)
    pop_sd <- sd(data) * sqrt((n - 1) / n)
    if (abs(mean(data)) > 0.01 || abs(pop_sd - 1) > 0.01) {
      stop("Data must be standardized (mean=0, sd=1)")
    }
    
    # Create histogram using numpy for consistency with Python/2D
    breaks <- seq(feature_grid_range[1], feature_grid_range[2], 
                  length.out = feature_grid_nbins + 1)
    np <- reticulate::import("numpy", convert = FALSE)
    hist_result <- np$histogram(data, bins = breaks, density = TRUE)
    # Use [[0]] to get the first element (counts) from the Python tuple
    features <- as.numeric(reticulate::py_to_r(hist_result[[0]]))
    
  } else if (ndim == 2) {
    # Check for standardization (using population SD to match Python)
    n <- nrow(data)
    pop_sd <- apply(data, 2, sd) * sqrt((n - 1) / n)
    if (any(abs(colMeans(data)) > 0.01) || any(abs(pop_sd - 1) > 0.01)) {
      stop("Data must be standardized (mean=0, sd=1 for each dimension)")
    }
    
    # Create 2D histogram
    breaks <- seq(feature_grid_range[1], feature_grid_range[2], 
                  length.out = feature_grid_nbins + 1)
    
    # Use reticulate to call numpy's histogram2d for consistency
    np <- reticulate::import("numpy", convert = FALSE)
    hist_result <- np$histogram2d(
      data[, 1], data[, 2],
      bins = breaks,
      density = TRUE
    )
    # Use [[0]] to get the first element (H) from the Python tuple
    features <- as.matrix(reticulate::py_to_r(hist_result[[0]]))
  }
  
  # Normalize to [0, 1]
  features <- (features - min(features)) / (max(features) - min(features))
  
  return(features)
}

#' Convert correlation to covariance matrix
#' @keywords internal
.correlation_to_covariance <- function(corr_matrix, std_vector) {
  std_matrix <- outer(std_vector, std_vector)
  cov_matrix <- corr_matrix * std_matrix
  return(cov_matrix)
}

#' Predict statistics of the underlying reference distribution from mixture distributions using RINet
#'
#' Automatically detects whether input data is 1D or 2D and calls the
#' appropriate prediction function. This is the main user-facing function.
#' It estimates the statistics of a "healthy" reference population from a
#' mixture of healthy and pathological measurements.
#'
#' @param data A numeric vector, matrix, or list. For 1D: vector or matrix with
#'   1 column. For 2D: matrix with 2 columns. Can also be a list of such objects.
#' @param feature_grid_range Numeric vector of length 2 specifying the range
#'   for histogram binning. Default is c(-4, 4).
#' @param feature_grid_nbins Integer specifying the number of histogram bins.
#'   Default is 100.
#' @param verbose Integer controlling verbosity (0 = silent). Default is 0.
#' @param log_scale Logical indicating whether to log-transform the data before
#'   prediction. If TRUE (default), returns log-scale statistics and calculates
#'   reference intervals in the original scale. Default is TRUE.
#' @param percentiles Numeric vector of length 2 specifying the lower and upper
#'   percentiles for the reference interval. Default is c(0.025, 0.975).
#' @param n_bootstrap Integer specifying the number of bootstrap resamples for
#'   confidence intervals. Default is 0 (no bootstrap). When > 0, confidence
#'   intervals are computed for all predicted statistics using batch inference.
#' @param confidence_level Numeric specifying the confidence level for bootstrap
#'   intervals. Default is 0.95.
#' @return A list of predictions. Each element contains:
#'   \item{mean}{Predicted mean(s) (log-scale if log_scale=TRUE)}
#'   \item{std}{Predicted standard deviation(s) (log-scale if log_scale=TRUE)}
#'   \item{covariance}{Predicted covariance matrix}
#'   \item{correlation}{Predicted correlation (NA for 1D)}
#'   \item{reference_fraction}{Predicted reference component fraction}
#'   \item{reference_interval}{Reference interval in original scale (if log_scale=TRUE)}
#'   \item{log_scale}{Logical indicating whether log-scaling was used}
#'   \item{bootstrap_ci}{List of bootstrap confidence intervals (if n_bootstrap > 0)}
#' @export
#' @examples
#' \dontrun{
#'   # 1D sample (using positive data for log-scale)
#'   sample_1d <- exp(rnorm(1000, mean = 2, sd = 0.5))
#'   result <- predict_rinet(sample_1d)
#' 
#'   # 2D sample (using positive data for log-scale)
#'   sample_2d <- exp(matrix(rnorm(2000, mean = 2, sd = 0.5), ncol = 2))
#'   result <- predict_rinet(sample_2d)
#' 
#'   # Multiple samples (automatically detected)
#'   samples <- list(exp(rnorm(1000, mean = 2, sd = 0.5)),
#'                   exp(rnorm(1000, mean = 2, sd = 0.5)))
#'   results <- predict_rinet(samples)
#' }
predict_rinet <- function(data, feature_grid_range = c(-4, 4),
                          feature_grid_nbins = 100, verbose = 0,
                          log_scale = TRUE, percentiles = c(0.025, 0.975),
                          n_bootstrap = 0, confidence_level = 0.95) {
  # Determine dimensionality from first sample
  if (is.list(data) && !is.data.frame(data)) {
    sample <- data[[1]]
  } else {
    sample <- data
  }
  
  # Convert to matrix if needed to check dimensions
  if (is.vector(sample) && !is.list(sample)) {
    ndim <- 1
  } else {
    sample_mat <- as.matrix(sample)
    ncols <- ncol(sample_mat)
    
    if (ncols == 1) {
      ndim <- 1
    } else if (ncols == 2) {
      ndim <- 2
    } else {
      stop("Data must be 1D (vector or single column) or 2D (2 columns). Found ", 
           ncols, " columns.")
    }
  }
  
  # Call appropriate function
  if (ndim == 1) {
    return(predict_rinet_1d(data, feature_grid_range, feature_grid_nbins, verbose,
                           log_scale, percentiles, n_bootstrap, confidence_level))
  } else {
    return(predict_rinet_2d(data, feature_grid_range, feature_grid_nbins, verbose,
                           log_scale, percentiles, n_bootstrap, confidence_level))
  }
}

#' Predict statistics of the underlying reference distribution from 1D mixture distributions using RINet
#'
#' Takes one or more 1D samples and predicts the underlying reference
#' population statistics (mean, std, reference fraction) from a mixture
#' of healthy and pathological measurements.
#'
#' @param data A numeric vector, matrix, or list of vectors. Each sample should
#'   contain observations from a 1D mixture distribution.
#' @param feature_grid_range Numeric vector of length 2 specifying the range
#'   for histogram binning. Default is c(-4, 4).
#' @param feature_grid_nbins Integer specifying the number of histogram bins.
#'   Default is 100.
#' @param verbose Integer controlling verbosity (0 = silent). Default is 0.
#' @param log_scale Logical indicating whether to log-transform the data before
#'   prediction. If TRUE (default), returns log-scale statistics and calculates
#'   reference intervals in the original scale. Default is TRUE.
#' @param percentiles Numeric vector of length 2 specifying the lower and upper
#'   percentiles for the reference interval. Default is c(0.025, 0.975).
#' @param n_bootstrap Integer specifying the number of bootstrap resamples for
#'   confidence intervals. Default is 0 (no bootstrap). When > 0, confidence
#'   intervals are computed for all predicted statistics.
#' @param confidence_level Numeric specifying the confidence level for bootstrap
#'   intervals. Default is 0.95.
#' @return A list of predictions. Each element contains:
#'   \item{mean}{Predicted mean (scalar, log-scale if log_scale=TRUE)}
#'   \item{std}{Predicted standard deviation (scalar, log-scale if log_scale=TRUE)}
#'   \item{covariance}{Covariance matrix (1x1 matrix)}
#'   \item{correlation}{Always NA for 1D}
#'   \item{reference_fraction}{Predicted reference component fraction}
#'   \item{reference_interval}{Reference interval in original scale (if log_scale=TRUE)}
#'   \item{log_scale}{Logical indicating whether log-scaling was used}
#'   \item{bootstrap_ci}{List of bootstrap confidence intervals (if n_bootstrap > 0):
#'     mean_ci, std_ci, reference_fraction_ci, reference_interval_lower_ci,
#'     reference_interval_upper_ci}
#' @export
#' @examples
#' \dontrun{
#'   # Single sample (using positive data for log-scale)
#'   sample1 <- exp(rnorm(1000, mean = 2, sd = 0.3))
#'   result <- predict_rinet_1d(sample1)
#'   print(result[[1]]$mean)
#' 
#'   # Multiple samples
#'   samples <- list(exp(rnorm(1000, 2, 0.3)), exp(rnorm(1000, 1.5, 0.4)))
#'   results <- predict_rinet_1d(samples)
#' }
predict_rinet_1d <- function(data, feature_grid_range = c(-4, 4),
                             feature_grid_nbins = 100, verbose = 0,
                             log_scale = TRUE, percentiles = c(0.025, 0.975),
                             n_bootstrap = 0, confidence_level = 0.95) {
  # Convert to list if single sample
  if (!is.list(data) || is.data.frame(data)) {
    data <- list(as.vector(as.matrix(data)))
  }

  # Apply log transformation if requested
  original_data <- data
  if (log_scale) {
    data <- lapply(data, function(x) {
      if (any(x <= 0)) {
        stop("Log-scaling requires all values to be positive. Found zero or negative values.")
      }
      log(x)
    })
  }

  # Load model and scaler
  model <- .load_model(1)
  scaler <- .load_scaler(1)
  np <- reticulate::import("numpy", convert = FALSE)

  # Helper function to compute population SD
  pop_sd <- function(x) {
    n <- length(x)
    sd(x) * sqrt((n - 1) / n)
  }

  # Helper function to process a single sample and return features + standardization params
  process_sample <- function(x) {
    m <- mean(x)
    s <- pop_sd(x)
    if (is.na(s) || s < 1e-10) {
      stop("Sample has zero or near-zero variance.\n",
           "RINet requires variability in the data for reference interval estimation.\n",
           "Constant or near-constant values are not valid inputs.")
    }
    x_std <- (x - m) / s
    feat <- .extract_features(x_std, ndim = 1,
                              feature_grid_range = feature_grid_range,
                              feature_grid_nbins = feature_grid_nbins)
    list(features = feat, mean = m, std = s)
  }

  # Process original samples
  processed <- lapply(data, process_sample)
  means <- lapply(processed, `[[`, "mean")
  stds <- lapply(processed, `[[`, "std")
  features <- lapply(processed, `[[`, "features")

  # If bootstrap requested, generate resampled features
  n_samples <- length(data)
  bootstrap_indices <- NULL  # Track which bootstrap belongs to which original sample

  if (n_bootstrap > 0) {
    bootstrap_features <- list()
    bootstrap_means <- list()
    bootstrap_stds <- list()
    bootstrap_indices <- integer(0)

    for (i in seq_len(n_samples)) {
      x <- data[[i]]
      n <- length(x)
      for (b in seq_len(n_bootstrap)) {
        # Resample with replacement
        x_boot <- sample(x, n, replace = TRUE)
        proc <- process_sample(x_boot)
        bootstrap_features <- c(bootstrap_features, list(proc$features))
        bootstrap_means <- c(bootstrap_means, list(proc$mean))
        bootstrap_stds <- c(bootstrap_stds, list(proc$std))
        bootstrap_indices <- c(bootstrap_indices, i)
      }
    }

    # Combine original and bootstrap features for batch prediction
    all_features <- c(features, bootstrap_features)
    all_means <- c(means, bootstrap_means)
    all_stds <- c(stds, bootstrap_stds)
  } else {
    all_features <- features
    all_means <- means
    all_stds <- stds
  }

  # Stack and add channel dimension for batch prediction
  features_array <- np$stack(all_features, axis = 0L)
  features_array <- np$expand_dims(features_array, axis = -1L)

  # Predict (single batch call for all samples + bootstraps)
  predictions <- NULL
  helpers <- tryCatch(get(".py_silence", envir = parent.env(environment()), inherits = TRUE), error = function(e) NULL)
  if (!is.null(helpers)) {
    predictions <- tryCatch(
      helpers$predict_silent(model, features_array, verbose = verbose),
      error = function(e) NULL
    )
  }
  if (is.null(predictions)) {
    predictions <- .silence_tf(model$predict(features_array, verbose = verbose))
  }

  # Inverse transform using scaler
  predictions <- reticulate::py_to_r(scaler$inverse_transform(predictions))

  # Helper to compute scaled results from raw predictions
  compute_result <- function(pred_row, m, s) {
    p_mean <- pred_row[1]
    p_std <- pred_row[2]
    scaled_mean <- (p_mean * s) + m
    scaled_std <- p_std * s
    ref_frac <- if (length(pred_row) > 2) pred_row[3] else NA

    ref_interval <- NULL
    if (log_scale) {
      lower <- exp(scaled_mean + qnorm(percentiles[1]) * scaled_std)
      upper <- exp(scaled_mean + qnorm(percentiles[2]) * scaled_std)
      ref_interval <- c(lower = lower, upper = upper)
    }

    list(mean = scaled_mean, std = scaled_std,
         reference_fraction = ref_frac, reference_interval = ref_interval)
  }

  # Process results
  results <- list()
  ci_alpha <- (1 - confidence_level) / 2
  ci_probs <- c(ci_alpha, 1 - ci_alpha)

  for (i in seq_len(n_samples)) {
    # Original sample result
    res <- compute_result(predictions[i, ], all_means[[i]], all_stds[[i]])

    cov <- .correlation_to_covariance(matrix(1, 1, 1), res$std)

    result_list <- list(
      mean = res$mean,
      std = res$std,
      covariance = cov,
      correlation = NA,
      reference_fraction = res$reference_fraction,
      reference_interval = res$reference_interval,
      log_scale = log_scale
    )

    # Add bootstrap CIs if requested
    if (n_bootstrap > 0) {
      # Find bootstrap predictions for this sample
      boot_idx <- which(bootstrap_indices == i) + n_samples
      boot_results <- lapply(boot_idx, function(j) {
        compute_result(predictions[j, ], all_means[[j]], all_stds[[j]])
      })

      boot_means <- sapply(boot_results, `[[`, "mean")
      boot_stds <- sapply(boot_results, `[[`, "std")
      boot_ref_fracs <- sapply(boot_results, `[[`, "reference_fraction")

      ci_list <- list(
        mean_ci = quantile(boot_means, probs = ci_probs, na.rm = TRUE),
        std_ci = quantile(boot_stds, probs = ci_probs, na.rm = TRUE),
        reference_fraction_ci = quantile(boot_ref_fracs, probs = ci_probs, na.rm = TRUE)
      )

      if (log_scale) {
        boot_ri <- do.call(rbind, lapply(boot_results, `[[`, "reference_interval"))
        ci_list$reference_interval_lower_ci <- quantile(boot_ri[, "lower"], probs = ci_probs, na.rm = TRUE)
        ci_list$reference_interval_upper_ci <- quantile(boot_ri[, "upper"], probs = ci_probs, na.rm = TRUE)
      }

      result_list$bootstrap_ci <- ci_list
      result_list$n_bootstrap <- n_bootstrap
    }

    results[[i]] <- result_list
  }

  return(results)
}

#' Predict statistics of the underlying reference distribution from 2D mixture distributions using RINet
#'
#' Takes one or more 2D samples and predicts the underlying reference
#' population statistics (means, stds, correlation, reference fraction)
#' from a mixture of healthy and pathological measurements.
#'
#' @param data A matrix or list of matrices. Each sample should be a matrix
#'   with 2 columns representing observations from a 2D mixture distribution.
#' @param feature_grid_range Numeric vector of length 2 specifying the range
#'   for histogram binning. Default is c(-4, 4).
#' @param feature_grid_nbins Integer specifying the number of histogram bins.
#'   Default is 100.
#' @param verbose Integer controlling verbosity (0 = silent). Default is 0.
#' @param log_scale Logical indicating whether to log-transform the data before
#'   prediction. If TRUE (default), returns log-scale statistics and calculates
#'   reference intervals in the original scale. Default is TRUE.
#' @param percentiles Numeric vector of length 2 specifying the lower and upper
#'   percentiles for the reference interval. Default is c(0.025, 0.975).
#' @param n_bootstrap Integer specifying the number of bootstrap resamples for
#'   confidence intervals. Default is 0 (no bootstrap). When > 0, confidence
#'   intervals are computed for all predicted statistics.
#' @param confidence_level Numeric specifying the confidence level for bootstrap
#'   intervals. Default is 0.95.
#' @return A list of predictions. Each element contains:
#'   \item{mean}{Predicted means (vector of length 2, log-scale if log_scale=TRUE)}
#'   \item{std}{Predicted standard deviations (vector of length 2, log-scale if log_scale=TRUE)}
#'   \item{covariance}{Predicted covariance matrix (2x2 matrix)}
#'   \item{correlation}{Predicted correlation coefficient (scalar)}
#'   \item{reference_fraction}{Predicted reference component fraction}
#'   \item{reference_interval}{Reference region ellipse vertices (100x2 matrix) in original scale (if log_scale=TRUE)}
#'   \item{log_scale}{Logical indicating whether log-scaling was used}
#'   \item{bootstrap_ci}{List of bootstrap confidence intervals (if n_bootstrap > 0):
#'     mean_ci (2x2 matrix), std_ci (2x2 matrix), correlation_ci, reference_fraction_ci}
#' @export
#' @examples
#' \dontrun{
#'   # Single 2D sample (using positive data for log-scale)
#'   sample1 <- exp(matrix(rnorm(2000, mean = 2, sd = 0.3), ncol = 2))
#'   result <- predict_rinet_2d(sample1)
#'   print(result[[1]]$mean)
#'   print(result[[1]]$covariance)
#' 
#'   # Multiple samples
#'   samples <- list(exp(matrix(rnorm(2000, mean = 2, sd = 0.3), ncol = 2)),
#'                   exp(matrix(rnorm(2000, mean = 2, sd = 0.3), ncol = 2)))
#'   results <- predict_rinet_2d(samples)
#' }
predict_rinet_2d <- function(data, feature_grid_range = c(-4, 4),
                             feature_grid_nbins = 100, verbose = 0,
                             log_scale = TRUE, percentiles = c(0.025, 0.975),
                             n_bootstrap = 0, confidence_level = 0.95) {
  # Convert to list if single sample
  if (!is.list(data) || is.data.frame(data)) {
    data <- list(as.matrix(data))
  }

  # Apply log transformation if requested
  original_data <- data
  if (log_scale) {
    data <- lapply(data, function(x) {
      if (any(x <= 0)) {
        stop("Log-scaling requires all values to be positive. Found zero or negative values.")
      }
      log(x)
    })
  }

  # Ensure all are matrices
  data <- lapply(data, as.matrix)

  # Check dimensions
  if (any(sapply(data, ncol) != 2)) {
    stop("All samples must have 2 columns for 2D prediction")
  }

  # Load model and scaler
  model <- .load_model(2)
  scaler <- .load_scaler(2)
  np <- reticulate::import("numpy", convert = FALSE)

  # Helper function to compute population SD
 pop_sd_vec <- function(x) {
    n <- nrow(x)
    apply(x, 2, sd) * sqrt((n - 1) / n)
  }

  # Helper function to process a single sample and return features + standardization params
  process_sample <- function(x) {
    m <- colMeans(x)
    s <- pop_sd_vec(x)
    zero_var <- which(is.na(s) | s < 1e-10)
    if (length(zero_var) > 0) {
      dims <- paste(zero_var, collapse = ", ")
      stop("Sample has zero or near-zero variance in dimension(s): ", dims, ".\n",
           "RINet requires variability in the data for reference interval estimation.\n",
           "Constant or near-constant values are not valid inputs.")
    }
    x_std <- sweep(sweep(x, 2, m, "-"), 2, s, "/")
    feat <- .extract_features(x_std, ndim = 2,
                              feature_grid_range = feature_grid_range,
                              feature_grid_nbins = feature_grid_nbins)
    list(features = feat, mean = m, std = s)
  }

  # Process original samples
  processed <- lapply(data, process_sample)
  means <- lapply(processed, `[[`, "mean")
  stds <- lapply(processed, `[[`, "std")
  features <- lapply(processed, `[[`, "features")

  # If bootstrap requested, generate resampled features
  n_samples <- length(data)
  bootstrap_indices <- NULL

  if (n_bootstrap > 0) {
    bootstrap_features <- list()
    bootstrap_means <- list()
    bootstrap_stds <- list()
    bootstrap_indices <- integer(0)

    for (i in seq_len(n_samples)) {
      x <- data[[i]]
      n <- nrow(x)
      for (b in seq_len(n_bootstrap)) {
        # Resample rows with replacement
        idx <- sample(n, n, replace = TRUE)
        x_boot <- x[idx, , drop = FALSE]
        proc <- process_sample(x_boot)
        bootstrap_features <- c(bootstrap_features, list(proc$features))
        bootstrap_means <- c(bootstrap_means, list(proc$mean))
        bootstrap_stds <- c(bootstrap_stds, list(proc$std))
        bootstrap_indices <- c(bootstrap_indices, i)
      }
    }

    # Combine original and bootstrap features for batch prediction
    all_features <- c(features, bootstrap_features)
    all_means <- c(means, bootstrap_means)
    all_stds <- c(stds, bootstrap_stds)
  } else {
    all_features <- features
    all_means <- means
    all_stds <- stds
  }

  # Stack and add channel dimension for batch prediction
  features_array <- np$stack(all_features, axis = 0L)
  features_array <- np$expand_dims(features_array, axis = -1L)

  # Predict (single batch call for all samples + bootstraps)
  predictions <- NULL
  helpers <- tryCatch(get(".py_silence", envir = parent.env(environment()), inherits = TRUE), error = function(e) NULL)
  if (!is.null(helpers)) {
    predictions <- tryCatch(
      helpers$predict_silent(model, features_array, verbose = verbose),
      error = function(e) NULL
    )
  }
  if (is.null(predictions)) {
    predictions <- .silence_tf(model$predict(features_array, verbose = verbose))
  }

  # Inverse transform using scaler
  predictions <- reticulate::py_to_r(scaler$inverse_transform(predictions))

  # Helper to compute scaled results from raw predictions
  compute_result <- function(pred_row, m, s) {
    p_mean <- pred_row[1:2]
    p_std <- pred_row[3:4]
    p_cor <- pred_row[5]

    scaled_mean <- (p_mean * s) + m
    scaled_std <- p_std * s
    ref_frac <- if (length(pred_row) > 5) pred_row[6] else NA

    corr_matrix <- matrix(c(1, p_cor, p_cor, 1), 2, 2)
    cov <- .correlation_to_covariance(corr_matrix, scaled_std)

    ref_interval <- NULL
    if (log_scale) {
      prob <- percentiles[2] - percentiles[1]
      chi2_val <- qchisq(prob, df = 2)
      eigen_decomp <- eigen(cov)
      eigenvalues <- eigen_decomp$values
      eigenvectors <- eigen_decomp$vectors
      theta <- seq(0, 2 * pi, length.out = 100)
      ellipse_log <- matrix(NA, nrow = 100, ncol = 2)
      for (j in 1:100) {
        point <- sqrt(chi2_val) * c(sqrt(eigenvalues[1]) * cos(theta[j]),
                                     sqrt(eigenvalues[2]) * sin(theta[j]))
        ellipse_log[j, ] <- eigenvectors %*% point + scaled_mean
      }
      ref_interval <- exp(ellipse_log)
      colnames(ref_interval) <- c("dim1", "dim2")
    }

    list(mean = scaled_mean, std = scaled_std, correlation = p_cor,
         covariance = cov, reference_fraction = ref_frac,
         reference_interval = ref_interval)
  }

  # Process results
  results <- list()
  ci_alpha <- (1 - confidence_level) / 2
  ci_probs <- c(ci_alpha, 1 - ci_alpha)

  for (i in seq_len(n_samples)) {
    # Original sample result
    res <- compute_result(predictions[i, ], all_means[[i]], all_stds[[i]])

    result_list <- list(
      mean = res$mean,
      std = res$std,
      covariance = res$covariance,
      correlation = res$correlation,
      reference_fraction = res$reference_fraction,
      reference_interval = res$reference_interval,
      log_scale = log_scale
    )

    # Add bootstrap CIs if requested
    if (n_bootstrap > 0) {
      boot_idx <- which(bootstrap_indices == i) + n_samples
      boot_results <- lapply(boot_idx, function(j) {
        compute_result(predictions[j, ], all_means[[j]], all_stds[[j]])
      })

      # Extract vectors for each dimension
      boot_mean1 <- sapply(boot_results, function(r) r$mean[1])
      boot_mean2 <- sapply(boot_results, function(r) r$mean[2])
      boot_std1 <- sapply(boot_results, function(r) r$std[1])
      boot_std2 <- sapply(boot_results, function(r) r$std[2])
      boot_cors <- sapply(boot_results, `[[`, "correlation")
      boot_ref_fracs <- sapply(boot_results, `[[`, "reference_fraction")

      ci_list <- list(
        mean_ci = rbind(
          quantile(boot_mean1, probs = ci_probs, na.rm = TRUE),
          quantile(boot_mean2, probs = ci_probs, na.rm = TRUE)
        ),
        std_ci = rbind(
          quantile(boot_std1, probs = ci_probs, na.rm = TRUE),
          quantile(boot_std2, probs = ci_probs, na.rm = TRUE)
        ),
        correlation_ci = quantile(boot_cors, probs = ci_probs, na.rm = TRUE),
        reference_fraction_ci = quantile(boot_ref_fracs, probs = ci_probs, na.rm = TRUE)
      )

      result_list$bootstrap_ci <- ci_list
      result_list$n_bootstrap <- n_bootstrap
    }

    results[[i]] <- result_list
  }

  return(results)
}

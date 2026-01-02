#' @importFrom stats sd
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
    
    message("Loading ", ndim, "D model...")
    if (!is.null(helpers) && !is.null(helpers$load_model_silent)) {
      .models[[model_key]] <- tryCatch(
        helpers$load_model_silent(model_path),
        error = function(e) NULL
      )
    }
    if (is.null(.models[[model_key]])) {
      .models[[model_key]] <- .silence_tf(keras::load_model_tf(model_path))
    }
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
#' @return A list of predictions. Each element contains:
#'   \item{mean}{Predicted mean(s)}
#'   \item{std}{Predicted standard deviation(s)}
#'   \item{covariance}{Predicted covariance matrix}
#'   \item{correlation}{Predicted correlation (NA for 1D)}
#'   \item{reference_fraction}{Predicted reference component fraction}
#' @export
#' @examples
#' \dontrun{
#'   # 1D sample
#'   sample_1d <- rnorm(1000)
#'   result <- predict_rinet(sample_1d)
#'   
#'   # 2D sample
#'   sample_2d <- matrix(rnorm(2000), ncol = 2)
#'   result <- predict_rinet(sample_2d)
#'   
#'   # Multiple samples (automatically detected)
#'   samples <- list(rnorm(1000), rnorm(1000))
#'   results <- predict_rinet(samples)
#' }
predict_rinet <- function(data, feature_grid_range = c(-4, 4),
                          feature_grid_nbins = 100, verbose = 0) {
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
    return(predict_rinet_1d(data, feature_grid_range, feature_grid_nbins, verbose))
  } else {
    return(predict_rinet_2d(data, feature_grid_range, feature_grid_nbins, verbose))
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
#' @return A list of predictions. Each element contains:
#'   \item{mean}{Predicted mean (scalar)}
#'   \item{std}{Predicted standard deviation (scalar)}
#'   \item{covariance}{Covariance matrix (1x1 matrix)}
#'   \item{correlation}{Always NA for 1D}
#'   \item{reference_fraction}{Predicted reference component fraction}
#' @export
#' @examples
#' \dontrun{
#'   # Single sample
#'   sample1 <- rnorm(1000, mean = 5, sd = 2)
#'   result <- predict_rinet_1d(sample1)
#'   print(result[[1]]$mean)
#'   
#'   # Multiple samples
#'   samples <- list(rnorm(1000, 3, 1.5), rnorm(1000, -2, 0.8))
#'   results <- predict_rinet_1d(samples)
#' }
predict_rinet_1d <- function(data, feature_grid_range = c(-4, 4), 
                             feature_grid_nbins = 100, verbose = 0) {
  # Convert to list if single sample
  if (!is.list(data) || is.data.frame(data)) {
    data <- list(as.vector(as.matrix(data)))
  }
  
  # Load model and scaler
  model <- .load_model(1)
  scaler <- .load_scaler(1)
  
  # Store original means and sds (using population SD to match Python/numpy ddof=0)
  means <- lapply(data, mean)
  stds <- lapply(data, function(x) {
    n <- length(x)
    sd(x) * sqrt((n - 1) / n)
  })
  
  # Check for zero or near-zero variance
  for (i in seq_along(stds)) {
    if (is.na(stds[[i]]) || stds[[i]] < 1e-10) {
      stop("Sample ", i, " has zero or near-zero variance.\n",
           "RINet requires variability in the data for reference interval estimation.\n",
           "Constant or near-constant values are not valid inputs.")
    }
  }
  
  # Standardize data
  data_std <- mapply(function(x, m, s) (x - m) / s, 
                     data, means, stds, SIMPLIFY = FALSE)
  
  # Extract features
  features <- lapply(data_std, .extract_features, ndim = 1,
                    feature_grid_range = feature_grid_range,
                    feature_grid_nbins = feature_grid_nbins)
  
  # Stack and add channel dimension (using np_array to ensure correct stacking order)
  np <- reticulate::import("numpy", convert = FALSE)
  features_array <- np$stack(features, axis = 0L)
  features_array <- np$expand_dims(features_array, axis = -1L)
  
  # Predict
  predictions <- NULL
  helpers <- tryCatch(get(".py_silence", envir = parent.env(environment()), inherits = TRUE), error = function(e) NULL)
  if (!is.null(helpers) && !is.null(helpers$predict_silent)) {
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
  
  # Process results
  results <- list()
  for (i in seq_len(nrow(predictions))) {
    p_mean <- predictions[i, 1]
    p_std <- predictions[i, 2]
    
    # Scale back to original
    scaled_mean <- (p_mean * stds[[i]]) + means[[i]]
    scaled_std <- p_std * stds[[i]]
    
    cov <- .correlation_to_covariance(matrix(1, 1, 1), scaled_std)
    
    results[[i]] <- list(
      mean = scaled_mean,
      std = scaled_std,
      covariance = cov,
      correlation = NA,
      reference_fraction = if (ncol(predictions) > 2) predictions[i, 3] else NA
    )
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
#' @return A list of predictions. Each element contains:
#'   \item{mean}{Predicted means (vector of length 2)}
#'   \item{std}{Predicted standard deviations (vector of length 2)}
#'   \item{covariance}{Predicted covariance matrix (2x2 matrix)}
#'   \item{correlation}{Predicted correlation coefficient (scalar)}
#'   \item{reference_fraction}{Predicted reference component fraction}
#' @export
#' @examples
#' \dontrun{
#'   # Single 2D sample
#'   sample1 <- matrix(rnorm(2000), ncol = 2)
#'   result <- predict_rinet_2d(sample1)
#'   print(result[[1]]$mean)
#'   print(result[[1]]$covariance)
#'   
#'   # Multiple samples
#'   samples <- list(matrix(rnorm(2000), ncol = 2), 
#'                   matrix(rnorm(2000), ncol = 2))
#'   results <- predict_rinet_2d(samples)
#' }
predict_rinet_2d <- function(data, feature_grid_range = c(-4, 4),
                             feature_grid_nbins = 100, verbose = 0) {
  # Convert to list if single sample
  if (!is.list(data) || is.data.frame(data)) {
    data <- list(as.matrix(data))
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
  
  # Store original means and sds (using population SD to match Python/numpy ddof=0)
  means <- lapply(data, colMeans)
  stds <- lapply(data, function(x) {
    n <- nrow(x)
    apply(x, 2, sd) * sqrt((n - 1) / n)
  })
  
  # Check for zero or near-zero variance
  for (i in seq_along(stds)) {
    zero_var <- which(is.na(stds[[i]]) | stds[[i]] < 1e-10)
    if (length(zero_var) > 0) {
      dims <- paste(zero_var, collapse = ", ")
      stop("Sample ", i, " has zero or near-zero variance in dimension(s): ", dims, ".\n",
           "RINet requires variability in the data for reference interval estimation.\n",
           "Constant or near-constant values are not valid inputs.")
    }
  }
  
  # Standardize data
  data_std <- mapply(function(x, m, s) {
    sweep(sweep(x, 2, m, "-"), 2, s, "/")
  }, data, means, stds, SIMPLIFY = FALSE)
  
  # Extract features
  features <- lapply(data_std, .extract_features, ndim = 2,
                    feature_grid_range = feature_grid_range,
                    feature_grid_nbins = feature_grid_nbins)
  
  # Stack and add channel dimension (using np_array to ensure correct stacking order)
  np <- reticulate::import("numpy", convert = FALSE)
  features_array <- np$stack(features, axis = 0L)
  features_array <- np$expand_dims(features_array, axis = -1L)
  
  # Predict
  predictions <- NULL
  helpers <- tryCatch(get(".py_silence", envir = parent.env(environment()), inherits = TRUE), error = function(e) NULL)
  if (!is.null(helpers) && !is.null(helpers$predict_silent)) {
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
  
  # Process results
  results <- list()
  for (i in seq_len(nrow(predictions))) {
    p_mean <- predictions[i, 1:2]
    p_std <- predictions[i, 3:4]
    p_cor <- predictions[i, 5]
    
    # Scale back to original
    scaled_mean <- (p_mean * stds[[i]]) + means[[i]]
    scaled_std <- p_std * stds[[i]]
    
    corr_matrix <- matrix(c(1, p_cor, p_cor, 1), 2, 2)
    cov <- .correlation_to_covariance(corr_matrix, scaled_std)
    
    results[[i]] <- list(
      mean = scaled_mean,
      std = scaled_std,
      covariance = cov,
      correlation = p_cor,
      reference_fraction = if (ncol(predictions) > 5) predictions[i, 6] else NA
    )
  }
  
  return(results)
}

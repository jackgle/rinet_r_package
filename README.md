# rinet

**RINet**: An R package for indirect estimation of clinical reference distributions using a neural network.

## What it does

RINet is designed to estimate the statistics of an underlying reference distribution (healthy population) from a mixture distribution of raw clinical measurements that include both healthy and pathological patients.

Given samples from 1D or 2D mixture distributions, RINet predicts:
- Mean(s) of the reference component
- Standard deviation(s) of the reference component
- Correlation (2D only)
- Covariance matrix
- Reference component fraction (the proportion of "healthy" samples in the mixture)

The package automatically detects whether your data is 1D or 2D and uses the appropriate model.

## Installation

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install rinet from GitHub
devtools::install_github("jackgle/rinet_r_package")
```

## Prerequisites

This package requires the `keras` R package with TensorFlow:

```r
# Install keras R package
install.packages(c("keras", "reticulate"))

# Install TensorFlow backend
library(keras)
install_keras()
```

## Usage

### 1D Example: Predicting Albumin (ALB) Reference Interval

```r
library(rinet)
library(reflimR)

# Load liver test data (mixture of reference and pathological samples)
data(livertests)

# Filter for females and extract ALB values
alb_mixture <- livertests[livertests$Sex == "f", "ALB"]

# Method 1: RINet prediction from mixture
rinet_result <- predict_rinet(alb_mixture)

cat("RINet Prediction (from mixture):\n")
cat(sprintf("  Reference Interval: [%.2f, %.2f]\n", 
    rinet_result[[1]]$reference_interval["lower"],
    rinet_result[[1]]$reference_interval["upper"]))
cat(sprintf("  Log-Mean: %.2f\n", rinet_result[[1]]$mean))
cat(sprintf("  Log-SD:   %.2f\n", rinet_result[[1]]$std))
cat(sprintf("  Reference Fraction: %.2f\n", rinet_result[[1]]$reference_fraction))

# Method 2: Direct calculation (requires knowing which samples are reference)
reference_only <- livertests[livertests$Sex == "f" & 
                             livertests$Category == "reference", "ALB"]

cat("\nDirect Method (reference samples only):\n")
log_ref <- log(reference_only)
log_mean <- mean(log_ref)
log_sd <- sd(log_ref)
cat(sprintf("  Reference Interval: [%.2f, %.2f]\n",
    exp(log_mean + qnorm(0.025) * log_sd),
    exp(log_mean + qnorm(0.975) * log_sd)))
cat(sprintf("  Log-Mean: %.2f\n", log_mean))
cat(sprintf("  Log-SD:   %.2f\n", log_sd))
cat(sprintf("  Reference Fraction: %.2f\n", 
    length(reference_only) / length(alb_mixture)))
```

### 2D Example: Predicting ALB and PROT Joint Distribution

```r
library(rinet)
library(reflimR)

# Load liver test data
data(livertests)

# Filter for males and extract ALB and PROT values
mixture_2d <- livertests[livertests$Sex == "m", c("ALB", "PROT")]
mixture_2d <- mixture_2d[complete.cases(mixture_2d), ]

# Method 1: RINet prediction from mixture
rinet_result <- predict_rinet(mixture_2d)

cat("RINet Prediction (from mixture):\n")
cat(sprintf("  Log-Mean ALB:  %.2f\n", rinet_result[[1]]$mean[1]))
cat(sprintf("  Log-Mean PROT: %.2f\n", rinet_result[[1]]$mean[2]))
cat(sprintf("  Log-SD ALB:    %.2f\n", rinet_result[[1]]$std[1]))
cat(sprintf("  Log-SD PROT:   %.2f\n", rinet_result[[1]]$std[2]))
cat(sprintf("  Correlation: %.3f\n", rinet_result[[1]]$correlation))
cat(sprintf("  Reference Fraction: %.2f\n", rinet_result[[1]]$reference_fraction))
cat(sprintf("  Reference Region: %d ellipse vertices (original scale)\n", 
    nrow(rinet_result[[1]]$reference_interval)))

# Optional: plot the reference ellipse
# plot(rinet_result[[1]]$reference_interval, type="l", 
#      xlab="ALB", ylab="PROT", main="95% Reference Region")

# Method 2: Direct calculation (requires knowing which samples are reference)
reference_only <- livertests[livertests$Sex == "m" & 
                             livertests$Category == "reference", 
                             c("ALB", "PROT")]
reference_only <- reference_only[complete.cases(reference_only), ]

cat("\nDirect Method (reference samples only):\n")
log_ref_2d <- log(reference_only)
log_mean_2d <- colMeans(log_ref_2d)
log_cov_2d <- cov(log_ref_2d)

# Create reference ellipse
prob <- 0.95
chi2_val <- qchisq(prob, df = 2)
eigen_decomp <- eigen(log_cov_2d)
theta <- seq(0, 2 * pi, length.out = 100)
ellipse_log <- matrix(NA, nrow = 100, ncol = 2)
for (j in 1:100) {
  point <- sqrt(chi2_val) * c(sqrt(eigen_decomp$values[1]) * cos(theta[j]),
                               sqrt(eigen_decomp$values[2]) * sin(theta[j]))
  ellipse_log[j, ] <- eigen_decomp$vectors %*% point + log_mean_2d
}
ellipse_original <- exp(ellipse_log)

cat(sprintf("  Log-Mean ALB:  %.2f\n", log_mean_2d[1]))
cat(sprintf("  Log-Mean PROT: %.2f\n", log_mean_2d[2]))
cat(sprintf("  Log-SD ALB:    %.2f\n", sqrt(log_cov_2d[1,1])))
cat(sprintf("  Log-SD PROT:   %.2f\n", sqrt(log_cov_2d[2,2])))
cat(sprintf("  Correlation: %.3f\n", cor(log_ref_2d$ALB, log_ref_2d$PROT)))
cat(sprintf("  Reference Fraction: %.2f\n", 
    nrow(reference_only) / nrow(mixture_2d)))
cat(sprintf("  Reference Region: %d ellipse vertices (original scale)\n", 
    nrow(ellipse_original)))
```

## How it Works

1. **Standardization**: Input data is standardized (mean=0, sd=1)
2. **Feature Extraction**: Histograms are computed and normalized
3. **CNN Prediction**: Features are fed through the trained CNN
4. **Inverse Transform**: Outputs are scaled back using the scaler
5. **Destandardization**: Statistics are transformed to original scale

Models are automatically loaded on first use and cached for efficiency.

## Package Structure

```
rinet/
├── DESCRIPTION
├── NAMESPACE
├── LICENSE
├── R/
│   └── model_functions.R      # predict_rinet_1d(), predict_rinet_2d()
├── man/
│   └── rinet-package.Rd
└── inst/
    └── models/
        ├── rinet_1d.keras     # 1D CNN model
        ├── rinet_2d.keras     # 2D CNN model
        ├── scaler_1d.pkl      # 1D scaler
        └── scaler_2d.pkl      # 2D scaler
```

## Development

Rebuild documentation:
```r
devtools::document()
```

Check package:
```r
devtools::check()
```

## License

MIT License - see LICENSE file for details

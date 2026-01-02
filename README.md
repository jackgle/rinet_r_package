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
cat(sprintf("  Mean: %.2f\n", rinet_result[[1]]$mean))
cat(sprintf("  SD:   %.2f\n", rinet_result[[1]]$std))
cat(sprintf("  Reference Fraction: %.2f\n", rinet_result[[1]]$reference_fraction))

# Method 2: Direct calculation (requires knowing which samples are reference)
reference_only <- livertests[livertests$Sex == "f" & 
                             livertests$Category == "reference", "ALB"]

cat("\nDirect Method (reference samples only):\n")
cat(sprintf("  Mean: %.2f\n", mean(reference_only)))
cat(sprintf("  SD:   %.2f\n", sd(reference_only)))
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
cat(sprintf("  Mean ALB:  %.2f\n", rinet_result[[1]]$mean[1]))
cat(sprintf("  Mean PROT: %.2f\n", rinet_result[[1]]$mean[2]))
cat(sprintf("  SD ALB:    %.2f\n", rinet_result[[1]]$std[1]))
cat(sprintf("  SD PROT:   %.2f\n", rinet_result[[1]]$std[2]))
cat(sprintf("  Correlation: %.3f\n", rinet_result[[1]]$correlation))
cat(sprintf("  Reference Fraction: %.2f\n", rinet_result[[1]]$reference_fraction))

# Method 2: Direct calculation (requires knowing which samples are reference)
reference_only <- livertests[livertests$Sex == "m" & 
                             livertests$Category == "reference", 
                             c("ALB", "PROT")]
reference_only <- reference_only[complete.cases(reference_only), ]

cat("\nDirect Method (reference samples only):\n")
cat(sprintf("  Mean ALB:  %.2f\n", mean(reference_only$ALB)))
cat(sprintf("  Mean PROT: %.2f\n", mean(reference_only$PROT)))
cat(sprintf("  SD ALB:    %.2f\n", sd(reference_only$ALB)))
cat(sprintf("  SD PROT:   %.2f\n", sd(reference_only$PROT)))
cat(sprintf("  Correlation: %.3f\n", cor(reference_only$ALB, reference_only$PROT)))
cat(sprintf("  Reference Fraction: %.2f\n", 
    nrow(reference_only) / nrow(mixture_2d)))
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

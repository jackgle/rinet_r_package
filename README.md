# rinet

**RINet**: An R package for predicting clinical reference intervals from mixture distributions using a neural network.

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

# Install rinet from local directory
devtools::install("/Users/jack/Code/rinet")
```

Or if you're in the package directory:

```r
devtools::install()
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

### 1D Example

```r
library(rinet)

# Create a 1D mixture distribution
# Reference component (healthy): 70% of samples, mean=100, sd=15
# Pathological component 1: 20% of samples, mean=150, sd=20
# Pathological component 2: 10% of samples, mean=60, sd=10

n <- 1000
n_ref <- 700
n_path1 <- 200
n_path2 <- 100

mixture_1d <- c(
  rnorm(n_ref, mean = 100, sd = 15),   # Reference (healthy)
  rnorm(n_path1, mean = 150, sd = 20), # Pathological high
  rnorm(n_path2, mean = 60, sd = 10)   # Pathological low
)

# Predict reference distribution statistics
result <- predict_rinet(mixture_1d)

# Access results
result[[1]]$mean                  # Should be close to 100
result[[1]]$std                   # Should be close to 15
result[[1]]$reference_fraction    # Should be close to 0.70
```

### 2D Example

```r
library(rinet)
library(MASS)  # for mvrnorm

# Create a 2D mixture distribution
# Reference component (healthy): 75% of samples
ref_mean <- c(120, 80)
ref_sigma <- matrix(c(100, 30, 30, 64), 2, 2)  # correlation ~0.375

# Pathological component: 25% of samples
path_mean <- c(90, 120)
path_sigma <- matrix(c(144, -40, -40, 100), 2, 2)  # correlation ~-0.33

n <- 1000
n_ref <- 750
n_path <- 250

mixture_2d <- rbind(
  mvrnorm(n_ref, ref_mean, ref_sigma),   # Reference (healthy)
  mvrnorm(n_path, path_mean, path_sigma) # Pathological
)

# Predict reference distribution statistics
result <- predict_rinet(mixture_2d)

# Access results
result[[1]]$mean                  # Should be close to c(120, 80)
result[[1]]$std                   # Should be close to c(10, 8)
result[[1]]$correlation           # Should be close to 0.375
result[[1]]$reference_fraction    # Should be close to 0.75
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

\name{rinet-package}
\alias{rinet-package}
\alias{rinet}
\docType{package}
\title{
  RINet: Predict Clinical Reference Intervals from Mixture Distributions
}
\description{
  Predict the statistics of an underlying reference distribution from a 
  mixture distribution of raw clinical measurements (healthy and pathological 
  patients). Uses pre-trained convolutional neural networks to estimate 
  means, standard deviations, correlations, and reference fractions from 
  1D or 2D sample data.
}
\details{
  The main functions in this package are:
  \itemize{
    \item \code{predict_rinet()}: Automatically detect and predict (recommended)
    \item \code{predict_rinet_1d()}: Predict statistics from 1D mixture samples
    \item \code{predict_rinet_2d()}: Predict statistics from 2D mixture samples
  }
  
  Models are automatically loaded on first use and cached for efficiency.
  Model files should be placed in \code{inst/models/} with names:
  \itemize{
    \item \code{rinet_1d.keras} - 1D model
    \item \code{rinet_2d.keras} - 2D model
    \item \code{scaler_1d.pkl} - 1D scaler
    \item \code{scaler_2d.pkl} - 2D scaler
  }
}
\author{
  Jack
  
  Maintainer: Jack <your.email@example.com>
}
\keyword{package}

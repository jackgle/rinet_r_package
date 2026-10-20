#' @section Setup:
#' The models run on the torch backend, which is downloaded separately from
#' the package. Run \code{torch::install_torch()} once after installing rinet.
#'
#' @section Model versions:
#' Two 1D models are included. \code{"v2.1"} (the default) was retrained on
#' revised simulations with pathological shifts and widths drawn relative to
#' the reference component, and \code{"v2.0"} is the model described in the
#' published paper. Select one with the \code{model_version} argument of
#' \code{\link{predict_rinet}} or \code{\link{predict_rinet_1d}}. The 2D model
#' is v2.0. See \code{news(package = "rinet")} for evaluation details. The
#' simulation and training code behind v2.1 is available at
#' \url{https://github.com/jackgle/RINet}.
#'
#' @keywords internal
"_PACKAGE"

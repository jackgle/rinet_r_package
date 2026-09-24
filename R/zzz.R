.onAttach <- function(libname, pkgname) {
  if (!torch::torch_is_installed()) {
    packageStartupMessage(
      "rinet needs a one-time setup. Run:\n",
      "  torch::install_torch()\n",
      "See ?predict_rinet for details."
    )
  }
}

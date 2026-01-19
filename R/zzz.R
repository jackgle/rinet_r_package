.onLoad <- function(libname, pkgname) {
  # Configure Python environment
  if (Sys.which("python") != "") {
    reticulate::use_python(Sys.which("python"), required = FALSE)
  }
  
  # Set TensorFlow logging to ERROR only (level 3 = most quiet)
  # These must be set BEFORE TensorFlow is imported
  Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3")
  Sys.setenv(GRPC_VERBOSITY = "ERROR")
  Sys.setenv(GLOG_minloglevel = "3")
  Sys.setenv(PYTHONWARNINGS = "ignore")
  Sys.setenv(TF_ENABLE_ONEDNN_OPTS = "0")
  Sys.setenv(TF_XLA_FLAGS = "--tf_xla_auto_jit=0")
  Sys.setenv(ABSL_LOG_SEVERITY = "3")

  # Pre-import TensorFlow with stderr redirected to devnull to swallow C++ logs
  # that bypass R's sink() (like absl and XLA initialization chatter)
  tryCatch({
    reticulate::py_run_string(
      paste(
        "import os, sys",
        "try:",
        "    # Redirect stderr (fd 2) to devnull",
        "    orig_stderr = os.dup(2)",
        "    devnull = os.open(os.devnull, os.O_WRONLY)",
        "    os.dup2(devnull, 2)",
        "    try:",
        "        import tensorflow as tf",
        "        _ = tf.constant(0)",
        "        # Silence absl logging if it was initialized",
        "        if 'absl.logging' in sys.modules:",
        "            import absl.logging",
        "            absl.logging.set_verbosity(absl.logging.ERROR)",
        "        ",
        "        # Define a helper for silent model loading used by R",
        "        def load_model_silent(path):",
        "            orig = os.dup(2)",
        "            dn = os.open(os.devnull, os.O_WRONLY)",
        "            os.dup2(dn, 2)",
        "            try:",
        "                return tf.keras.models.load_model(path)",
        "            finally:",
        "                os.dup2(orig, 2)",
        "                os.close(orig)",
        "                os.close(dn)",
        "        ",
        "        # Define a helper for silent prediction used by R",
        "        def predict_silent(model, x, verbose=0):",
        "            orig = os.dup(2)",
        "            dn = os.open(os.devnull, os.O_WRONLY)",
        "            os.dup2(dn, 2)",
        "            try:",
        "                return model.predict(x, verbose=verbose)",
        "            finally:",
        "                os.dup2(orig, 2)",
        "                os.close(orig)",
        "                os.close(dn)",
        "    finally:",
        "        # Restore stderr",
        "        os.dup2(orig_stderr, 2)",
        "        os.close(orig_stderr)",
        "        os.close(devnull)",
        "except Exception:",
        "    pass",
        sep = "\n"
      )
    )
    # Export the Python helpers to the package namespace for use in model_functions.R
    assign(".py_silence", reticulate::import_main(), envir = asNamespace(pkgname))
  }, error = function(e) {})
  
  invisible()
}

.onAttach <- function(libname, pkgname) {
  # Check for all required Python packages
  missing_packages <- c()
  if (!reticulate::py_module_available("tensorflow")) {
    missing_packages <- c(missing_packages, "tensorflow")
  }
  if (!reticulate::py_module_available("keras")) {
    missing_packages <- c(missing_packages, "keras")
  }
  if (!reticulate::py_module_available("sklearn")) {
    missing_packages <- c(missing_packages, "scikit-learn")
  }
  
  if (length(missing_packages) > 0) {
    packageStartupMessage(
      "rinet requires Python packages: ", paste(missing_packages, collapse = ", "), "\n",
      "Install them with:\n",
      "  library(reticulate)\n",
      "  py_install(c(\"tensorflow\", \"keras\", \"scikit-learn\"))\n",
      "See README for details."
    )
  }
}
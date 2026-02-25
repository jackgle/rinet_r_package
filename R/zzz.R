.onLoad <- function(libname, pkgname) {
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
  # TF silencing block - optional, failure is fine
  tryCatch({
    reticulate::py_run_string(
      paste(
        "import os, sys",
        "try:",
        "    orig_stderr = os.dup(2)",
        "    devnull = os.open(os.devnull, os.O_WRONLY)",
        "    os.dup2(devnull, 2)",
        "    try:",
        "        import tensorflow as tf",
        "        _ = tf.constant(0)",
        "        if 'absl.logging' in sys.modules:",
        "            import absl.logging",
        "            absl.logging.set_verbosity(absl.logging.ERROR)",
        "    finally:",
        "        os.dup2(orig_stderr, 2)",
        "        os.close(orig_stderr)",
        "        os.close(devnull)",
        "except Exception:",
        "    pass",
        sep = "\n"
      )
    )
  }, error = function(e) {})

  # Define helpers independently - only requires os, no TF at definition time
  tryCatch({
    reticulate::py_run_string(
      paste(
        "import os as _os",
        "",
        "def load_model_silent(path):",
        "    import tensorflow as _tf",
        "    orig = _os.dup(2)",
        "    dn = _os.open(_os.devnull, _os.O_WRONLY)",
        "    _os.dup2(dn, 2)",
        "    try:",
        "        return _tf.keras.models.load_model(path)",
        "    finally:",
        "        _os.dup2(orig, 2)",
        "        _os.close(orig)",
        "        _os.close(dn)",
        "",
        "def predict_silent(model, x, verbose=0):",
        "    orig = _os.dup(2)",
        "    dn = _os.open(_os.devnull, _os.O_WRONLY)",
        "    _os.dup2(dn, 2)",
        "    try:",
        "        return model.predict(x, verbose=verbose)",
        "    finally:",
        "        _os.dup2(orig, 2)",
        "        _os.close(orig)",
        "        _os.close(dn)",
        sep = "\n"
      )
    )
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
      "Then restart R and reload rinet."
    )
  }
}
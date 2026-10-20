## Submission summary

This is an update from rinet 0.1.1. It replaces the package's Python
dependency (reticulate/TensorFlow) with the torch R package, and adds a
retrained default 1D model while keeping the published model available.
See NEWS.md for details.

## Maintainer email change

The maintainer email changes from jackgl4124@gmail.com to
jlebien@abartyshealth.com (same maintainer). I will confirm the change from
the previous address.

## Package size

The installed package is about 7.6 MB, of which 7.4 MB are the pretrained
model weights in inst/models. These are required for the package's only
functionality (prediction) and cannot be downloaded separately.

## Examples

Examples are wrapped in \dontrun{} because they require the libtorch
backend, which the torch package downloads at runtime with
torch::install_torch() and which is not available on CRAN check machines.

## Test environments

* local macOS (aarch64), R 4.6.1
* win-builder (R-release, 4.6.1 ucrt)
* win-builder (R-devel): TODO

## R CMD check results

0 errors | 0 warnings | 1 note

* New maintainer (email change, see above).

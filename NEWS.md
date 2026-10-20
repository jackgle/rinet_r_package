# rinet 0.2.0

## Backend: torch replaces Python/TensorFlow

* The package no longer requires Python, TensorFlow, Keras, or scikit-learn.
  Models now run through the `torch` R package. After installing rinet, run
  `torch::install_torch()` once to download the backend.
* This is an infrastructure change only: predictions from the v2.0 model are
  identical to rinet 0.1.x up to floating-point precision.

## New default 1D model (v2.1)

* The default 1D model is now **v2.1**. It has the same inputs and outputs as
  v2.0 but was retrained on revised simulated data: pathological component
  shifts and widths are drawn relative to the reference component (instead of
  in absolute units with clamping), and the training set is five times larger.
  This gives more even coverage of mixtures in which the pathological
  component is far from, or narrower than, the reference.
* The model described in the published paper is still available as
  `model_version = "v2.0"` in `predict_rinet()` and `predict_rinet_1d()`.
* Each prediction now includes a `model_version` element.
* The 2D model is unchanged (v2.0).
* The updated simulation and training code is available at
  <https://github.com/jackgle/RINet>.

### Evaluation of v2.1 against v2.0

* RIbench benchmark (Ammer et al., 2022; 5,760 simulated test sets, 10
  analytes; absolute z-score deviation, lower is better): overall score
  improved from 0.395 to 0.356.
  * Largest gains at high pathological fractions, e.g. normally distributed
    analytes improved from 0.523 to 0.379.
  * Some loss at low pathological fractions for skewed analytes: skewed
    0.364 to 0.371, heavily skewed 0.249 to 0.267.
* Liver function test dataset (10 analytes x 2 sexes, normalized error
  against the direct method on reference-labelled patients): mean error
  0.061 (v2.0) vs 0.067 (v2.1), within the variation expected for 20 cases.

# rinet 0.1.1

* More robust Python environment setup and clearer installation message.

# rinet 0.1.0

* Initial CRAN release.

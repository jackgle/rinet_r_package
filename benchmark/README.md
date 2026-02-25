# RIbench Benchmark Evaluation for RINet

Evaluates RINet against the [RIbench](https://cran.r-project.org/web/packages/RIbench/) benchmark suite (Ammer et al., Clin Chem 2022), which contains 5,760 simulated test sets across 10 biomarkers.

## Setup

```r
install.packages("RIbench")
```

The `rinet` package must also be installed (or loaded via `devtools::load_all()`).

## Usage

From the package root directory:

```bash
# Quick test (3 test sets per biomarker)
Rscript benchmark/run_benchmark.R

# Full benchmark (5,760 test sets — takes hours)
Rscript benchmark/run_benchmark.R all

# Single biomarker
Rscript benchmark/run_benchmark.R Ca

# Single distribution type: normal, skewed, heavilySkewed, shifted
Rscript benchmark/run_benchmark.R skewed
```

## Files

- `rinet_wrapper.R` — Wrapper function adapting `predict_rinet_1d()` to the RIbench interface (`requirePercentiles = TRUE`)
- `run_benchmark.R` — Main script that generates test sets, runs RINet, and evaluates results
- `results/` — Generated output directory containing data, results, and evaluation plots

## Notes

- RINet assumes **lognormal** distributions (`log_scale = TRUE`) for all biomarkers. This may underperform on the normal-distribution biomarkers (Hb, Ca, FT4) in the benchmark.
- The benchmark score uses absolute z-score deviation as the default error metric, with a cutoff of 5 for implausible results.

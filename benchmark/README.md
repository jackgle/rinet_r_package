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
Rscript benchmark/run_benchmark_rinet.R test

# Full benchmark (5,760 test sets)
Rscript benchmark/run_benchmark_rinet.R all

# Single biomarker
Rscript benchmark/run_benchmark_rinet.R Ca

# Single distribution type: normal, skewed, heavilySkewed, shifted
Rscript benchmark/run_benchmark_rinet.R skewed

# Disable log-scaling for normal-distributed analytes (Hb, Ca, FT4)
Rscript benchmark/run_benchmark_rinet.R --no-log-normal

# With outlier removal (Tukey 1.5x IQR on log-scale)
Rscript benchmark/run_benchmark_rinet.R --outliers

# Custom output directory
Rscript benchmark/run_benchmark_rinet.R --outdir benchmark/results_no_log_normal --no-log-normal

# Custom working directory (where RIbench Data/ lives)
Rscript benchmark/run_benchmark_rinet.R --workdir ../rinet_v1/data/RIbench --outdir benchmark/results
```

### eval_refiner.R

Evaluates pre-computed refineR prediction CSVs (no model run).

```bash
Rscript benchmark/eval_refiner.R <predictions_dir> [--outdir PATH] [subset]

# Example
Rscript benchmark/eval_refiner.R ../rinet_v1/evaluation/ribench/refineR_predictions \
  --outdir benchmark/results_refiner
```

## Files

- `run_benchmark_rinet.R` — Generates test sets, runs RINet, and evaluates results
- `eval_refiner.R` — Evaluates pre-computed refineR prediction CSVs using the same RIbench pipeline
- `results/` — Generated output directory containing data, results, and evaluation plots

## Notes

- RINet uses `log_scale = TRUE` by default for all biomarkers. Use `--no-log-normal` to disable this for normal-distributed analytes (Hb, Ca, FT4).
- The benchmark score uses absolute z-score deviation as the default error metric, with a cutoff of 5 for implausible results.

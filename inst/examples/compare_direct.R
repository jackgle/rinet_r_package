library(rinet)
library(reflimR)

# Load liver test data (mixture of reference and pathological samples)
data(livertests)

# Filter for females and extract ALB values
alb_mixture <- livertests[livertests$Sex == "f", "ALB"]

# Method 1: RINet prediction from mixture
rinet_result <- predict_rinet(alb_mixture)
# optionally with bootstrapping for confidence intervals
# rinet_result <- predict_rinet(alb_mixture, n_bootstrap = 100, confidence_level = 0.95)

cat("RINet Prediction (from mixture):\n")
cat(sprintf("  Reference Interval: [%.2f, %.2f]\n", 
    rinet_result[[1]]$reference_interval["lower"],
    rinet_result[[1]]$reference_interval["upper"]))
cat(sprintf("  Log-Mean: %.2f\n", rinet_result[[1]]$mean))
cat(sprintf("  Log-SD:   %.2f\n", rinet_result[[1]]$std))
cat(sprintf("  Reference Fraction: %.2f\n", rinet_result[[1]]$reference_fraction))

# Method 2: Direct calculation (requires knowing which samples are reference)
reference_only <- livertests[livertests$Sex == "f" & 
                             livertests$Category == "reference", "ALB"]

cat("\nDirect Method (reference samples only):\n")
log_ref <- log(reference_only)
log_mean <- mean(log_ref)
log_sd <- sd(log_ref)
cat(sprintf("  Reference Interval: [%.2f, %.2f]\n",
    exp(log_mean + qnorm(0.025) * log_sd),
    exp(log_mean + qnorm(0.975) * log_sd)))
cat(sprintf("  Log-Mean: %.2f\n", log_mean))
cat(sprintf("  Log-SD:   %.2f\n", log_sd))
cat(sprintf("  Reference Fraction: %.2f\n", 
    length(reference_only) / length(alb_mixture)))
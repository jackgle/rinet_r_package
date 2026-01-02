library(rinet)
library(reflimR)

# List of analytes
analytes <- c("ALB", "ALT", "AST", "BIL", "CHE", "CREA", "GGT", "PROT")
sexes <- unique(livertests$Sex)

results_list <- list()

for (s in sexes) {
  for (a in analytes) {
    # Filter data
    data_subset <- livertests[livertests$Sex == s, a]
    
    # Remove NAs if any
    data_subset <- data_subset[!is.na(data_subset)]
    
    if (length(data_subset) > 10) { # Minimum sample size check
      message(sprintf("Predicting for Sex: %s, Analyte: %s", s, a))
      
      # predict_rinet returns a list of results (one for each sample)
      # Since we pass one vector, we take the first element
      pred <- predict_rinet(data_subset)[[1]]
      
      results_list[[length(results_list) + 1]] <- data.frame(
        Sex = s,
        Analyte = a,
        Mean = pred$mean,
        Std = pred$std,
        RefFraction = pred$reference_fraction
      )
    }
  }
}

# Combine results
final_results <- do.call(rbind, results_list)
print(final_results)

# --- 2D Tests ---
message("\n--- Running 2D Tests ---")
pairs <- list(c("ALB", "PROT"), c("ALT", "AST"))
results_2d <- list()

for (s in sexes) {
  for (p in pairs) {
    # Filter data for both analytes and remove rows with any NAs
    data_subset <- livertests[livertests$Sex == s, p]
    data_subset <- data_subset[complete.cases(data_subset), ]
    
    if (nrow(data_subset) > 10) {
      message(sprintf("Predicting 2D for Sex: %s, Analytes: %s & %s", s, p[1], p[2]))
      
      pred <- predict_rinet(data_subset)[[1]]
      
      results_2d[[length(results_2d) + 1]] <- data.frame(
        Sex = s,
        Analyte1 = p[1],
        Analyte2 = p[2],
        Mean1 = pred$mean[1],
        Mean2 = pred$mean[2],
        Std1 = pred$std[1],
        Std2 = pred$std[2],
        Corr = pred$correlation,
        RefFraction = pred$reference_fraction
      )
    }
  }
}

final_results_2d <- do.call(rbind, results_2d)
print(final_results_2d)

#' DSEM Sensitivity Analysis Workflow
#'
#' This script automates a stress-test of the DSEM methodology by iteratively 
#' degrading data quality (missingness), increasing noise (variance), and 
#' shrinking sample sizes (n_steps). It records the Negative Log-Likelihood (NLL) 
#' percentile to identify the exact thresholds where the structural model breaks down.

# ==========================================
# 0. SETUP & SCRIPT SOURCING
# ==========================================
library(ggplot2)
library(dplyr)

# Source your established pipeline functions
# source("simulate_dsem_data.R")
# source("fit_generalized_dsem.R")
# source("randomize_dsem_graph.R")
# source("count_dsem_permutations.R")
# source("evaluate_random_graphs.R")
# source("fit_random_graphs.R")
# source("compare_dsem_distributions.R")

# Base Inputs
adj_file <-here::here('data-raw','example_unknown_dag.csv')
lags_file <- here::here('data-raw', 'example_lags_dag.csv')

# Define Node Classes based on your graph
node_classes <- list(
  source = c("bot.temp", "diatom"),
  intermediate = c("juv.bio", "adult.bio"),
  sink = c("catch"),
  measurement = c("growth", "stratification", "diatom.prop")
  # Note: 'food.quality' is Latent and already 100% NA, so we skip injecting NAs into it.
)

latent_dict <- list('food.quality' = c('growth', 'stratification', 'diatom.prop'))

# Helper Wrapper: Runs a single sensitivity iteration
run_sensitivity_iteration <- function(sim_data, ref_adj_file, ref_lags_file, iter_name, N_rand = 30) {
  
  # 1. Save temporary data
  tmp_data_file <- tempfile(fileext = ".rds")
  saveRDS(sim_data, tmp_data_file)
  
  # 2. Fit Reference
  ref_fit <- tryCatch({
    fit_generalized_dsem(ref_adj_file, ref_lags_file, tmp_data_file, 
                         latent_dict = latent_dict, detrend = TRUE, standardize = TRUE)
  }, error = function(e) return(NULL))
  
  if (is.null(ref_fit)) return(data.frame(Scenario = iter_name, NLL_Percentile = NA, Status = "Ref_Failed"))
  
  # 3. Generate Randoms & Fit Batch
  ref_bin <- matrix(as.numeric(read.csv(ref_adj_file, row.names = 1) != 0), nrow = 9)
  dimnames(ref_bin) <- dimnames(read.csv(ref_adj_file, row.names = 1))
  
  batch_eval <- evaluate_random_graphs(ref_bin, read.csv(ref_lags_file, row.names = 1), 
                                       N = N_rand, latent_dict = latent_dict)
  
  batch_fits <- fit_random_graphs(batch_eval, tmp_data_file, cores = 12,#parallel::detectCores() - 1, 
                                  latent_dict = latent_dict)
  
  # 4. Compare Distributions
  dist_results <- compare_dsem_distributions(ref_fit, batch_fits, plot = FALSE)
  
  nll_pct <- dist_results$model_fit_comparison$Percentile[1]
  
  return(data.frame(Scenario = iter_name, NLL_Percentile = nll_pct, Status = "Success"))
}

# Base Simulation Parameters
base_mu <- c('bot.temp'=5, 'diatom'=1E6, 'juv.bio'=1E5, 'adult.bio'=1E4, 'catch'=1E3, 'food.quality'=100, 'growth'=10, 'stratification'=10, 'diatom.prop'=0.7)
base_sd <- c('bot.temp'=1, 'diatom'=1E3, 'juv.bio'=1E3, 'adult.bio'=1E3, 'catch'=1E2, 'food.quality'=10, 'growth'=3, 'stratification'=3, 'diatom.prop'=0.05)
base_slope <- c('bot.temp'=0.025, 'diatom'=0, 'juv.bio'=-0.01, 'adult.bio'=-0.05, 'catch'=0, 'food.quality'=0.025, 'growth'=0.05, 'stratification'=0.01, 'diatom.prop'=0)

# Initialize Results Storage
sensitivity_results <- data.frame()


# ==========================================
# 0.5 PRE-FLIGHT CHECK: BASELINE ROBUSTNESS
# ==========================================
message("--- Pre-Flight Check: Validating Baseline Robustness ---")

# 1) Set a seed for reproducibility
set.seed(42) 

# 2) Define run parameters for the main reference model
M_preflight <- 30 
max_ts_length <- 400

preflight_params <- list()
preflight_nll <- numeric()
failed_iters <- 0

# 3) Run data generation and DSEM fitting M times
message(sprintf("Running baseline data generation and DSEM fitting %d times...", M_preflight))
for (i in 1:M_preflight) {
  
  sim_data <- simulate_dsem_data(adj_file, lags_file, n_steps = max_ts_length, 
                                 var.mu = base_mu, var.sd = base_sd, var.slope = base_slope, 
                                 latent_vars = 'food.quality')
  
  tmp_preflight_file <- tempfile(fileext = ".rds")
  saveRDS(sim_data, tmp_preflight_file)
  
  fit <- tryCatch({
    fit_generalized_dsem(adj_file, lags_file, tmp_preflight_file, 
                         latent_dict = latent_dict, detrend = TRUE, standardize = TRUE)
  }, error = function(e) return(NULL))
  
  # 4) Track failures, fit stats, and parameter estimates
  if (is.null(fit)) {
    failed_iters <- failed_iters + 1
  } else {
    preflight_nll <- c(preflight_nll, fit$fit$opt$objective)
    
    est <- fit$estimates
    est <- est[est$Source != "Time", ] # Remove Time covariate drift
    est$Path <- paste0(est$Source, " -> ", est$Target, " (Lag ", est$Lag, ")")
    est$Iter <- i
    preflight_params[[length(preflight_params) + 1]] <- est[, c("Iter", "Path", "Estimate", "Std_Error", "p_value")]
  }
}

success_rate <- (M_preflight - failed_iters) / M_preflight
message(sprintf("Pre-Flight Results: %d out of %d iterations succeeded (%.1f%% Success Rate).", 
                (M_preflight - failed_iters), M_preflight, success_rate * 100))

if (success_rate < 0.5) {
  stop("\nPRE-FLIGHT FAILED: The baseline model failed to solve on >50% of the simulated datasets. The model is too unstable to proceed to sensitivity testing. Consider boosting signal-to-noise ratio or increasing time series length.")
}

preflight_df <- do.call(rbind, preflight_params)

# 6) Return a summary table of the robustness across iterations
preflight_summary <- preflight_df %>%
  dplyr::group_by(Path) %>%
  dplyr::summarize(
    Mean_Estimate = round(mean(Estimate, na.rm = TRUE), 4),
    SD_Estimate = round(sd(Estimate, na.rm = TRUE), 4),
    Mean_SE = round(mean(Std_Error, na.rm = TRUE), 4),
    .groups = 'drop'
  )

message("\nPre-Flight Parameter Robustness Summary:")
print(preflight_summary)
message(sprintf("\nMean Baseline NLL: %.2f (SD: %.2f)", mean(preflight_nll), sd(preflight_nll)))

# 5) Plot a histogram of the parameter estimates
p_preflight <- ggplot(preflight_df, aes(x = Estimate)) +
  geom_histogram(bins = 20, fill = "seagreen", color = "darkgreen", alpha = 0.8) +
  facet_wrap(~Path, scales = "free") +
  theme_minimal() +
  labs(title = "Pre-Flight Check: Baseline Parameter Stability",
       subtitle = sprintf("Distribution of estimates across %d baseline simulations", M_preflight - failed_iters),
       x = "Parameter Estimate", y = "Frequency")
print(p_preflight)

message("\nSUCCESS: Pre-flight robustness check passed! Generating final Master Baseline Dataset...")

# Generate ONE master perfect dataset to ensure apples-to-apples comparison 
# for the missing data and time-series length experiments.
sim_data_master <- simulate_dsem_data(adj_file, lags_file, n_steps = max_ts_length, 
                                      var.mu = base_mu, var.sd = base_sd, var.slope = base_slope, 
                                      latent_vars = 'food.quality')


# ==========================================
# EXPERIMENT 1: MISSING DATA TOLERANCE
# ==========================================
message("Starting Exp 1: Missing Data...")
missing_fractions <- c(0, 0.1, 0.25, 0.50) # 0%, 10%, 25%, 50% missing data

for (class_name in names(node_classes)) {
  for (frac in missing_fractions) {
    
    # Copy the master data
    sim_data <- sim_data_master
    
    # Inject NAs into target class
    targets <- node_classes[[class_name]]
    for (node in targets) {
      if (frac > 0) {
        na_indices <- sample(1:max_ts_length, size = floor(max_ts_length * frac))
        sim_data[na_indices, node] <- NA
      }
    }
    
    iter_label <- paste0("Missing_", class_name, "_", frac*100, "pct")
    message("  Running: ", iter_label)
    
    res <- run_sensitivity_iteration(sim_data, adj_file, lags_file, iter_label)
    res$Experiment <- "1_Missing_Data"
    res$Class <- class_name
    res$Level <- frac
    
    sensitivity_results <- rbind(sensitivity_results, res)
  }
}

# ==========================================
# EXPERIMENT 2: TIME-SERIES VARIANCE (NOISE)
# ==========================================
message("Starting Exp 2: Variance Overload...")
variance_multipliers <- c(1, 2, 5, 10) # 1x (Baseline), 2x, 5x, 10x noise

for (class_name in names(node_classes)) {
  for (mult in variance_multipliers) {
    
    # Because we are changing the underlying generation math, we MUST simulate new data here
    mod_sd <- base_sd
    targets <- node_classes[[class_name]]
    mod_sd[targets] <- mod_sd[targets] * mult
    
    sim_data <- simulate_dsem_data(adj_file, lags_file, n_steps = max_ts_length, 
                                   var.mu = base_mu, var.sd = mod_sd, var.slope = base_slope, 
                                   latent_vars = 'food.quality')
    
    iter_label <- paste0("Noise_", class_name, "_", mult, "x")
    message("  Running: ", iter_label)
    
    res <- run_sensitivity_iteration(sim_data, adj_file, lags_file, iter_label)
    res$Experiment <- "2_Variance"
    res$Class <- class_name
    res$Level <- mult
    
    sensitivity_results <- rbind(sensitivity_results, res)
  }
}

# ==========================================
# EXPERIMENT 3: EDGE DENSITY TO TS LENGTH
# ==========================================
message("Starting Exp 3: Time Series Length...")
ts_lengths <- c(200, 100, 50, 30)

for (n_len in ts_lengths) {
  
  # Simply subset the master dataset
  sim_data <- sim_data_master[1:n_len, ]
  
  iter_label <- paste0("Length_N", n_len)
  message("  Running: ", iter_label)
  
  res <- run_sensitivity_iteration(sim_data, adj_file, lags_file, iter_label)
  res$Experiment <- "3_TS_Length"
  res$Class <- "All"
  res$Level <- n_len
  
  sensitivity_results <- rbind(sensitivity_results, res)
}

saveRDS(sensitivity_results, 'Z:/atlantiseof/data/sensitivity_results.rds')
# ==========================================
# VIEW RESULTS & PLOT
# ==========================================
print(sensitivity_results)

ggplot(sensitivity_results %>% filter(!is.na(NLL_Percentile)), aes(x = Level, y = NLL_Percentile, color = Class)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  facet_wrap(~Experiment, scales = "free_x") +
  theme_minimal() +
  labs(title = "DSEM Sensitivity Analysis",
       subtitle = "Red line indicates p=0.05 (significance threshold)",
       y = "Model Fit (NLL) Percentile vs Random Null",
       x = "Degradation Level (Fraction NA / Noise Multiplier / TS Length)")
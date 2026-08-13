#' Run a Single DSEM Sensitivity Iteration (Helper)
#' 
#' @export
run_sensitivity_iteration <- function(sim_data, ref_adj_file, ref_lags_file, iter_name, 
                                      latent_dict, N_rand, cores) {
  
  tmp_data_file <- tempfile(fileext = ".rds")
  saveRDS(sim_data, tmp_data_file)
  
  message("    -> Fitting reference model...")
  ref_fit <- tryCatch({
    fit_generalized_dsem(ref_adj_file, ref_lags_file, tmp_data_file, 
                         latent_dict = latent_dict, detrend = TRUE, standardize = TRUE)
  }, error = function(e) return(NULL))
  
  if (is.null(ref_fit)) {
    message("    -> [FAILED] Reference model could not solve this dataset.")
    return(data.frame(Scenario = iter_name, NLL_Percentile = NA, Status = "Ref_Failed"))
  }
  
  ref_adj_raw <- read.csv(ref_adj_file, row.names = 1)
  ref_bin <- matrix(as.numeric(ref_adj_raw != 0), nrow = nrow(ref_adj_raw))
  dimnames(ref_bin) <- dimnames(ref_adj_raw)
  
  message(sprintf("    -> Generating %d random topologies...", N_rand))
  batch_eval <- evaluate_random_graphs(ref_bin, read.csv(ref_lags_file, row.names = 1), 
                                       N = N_rand, latent_dict = latent_dict)
  
  message(sprintf("    -> Fitting %d random models in parallel (cores = %d)...", N_rand, cores))
  batch_fits <- fit_random_graphs(batch_eval, tmp_data_file, cores = cores, latent_dict = latent_dict)
  
  dist_results <- compare_dsem_distributions(ref_fit, batch_fits, plot = FALSE)
  nll_pct <- dist_results$model_fit_comparison$Percentile[1]
  
  message(sprintf("    -> [COMPLETE] NLL Percentile: %.4f", nll_pct))
  
  # Aggressively clear large objects and force RAM garbage collection
  rm(batch_fits, batch_eval, ref_fit)
  gc()
  
  return(data.frame(Scenario = iter_name, NLL_Percentile = nll_pct, Status = "Success"))
}

#' Run DSEM Sensitivity Analysis Workflow
#'
#' This function automates a stress-test of the DSEM methodology by iteratively 
#' degrading data quality (missingness), increasing noise (variance), and 
#' shrinking sample sizes (n_steps). It records the Negative Log-Likelihood (NLL) 
#' percentile to identify the exact thresholds where the structural model breaks down.
#'
#' @param adj_file Character. Path to the reference DAG adjacency matrix CSV.
#' @param lags_file Character. Path to the reference DAG temporal lags matrix CSV.
#' @param node_classes List. Named list mapping node classes to variable names 
#'        (e.g., list(source = c("A"), intermediate = c("B"), sink = c("C"), measurement = c("D"))).
#' @param latent_dict List. Optional named list mapping latent variables to their indicators. Defaults to NULL.
#' @param out_dir Character. Directory path where the results RDS and plots will be saved.
#' @param base_mu Numeric vector. A named vector of baseline means for simulation.
#' @param base_sd Numeric vector. A named vector of baseline standard deviations for simulation.
#' @param base_slope Numeric vector. A named vector of baseline time-trend slopes for simulation.
#' @param max_ts_length Integer. The maximum time series length for the baseline dataset. Defaults to 400.
#' @param M_preflight Integer. Number of iterations for the baseline robustness pre-flight check. Defaults to 30.
#' @param N_rand Integer. Number of random topologies to evaluate per sensitivity step. Defaults to 30.
#' @param cores Integer. Number of CPU cores to use for parallel processing. Defaults to 1.
#'
#' @return A list containing:
#'         \itemize{
#'           \item \code{results_df}: The full sensitivity analysis data frame.
#'           \item \code{preflight_summary}: The baseline robustness summary table.
#'           \item \code{plot_preflight}: The ggplot object for pre-flight baseline stability.
#'           \item \code{plot_sensitivity}: The ggplot object for the sensitivity degradation curves.
#'         }
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline facet_wrap theme_minimal labs ggsave geom_histogram
#' @importFrom dplyr group_by summarize filter %>%
#' @export
run_dsem_sensitivity <- function(adj_file, lags_file, node_classes, latent_dict = NULL, out_dir,
                                 base_mu, base_sd, base_slope, 
                                 max_ts_length = 400, M_preflight = 30, N_rand = 30, cores = 1) {
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  sensitivity_results <- data.frame()
  
  # ==========================================
  # 0.5 PRE-FLIGHT CHECK: BASELINE ROBUSTNESS
  # ==========================================
  message("\n==========================================")
  message("--- Pre-Flight Check: Validating Baseline Robustness ---")
  message("==========================================")
  
  preflight_params <- list()
  preflight_nll <- numeric()
  failed_iters <- 0
  
  message(sprintf("Running baseline data generation and DSEM fitting %d times...", M_preflight))
  
  for (i in 1:M_preflight) {
    
    # Progress print for the pre-flight loop (prints every 5 iterations to avoid spam)
    if (i %% 5 == 0 || i == 1) {
      message(sprintf("  -> Pre-flight iteration %d of %d...", i, M_preflight))
    }
    
    sim_data <- simulate_dsem_data(adj_file, lags_file, n_steps = max_ts_length, 
                                   var.mu = base_mu, var.sd = base_sd, var.slope = base_slope, 
                                   latent_vars = names(latent_dict))
    
    tmp_preflight_file <- tempfile(fileext = ".rds")
    saveRDS(sim_data, tmp_preflight_file)
    
    fit <- tryCatch({
      fit_generalized_dsem(adj_file, lags_file, tmp_preflight_file, 
                           latent_dict = latent_dict, detrend = TRUE, standardize = TRUE)
    }, error = function(e) return(NULL))
    
    if (is.null(fit)) {
      failed_iters <- failed_iters + 1
    } else {
      preflight_nll <- c(preflight_nll, fit$fit$opt$objective)
      
      est <- fit$estimates
      est <- est[est$Source != "Time", ] 
      est$Path <- paste0(est$Source, " -> ", est$Target, " (Lag ", est$Lag, ")")
      est$Iter <- i
      preflight_params[[length(preflight_params) + 1]] <- est[, c("Iter", "Path", "Estimate", "Std_Error", "p_value")]
    }
  }
  
  success_rate <- (M_preflight - failed_iters) / M_preflight
  message(sprintf("\nPre-Flight Results: %d out of %d iterations succeeded (%.1f%% Success Rate).", 
                  (M_preflight - failed_iters), M_preflight, success_rate * 100))
  
  if (success_rate < 0.5) {
    stop("\nPRE-FLIGHT FAILED: The baseline model failed to solve on >50% of datasets. The model is too unstable to proceed.")
  }
  
  preflight_df <- do.call(rbind, preflight_params)
  
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
  
  p_preflight <- ggplot2::ggplot(preflight_df, ggplot2::aes(x = Estimate)) +
    ggplot2::geom_histogram(bins = 20, fill = "seagreen", color = "darkgreen", alpha = 0.8) +
    ggplot2::facet_wrap(~Path, scales = "free") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Pre-Flight Check: Baseline Parameter Stability",
                  subtitle = sprintf("Distribution of estimates across %d baseline simulations", M_preflight - failed_iters),
                  x = "Parameter Estimate", y = "Frequency")
  
  ggplot2::ggsave(filename = file.path(out_dir, "preflight_stability.png"), plot = p_preflight, width = 10, height = 8)
  
  message("\nSUCCESS: Pre-flight check passed! Generating final Master Baseline Dataset...")
  
  
  sim_data_master <- simulate_dsem_data(adj_file, lags_file, n_steps = max_ts_length, 
                                        var.mu = base_mu, var.sd = base_sd, var.slope = base_slope, 
                                        latent_vars = names(latent_dict))
  
  
  # ==========================================
  # EXPERIMENT 1: MISSING DATA TOLERANCE
  # ==========================================
  message("\n==========================================")
  message("Starting Exp 1: Missing Data...")
  message("==========================================")
  missing_fractions <- c(0, 0.1, 0.25, 0.50)
  
  total_exp1 <- length(names(node_classes)) * length(missing_fractions)
  counter <- 1
  
  for (class_name in names(node_classes)) {
    for (frac in missing_fractions) {
      sim_data <- sim_data_master
      targets <- node_classes[[class_name]]
      
      for (node in targets) {
        if (frac > 0) {
          na_indices <- sample(1:max_ts_length, size = floor(max_ts_length * frac))
          sim_data[na_indices, node] <- NA
        }
      }
      
      iter_label <- paste0("Missing_", class_name, "_", frac*100, "pct")
      message(sprintf("\n[Exp 1: %d/%d] Running: %s", counter, total_exp1, iter_label))
      
      res <- run_sensitivity_iteration(sim_data, adj_file, lags_file, iter_label, latent_dict, N_rand, cores)
      res$Experiment <- "1_Missing_Data"
      res$Class <- class_name
      res$Level <- frac
      
      sensitivity_results <- rbind(sensitivity_results, res)
      counter <- counter + 1
    }
  }
  saveRDS(sensitivity_results,file.path(out_dir,'sensitivity_results_exp1.rds'))
  
  # ==========================================
  # EXPERIMENT 2: TIME-SERIES VARIANCE (NOISE)
  # ==========================================
  message("\n==========================================")
  message("Starting Exp 2: Variance Overload...")
  message("==========================================")
  variance_multipliers <- c(1, 2, 5, 10)
  
  total_exp2 <- length(names(node_classes)) * length(variance_multipliers)
  counter <- 1
  
  for (class_name in names(node_classes)) {
    for (mult in variance_multipliers) {
      mod_sd <- base_sd
      targets <- node_classes[[class_name]]
      mod_sd[targets] <- mod_sd[targets] * mult
      
      sim_data <- simulate_dsem_data(adj_file, lags_file, n_steps = max_ts_length, 
                                     var.mu = base_mu, var.sd = mod_sd, var.slope = base_slope, 
                                     latent_vars = names(latent_dict))
      
      iter_label <- paste0("Noise_", class_name, "_", mult, "x")
      message(sprintf("\n[Exp 2: %d/%d] Running: %s", counter, total_exp2, iter_label))
      
      res <- run_sensitivity_iteration(sim_data, adj_file, lags_file, iter_label, latent_dict, N_rand, cores)
      res$Experiment <- "2_Variance"
      res$Class <- class_name
      res$Level <- mult
      
      sensitivity_results <- rbind(sensitivity_results, res)
      counter <- counter + 1
    }
  }
  saveRDS(sensitivity_results,file.path(out_dir,'sensitivity_results_exp2.rds'))
  
  
  # ==========================================
  # EXPERIMENT 3: EDGE DENSITY TO TS LENGTH
  # ==========================================
  message("\n==========================================")
  message("Starting Exp 3: Time Series Length...")
  message("==========================================")
  ts_lengths <- c(max_ts_length, floor(max_ts_length/2), floor(max_ts_length/4), 30)
  
  total_exp3 <- length(ts_lengths)
  counter <- 1
  for (n_len in ts_lengths) {
    sim_data <- sim_data_master[1:n_len, ]
    
    iter_label <- paste0("Length_N", n_len)
    message(sprintf("\n[Exp 3: %d/%d] Running: %s", counter, total_exp3, iter_label))
    
    res <- run_sensitivity_iteration(sim_data, adj_file, lags_file, iter_label, latent_dict, N_rand, cores)
    res$Experiment <- "3_TS_Length"
    res$Class <- "All"
    res$Level <- n_len
    
    sensitivity_results <- rbind(sensitivity_results, res)
    counter <- counter + 1
  }
  saveRDS(sensitivity_results,file.path(out_dir,'sensitivity_results_exp3.rds'))
  
  # ==========================================
  # SAVE RESULTS & PLOT
  # ==========================================
  message("\n==========================================")
  message("Saving final results and generating plots...")
  message("==========================================")
  
  
  
  p_sens <- ggplot2::ggplot(sensitivity_results %>% dplyr::filter(!is.na(NLL_Percentile)), 
                            ggplot2::aes(x = Level, y = NLL_Percentile, color = Class)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    ggplot2::facet_wrap(~Experiment, scales = "free_x") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "DSEM Sensitivity Analysis",
                  subtitle = "Red line indicates p=0.05 (significance threshold)",
                  y = "Model Fit (NLL) Percentile vs Random Null",
                  x = "Degradation Level (Fraction NA / Noise Multiplier / TS Length)")
  
  ggplot2::ggsave(filename = file.path(out_dir, "sensitivity_plot.png"), plot = p_sens, width = 12, height = 6)
   out.ls = list(
     results_df = sensitivity_results,
     preflight_summary = preflight_summary,
     plot_preflight = p_preflight,
     plot_sensitivity = p_sens
   )
  saveRDS(out.ls, file.path(out_dir, "sensitivity_results_all.rds"))
  
  message("Workflow complete!")
  
  return(list(
    results_df = sensitivity_results,
    preflight_summary = preflight_summary,
    plot_preflight = p_preflight,
    plot_sensitivity = p_sens
  ))
}
#' Run a Single DSEM Sensitivity Iteration (Helper)
#' 
#' @export
run_sensitivity_iteration <- function(sim_data, ref_adj_file, ref_lags_file, iter_name, 
                                      latent_dict, N_rand, cores, detrend, standardize, control_time_drift) {
  
  tmp_data_file <- tempfile(fileext = ".rds")
  saveRDS(sim_data, tmp_data_file)
  
  message("    -> Fitting reference model...")
  ref_fit <- tryCatch({
    fit_generalized_dsem(ref_adj_file, ref_lags_file, tmp_data_file, 
                         latent_dict = latent_dict,
                         detrend = detrend, 
                         standardize = standardize, 
                         control_time_drift = control_time_drift)
  }, error = function(e) return(NULL))
  
  if (is.null(ref_fit)) {
    message("    -> [FAILED] Reference model could not solve this dataset.")
    unlink(tmp_data_file)
    return(data.frame(Scenario = iter_name, NLL_Percentile = NA, Status = "Ref_Failed"))
  }
  
  ref_adj_raw <- read.csv(ref_adj_file, row.names = 1)
  ref_bin <- matrix(as.numeric(ref_adj_raw != 0), nrow = nrow(ref_adj_raw))
  dimnames(ref_bin) <- dimnames(ref_adj_raw)
  
  message(sprintf("    -> Generating %d random topologies...", N_rand))
  batch_eval <- evaluate_random_graphs(ref_bin, read.csv(ref_lags_file, row.names = 1), 
                                       N = N_rand, latent_dict = latent_dict)
  
  message(sprintf("    -> Fitting %d random models in parallel (cores = %d)...", N_rand, cores))
  batch_fits <- fit_random_graphs(batch_eval, tmp_data_file, cores = cores, latent_dict = latent_dict,
                                  detrend = detrend, 
                                  standardize = standardize, 
                                  control_time_drift = control_time_drift)
  
  dist_results <- compare_dsem_distributions(ref_fit, batch_fits, plot = FALSE)
  nll_pct <- dist_results$model_fit_comparison$Percentile[1]
  
  message(sprintf("    -> [COMPLETE] NLL Percentile: %.4f", nll_pct))
  
  # Aggressively clear large objects and force RAM garbage collection
  unlink(tmp_data_file)
  rm(batch_fits, batch_eval, ref_fit)
  gc(verbose = FALSE, reset = TRUE) # Deep garbage collection
  
  return(data.frame(Scenario = iter_name, NLL_Percentile = nll_pct, Status = "Success"))
}

#' Run DSEM Sensitivity Analysis Workflow
#'
#' @export
run_dsem_sensitivity <- function(weights_file, adj_file, lags_file, node_classes, latent_dict = NULL, out_dir,
                                 base_mu, base_sd, base_slope, 
                                 max_ts_length = 400, M_preflight = 30, N_rand = 30, cores = 1,
                                 missing_fractions = c(0, 0.1, 0.25, 0.50),
                                 variance_multipliers = c(1, 2, 5, 10),
                                 ts_length_fractions = c(1, 0.5, 0.25, 0.1),
                                 detrend = detrend, 
                                 standardize = standardize, 
                                 control_time_drift = control_time_drift,
                                 jumpstart = FALSE) {
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Create intermediates directory for memory-safe disk writing
  int_dir <- file.path(out_dir, "intermediates")
  if (!dir.exists(int_dir)) dir.create(int_dir, recursive = TRUE)
  
  # Paths for persistent pre-flight objects
  pf_summary_file <- file.path(int_dir, "preflight_summary.rds")
  pf_plot_file <- file.path(int_dir, "preflight_plot.rds")
  master_data_file <- file.path(int_dir, "sim_data_master.rds")
  
  # ==========================================
  # 0.5 PRE-FLIGHT CHECK: BASELINE ROBUSTNESS
  # ==========================================
  message("\n==========================================")
  message("--- Pre-Flight Check: Validating Baseline Robustness ---")
  message("==========================================")
  
  if (jumpstart && file.exists(pf_summary_file) && file.exists(master_data_file) && file.exists(pf_plot_file)) {
    message("Jumpstart enabled: Pre-flight and master dataset found. Skipping pre-flight loop...")
    preflight_summary <- readRDS(pf_summary_file)
    p_preflight <- readRDS(pf_plot_file)
    sim_data_master <- readRDS(master_data_file)
  } else {
    message(sprintf("Running baseline data generation and DSEM fitting %d times in parallel (cores = %d)...", M_preflight, cores))
    
    cl <- parallel::makeCluster(cores)
    parallel::clusterExport(cl, varlist = c(
      "weights_file", "adj_file", "lags_file", "max_ts_length", 
      "base_mu", "base_sd", "base_slope", "latent_dict",
      "simulate_dsem_data", "fit_generalized_dsem"
    ), envir = environment())
    
    preflight_results_list <- parallel::parLapply(cl, 1:M_preflight, function(i) {
      sim_data <- simulate_dsem_data(weights_file, lags_file, n_steps = max_ts_length, 
                                     var.mu = base_mu, var.sd = base_sd, var.slope = base_slope, 
                                     latent_vars = names(latent_dict),diagnostics = F)
      
      tmp_preflight_file <- tempfile(fileext = ".rds")
      saveRDS(sim_data, tmp_preflight_file)
      
      fit <- tryCatch({
        fit_generalized_dsem(adj_file, lags_file, tmp_preflight_file, 
                             latent_dict = latent_dict,
                             detrend = detrend, 
                             standardize = standardize, 
                             control_time_drift = control_time_drift)
      }, error = function(e) return(NULL))
      
      unlink(tmp_preflight_file) 
      
      if (is.null(fit)) {
        return(list(failed = TRUE, nll = NA, est = NULL))
      } else {
        est <- fit$estimates
        est <- est[est$Source != "Time", ] 
        est$Path <- paste0(est$Source, " -> ", est$Target, " (Lag ", est$Lag, ")")
        est$Iter <- i
        return(list(failed = FALSE, 
                    nll = fit$fit$opt$objective, 
                    est = est[, c("Iter", "Path", "Estimate", "Std_Error", "p_value")]))
      }
    })
    
    parallel::stopCluster(cl)
    gc(verbose = FALSE, reset = TRUE)
    
    failed_iters <- sum(sapply(preflight_results_list, function(x) x$failed))
    preflight_nll <- sapply(preflight_results_list, function(x) x$nll)
    preflight_nll <- preflight_nll[!is.na(preflight_nll)]
    
    preflight_params <- lapply(preflight_results_list, function(x) x$est)
    preflight_params <- preflight_params[!sapply(preflight_params, is.null)]
    
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
    saveRDS(preflight_summary, pf_summary_file)
    message(sprintf("\nMean Baseline NLL: %.2f (SD: %.2f)", mean(preflight_nll), sd(preflight_nll)))
    
    p_preflight <- ggplot2::ggplot(preflight_df, ggplot2::aes(x = Estimate)) +
      ggplot2::geom_histogram(bins = 20, fill = "seagreen", color = "darkgreen", alpha = 0.8) +
      ggplot2::facet_wrap(~Path, scales = "free") +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Pre-Flight Check: Baseline Parameter Stability",
                    subtitle = sprintf("Distribution of estimates across %d baseline simulations", M_preflight - failed_iters),
                    x = "Parameter Estimate", y = "Frequency")
    
    ggplot2::ggsave(filename = file.path(out_dir, "preflight_stability.png"), plot = p_preflight, width = 10, height = 8)
    saveRDS(p_preflight, pf_plot_file)
    
    message("\nSUCCESS: Pre-flight check passed! Generating final Master Baseline Dataset...")
    sim_data_master <- simulate_dsem_data(weights.file = weights_file, lags.file = lags_file, n_steps = max_ts_length, 
                                          var.mu = base_mu, var.sd = base_sd, var.slope = base_slope, 
                                          latent_vars = names(latent_dict))
    saveRDS(sim_data_master, master_data_file)
  }
  
  # ==========================================
  # EXPERIMENT 1: MISSING DATA TOLERANCE
  # ==========================================
  message("\n==========================================")
  message("Starting Exp 1: Missing Data...")
  message("==========================================")
  
  total_exp1 <- length(names(node_classes)) * length(missing_fractions)
  counter <- 1
  
  for (class_name in names(node_classes)) {
    for (frac in missing_fractions) {
      iter_label <- paste0("Missing_", class_name, "_", round(frac*100, 1), "pct")
      out_file <- file.path(int_dir, sprintf("exp1_iter_%03d.rds", counter))
      
      if (jumpstart && file.exists(out_file)) {
        message(sprintf("\n[Exp 1: %d/%d] Jumpstart: Skipping %s (Already completed)", counter, total_exp1, iter_label))
        counter <- counter + 1
        next
      }
      
      message(sprintf("\n[Exp 1: %d/%d] Running: %s", counter, total_exp1, iter_label))
      
      sim_data <- sim_data_master
      targets <- node_classes[[class_name]]
      
      for (node in targets) {
        if (frac > 0) {
          na_indices <- sample(1:max_ts_length, size = floor(max_ts_length * frac))
          sim_data[na_indices, node] <- NA
        }
      }
      
      res <- run_sensitivity_iteration(sim_data, adj_file, lags_file, iter_label, latent_dict, N_rand, cores,
                                       detrend = detrend, 
                                       standardize = standardize, 
                                       control_time_drift = control_time_drift)
      res$Experiment <- "1_Missing_Data"
      res$Class <- class_name
      res$Level <- frac
      
      # WRITE TO DISK AND PURGE MEMORY
      saveRDS(res, out_file)
      rm(res, sim_data)
      gc(verbose = FALSE, reset = TRUE)
      
      counter <- counter + 1
    }
  }
  
  # ==========================================
  # EXPERIMENT 2: TIME-SERIES VARIANCE (NOISE)
  # ==========================================
  message("\n==========================================")
  message("Starting Exp 2: Variance Overload...")
  message("==========================================")
  
  total_exp2 <- length(names(node_classes)) * length(variance_multipliers)
  counter <- 1
  
  for (class_name in names(node_classes)) {
    for (mult in variance_multipliers) {
      iter_label <- paste0("Noise_", class_name, "_", round(mult, 2), "x")
      out_file <- file.path(int_dir, sprintf("exp2_iter_%03d.rds", counter))
      
      if (jumpstart && file.exists(out_file)) {
        message(sprintf("\n[Exp 2: %d/%d] Jumpstart: Skipping %s (Already completed)", counter, total_exp2, iter_label))
        counter <- counter + 1
        next
      }
      
      message(sprintf("\n[Exp 2: %d/%d] Running: %s", counter, total_exp2, iter_label))
      
      mod_sd <- base_sd
      targets <- node_classes[[class_name]]
      mod_sd[targets] <- mod_sd[targets] * mult
      
      sim_data <- simulate_dsem_data(weights.file = weights_file, lags_file, n_steps = max_ts_length, 
                                     var.mu = base_mu, var.sd = mod_sd, var.slope = base_slope, 
                                     latent_vars = names(latent_dict))
      
      res <- run_sensitivity_iteration(sim_data, adj_file, lags_file, iter_label, latent_dict, N_rand, cores,
                                       detrend = detrend, 
                                       standardize = standardize, 
                                       control_time_drift = control_time_drift)
      res$Experiment <- "2_Variance"
      res$Class <- class_name
      res$Level <- mult
      
      # WRITE TO DISK AND PURGE MEMORY
      saveRDS(res, out_file)
      rm(res, sim_data, mod_sd)
      gc(verbose = FALSE, reset = TRUE)
      
      counter <- counter + 1
    }
  }
  
  # ==========================================
  # EXPERIMENT 3: EDGE DENSITY TO TS LENGTH
  # ==========================================
  message("\n==========================================")
  message("Starting Exp 3: Time Series Length...")
  message("==========================================")
  ts_lengths <- floor(ts_length_fractions * max_ts_length)
  
  total_exp3 <- length(ts_lengths)
  counter <- 1
  for (n_len in ts_lengths) {
    iter_label <- paste0("Length_N", n_len)
    out_file <- file.path(int_dir, sprintf("exp3_iter_%03d.rds", counter))
    
    if (jumpstart && file.exists(out_file)) {
      message(sprintf("\n[Exp 3: %d/%d] Jumpstart: Skipping %s (Already completed)", counter, total_exp3, iter_label))
      counter <- counter + 1
      next
    }
    
    message(sprintf("\n[Exp 3: %d/%d] Running: %s", counter, total_exp3, iter_label))
    
    sim_data <- sim_data_master[1:n_len, ]
    
    res <- run_sensitivity_iteration(sim_data, adj_file, lags_file, iter_label, latent_dict, N_rand, cores,
                                     detrend = detrend, 
                                     standardize = standardize, 
                                     control_time_drift = control_time_drift)
    res$Experiment <- "3_TS_Length"
    res$Class <- "All"
    res$Level <- n_len
    
    # WRITE TO DISK AND PURGE MEMORY
    saveRDS(res, out_file)
    rm(res, sim_data)
    gc(verbose = FALSE, reset = TRUE)
    
    counter <- counter + 1
  }
  
  # ==========================================
  # AGGREGATE RESULTS & PLOT
  # ==========================================
  message("\n==========================================")
  message("Aggregating intermediate files and generating plots...")
  message("==========================================")
  
  # Read only the experiment iteration files back into memory
  int_files <- list.files(int_dir, pattern = "^exp[123]_iter_[0-9]{3}\\.rds$", full.names = TRUE)
  if(length(int_files) > 0) {
    sensitivity_results <- do.call(rbind, lapply(int_files, readRDS))
  } else {
    stop("No intermediate experiment results found.")
  }
  
  # Save the finalized dataframe
  saveRDS(sensitivity_results, file.path(out_dir, "sensitivity_results_all_experiments.rds"))
  
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
  
  out.ls <- list(
    results_df = sensitivity_results,
    preflight_summary = preflight_summary,
    plot_preflight = p_preflight,
    plot_sensitivity = p_sens
  )
  saveRDS(out.ls, file.path(out_dir, "sensitivity_results_all.rds"))
  
  message("Workflow complete! Memory successfully managed.")
  
  return(out.ls)
}
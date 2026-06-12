#' Compare Reference DSEM against Randomized Graph Distributions
#'
#' Evaluates the uniqueness and fit of a reference DSEM model by comparing its 
#' parameter estimates and overall model fit (Negative Log-Likelihood) against 
#' an empirical null distribution generated from a batch of randomized graphs.
#'
#' @param ref_fit List. The fitted reference model (output of `fit_generalized_dsem`).
#' @param batch_fits List. The batch of randomized fits (output of `fit_random_graphs`).
#' @param plot Logical. If TRUE, generates and returns ggplot2 visualizations.
#'
#' @return A list containing:
#'         \itemize{
#'           \item \code{parameter_comparison}: Data frame of path percentiles.
#'           \item \code{model_fit_comparison}: Data frame of overall fit percentiles.
#'           \item \code{plots}: A list of ggplot objects (if plot = TRUE).
#'         }
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming you have a reference fit and batch fits:
#' dist_results <- compare_dsem_distributions(ref_fit, batch_fits, plot = TRUE)
#' 
#' # Print the tables
#' print(dist_results$model_fit_comparison)
#' print(dist_results$parameter_comparison)
#' 
#' # View the plots
#' dist_results$plots$fit_plot
#' dist_results$plots$param_plot
#' }

compare_dsem_distributions <- function(ref_fit, batch_fits, plot = FALSE) {
  
  # 1. Extract Reference Data
  ref_est <- ref_fit$estimates
  ref_est <- ref_est[ref_est$Source != "Time", ] # Remove Time covariate drift
  ref_est$Path <- paste0(ref_est$Source, " -> ", ref_est$Target, " (Lag ", ref_est$Lag, ")")
  
  ref_nll <- ref_fit$fit$opt$objective
  
  # 2. Extract Batch Data
  valid_iters <- which(!sapply(batch_fits, is.null))
  
  batch_nll <- sapply(batch_fits[valid_iters], function(x) x$fit$opt$objective)
  
  batch_param_list <- lapply(valid_iters, function(i) {
    b_fit <- batch_fits[[i]]
    est <- b_fit$estimates
    est <- est[est$Source != "Time", ]
    est$Path <- paste0(est$Source, " -> ", est$Target, " (Lag ", est$Lag, ")")
    data.frame(Iter = i, Path = est$Path, Estimate = est$Estimate, stringsAsFactors = FALSE)
  })
  batch_params <- do.call(rbind, batch_param_list)
  
  # 3. Parameter Estimation Comparison
  # We must compare the reference paths against the random graphs. 
  # If a random graph did NOT include a specific reference path, its structural effect is 0.
  param_grid <- expand.grid(Iter = valid_iters, Path = ref_est$Path, stringsAsFactors = FALSE)
  merged_params <- merge(param_grid, batch_params, by = c("Iter", "Path"), all.x = TRUE)
  merged_params$Estimate[is.na(merged_params$Estimate)] <- 0 
  
  param_comp_list <- lapply(ref_est$Path, function(p) {
    ref_val <- ref_est$Estimate[ref_est$Path == p]
    rand_vals <- merged_params$Estimate[merged_params$Path == p]
    
    # Calculate the percentile of the reference value within the empirical null
    pctile <- mean(rand_vals <= ref_val)
    
    data.frame(
      Path = p,
      Ref_Estimate = round(ref_val, 4),
      Rand_Mean = round(mean(rand_vals), 4),
      Rand_SD = round(sd(rand_vals), 4),
      Percentile = round(pctile, 4),
      stringsAsFactors = FALSE
    )
  })
  param_comp <- do.call(rbind, param_comp_list)
  
  # 4. Overall Model Fit Comparison (Negative Log-Likelihood)
  # For NLL, lower is better. A low percentile means the reference model 
  # had a smaller NLL than most random graphs (indicating superior fit).
  fit_pctile <- mean(batch_nll <= ref_nll)
  
  fit_comp <- data.frame(
    Metric = "Negative Log-Likelihood (NLL)",
    Ref_Value = round(ref_nll, 4),
    Rand_Mean = round(mean(batch_nll), 4),
    Rand_SD = round(sd(batch_nll), 4),
    Percentile = round(fit_pctile, 4),
    stringsAsFactors = FALSE
  )
  
  # 5. Optional Plotting
  out_plots <- list()
  if (plot) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("The 'ggplot2' package is required for plotting. Returning tables without plots.")
    } else {
      
      # Parameter Distributions Facet Plot
      p_param <- ggplot2::ggplot(merged_params, ggplot2::aes(x = Estimate)) +
        ggplot2::geom_histogram(bins = 30, fill = "gray75", color = "gray30") +
        ggplot2::geom_vline(data = ref_est, ggplot2::aes(xintercept = Estimate), 
                            color = "firebrick", linetype = "dashed", linewidth = 1) +
        ggplot2::facet_wrap(~Path, scales = "free") +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = "Structural Parameters: Reference vs. Random Topologies",
          subtitle = "Red dashed line indicates the parameter estimate from the reference model.",
          x = "Parameter Estimate", y = "Frequency",
          caption = "Note: A value of 0 indicates the path was absent in a randomized topology."
        )
      
      # Model Fit (NLL) Plot
      nll_df <- data.frame(NLL = batch_nll)
      p_fit <- ggplot2::ggplot(nll_df, ggplot2::aes(x = NLL)) +
        ggplot2::geom_histogram(bins = 30, fill = "steelblue", color = "midnightblue") +
        ggplot2::geom_vline(xintercept = ref_nll, color = "firebrick", 
                            linetype = "dashed", linewidth = 1.2) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = "Overall Model Fit: Reference vs. Random Topologies",
          subtitle = "Lower Negative Log-Likelihood (NLL) indicates a better fit to the data.",
          x = "Negative Log-Likelihood (NLL)", y = "Frequency",
          caption = "Red dashed line indicates the Reference Model NLL."
        )
      
      out_plots$param_plot <- p_param
      out_plots$fit_plot <- p_fit
    }
  }
  
  return(list(
    parameter_comparison = param_comp,
    model_fit_comparison = fit_comp,
    plots = out_plots
  ))
}
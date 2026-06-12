#' Compare DSEM Estimates against True Weights Matrix
#'
#' Compares the structural path coefficients estimated by \code{fit_generalized_dsem} 
#' against the known data-generating weights matrix. It handles structural mismatches 
#' by identifying successfully recovered paths, missing paths (False Negatives), 
#' and incorrectly specified paths (False Positives).
#'
#' @param dsem_output List. The output object returned by \code{fit_generalized_dsem}.
#' @param true_weights_file Character. Path to the CSV file containing the true 
#'        weighted adjacency matrix used to generate the data.
#'
#' @return A list containing:
#'         \itemize{
#'           \item \code{comparison}: A data frame aligned by path, showing True vs. 
#'                 Estimated weights, standard errors, p-values, and estimation Bias.
#'           \item \code{metrics}: A named vector of summary statistics including Mean 
#'                 Absolute Error (MAE), Root Mean Squared Error (RMSE), False Positives, 
#'                 and False Negatives.
#'         }
#' @export
#'
#' @examples
#' \dontrun{
#' # Fit a structural hypothesis (could be the true DAG or a modified test DAG)
#' dsem_results <- fit_generalized_dsem("test_adj.csv", "lags.csv", "sim_data.rds")
#' 
#' # Compare results back to the original ground-truth weights
#' evaluation <- compare_dsem_results(dsem_results, "original_weights.csv")
#' 
#' # Inspect structural accuracy metrics
#' print(evaluation$metrics)
#' 
#' # Inspect parameter-level bias
#' print(evaluation$comparison)
#' }



compare_dsem_results <- function(dsem_output, true_weights_file) {
  
    # 1. Ingest Ground-Truth Matrix
    true_mat <- as.matrix(read.csv(true_weights_file, row.names = 1))
    
    # Convert the true matrix into a long data frame of valid edges (Weight != 0)
    true_edges <- as.data.frame(as.table(true_mat), stringsAsFactors = FALSE)
    colnames(true_edges) <- c("Source", "Target", "True_Weight")
    true_edges <- true_edges[true_edges$True_Weight != 0, ]
    
    # 2. Extract and Clean Estimated Edges
    est_edges <- dsem_output$estimates
    
    # Handle potential capitalization differences from raw dsem summary output
    if ("p_value" %in% colnames(est_edges)) colnames(est_edges)[colnames(est_edges) == "p_value"] <- "P_Value"
    
    # Strip out the 'Time' control covariate, as it isn't part of the causal matrix
    est_edges <- est_edges[est_edges$Source != "Time", ]
    
    # Rename Lag to Estimated_Lag for final reporting
    colnames(est_edges)[colnames(est_edges) == "Lag"] <- "Estimated_Lag"
    
    # 3. Align True and Estimated Structures via Full Outer Join
    # (dsem already provides clean Source names and integer Lags, so no string parsing is needed)
    comp_df <- merge(true_edges, 
                     est_edges[, c("Source", "Target", "Estimated_Lag", "Estimate", "Std_Error", "P_Value")], 
                     by = c("Source", "Target"), 
                     all = TRUE)
    
    # Clean up missing values resulting from structural mismatches
    comp_df$True_Weight[is.na(comp_df$True_Weight)] <- 0
    comp_df$Status <- "Matched"
    
    # Identify Structural Discrepancies
    comp_df$Status[comp_df$True_Weight != 0 & is.na(comp_df$Estimate)] <- "False Negative (Missed Path)"
    comp_df$Status[comp_df$True_Weight == 0 & !is.na(comp_df$Estimate)] <- "False Positive (Extra Path)"
    
    # Calculate parameter estimation bias for matched or extra paths
    comp_df$Bias <- NA
    valid_est <- !is.na(comp_df$Estimate)
    comp_df$Bias[valid_est] <- comp_df$Estimate[valid_est] - comp_df$True_Weight[valid_est]
    
    # Organize columns for clear reading
    comp_df <- comp_df[, c("Source", "Target", "Status", "True_Weight", "Estimated_Lag", "Estimate", "Bias", "Std_Error", "P_Value")]
    rownames(comp_df) <- NULL
    
    
    
  
  # 4. Compute Global Network Validation Metrics
  matched_paths <- comp_df[comp_df$Status == "Matched", ]
  
  mae <- if(nrow(matched_paths) > 0) mean(abs(matched_paths$Bias)) else NA
  rmse <- if(nrow(matched_paths) > 0) sqrt(mean(matched_paths$Bias^2)) else NA
  false_positives <- sum(comp_df$Status == "False Positive (Extra Path)")
  false_negatives <- sum(comp_df$Status == "False Negative (Missed Path)")
  
  global_metrics <- c(
    Matched_Paths_MAE = round(mae, 4),
    Matched_Paths_RMSE = round(rmse, 4),
    False_Positives = false_positives,
    False_Negatives = false_negatives
  )
  
  return(list(
    comparison = comp_df,
    metrics = global_metrics
  ))
}
#' Batch Fit DSEM on Randomized Graphs
#'
#' Takes the output from \code{evaluate_random_graphs} and fits a generalized 
#' DSEM model to each generated graph structure using the provided simulated 
#' data and temporal lags.
#'
#' @param batch_results List. The output object from \code{evaluate_random_graphs} 
#'        containing \code{$graphs}, \code{$lags}, and \code{$diagnostics}.
#' @param data.file Character. Path to the simulated time series data. Supports `.csv` or `.rds`.
#' @param cores Integer. Number of CPU cores to use for parallel processing. Defaults to 1 (sequential).
#'
#' @return A named list of fitted DSEM objects, corresponding to each graph in the batch.
#'         If a specific topology fails to converge, its list element will be NULL.
#' @importFrom parallel makeCluster stopCluster parLapply clusterExport clusterEvalQ
#' @export
#'
#' @examples
#' \dontrun{
#' # 1. Generate unique random graphs
#' batch_results <- evaluate_random_graphs(true_adj, true_lags, N = 100, deduplicate = TRUE)
#' 
#' # 2. Fit them all to the simulated data in parallel (e.g., using 4 cores)
#' batch_fits <- fit_random_graphs(batch_results, "sim_data.rds", cores = 4)
#' 
#' # Inspect the first model's estimates
#' print(batch_fits[[1]]$estimates)
#' }

fit_random_graphs <- function(batch_results, data.file, cores = 3) {
  
  graphs_list <- batch_results$graphs
  lags_list <- batch_results$lags
  diag_df <- batch_results$diagnostics
  
  num_graphs <- length(graphs_list)
  
  message(sprintf("Starting batch fit for %d graphs using %d core(s)...", num_graphs, cores))
  
  # Define the core fitting logic as a standalone function for lapply/parLapply
  fit_single_graph <- function(i) {
    current_adj <- graphs_list[[i]]
    current_lags <- lags_list[[i]]
    
    # 1. Write the current adjacency and lag matrices to temporary CSV files
    tmp_adj_file <- tempfile(fileext = ".csv")
    tmp_lags_file <- tempfile(fileext = ".csv")
    
    write.csv(current_adj, file = tmp_adj_file, row.names = TRUE)
    write.csv(current_lags, file = tmp_lags_file, row.names = TRUE)
    
    # 2. Fit the model using a tryCatch block
    fit_res <- tryCatch({
      fit_generalized_dsem(adj.file = tmp_adj_file, 
                           lags.file = tmp_lags_file, 
                           data.file = data.file)
    }, error = function(e) {
      warning(sprintf("Fitting failed for graph index %d (Iter %d): %s", 
                      i, diag_df$Iter[i], e$message), call. = FALSE)
      return(NULL)
    })
    
    # 3. Cleanup temp files
    unlink(tmp_adj_file)
    unlink(tmp_lags_file)
    
    return(fit_res)
  }
  
  # Execute via Parallel Cluster or Sequential Loop
  if (cores > 1) {
    
    cl <- parallel::makeCluster(cores)
    
    # Export necessary variables and the fitting function to the worker nodes
    parallel::clusterExport(cl, 
                            varlist = c("graphs_list", "lags_list", "diag_df", 
                                        "data.file", "fit_generalized_dsem"), 
                            envir = environment())
    
    # Ensure workers have the required package loaded
    parallel::clusterEvalQ(cl, { library(dsem) })
    
    # Run the batch
    results_list <- parallel::parLapply(cl, 1:num_graphs, fit_single_graph)
    
    parallel::stopCluster(cl)
    
  } else {
    
    # Sequential execution with progress tracking
    results_list <- vector("list", num_graphs)
    for (i in 1:num_graphs) {
      results_list[[i]] <- fit_single_graph(i)
      if (i %% 10 == 0 || i == num_graphs) {
        message(sprintf("Processed %d / %d graphs...", i, num_graphs))
      }
    }
    
  }
  
  # Name the output list based on the Iteration (and Exact_Match_ID if available)
  if ("Exact_Match_ID" %in% colnames(diag_df)) {
    names(results_list) <- paste0("Iter_", diag_df$Iter, "_MatchID_", diag_df$Exact_Match_ID)
  } else {
    names(results_list) <- paste0("Iter_", diag_df$Iter)
  }
  
  message("Batch fitting complete.")
  return(results_list)
}
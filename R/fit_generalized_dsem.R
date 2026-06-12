#' Perform Generalized DSEM Analysis via the `dsem` Package
#'
#' Ingests a binary directed acyclic graph (DAG) and temporal lags, translates them
#' into the "arrow-and-lag" syntax required by the `dsem` package, and fits the 
#' dynamic structural equation model. `dsem` handles temporal alignments, missing 
#' data, and simultaneous estimation natively via TMB.
#'
#' @param adj.file Character. Path to the CSV file containing the binary adjacency matrix.
#'        1 indicates an assumed causal edge, 0 indicates no edge.
#' @param lags.file Character. Path to the CSV file containing the temporal lag matrix.
#' @param data.file Character. Path to the simulated time series data. Supports `.csv` or `.rds`.
#'
#' @return A list containing:
#'         \itemize{
#'           \item \code{syntax}: Character string of the generated `dsem` model syntax.
#'           \item \code{fit}: The fitted `dsem` object.
#'           \item \code{estimates}: A clean data frame of the estimated structural path coefficients.
#'         }
#' @importFrom dsem dsem dsem_control
#' @importFrom stats ts
#' @export
#'
#' @examples
#' \dontrun{
#' dsem_results <- fit_generalized_dsem("binary_adj.csv", "lags.csv", "sim_data.rds")
#' 
#' # View the generated dsem arrow-and-lag syntax
#' cat(dsem_results$syntax)
#' 
#' # View the formal SEM estimates
#' print(dsem_results$estimates)
#' }

fit_generalized_dsem <- function(adj.file, lags.file, data.file) {
  
  # 1. Ingest Matrices
  adj <- as.matrix(read.csv(adj.file, row.names = 1))
  lags <- as.matrix(read.csv(lags.file, row.names = 1))
  
  if (!all(dim(adj) == dim(lags))) {
    stop("Dimensions of adjacency and lag matrices must match exactly.")
  }
  
  # 2. Ingest Simulated Data
  if (grepl("\\.csv$", data.file, ignore.case = TRUE)) {
    ts_data <- read.csv(data.file)
  } else if (grepl("\\.rds$", data.file, ignore.case = TRUE)) {
    ts_data <- readRDS(data.file)
  } else {
    stop("data.file must be a .csv or .rds file.")
  }
  
  n_steps <- nrow(ts_data)
  
  # Add a Time covariate to control for the baseline drift generated during simulation
  ts_data$Time <- 1:n_steps
  
  # Identify numeric columns (excluding the Time covariate we added)
  vars_to_scale <- setdiff(colnames(ts_data), "Time")
  
  # Apply standard scaling (mean = 0, sd = 1) to those columns to prevent Hessian errors
  ts_data[vars_to_scale] <- scale(ts_data[vars_to_scale])
  
  # Convert data to a formal time-series object required by dsem
  ts_data_matrix <- stats::ts(as.matrix(ts_data))
  
  # 3. Build dsem Syntax and Parameter Mapping
  model_syntax <- c()
  
  # Create a tracking data frame to easily map the raw dsem output back to readable relationships
  path_map <- data.frame(Source = character(), Target = character(), 
                         Lag = numeric(), Param = character(), stringsAsFactors = FALSE)
  
  for (target in colnames(adj)) {
    
    # Identify parents (causes) based on the binary test structure
    parents <- rownames(adj)[adj[, target] == 1]
    
    # Skip if exogenous (no incoming edges)
    if (length(parents) == 0) next
    
    # Formulate the Time trend control path (Time -> Target, lag = 0)
    time_param <- paste0("drift_", target)
    model_syntax <- c(model_syntax, paste("Time ->", target, ", 0 ,", time_param))
    path_map <- rbind(path_map, data.frame(Source = "Time", Target = target, Lag = 0, Param = time_param))
    
    for (p in parents) {
      lag_val <- lags[p, target]
      
      # Name the parameter uniquely so we can extract it cleanly later
      param_name <- paste0("path_", p, "_to_", target, "_lag", lag_val)
      
      # dsem syntax: source -> target, lag, parameter_name
      path_eq <- paste(p, "->", target, ",", lag_val, ",", param_name)
      
      model_syntax <- c(model_syntax, path_eq)
      path_map <- rbind(path_map, data.frame(Source = p, Target = target, Lag = lag_val, Param = param_name))
    }
  }
  path_map$parameter.id = 1:nrow(path_map)
  
  # Collapse the syntax vector into a single string separated by newlines
  final_syntax <- paste(model_syntax, collapse = "\n")
  
  # 4. Fit the SEM using dsem
  # We suppress output during fitting to keep the console clean in batch tests
  fit <- dsem::dsem(sem = final_syntax, 
                    tsdata = ts_data_matrix,
                    control = dsem::dsem_control(quiet = TRUE))
  
  # 5. Extract Parameter Estimates
  est_raw <- as.data.frame(summary(fit))
  est_raw = dplyr::rename(est_raw,parameter.id = 'parameter')
  
  # Merge the raw TMB estimates with our path map to filter out variances/intercepts
  estimates_df <- merge(path_map, est_raw, by = "parameter.id", all.x = TRUE)
  
  # Rename standard TMB output columns to match the previous reporting structure
  colnames(estimates_df)[colnames(estimates_df) == "Std. Error"] <- "Std_Error"
  colnames(estimates_df)[colnames(estimates_df) == "z value"] <- "z_Value"
  colnames(estimates_df)[colnames(estimates_df) == "Pr(>|z|)"] <- "p_Value"
  
  # Order and clean columns
  estimates_df <- estimates_df[, c("Target", "Source", "Lag", "Estimate", "Std_Error", "z_value", "p_value")]
  
  # Sort so Time drift components appear first for each target, followed by structural edges
  estimates_df <- estimates_df[order(estimates_df$Target, estimates_df$Source != "Time"), ]
  rownames(estimates_df) <- NULL
  
  return(list(
    syntax = final_syntax,
    fit = fit,
    estimates = estimates_df
  ))
}
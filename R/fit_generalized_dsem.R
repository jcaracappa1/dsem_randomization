#' Perform Generalized DSEM Analysis via the `dsem` Package
#'
#' Ingests a binary directed acyclic graph (DAG) and temporal lags, translates them
#' into the "arrow-and-lag" syntax required by the `dsem` package, and fits the 
#' dynamic structural equation model. Includes native handling of "Latent Complexes" 
#' (measurement models) by automatically fixing anchor parameters for identifiability.
#'
#' @param adj.file Character. Path to the CSV file containing the binary adjacency matrix.
#'        1 indicates an assumed causal edge, 0 indicates no edge.
#' @param lags.file Character. Path to the CSV file containing the temporal lag matrix.
#' @param data.file Character. Path to the simulated time series data. Supports `.csv` or `.rds`.
#' @param latent_dict List. Optional named list mapping latent variables to their 
#'        indicators. e.g., list("Phytoplankton" = c("diatom", "dino")).
#' @param detrend Logical. If TRUE, removes linear trends from all observed variables via `lm()` 
#'        prior to SEM fitting (the "Pre-Detrending" gold standard). Defaults to FALSE.
#' @param standardize Logical. If TRUE, applies `scale()` (Mean = 0, Variance = 1) to 
#'        observed variables to prevent computationally singular Hessian matrices. Defaults to TRUE.
#' @param control_time_drift Logical. If TRUE, injects a Time covariate as a direct structural cause 
#'        to partial out drift. Defaults to FALSE.
#'
#' @return A list containing:
#'         \itemize{
#'           \item \code{syntax}: Character string of the generated `dsem` model syntax.
#'           \item \code{fit}: The fitted `dsem` object.
#'           \item \code{estimates}: A clean data frame of the estimated structural path coefficients.
#'         }
#' @importFrom dsem dsem dsem_control
#' @importFrom stats ts lm resid na.exclude
#' @export

fit_generalized_dsem <- function(adj.file, lags.file, data.file, latent_dict = NULL, 
                                 detrend = FALSE, standardize = TRUE, control_time_drift = FALSE) {
  
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
  
  # Always generate Time as a covariate for detrending or drift control
  ts_data$Time <- 1:n_steps
  
  # CRITICAL: Identify observed variables vs pure latent variables (which are 100% NA)
  obs_vars <- colnames(ts_data)[sapply(ts_data, function(x) !all(is.na(x)))]
  obs_vars_no_time <- setdiff(obs_vars, "Time")
  
  # --- PRE-DETRENDING ---
  if (detrend) {
    for (v in obs_vars_no_time) {
      if (var(ts_data[[v]], na.rm = TRUE) > 0) {
        fit_lm <- stats::lm(ts_data[[v]] ~ ts_data$Time, na.action = stats::na.exclude)
        ts_data[[v]] <- stats::resid(fit_lm)
      }
    }
  }
  
  # --- STANDARDIZATION / Z-SCORING ---
  if (standardize) {
    vars_to_scale <- obs_vars_no_time
    if (control_time_drift) vars_to_scale <- c(vars_to_scale, "Time")
    
    for (v in vars_to_scale) {
      if (var(ts_data[[v]], na.rm = TRUE) > 0) {
        ts_data[[v]] <- scale(ts_data[[v]])[, 1]
      }
    }
  }
  
  # Convert data to a formal time-series object required by dsem
  ts_data_matrix <- stats::ts(as.matrix(ts_data))
  
  # 3. Build dsem Syntax and Parameter Mapping
  model_syntax <- c()
  
  path_map <- data.frame(Source = character(), Target = character(), 
                         Lag = numeric(), Param = character(), stringsAsFactors = FALSE)
  
  for (target in colnames(adj)) {
    
    parents <- rownames(adj)[adj[, target] == 1]
    if (length(parents) == 0) next
    
    is_indicator <- FALSE
    if (!is.null(latent_dict)) {
      if (target %in% unname(unlist(latent_dict))) is_indicator <- TRUE
    }
    
    # Formulate the Time trend control path ONLY if requested
    if (control_time_drift && !is_indicator) {
      time_param <- paste0("drift_", target)
      model_syntax <- c(model_syntax, paste("Time ->", target, ", 0 ,", time_param))
      path_map <- rbind(path_map, data.frame(Source = "Time", Target = target, Lag = 0, Param = time_param))
    }
    
    for (p in parents) {
      lag_val <- lags[p, target]
      param_name <- paste0("path_", p, "_to_", target, "_lag", lag_val)
      
      # IDENTIFIABILITY FIX 1: Handle Latent Measurement Model Anchor
      if (!is.null(latent_dict) && p %in% names(latent_dict) && target %in% latent_dict[[p]]) {
        if (target == latent_dict[[p]][1]) {
          # Fix the FIRST indicator's path mathematically to 1 to set Latent variable scale
          param_name <- "1" 
        } else {
          param_name <- paste0("loading_", p, "_to_", target)
        }
      }
      
      path_eq <- paste(p, "->", target, ",", lag_val, ",", param_name)
      model_syntax <- c(model_syntax, path_eq)
      path_map <- rbind(path_map, data.frame(Source = p, Target = target, Lag = lag_val, Param = param_name))
    }
    
    # IDENTIFIABILITY FIX 2: 1- and 2-Indicator Variance Traps
    if (!is.null(latent_dict)) {
      for (l_node in names(latent_dict)) {
        inds <- latent_dict[[l_node]]
        
        # If 1 or 2 indicators, the math is technically underidentified.
        # Safest generalized fallback: lock measurement errors to a small constant (0.01)
        # to force convergence without causing empirical loading contradictions.
        if (length(inds) <= 2 && target %in% inds) {
          model_syntax <- c(model_syntax, paste(target, "<->", target, ", 0 , 0.01"))
        }
      }
    }
  }
  
  # Safely assign parameter IDs tracking only estimated variables (skipping fixed '1' values)
  path_map$parameter.id <- NA
  est_counter <- 1
  for (i in seq_len(nrow(path_map))) {
    if (path_map$Param[i] != "1") {
      path_map$parameter.id[i] <- est_counter
      est_counter <- est_counter + 1
    }
  }
  
  final_syntax <- paste(model_syntax, collapse = "\n")
  
  # 4. Fit the SEM using dsem
  fit <- dsem::dsem(sem = final_syntax, 
                    tsdata = ts_data_matrix,
                    estimate_mu = obs_vars,
                    control = dsem::dsem_control(quiet = TRUE, newton_loops = 2))
  
  # 5. Extract Parameter Estimates
  est_raw <- as.data.frame(summary(fit))
  
  if ("parameter" %in% colnames(est_raw)) {
    colnames(est_raw)[colnames(est_raw) == "parameter"] <- "parameter.id"
  } else {
    est_raw$parameter.id <- 1:nrow(est_raw) 
  }
  
  estimates_df <- merge(path_map, est_raw, by = "parameter.id", all.x = TRUE)
  
  # Cleanly inject data for any mathematically fixed paths
  fixed_idx <- which(estimates_df$Param == "1")
  if (length(fixed_idx) > 0) {
    estimates_df$Estimate[fixed_idx] <- 1
    if ("Std. Error" %in% colnames(estimates_df)) estimates_df$"Std. Error"[fixed_idx] <- 0
  }
  
  colnames(estimates_df)[colnames(estimates_df) == "Std. Error"] <- "Std_Error"
  colnames(estimates_df)[tolower(colnames(estimates_df)) == "z value"] <- "z_value"
  colnames(estimates_df)[tolower(colnames(estimates_df)) == "pr(>|z|)"] <- "p_value"
  
  if (!"z_value" %in% colnames(estimates_df)) estimates_df$z_value <- NA
  if (!"p_value" %in% colnames(estimates_df)) estimates_df$p_value <- NA
  
  estimates_df <- estimates_df[, c("Target", "Source", "Lag", "Estimate", "Std_Error", "z_value", "p_value")]
  
  if (control_time_drift) {
    estimates_df <- estimates_df[order(estimates_df$Target, estimates_df$Source != "Time"), ]
  } else {
    estimates_df <- estimates_df[order(estimates_df$Target, estimates_df$Source), ]
  }
  
  rownames(estimates_df) <- NULL
  
  return(list(
    syntax = final_syntax,
    fit = fit,
    estimates = estimates_df
  ))
}
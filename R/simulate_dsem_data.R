#' Generate Time Series Data from Weighted and Lagged Adjacency Matrices
#'
#' Reads two CSV files containing named matrices for causal weights and temporal lags. 
#' Organizes the instantaneous effects into a DAG, and generates simulated time-series 
#' data respecting both immediate and lagged relationships. Causal effects are 
#' propagated as standardized anomalies internally to prevent variance cascading, 
#' but the final output retains its raw, realistic ecological scales and trends.
#'
#' @param weights.file Character. Path to the CSV file containing the weighted adjacency matrix. 
#' @param lags.file Character. Path to the CSV file containing the lag matrix. 
#' @param n_steps Integer. The number of time steps (rows) to generate for the time series.
#' @param var.mu Numeric vector. A named vector of baseline means (intercepts).
#' @param var.sd Numeric vector. A named vector of standard deviations for the noise term.
#' @param var.slope Numeric vector. A named vector of time-trend slopes.
#' @param latent_vars Character vector. Optional. Names of nodes to treat as latent. 
#'        Their columns will be set to NA in the output data to force state-space estimation.
#' @param diagnostics Logical. If TRUE, plots the DAG, creates a faceted plot of the data, 
#'        and prints a diagnostic table of expected weights vs. empirical correlations.
#'
#' @return If `diagnostics = FALSE`, returns a data frame of the generated time series. 
#'         If `diagnostics = TRUE`, returns a list containing `$data` (the time series with 
#'         latents masked), `$true_states` (the Oracle dataset with true latent values), 
#'         and `$diagnostic_summary`.
#' @importFrom igraph graph_from_adjacency_matrix topo_sort is_dag E plot.igraph get.edgelist
#' @importFrom ggplot2 ggplot aes geom_line facet_wrap theme_minimal labs
#' @importFrom tidyr pivot_longer
#' @importFrom stats rnorm cor
#' @export

simulate_dsem_data <- function(weights.file, lags.file, n_steps, var.mu, var.sd, var.slope, latent_vars = NULL, diagnostics = FALSE) {
  
  # 1. Ingest and Validate
  adj_matrix <- as.matrix(read.csv(weights.file, row.names = 1))
  lag_matrix <- as.matrix(read.csv(lags.file, row.names = 1))
  
  # Guard against NAs caused by empty cells in CSV files
  adj_matrix[is.na(adj_matrix)] <- 0
  lag_matrix[is.na(lag_matrix)] <- 0
  
  node_names <- colnames(adj_matrix)
  
  # Safeguard scalar arguments: expand single numbers to fully named vectors
  if (length(var.mu) == 1 && is.null(names(var.mu))) var.mu <- setNames(rep(var.mu, length(node_names)), node_names)
  if (length(var.sd) == 1 && is.null(names(var.sd))) var.sd <- setNames(rep(var.sd, length(node_names)), node_names)
  if (length(var.slope) == 1 && is.null(names(var.slope))) var.slope <- setNames(rep(var.slope, length(node_names)), node_names)
  
  inst_adj_matrix <- adj_matrix
  inst_adj_matrix[lag_matrix != 0] <- 0 
  g_inst <- igraph::graph_from_adjacency_matrix(inst_adj_matrix, mode = "directed", weighted = TRUE)
  
  if (!igraph::is_dag(g_inst)) {
    stop("The INSTANTANEOUS relationships (lag = 0) contain cycles.")
  }
  
  topological_order <- names(igraph::topo_sort(g_inst))
  
  # Initialize the output matrix
  sim_data <- matrix(0, nrow = n_steps, ncol = length(node_names))
  colnames(sim_data) <- node_names
  
  # 2. Iterative Time-Step Simulation Loop
  for (t in 1:n_steps) {
    for (node in topological_order) {
      
      # Generate unique native baseline (intercept + slope + residual noise)
      baseline <- var.mu[node] + (var.slope[node] * t) + rnorm(1, mean = 0, sd = var.sd[node])
      
      causal_effect <- 0
      parent_nodes <- rownames(adj_matrix)[adj_matrix[, node] != 0]
      
      for (parent in parent_nodes) {
        weight <- adj_matrix[parent, node]
        lag <- lag_matrix[parent, node]
        
        t_parent <- t - lag
        
        # If history doesn't exist yet (t <= lag), skip to assume 0 causal input
        if (t_parent <= 0) next
        
        parent_val <- sim_data[t_parent, parent]
        
        # Internally calculate anomaly to prevent massive scales from cascading downstream
        parent_expected_baseline <- var.mu[parent] + (var.slope[parent] * t_parent)
        parent_anomaly <- parent_val - parent_expected_baseline
        
        p_sd <- if (is.na(var.sd[parent]) || var.sd[parent] <= 0) 1e-6 else var.sd[parent]
        c_sd <- if (is.na(var.sd[node]) || var.sd[node] <= 0) 1e-6 else var.sd[node]
        
        parent_z <- parent_anomaly / p_sd
        
        # Multiply weight as a standardized correlation, scaled to child's native variance
        scaled_effect <- weight * parent_z * c_sd
        causal_effect <- causal_effect + scaled_effect
      }
      
      sim_data[t, node] <- baseline + causal_effect
    }
  }
  
  # 3. Handle Latent Variable Masking (Retaining Raw Scales)
  df_raw <- as.data.frame(sim_data) 
  df_true <- df_raw # "Oracle" dataset retains its raw, realistic scales
  df_sim <- df_true
  
  if (!is.null(latent_vars)) {
    for (lv in latent_vars) {
      if (lv %in% colnames(df_sim)) {
        df_sim[, lv] <- NA
      } else {
        warning(sprintf("Latent variable '%s' not found in the adjacency matrix.", lv))
      }
    }
  }
  
  # 4. Diagnostics Suite
  if (diagnostics) {
    message("--- Running DSEM Diagnostics ---")
    
    g_full <- igraph::graph_from_adjacency_matrix(adj_matrix != 0, mode = "directed")
    edge_list <- igraph::get.edgelist(g_full)
    
    edge_labels <- apply(edge_list, 1, function(e) {
      w <- adj_matrix[e[1], e[2]]
      l <- lag_matrix[e[1], e[2]]
      paste0("W: ", w, "\nL: ", l)
    })
    
    igraph::E(g_full)$label <- edge_labels
    graphics::plot(g_full, main = "DSEM Causal DAG\n(W = Weight, L = Lag)", 
                   edge.label.cex = 0.8, vertex.color = "lightblue", vertex.size = 30)
    
    df_long <- df_raw |> 
      dplyr::mutate(Time = 1:n_steps) |> 
      tidyr::pivot_longer(cols =-Time, names_to = "Variable", values_to = "Value")
    p <- ggplot2::ggplot(df_long, ggplot2::aes(x = Time, y = Value, color = Variable)) + 
      ggplot2::geom_line() + 
      ggplot2::facet_wrap(~Variable, scales = "free_y", ncol = 1) +
      ggplot2::theme_minimal() + 
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(title = "Simulated Time Series by Node (Raw Scales)")
    print(p)
    
    diag_results <- data.frame(Parent = character(), Child = character(), 
                               Input_Weight = numeric(), Input_Lag = numeric(), 
                               Empirical_Cor = numeric(), stringsAsFactors = FALSE)
    
    for (i in 1:nrow(edge_list)) {
      p_node <- edge_list[i, 1]
      c_node <- edge_list[i, 2]
      w <- adj_matrix[p_node, c_node]
      l <- lag_matrix[p_node, c_node]
      
      if (l == 0) {
        p_vec <- df_raw[, p_node]
        c_vec <- df_raw[, c_node]
      } else {
        p_vec <- df_raw[1:(n_steps - l), p_node]
        c_vec <- df_raw[(l + 1):n_steps, c_node]
      }
      
      emp_cor <- suppressWarnings(stats::cor(p_vec, c_vec))
      
      diag_results <- rbind(diag_results, data.frame(
        Parent = p_node, Child = c_node, 
        Input_Weight = w, Input_Lag = l, Empirical_Cor = round(emp_cor, 3)
      ))
    }
    
    print(diag_results)
    
    return(list(data = df_sim, true_states = df_true, diagnostic_summary = diag_results))
  }
  
  return(df_sim)
}
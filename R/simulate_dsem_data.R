#' Generate Time Series Data from Weighted and Lagged Adjacency Matrices
#'
#' Reads two CSV files containing named matrices for causal weights and temporal lags. 
#' Organizes the instantaneous effects into a DAG, and generates simulated time-series 
#' data respecting both immediate and lagged relationships. Includes an optional 
#' diagnostics suite to visualize the DAG, plot the series, and report empirical correlations.
#'
#' @param weights.file Character. Path to the CSV file containing the weighted adjacency matrix. 
#' @param lags.file Character. Path to the CSV file containing the lag matrix. 
#' @param n_steps Integer. The number of time steps (rows) to generate for the time series.
#' @param var.mu Numeric vector. A named vector of baseline means (intercepts).
#' @param var.sd Numeric vector. A named vector of standard deviations for the noise term.
#' @param var.slope Numeric vector. A named vector of time-trend slopes.
#' @param diagnostics Logical. If TRUE, plots the DAG, creates a faceted plot of the data, 
#'        and prints a diagnostic table of expected weights vs. empirical correlations.
#'
#' @return If `diagnostics = FALSE`, returns a data frame of the generated time series. 
#'         If `diagnostics = TRUE`, returns a list containing `$data` (the time series) 
#'         and `$diagnostic_summary` (a data frame of edge checks).
#' @importFrom igraph graph_from_adjacency_matrix topo_sort is_dag E plot.igraph get.edgelist
#' @importFrom ggplot2 ggplot aes geom_line facet_wrap theme_minimal labs
#' @importFrom tidyr pivot_longer
#' @importFrom stats rnorm cor
#' @export
#'
#' @examples
#' \dontrun{
#' sim_results <- simulate_dsem_data("weights.csv", "lags.csv", n_steps = 100,
#'                                   var.mu = mu, var.sd = sd, var.slope = slopes,
#'                                   diagnostics = TRUE)
#' }


simulate_dsem_data <- function(weights.file, lags.file, n_steps, var.mu, var.sd, var.slope, diagnostics = FALSE) {
  
  # 1. Ingest and Validate
  adj_matrix <- as.matrix(read.csv(weights.file, row.names = 1))
  lag_matrix <- as.matrix(read.csv(lags.file, row.names = 1))
  node_names <- colnames(adj_matrix)
  
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
      
      # Generate unique native baseline for this specific time step
      baseline <- var.mu[node] + (var.slope[node] * t) + rnorm(1, mean = 0, sd = var.sd[node])
      
      causal_effect <- 0
      parent_nodes <- names(adj_matrix[, node][adj_matrix[, node] != 0])
      
      for (parent in parent_nodes) {
        weight <- adj_matrix[parent, node]
        lag <- lag_matrix[parent, node]
        
        t_parent <- t - lag
        
        # If history doesn't exist yet (t <= lag), skip to assume 0 causal input
        if (t_parent <= 0) next
        
        parent_val <- sim_data[t_parent, parent]
        
        # Standardize the parent value based on its expected baseline at that historical time
        parent_baseline_at_t <- var.mu[parent] + (var.slope[parent] * t_parent)
        parent_standardized <- (parent_val - parent_baseline_at_t) / var.sd[parent]
        
        # Scale by the child's standard deviation so weight acts as a correlation
        scaled_effect <- weight * parent_standardized * var.sd[node]
        causal_effect <- causal_effect + scaled_effect
      }
      
      sim_data[t, node] <- baseline + causal_effect
    }
  }
  
  df_sim <- as.data.frame(sim_data)
  
  # 4. Diagnostics Suite
  if (diagnostics) {
    message("--- Running DSEM Diagnostics ---")
    
    # Diagnostic 1: Visualize Full DAG with Lags
    g_full <- igraph::graph_from_adjacency_matrix(adj_matrix != 0, mode = "directed")
    edge_list <- igraph::get.edgelist(g_full)
    
    # Create labels showing Weight and Lag for each edge
    edge_labels <- apply(edge_list, 1, function(e) {
      w <- adj_matrix[e[1], e[2]]
      l <- lag_matrix[e[1], e[2]]
      paste0("W: ", w, "\nL: ", l)
    })
    
    igraph::E(g_full)$label <- edge_labels
    graphics::plot(g_full, main = "DSEM Causal DAG\n(W = Weight, L = Lag)", 
                   edge.label.cex = 0.8, vertex.color = "lightblue", vertex.size = 30)
    
    # Diagnostic 2: Faceted Plot of Time Series
    df_long <- df_sim |> 
      dplyr::mutate(Time = 1:n_steps) |> 
      tidyr::pivot_longer(cols =-Time,
                          names_to = "Variable",
                          values_to = "Value")
    p <- ggplot2::ggplot(df_long, ggplot2::aes(x = Time, y = Value, color = Variable)) + 
      ggplot2::geom_line() + 
      ggplot2::facet_wrap(~Variable, scales = "free_y", ncol = 1) +
      ggplot2::theme_minimal() + 
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(title = "Simulated Time Series by Node")
    print(p)
    
    # Diagnostic 3: Check correlations
    diag_results <- data.frame(Parent = character(), Child = character(), 
                               Input_Weight = numeric(), Input_Lag = numeric(), 
                               Empirical_Cor = numeric(), stringsAsFactors = FALSE)
    
    for (i in 1:nrow(edge_list)) {
      p_node <- edge_list[i, 1]
      c_node <- edge_list[i, 2]
      w <- adj_matrix[p_node, c_node]
      l <- lag_matrix[p_node, c_node]
      
      # Shift arrays to align the causal lag before calculating correlation
      if (l == 0) {
        p_vec <- df_sim[, p_node]
        c_vec <- df_sim[, c_node]
      } else {
        p_vec <- df_sim[1:(n_steps - l), p_node]
        c_vec <- df_sim[(l + 1):n_steps, c_node]
      }
      
      emp_cor <- stats::cor(p_vec, c_vec)
      
      diag_results <- rbind(diag_results, data.frame(
        Parent = p_node, Child = c_node, 
        Input_Weight = w, Input_Lag = l, Empirical_Cor = round(emp_cor, 3)
      ))
    }
    
    print(diag_results)
    
    return(list(data = df_sim, diagnostic_summary = diag_results))
  }
  
  return(df_sim)
}
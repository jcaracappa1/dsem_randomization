#' Generate Randomized Directed Graph Preserving Source/Sink Constraints
#' @param ref_adj Matrix. The reference directed adjacency matrix (weights or binary).
#'        Rows are sources, columns are targets.
#' @param lags_file Character. Path to the CSV file containing the reference lag matrix.
#'
#' @return A list containing:
#'         \itemize{
#'           \item \code{adj}: A binary adjacency matrix representing the randomized DAG structure.
#'           \item \code{lags}: A matched temporal lag matrix ensuring mathematical solvability.
#'         }
#' @export
randomize_dsem_graph <- function(ref_adj, lags_file) {
  
  # Helper to prevent R's notorious sample(numeric, 1) sequence-expansion bug
  safe_sample <- function(x, size = 1) {
    if (length(x) <= 1) return(x)
    sample(x, size, replace = FALSE)
  }
  
  # 1. Identify Node Classifications
  node_names <- colnames(ref_adj)
  if (is.null(node_names)) stop("Reference matrix must have column/row names.")
  
  # Extract lag proportions and max lag from the reference file
  ref_lags <- as.matrix(read.csv(lags_file, row.names = 1))
  if (!all(dim(ref_adj) == dim(ref_lags))) stop("Adjacency and Lag matrix dimensions must match.")
  
  ref_edges <- ref_adj != 0
  ref_lag_vals <- ref_lags[ref_edges]
  
  max_lag <- max(c(0, ref_lag_vals), na.rm = TRUE)
  if (max_lag == 0) max_lag <- 1 # Ensure we have at least Lag 1 available to break loops
  
  prob_lag0 <- sum(ref_lag_vals == 0) / length(ref_lag_vals)
  if (is.nan(prob_lag0)) prob_lag0 <- 1
  
  col_s <- colSums(ref_adj != 0)
  row_s <- rowSums(ref_adj != 0)
  
  n_nodes <- node_names[col_s == 0 & row_s > 0]
  m_nodes <- node_names[row_s == 0 & col_s > 0]
  k_nodes <- setdiff(node_names, c(n_nodes, m_nodes))
  
  E <- sum(ref_adj != 0)
  if (E == 0) return(matrix(0, length(node_names), length(node_names), dimnames = list(node_names, node_names)))
  
  # 2. Assign Ranks
  ranks <- setNames(numeric(length(node_names)), node_names)
  ranks[n_nodes] <- 0
  
  if (length(k_nodes) > 1) {
    ranks[k_nodes] <- sample(1:(length(k_nodes) - 1), length(k_nodes), replace = TRUE)
  } else if (length(k_nodes) == 1) {
    ranks[k_nodes] <- 1
  }
  
  max_k_rank <- if (length(k_nodes) > 0) max(ranks[k_nodes]) else 0
  ranks[m_nodes] <- max_k_rank + 1
  
  # 3. Determine Active 'k' Nodes mathematically capable of holding E edges
  valid_S <- integer(0)
  N_n <- length(n_nodes)
  N_m <- length(m_nodes)
  
  for (S in 0:length(k_nodes)) {
    if (S <= E - 1) {
      # Max capacity formula = (sources to active_k/m) + (active_k to active_k/m)
      max_capacity <- (N_n * (S + N_m)) + (S * (max(0, S - 1) + N_m))
      if (max_capacity >= E) {
        valid_S <- c(valid_S, S)
      }
    }
  }
  
  if (length(valid_S) == 0) {
    stop("The reference graph has more edges than can legally be placed under the n, m, k, and no-self-loop constraints.")
  }
  
  # Use safe_sample to prevent picking a smaller subset than mathematically required
  target_S <- safe_sample(valid_S, 1)
  K_act <- if (target_S > 0) safe_sample(k_nodes, target_S) else character(0)
  
  adj <- matrix(0, nrow = length(node_names), ncol = length(node_names), 
                dimnames = list(node_names, node_names))
  
  # 4. Build the Hierarchical Backbone
  for (u in K_act) {
    valid_targets <- c(K_act, m_nodes)
    valid_targets <- valid_targets[ranks[valid_targets] > ranks[u]]
    
    if (length(valid_targets) == 0) valid_targets <- m_nodes
    
    v <- safe_sample(valid_targets, 1)
    adj[u, v] <- 1
  }
  
  s <- safe_sample(n_nodes, 1)
  valid_n_targets <- c(K_act, m_nodes)
  v <- safe_sample(valid_n_targets, 1)
  adj[s, v] <- 1
  
  # 5. Fill Remaining Edges Randomly
  edges_needed <- E - sum(adj)
  
  if (edges_needed > 0) {
    allowed_nodes <- c(n_nodes, K_act, m_nodes)
    grid <- expand.grid(source = allowed_nodes, target = allowed_nodes, stringsAsFactors = FALSE)
    
    # Apply Constraints
    grid <- grid[grid$source != grid$target, ]      
    grid <- grid[!(grid$target %in% n_nodes), ]     
    grid <- grid[!(grid$source %in% m_nodes), ]     
    
    # Remove edges already in the backbone
    already_placed <- apply(grid, 1, function(row) adj[row[1], row[2]] == 1)
    grid <- grid[!already_placed, , drop = FALSE]
    
    # seq_len safely iterates rows without triggering sample expansion errors
    chosen_indices <- safe_sample(seq_len(nrow(grid)), edges_needed)
    
    for (idx in chosen_indices) {
      adj[grid$source[idx], grid$target[idx]] <- 1
    }
  }
  
  # 6. Apply Hybrid Lag Strategy
  rand_lags <- matrix(0, nrow = length(node_names), ncol = length(node_names), 
                      dimnames = list(node_names, node_names))
  
  # Identify cycles in the raw Lag 0 skeleton
  g_adj <- igraph::graph_from_adjacency_matrix(adj, mode = "directed")
  fas <- igraph::feedback_arc_set(g_adj)
  fas_edges <- igraph::as_edgelist(g_adj)[fas, , drop = FALSE]
  
  # Distribute lags
  edge_ind <- which(adj != 0, arr.ind = TRUE)
  for (i in seq_len(nrow(edge_ind))) {
    u <- edge_ind[i, 1]
    v <- edge_ind[i, 2]
    u_name <- rownames(adj)[u]
    v_name <- colnames(adj)[v]
    
    is_fas <- FALSE
    if (nrow(fas_edges) > 0) {
      is_fas <- any(fas_edges[, 1] == u_name & fas_edges[, 2] == v_name)
    }
    
    if (is_fas) {
      # Force temporal delay to unroll the cycle and prevent DSEM crashing
      rand_lags[u, v] <- safe_sample(1:max_lag, 1)
    } else {
      # Assign based on historical probability
      if (runif(1) <= prob_lag0) {
        rand_lags[u, v] <- 0
      } else {
        rand_lags[u, v] <- safe_sample(1:max_lag, 1)
      }
    }
  }
  
  return(list(adj = adj, lags = rand_lags))
}
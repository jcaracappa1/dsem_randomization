generate_network <- function(n, m, k, p_random_k = 0.3, p_n_out = 0.3) {
  N <- n + m + k
  
  sources <- if(n > 0) 1:n else integer(0)
  k_nodes <- if(k > 0) (n + 1):(n + k) else integer(0)
  sinks   <- if(m > 0) (n + k + 1):N else integer(0)
  
  A <- matrix(0, nrow = N, ncol = N)
  node_levels <- numeric(N)
  
  # 1. Assign Layout Ranks/Levels
  if (m > 0) node_levels[sinks] <- 0
  if (n > 0) node_levels[sources] <- k + 1
  if (k > 0) {
    k_ranks <- sample(1:k, size = k, replace = TRUE)
    node_levels[k_nodes] <- k_ranks
  }
  
  # 2. 'n' out-edges (Can have any or no connections)
  if (n > 0) {
    valid_n_targets <- c(k_nodes, sinks)
    if (length(valid_n_targets) > 0) {
      for (s in sources) {
        for (t in valid_n_targets) {
          if (runif(1) < p_n_out) A[s, t] <- 1
        }
      }
    }
  }
  
  # 3. Ensure all 'm' have an in_degree of at least 1
  if (m > 0) {
    valid_m_sources <- c(sources, k_nodes)
    if (length(valid_m_sources) > 0) {
      for (sink in sinks) {
        if (sum(A[valid_m_sources, sink]) == 0) {
          src <- if (length(valid_m_sources) == 1) valid_m_sources else sample(valid_m_sources, 1)
          A[src, sink] <- 1
        }
      }
    }
  }
  
  # 4. Construct 'k' backbone using the new ranks
  if (k > 0 && m > 0) {
    for (i in seq_along(k_nodes)) {
      u <- k_nodes[i]
      r <- k_ranks[i]
      
      # Backbone must step DOWN to a strictly lower rank, or directly to 'm'
      valid_lower_k <- k_nodes[k_ranks < r]
      possible_targets <- c(valid_lower_k, sinks)
      
      target <- if(length(possible_targets) == 1) possible_targets else sample(possible_targets, 1)
      A[u, target] <- 1
    }
  }
  
  # 5. Add random edges among 'k' (Zero Reciprocity Enforced)
  if (k > 1) {
    for (u in k_nodes) {
      for (v in k_nodes) {
        # NEW RULE: Check that A[v, u] == 0 to ensure the reverse edge does not exist
        if (u != v && A[v, u] == 0 && runif(1) < p_random_k) {
          A[u, v] <- 1
        }
      }
    }
  }
  
  # Strictly no self-loops
  diag(A) <- 0
  
  # Names
  node_names <- c()
  if(n > 0) node_names <- c(node_names, paste0("n", 1:n))
  if(k > 0) node_names <- c(node_names, paste0("k", 1:k))
  if(m > 0) node_names <- c(node_names, paste0("m", 1:m))
  rownames(A) <- node_names
  colnames(A) <- node_names
  
  return(list(matrix = A, levels = node_levels))
}

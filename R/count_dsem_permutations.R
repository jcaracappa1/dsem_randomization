#' Calculate Total Possible DSEM Graph Permutations
#'
#' Evaluates the total number of valid directed graph permutations that can be 
#' generated from a reference graph while respecting specific constraints 
#' (preserved sources/sinks, active reachability, no self-loops, conserved edges).
#' Uses exact combinatorics for small search spaces and Monte Carlo estimation 
#' for large graph spaces.
#'
#' @param ref_adj Matrix. The reference directed adjacency matrix.
#' @param method Character. "auto" (default) chooses automatically based on the size 
#'        of the permutation space. "exact" forces brute-force combinatorics. 
#'        "estimate" forces Monte Carlo simulation.
#' @param max_iter Numeric. The maximum number of permutations to evaluate. If the 
#'        total possible combinations exceed this, the function defaults to estimation.
#' @param latent_dict List. Optional named list mapping latent variables to their indicators.
#'
#' @return A list containing the total unconstrained combinations, the estimated 
#'         (or exact) number of valid permutations, the constraint pass rate, 
#'         and the method used.
#' @importFrom igraph graph_from_adjacency_matrix distances
#' @importFrom utils combn
#' @export
count_dsem_permutations <- function(ref_adj, method = c("auto", "exact", "estimate"), max_iter = 100000, latent_dict = NULL) {
  
  method <- match.arg(method)
  
  # --- LATENT COMPLEX FILTER ---
  # If a measurement model exists, it is fixed. We only calculate permutations 
  # for the purely structural nodes to avoid combinatorial explosion.
  if (!is.null(latent_dict)) {
    indicators <- unname(unlist(latent_dict))
    struct_nodes <- setdiff(colnames(ref_adj), indicators)
    ref_adj <- ref_adj[struct_nodes, struct_nodes, drop = FALSE]
  }
  
  # 1. Extract Node Classifications
  ref_bin <- matrix(as.numeric(ref_adj != 0), nrow = nrow(ref_adj), dimnames = dimnames(ref_adj))
  node_names <- colnames(ref_bin)
  
  col_s <- colSums(ref_bin)
  row_s <- rowSums(ref_bin)
  
  ref_n <- sort(node_names[col_s == 0 & row_s > 0])
  ref_m <- sort(node_names[row_s == 0 & col_s > 0])
  ref_k <- sort(setdiff(node_names, c(ref_n, ref_m)))
  ref_E <- sum(ref_bin)
  
  if (ref_E == 0) return(list(Total_Unconstrained = 1, Valid_Permutations = 1, Method = "Exact", Pass_Rate = 1))
  
  # 2. Define the Unconstrained Edge Pool
  # We natively filter out edges that violate universal constraints (no self-loops, n->in, m->out)
  grid <- expand.grid(source = node_names, target = node_names, stringsAsFactors = FALSE)
  grid <- grid[grid$source != grid$target, ]
  grid <- grid[!(grid$target %in% ref_n), ]
  grid <- grid[!(grid$source %in% ref_m), ]
  
  E_pool <- nrow(grid)
  
  # Calculate absolute maximum combinations: "E_pool choose ref_E"
  total_combinations <- choose(E_pool, ref_E)
  
  if (method == "auto") {
    method <- if (total_combinations <= max_iter) "exact" else "estimate"
  }
  
  if (method == "exact" && total_combinations > 2e6) {
    warning("Running exact combinatorics on > 2 million permutations may take a significant amount of time.")
  }
  
  # 3. Constraint Checking Function
  check_constraints <- function(adj) {
    rand_col_s <- colSums(adj)
    rand_row_s <- rowSums(adj)
    
    # C1 & C2: Source and Sink definitions are respected
    if (!all(rand_col_s[ref_n] == 0) || !any(rand_row_s[ref_n] > 0)) return(FALSE)
    if (!all(rand_row_s[ref_m] == 0) || !any(rand_col_s[ref_m] > 0)) return(FALSE)
    
    # No imposters
    rand_n <- node_names[rand_col_s == 0 & rand_row_s > 0]
    rand_m <- node_names[rand_row_s == 0 & rand_col_s > 0]
    if (!all(rand_n %in% ref_n)) return(FALSE)
    if (!all(rand_m %in% ref_m)) return(FALSE)
    
    g_rand <- igraph::graph_from_adjacency_matrix(adj, mode = "directed")
    
    # C3: Active K reach M
    active_k <- ref_k[rand_row_s[ref_k] > 0 | rand_col_s[ref_k] > 0]
    if (length(active_k) > 0) {
      dists_k <- suppressWarnings(igraph::distances(g_rand, v = active_k, to = ref_m, mode = "out"))
      if (!all(apply(dists_k, 1, function(row) any(is.finite(row))))) return(FALSE)
    }
    
    # C4: At least one N reaches an M
    dists_n <- suppressWarnings(igraph::distances(g_rand, v = ref_n, to = ref_m, mode = "out"))
    if (!any(is.finite(dists_n))) return(FALSE)
    
    return(TRUE)
  }
  
  base_adj <- matrix(0, nrow = length(node_names), ncol = length(node_names), dimnames = list(node_names, node_names))
  
  # 4. Execute Combinatorics or Estimation
  if (method == "exact") {
    
    # Iterate through permutations and cleanly return logicals (TRUE/FALSE)
    res <- utils::combn(E_pool, ref_E, FUN = function(idx) {
      adj <- base_adj
      adj[cbind(grid$source[idx], grid$target[idx])] <- 1
      return(check_constraints(adj))
    }, simplify = FALSE)
    
    # Summing the TRUEs gives us the exact count without scoping issues
    pass_count <- sum(unlist(res))
    valid_perms <- pass_count
    pass_rate <- pass_count / total_combinations
    
  } else {
    
    # Monte Carlo Estimation
    pass_count <- 0
    for (i in seq_len(max_iter)) {
      idx <- sample(seq_len(E_pool), ref_E, replace = FALSE)
      adj <- base_adj
      adj[cbind(grid$source[idx], grid$target[idx])] <- 1
      if (check_constraints(adj)) pass_count <- pass_count + 1
    }
    
    pass_rate <- pass_count / max_iter
    valid_perms <- round(total_combinations * pass_rate)
  }
  
  return(list(
    Total_Unconstrained_Combinations = total_combinations,
    Constraint_Pass_Rate = round(pass_rate, 4),
    Valid_Permutations = valid_perms,
    Method = method
  ))
}
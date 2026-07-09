#' Generate and Evaluate N Randomized DSEM Graphs
#'
#' Runs the `randomize_dsem_graph` function N times to create a batch of 
#' randomized adjacency matrices along with their matched temporal lags. 
#' It computes a comprehensive suite of topological diagnostics to ensure 
#' the structural constraints were met, tests for graph isomorphism and 
#' Markov equivalence against the reference graph, and identifies identical 
#' graphs (both structure and lags) generated within the batch.
#'
#' @param ref_adj Matrix. The reference directed adjacency matrix.
#' @param ref_lags Matrix. The reference temporal lag matrix.
#' @param N Integer. The number of randomized graphs to generate and evaluate.
#' @param deduplicate Logical. If \code{TRUE}, filters the results to return only unique graph structures.
#' @param latent_dict List. Optional named list mapping latent variables to their indicators.
#'
#' @return A list containing:
#'         \itemize{
#'           \item \code{graphs}: A list of N binary adjacency matrices.
#'           \item \code{lags}: A list of N corresponding temporal lag matrices.
#'           \item \code{diagnostics}: A data frame with N rows containing boolean 
#'                 evaluations for constraints, isomorphism, Markov equivalence, 
#'                 and an Exact_Match_ID to group identical structures.
#'         }
#' @importFrom igraph graph_from_adjacency_matrix distances isomorphic is_dag
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate 100 randomized versions of a true graph
#' batch_results <- evaluate_random_graphs(true_adj, true_lags, N = 100)
#' 
#' # View the diagnostic results
#' head(batch_results$diagnostics)
#' 
#' # Check how many UNIQUE graphs were produced in the batch of 100
#' length(unique(batch_results$diagnostics$Exact_Match_ID))
#' }

evaluate_random_graphs <- function(ref_adj, ref_lags, N = 10, deduplicate = FALSE, latent_dict = NULL) {
  
  # Ensure input is binary for topological checking
  ref_bin <- matrix(as.numeric(ref_adj != 0), nrow = nrow(ref_adj), dimnames = dimnames(ref_adj))
  node_names <- colnames(ref_bin)
  
  # 1. Extract Reference Properties
  col_s <- colSums(ref_bin)
  row_s <- rowSums(ref_bin)
  
  ref_n <- sort(node_names[col_s == 0 & row_s > 0])
  ref_m <- sort(node_names[row_s == 0 & col_s > 0])
  ref_k <- sort(setdiff(node_names, c(ref_n, ref_m)))
  ref_E <- sum(ref_bin)
  
  g_ref <- igraph::graph_from_adjacency_matrix(ref_bin, mode = "directed")
  
  # Helper to find v-structures (unshielded colliders: A -> B <- C without A-C edge)
  get_v_structures <- function(adj) {
    v_structs <- c()
    for (node in colnames(adj)) {
      parents <- colnames(adj)[adj[, node] == 1]
      if (length(parents) >= 2) {
        combos <- combn(parents, 2, simplify = FALSE)
        for (pair in combos) {
          # Unshielded check
          if (adj[pair[1], pair[2]] == 0 && adj[pair[2], pair[1]] == 0) {
            v_structs <- c(v_structs, paste(sort(pair), "->", node, collapse = "_"))
          }
        }
      }
    }
    return(sort(v_structs))
  }
  
  ref_skel <- ref_bin | t(ref_bin) # Undirected skeleton
  ref_v_structs <- get_v_structures(ref_bin)
  
  # 2. Initialize Outputs
  graphs_list <- vector("list", N)
  lags_list <- vector("list", N)
  graph_hashes <- character(N) # Used to track identical matrix structures
  
  diag_df <- data.frame(
    Iter = 1:N,
    C1_Preserved_Sources = logical(N),
    C2_Preserved_Sinks = logical(N),
    C3_Active_K_Reach_M = logical(N),
    C4_Path_N_to_M = logical(N),
    C5_No_Self_Loops = logical(N),
    C6_Edges_Conserved = logical(N),
    Is_DAG = logical(N),
    Isomorphic = logical(N),
    Markov_Equivalent = logical(N),
    Pass_All_Constraints = logical(N)
  )
  
  # 3. Batch Generation and Evaluation Loop
  for (i in 1:N) {
    
    # Generate graph and matched lags, passing the latent dictionary to preserve the complex
    rand_res <- randomize_dsem_graph(ref_bin, ref_lags, latent_dict = latent_dict)
    rand_adj <- rand_res$adj
    rand_lags <- rand_res$lags
    
    graphs_list[[i]] <- rand_adj
    lags_list[[i]] <- rand_lags
    
    # Create a unique string hash of both the matrix and lags to identify exact duplicates later
    adj_str <- paste(as.vector(rand_adj), collapse = "")
    lag_str <- paste(as.vector(rand_lags), collapse = "")
    graph_hashes[i] <- paste(adj_str, lag_str, sep = "_")
    
    g_rand <- igraph::graph_from_adjacency_matrix(rand_adj, mode = "directed")
    
    # Extract structural summaries of the new randomized graph
    rand_col_s <- colSums(rand_adj)
    rand_row_s <- rowSums(rand_adj)
    rand_n <- sort(node_names[rand_col_s == 0 & rand_row_s > 0])
    rand_m <- sort(node_names[rand_row_s == 0 & rand_col_s > 0])
    
    # C1: Preserved Sources
    # Original sources must never receive edges, at least one must emit an edge
    diag_df$C1_Preserved_Sources[i] <- all(rand_col_s[ref_n] == 0) && 
      any(rand_row_s[ref_n] > 0)
    
    # C2: Preserved Sinks
    # Original sinks must never emit edges, at least one must receive an edge,
    # and no non-sink nodes can accidentally act as absolute sinks.
    diag_df$C2_Preserved_Sinks[i] <- all(rand_row_s[ref_m] == 0) && 
      any(rand_col_s[ref_m] > 0) && 
      all(rand_m %in% ref_m)
    
    # C3: Active Intermediates reach Sink
    # An active 'k' is any 'k' node with at least one in or out edge
    active_k <- ref_k[rand_row_s[ref_k] > 0 | rand_col_s[ref_k] > 0]
    if (length(active_k) > 0) {
      dists_k_m <- suppressWarnings(igraph::distances(g_rand, v = active_k, to = ref_m, mode = "out"))
      # Check if every active k has at least one finite distance to an m
      diag_df$C3_Active_K_Reach_M[i] <- all(apply(dists_k_m, 1, function(row) any(is.finite(row))))
    } else {
      diag_df$C3_Active_K_Reach_M[i] <- TRUE # Trivially true if no k are active
    }
    
    # C4: At least one Source reaches a Sink
    dists_n_m <- suppressWarnings(igraph::distances(g_rand, v = ref_n, to = ref_m, mode = "out"))
    diag_df$C4_Path_N_to_M[i] <- any(is.finite(dists_n_m))
    
    # C5: No Self Loops
    diag_df$C5_No_Self_Loops[i] <- (sum(diag(rand_adj)) == 0)
    
    # C6: Edge Count Conservation
    diag_df$C6_Edges_Conserved[i] <- (sum(rand_adj) == ref_E)
    
    # Advanced: DAG Check
    is_dag <- igraph::is_dag(g_rand)
    diag_df$Is_DAG[i] <- is_dag
    
    # Advanced: Isomorphism
    diag_df$Isomorphic[i] <- igraph::isomorphic(g_ref, g_rand)
    
    # Advanced: Markov Equivalence (Requires DAG, Identical Skeleton, Identical V-Structures)
    if (is_dag) {
      rand_skel <- rand_adj | t(rand_adj)
      rand_v_structs <- get_v_structures(rand_adj)
      
      skel_match <- identical(ref_skel, rand_skel)
      v_match <- identical(ref_v_structs, rand_v_structs)
      
      diag_df$Markov_Equivalent[i] <- skel_match && v_match
    } else {
      # Cyclic graphs cannot be evaluated for standard DAG Markov Equivalence
      diag_df$Markov_Equivalent[i] <- FALSE 
    }
    
    # Meta Constraint Boolean
    diag_df$Pass_All_Constraints[i] <- all(
      diag_df$C1_Preserved_Sources[i],
      diag_df$C2_Preserved_Sinks[i],
      diag_df$C3_Active_K_Reach_M[i],
      diag_df$C4_Path_N_to_M[i],
      diag_df$C5_No_Self_Loops[i],
      diag_df$C6_Edges_Conserved[i]
    )
  }
  
  # 4. Group Identical Graphs using Hashes
  # Converts the unique string hashes into clean integer IDs
  diag_df$Exact_Match_ID <- as.numeric(factor(graph_hashes, levels = unique(graph_hashes)))
  
  # Reorder columns to put Exact_Match_ID right next to Iteration
  cols <- c("Iter", "Exact_Match_ID", setdiff(colnames(diag_df), c("Iter", "Exact_Match_ID")))
  diag_df <- diag_df[, cols]
  
  # 5. Deduplicate Results
  if (deduplicate) {
    # Keep only the first instance of each exact match ID
    unique_idx <- !duplicated(diag_df$Exact_Match_ID)
    
    diag_df <- diag_df[unique_idx, ]
    graphs_list <- graphs_list[unique_idx]
    lags_list <- lags_list[unique_idx]
    
    # Reset row names for clean reporting
    rownames(diag_df) <- NULL
  }
  
  return(list(
    graphs = graphs_list,
    lags = lags_list,
    diagnostics = diag_df
  ))
}
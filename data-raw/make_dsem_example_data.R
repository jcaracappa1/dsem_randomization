# 1. Generate the random matrix
ref_adj =  read.csv(here::here('data-raw','example_weight_dag.csv'),row.names = 1) |> as.matrix()
ref_adj[which(ref_adj!=0)] =1

p <- igraph::graph_from_adjacency_matrix(ref_adj, mode = "directed")

# 3. Plot using a hierarchical layout
igraph::plot.igraph(p, 
                    layout = igraph::layout_with_sugiyama(p)$layout, 
                    vertex.color = "lightblue", 
                    vertex.size = 30, 
                    edge.arrow.size = 0.6,
                    main = "Reference Hierarchical DAG")

dsem.data = DSEMrand::simulate_dsem_data(weights.file = here::here('data-raw','example_weight_dag.csv'),
                   lags.file = here::here('data-raw','example_lags_dag.csv'),
                   n_steps = 100,
                   var.mu = c('bot.temp' = 5,'diatom' = 1E6,'juv.bio' = 1E5,'adult.bio' = 1E4,'catch' = 1E3),
                   var.sd = c('bot.temp' = 1,'diatom' = 1E3,'juv.bio' = 1E5,'adult.bio' = 1E3,'catch' = 1E2),
                   var.slope = c('bot.temp' = 0.05,'diatom' = 0,'juv.bio' = -0.01,'adult.bio' = -0.05,'catch' = 0),
                   diagnostics = T
)

saveRDS(dsem.data$data,here::here('data-raw','dsem_example_data.rds'))

dsem.data.test = DSEMrand::fit_generalized_dsem(adj.file = here::here('data-raw','example_unknown_dag.csv'),
                     lags.file = here::here('data-raw','example_lags_dag.csv'),
                     data.file = here::here('data-raw','dsem_example_data.rds')
)

saveRDS(dsem.data.test, here::here('data-raw','dsem_example_results.rds'))

DSEMrand::compare_dsem_results(true_weights_file = here::here('data-raw','example_weight_dag.csv'),
                     dsem_output = readRDS(here::here('data-raw','dsem_example_results.rds'))
)

# 1. Generate the random matrix
ref_adj =  read.csv(here::here('data-raw','example_weight_dag.csv'),row.names = 1) |> as.matrix() |>  sign()
rand_adj <- DSEMrand::randomize_dsem_graph(ref_adj)

# 2. Convert it to an igraph object
g <- igraph::graph_from_adjacency_matrix(rand_adj, mode = "directed")

# 3. Plot using a hierarchical layout
igraph::plot.igraph(g, 
             layout = igraph::layout_with_sugiyama(g)$layout, 
             vertex.color = "lightblue", 
             vertex.size = 30, 
             edge.arrow.size = 0.6,
             main = "Randomized Hierarchical DAG")

test.rand = DSEMrand::evaluate_random_graphs(ref_adj,
                                             ref_lags = here::here('data-raw','example_lags_dag.csv'),
                                             N = 500,deduplicate =T)


batch_results = test.rand
data.file = here::here('data-raw','dsem_example_data.rds')

# tictoc::tic()
# batch_fits = DSEMrand::fit_random_graphs(batch_results = test.rand,
#                                       data.file = here::here('data-raw','dsem_example_data.rds'),
#                                       cores = 1
# )
# tictoc::toc()

tictoc::tic()
batch_fits = DSEMrand::fit_random_graphs(batch_results = test.rand,
                                         data.file = here::here('data-raw','dsem_example_data.rds'),
                                         cores = 8
)
tictoc::toc()

results = DSEMrand::compare_dsem_distributions(ref_fit = dsem.data.test, batch_fits, plot =T)
results$plot
results$model_fit_comparison
results$parameter_comparison


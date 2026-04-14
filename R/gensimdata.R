# Function that creates data from an generic adj matrix
# adj[i,j] = 1; i -> j

#' @export
gensimdata <- function (adj,
                        N      = 500,
                        b0     = 0,
                        ss     = 1,
                        s      = 1,
                        T      = 1,
                        ss_lag = 0.5) {

  if (!is.matrix(adj)) {
    stop ("'adj' must be a matrix.")
  }
  if (nrow(adj) != ncol(adj)) {
    stop ("'adj' must be a square matrix.")
  }
  if (!all(adj %in% c(0, 1))) {
    stop ("'adj' must contain only 0 and 1.")
  }

  # Name the rows and columns
  p <- nrow(adj)
  if (is.null(colnames(adj))) {
    colnames(adj) <- rownames(adj) <- paste0("V", seq_len(p))
  }
  node_names <- colnames(adj)

  # check if its a valid graph
  g <- igraph::graph_from_adjacency_matrix(adj)
  if (!igraph::is_dag(g)) {
    stop ("'adj' contains a cycle. The graph must be a DAG.")
  }

  # topologically order the graph nodes
  topo_idx <- as.integer(igraph::topo_sort(g))

  # helper function to simulate one slice at a time
  # Only simulates normal data
  sim_slice <- function (prev_slice = NULL) {

    # Initialize an empty data matrix
    X <- matrix(0, nrow = N, ncol = p)
    colnames(X) <- node_names

    # Loop over all Nodes. Has two If statments for the time steps
    for (j in topo_idx) {

      # First identify the parents of node_j. Parents are the rows
      intra_parents <- which(adj[, j] == 1)

      # T = 1
      if (is.null(prev_slice)) {

        if (length(intra_parents) == 0L) {

          X[, j] <- sNorm(N  = N,
                          b0 = b0,
                          s  = s)

        } else {

          X[, j] <- cNorm(N          = N,
                          parentData = lapply(intra_parents, # returns a list of node data (cols)
                                             function (i) X[, i]),
                          b0         = b0,
                          b1         = rep(ss, length(intra_parents)), # same ss for all parents
                          s          = s)

        }

        # T > 1
      } else {

        # Data from the T - 1 node
        lag_col <- prev_slice[, j]

        # root node in a slice
        if (length(intra_parents) == 0L) {

          X[, j] <- cNorm(N          = N,
                          # only input data is the previous T- 1 node data
                          parentData = list(lag_col), 
                          b0         = b0,
                          b1         = ss_lag,
                          s          = s)

        } else {

          X[, j] <- cNorm(N          = N,
                          parentData = c(lapply(intra_parents,
                                                function (i) X[, i]),
                                         list(lag_col)),
                          b0         = b0,
                          b1         = c(rep(ss, length(intra_parents)),
                                         ss_lag),
                          s          = s)

        }

      }

    }

    return (X)

  }

  # First Slice
  slice1 <- sim_slice(prev_slice = NULL)

  if (T == 1L) {
    return (slice1)
  }

  data_list = vector("list", T) # pre-allocation
  names(data_list) = paste0("T_", seq_len(T))
  data_list[[1L]] = slice1

  # Filling in the simulation data for remaining time points
  for (t in seq(2L, T)) {
    data_list[[t]] <- sim_slice(prev_slice = data_list[[t - 1L]])
  }

  return (data_list)

}

# ---------------
adj_gn4 <- matrix(c(0, 1, 1, 0,   # T1 --> T2, T1 --> T3
                    0, 0, 0, 1,   # T2 --> T4
                    0, 0, 0, 0,   # T3 (sink)
                    0, 0, 1, 0),  # T4 --> T3
                  nrow = 4, byrow = TRUE)

colnames(adj_gn4) <- rownames(adj_gn4) <- paste0("V", seq_len(nrow(adj_gn4)))
g <- igraph::graph_from_adjacency_matrix(adj_gn4)
igraph::topo_sort(g)

paste0("V_", seq_len(4))

#' @importFrom igraph components graph_from_adjacency_matrix
#' @importFrom expm expm

# build_adj
#
# Converts a vector of edge states ("currentES") into an adjacency matrix
# representation of a directed graph.
#
# @param coord A matrix specifying the row and column indices for each possible
# edge. Must contain rows named "Rows" and "Columns", where each column
# corresponds to an edge.
#
# @param currentES A vector of edge states for the current graph. Each entry
# indicates the direction or absence of an edge:
#   0 = forward edge (i → j),
#   1 = reverse edge (j → i),
#   2 = no edge.
#
# @param n An integer specifying the number of nodes in the graph.
#
# @return A numeric adjacency matrix of dimension n × n, where a value of 1
# indicates a directed edge and 0 indicates no edge.
build_adj <- function(coord, currentES, n) {
  A <- matrix(0, nrow = n, ncol = n)
  rows <- coord["Rows", ]
  cols <- coord["Columns", ]
  
  for (k in seq_along(currentES)) {
    type <- currentES[k]
    i <- rows[k]
    j <- cols[k]
    
    if (type == 0) {
      # forward edge i → j
      A[i, j] <- 1
      
    } else if (type == 1) {
      A[j, i] <- 1
      
    } else if (type == 2) {
      # no edge
    }
  }
  return(A)
}

# adj_to_currentES
#
# Converts an adjacency matrix representation of a directed graph into a
# vector of edge states ("currentES").
#
# @param adj A numeric adjacency matrix representing a directed graph, where
# a value of 1 indicates a directed edge and 0 indicates no edge.
#
# @param coord A matrix specifying the row and column indices for each possible
# edge. Must contain rows named "Rows" and "Columns", where each column
# corresponds to an edge.
#
# @return A vector of edge states for the current graph. Each entry encodes
# the direction or absence of an edge:
#   0 = forward edge (i → j),
#   1 = reverse edge (j → i),
#   2 = no edge.
adj_to_currentES <- function(adj, coord) {
  rows <- coord["Rows", ]
  cols <- coord["Columns", ]
  
  currentES <- integer(length(rows))
  
  for (k in seq_along(rows)) {
    i <- rows[k]
    j <- cols[k]
    
    if (adj[i, j] == 1) {
      currentES[k] <- 0     # forward edge
    } else if (adj[j, i] == 1) {
      currentES[k] <- 1     # backward edge
    } else {
      currentES[k] <- 2     # no edge
    }
  }
  
  return(currentES)
}

# dag_constraint
#
# Computes a smooth acyclicity constraint for a directed graph represented
# by an adjacency matrix. The constraint is zero if and only if the graph
# is a Directed Acyclic Graph (DAG).
#
# @param A A square numeric adjacency matrix representing a directed graph,
# where A[i, j] = 1 indicates a directed edge from node i to node j.
#
# @return A numeric scalar giving the value of the DAG constraint. A value of
# zero indicates that the graph is acyclic, while positive values indicate
# the presence of one or more directed cycles.
dag_constraint <- function(A) {
  B <- A * A
  expB <- expm(B)
  m <- nrow(A)
  return(sum(diag(expB)) - m)
}

# make_acyclic
#
# Iteratively removes directed cycles from a graph by randomly selecting
# edges that participate in a cycle and optionally reversing their direction
# according to a prior distribution. The procedure terminates once the graph
# is acyclic.
#
# @param adj A square numeric adjacency matrix representing a directed graph,
# where adj[i, j] = 1 indicates a directed edge from node i to node j.
#
# @param prior A numeric vector of prior probabilities for edge states.
# The entries correspond to:
#   prior[1] = probability of reverse edge (j → i),
#   prior[2] = probability of forward edge (i → j),
#   prior[3] = probability of no edge.
#
# @return A numeric adjacency matrix with all directed cycles removed,
# representing a Directed Acyclic Graph (DAG).

make_acyclic <- function(adj, prior, coord, edgeType, nCPh, pmr) {
  
  # Repeat until the adjacency matrix contains no directed cycles
  repeat {
    
    # Identify all edges that are part of at least one cycle.
    # Each element of `ces` is a length-2 vector (i, j) indicating
    # a directed edge i -> j that participates in a cycle.
    ces <- find_cycle_edges(adj)
    
    # If no cycle edges are found, the graph is acyclic; exit loop.
    if (is.null(ces)) break
    
    # Randomly select ONE edge from the set of cycle-forming edges.
    # This edge will be modified to help break the cycle.
    e <- ces[[sample.int(length(ces), 1)]]
    
    # The adjacency matrix encodes direction via (row, col),
    # but `coord` stores edges as unordered node pairs (min, max).
    # We reorder the selected edge indices to match `coord`.
    eEdgeType <- e
    if (e[1] > e[2]) {
      eEdgeType <- c(e[2], e[1])
    }
    
    # Find the column index `i` in `coord` that corresponds to
    # this unordered edge pair. This index is used to:
    #   1) look up edgeType
    #   2) apply correct prior constraints
    if (any(edgeType == 1) && (pmr || nCPh >= 1)) {
      for (i in 1:ncol(coord)) {
        if (all(as.vector(coord[, i]) == eEdgeType)) {
          break
        }
      }
    } else { # If no edgetypes, then all in edgetype are 0 and it doesn't matter what i is
      i = 1
    }
    
    # Determine the CURRENT edge state based on adjacency direction:
    #
    # If e[1] < e[2]:
    #   adj[e[1], e[2]] == 1  -> state 0 (forward)
    #
    # If e[1] > e[2]:
    #   adj[e[2], e[1]] == 1  -> state 1 (backward)
    #
    # The edge can only transition to the two OTHER states.
    if (e[1] < e[2]) { # Current state = 0 (forward)
      
      # Allowed transitions: 1 (backward) or 2 (no edge)
      states <- c(1, 2)
      
      # Corresponding indices in the prior vector:
      #   prior[2] → state 1
      #   prior[3] → state 2
      prPos <- c(2, 3)
      
    } else { # Current state = 1 (backward)
      
      # Allowed transitions: 0 (forward) or 2 (no edge)
      states <- c(0, 2)
      
      # Corresponding indices in the prior vector:
      #   prior[1] -> state 0
      #   prior[3] -> state 2
      prPos <- c(1, 3)
    }
    
    # Compute the conditional prior probabilities for the two
    # allowed transitions, taking into account:
    #   - edge type (gv-ge, gv-gv, etc.)
    #   - PMR constraints
    #   - presence of clinical phenotypes
    probability <- cPrior(
      prPos    = prPos,
      edgeType = edgeType[[i]],
      nCPh     = nCPh,
      pmr      = pmr,
      prior    = prior
    )
    
    # Sample the new edge state (0, 1, or 2) according to the
    # conditional prior probabilities returned by cPrior().
    newState <- sample(x = states, size = 1, prob = probability)
    
    # Remove the current edge entirely (both directions) so that
    # the new state can be cleanly applied.
    adj[e[1], e[2]] <- 0
    adj[e[2], e[1]] <- 0
    
    # Re-insert the edge based on the sampled new state:
    #
    # state 0 -> forward edge  (min → max)
    # state 1 -> backward edge (max → min)
    # state 2 -> no edge (do nothing, already done)
    if (newState == 0) {
      adj[eEdgeType[1], eEdgeType[2]] <- 1
    } else if (newState == 1) {
      adj[eEdgeType[2], eEdgeType[1]] <- 1
    }
  }
  
  # Return the modified adjacency matrix, guaranteed to be acyclic
  adj
}

# find_cycle_edges
#
# Identifies edges that participate in directed cycles within a graph by
# detecting strongly connected components (SCCs).
#
# @param adj A square numeric adjacency matrix representing a directed graph,
# where adj[i, j] = 1 indicates a directed edge from node i to node j.
#
# @return A list of integer vectors of length two, where each vector
# (u, v) represents a directed edge u → v that belongs to a cycle.
# Returns NULL if no cycles are present in the graph.
find_cycle_edges <- function(adj) {
  g <- graph_from_adjacency_matrix(adj, mode = "directed")
  scc <- components(g, mode = "strong")
  
  cycle_edges <- list()
  
  # Look at each SCC with >= 2 nodes (a cycle must have at least 2)
  for (comp_id in unique(scc$membership)) {
    nodes <- which(scc$membership == comp_id)
    if (length(nodes) < 2) next
    
    # extract edges whose endpoints are both inside this SCC
    for (u in nodes) {
      for (v in which(adj[u,] == 1)) {
        if (v %in% nodes) {
          cycle_edges <- append(cycle_edges, list(c(u, v)))
        }
      }
    }
  }
  
  if (length(cycle_edges) > 0) cycle_edges else NULL
}

# cycleRmvr
#
# Changes the edge direction of the current graph until there are no directed
# cycles.
#
# @param coord A matrix specifying the row and column indices for each possible
# edge. Must contain rows named "Rows" and "Columns", where each column
# corresponds to an edge.
#
# @param nNodes The number of nodes in the graph.
#
# @param currentES A vector of edge states for the current graph.
#
# @param prior A vector containing the prior probability of seeing each edge
# direction.
#
# @return A vector of the DNA of the individual with the cycles removed.
#
cycleRmvr <- function(currentES,
                      nNodes,
                      coord,
                      prior,
                      edgeType,
                      nCPh,
                      pmr) {
  
  # Convert currentES ito adj matrix
  A <- build_adj(coord, currentES, nNodes)
  
  # Check if it is not a DAG
  #if (dag_constraint(A) > 0) {
  A <- make_acyclic(A, prior, coord, edgeType, nCPh, pmr)
  currentES <- adj_to_currentES(A, coord)
  #}
  
  currentES
}
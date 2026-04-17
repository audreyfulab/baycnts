# build_master_data
#
# Converts a list of T data matrices into a single master data matrix with
# type-major column ordering so that cll() dispatch thresholds remain valid
#
# @param data_list A list of T matrices, each of dimension N x K.
# @param nGV       Number of genetic variant nodes (first nGV columns of each
#   slice).
# @param nCPh      Number of clinical phenotype nodes (last nCPh columns of
#   each slice).
#
# @return A single N x (K * T) matrix with type-major column ordering.
#
#' @export
build_master_data <- function(data_list, nGV, nCPh) {
  
  if (nGV > 0L || nCPh > 0) {
    stop("nGV, nCPh not supported yet")
  }
  
  ms <- do.call(cbind, data_list)
  ms

}

# build_master_AM
#
# Constructs the (K*T) x (K*T) master adjacency matrix for a DBN with T time
# steps and K nodes per slice.
#
# @param intra_AM A K x K adjacency matrix for intra-slice edges.
# @param inter_AM A K x K adjacency matrix for inter-slice (lag-1) edges.
# @param T_steps  Number of time steps.
# @param nGV      Number of genetic variant nodes.
# @param nCPh     Number of clinical phenotype nodes.
# @param K        Number of nodes per time slice.
#
# @return A (K*T) x (K*T) integer adjacency matrix.
#
#' @export
build_master_AM <- function(intra_AM, inter_AM, T_steps, nGV, nCPh, K) {

  # pre-allocate
  nM  <- K * T_steps
  mAM <- matrix(0L, nrow = nM, ncol = nM)

  for (t in seq_len(T_steps)) {

    # generate vector of idx to replace
    idx_t <- vapply(seq_len(K),
                    function(j) master_node_idx(j, t, K, T_steps, nGV, nCPh),
                    integer(1L))

    mAM[idx_t, idx_t] <- intra_AM

    if (t > 1L) {

      idx_prev <- vapply(seq_len(K),
                         function(j) master_node_idx(j, t - 1L, K, T_steps, nGV, nCPh),
                         integer(1L))
      # stack the lag block on top of the inter slice adj
      mAM[idx_prev, idx_t] <- inter_AM

    }

  }

  mAM

}

# master_node_idx
#
# @param j        Original node index within a single time slice (1-based).
# @param t        Time step
# @param K        Number of nodes per time slice.
# @param T_steps  Number of time steps.
# @param nGV      Number of genetic variant nodes.
# @param nCPh     Number of clinical phenotype nodes.
#
# @return Integer master column index.
#
master_node_idx <- function(j, t, K, T_steps, nGV = 0, nCPh = 0) {

  if (nGV > 0L || nCPh > 0L) {
    stop('nGV and nCPh not supported yed')
  }

  idx = (t - 1L) * K + j
  idx
  
}

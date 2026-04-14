#' mhEdge
#'
#' The main function for the Metropolis-Hastings algorithm. It returns the
#' posterior distribution of the edge directions.
#'
#' @param data A matrix with the variables across the columns and the
#' observations down the rows. If there are genetic variants in the data these
#' variables must come before the remaining variables. If there are clinical
#' phenotypes in the data these variables must come after all other variables.
#' For example, if there is a data set with one genetic variant variable, three
#' gene expression variables, and one clinical phenotype variable the first
#' column in the data matrix must contain the genetic variant data, the next
#' three columns will contain the gene expression data, and the last column will
#' contain the clinical phenotype data.
#'
#' @param adjMatrix An adjacency matrix indicating the edges that will be
#' considered by the Metropolis-Hastings algorithm. An adjacency matrix is a
#' matrix of zeros and ones. The ones represent an edge and its direction
#' between two nodes.
#'
#' @param prior A vector containing the prior probability for the three edge
#' states.
#'
#' @param nCPh The number of clinical phenotypes in the graph.
#'
#' @param nGV The number of genetic variants in the graph.
#'
#' @param pmr Logical. If true the Metropolis-Hastings algorithm will use the
#' Principle of Mendelian Randomization (PMR). This prevents the direction of an
#' edge pointing from a gene expression or a clinical phenotype node to a
#' genetic variant node.
#'
#' @param burnIn A number between 0 and 1 indicating the percentage of the
#' sample that will be discarded.
#'
#' @param iterations An integer for the number of iterations to run the MH
#' algorithm.
#'
#' @param thinTo An integer indicating the number of observations the chain
#' should be thinned to.
#'
#' @param progress Logical. If TRUE the runtime in seconds for the cycle finder
#' and log likelihood functions and a progress bar will be printed.
#'
#' @param threads An integer specifying the number of threads to use for finding
#' cycles. Default is 1
#'
#' @param maxParents Maximum number of parents allowed for any node. If NULL
#' (default), no limit is imposed. This helps prevent overfitting and numerical
#' issues with compositional data.
#' 
#' @param interAm A square matrix with the same dimensions as adjMatrix 
#' specifying the fixed inter-slice (lag-1) connections. Default is NULL.
#'
#' @return An object of class baycn containing 9 elements:
#'
#' \itemize{
#'
#' \item burnIn -- The percentage of MCMC iterations that will be discarded from
#' the beginning of the chain.
#'
#' \item chain -- A matrix where each row contains the vector of edge states for
#' the accepted graph.
#'
#' \item decimal -- A vector of decimal numbers. Each element in the vector is
#' the decimal of the accepted graph.
#'
#' \item iterations -- The number of iterations for which the
#' Metropolis-Hastings algorithm is run.
#'
#' \item posteriorES -- A matrix of posterior probabilities for all three edge
#' states for each edge in the network.
#'
#' \item posteriorPM -- A posterior probability adjacency matrix.
#'
#' \item likelihood -- A vector of log likelihood values. Each element in the
#' vector is the log likelihood of the accepted graph.
#'
#' \item stepSize -- The number of MCMC iterations discarded between each
#' accepted graph.
#'
#' \item time -- The runtime of the Metropolis-Hastings algorithm in seconds.
#'
#' }
#'
#' @references Martin, E. A. and Fu, A. Q. (2019). A Bayesian approach to
#' directed acyclic graphs with a candidate graph. \emph{arXiv preprint
#' arXiv:1909.10678}.
#'
#' @importFrom methods new
#'
#' @importFrom stats dbinom
#'
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @examples
#' # Generate data under topology m1_gv.
#' # Use ?simdata for a description and graph of m1_gv.
#' data_m1 <- simdata(b0 = 0,
#'                    N = 200,
#'                    s = 1,
#'                    graph = 'm1_gv',
#'                    ss = 1,
#'                    q = 0.27)
#'
#' # Create an adjacency matrix with the true edges.
#' am_m1 <- matrix(c(0, 1, 0,
#'                   0, 0, 1,
#'                   0, 0, 0),
#'                 byrow = TRUE,
#'                 nrow = 3)
#'
#' # Run the Metropolis-Hastings algorithm on the data from m1_gv using the
#' # Principle of Mendelian Randomization (PMR) and the true edges as the input.
#' mh_m1_pmr <- mhEdge(data = data_m1,
#'                     adjMatrix = am_m1,
#'                     prior = c(0.05,
#'                               0.05,
#'                               0.9),
#'                     nCPh = 0,
#'                     nGV = 1,
#'                     pmr = TRUE,
#'                     burnIn = 0.2,
#'                     iterations = 1000,
#'                     thinTo = 200,
#'                     progress = FALSE)
#'
#' summary(mh_m1_pmr)
#'
#' # Generate data under topology gn4.
#' # Use ?simdata for a description and graph of gn4.
#' data_gn4 <- simdata(b0 = 0,
#'                     N = 200,
#'                     s = 1,
#'                     graph = 'gn4',
#'                     ss = 1)
#'
#' # Create an adjacency matrix with the true edges.
#' am_gn4 <- matrix(c(0, 1, 1, 0,
#'                    0, 0, 0, 1,
#'                    0, 0, 0, 0,
#'                    0, 0, 1, 0),
#'                  byrow = TRUE,
#'                  nrow = 4)
#'
#' # Run the Metropolis-Hastings algorithm on the data from gn4 with the true
#' # edges as the input.
#' mh_gn4 <- mhEdge(data = data_gn4,
#'                  adjMatrix = am_gn4,
#'                  prior = c(0.05,
#'                            0.05,
#'                            0.9),
#'                  nCPh = 0,
#'                  nGV = 0,
#'                  pmr = FALSE,
#'                  burnIn = 0.2,
#'                  iterations = 1000,
#'                  thinTo = 200,
#'                  progress = FALSE)
#'
#' summary(mh_gn4)
#' 
#' @export
#'
mhEdge <- function (data,
                    adjMatrix,
                    prior = c(0.05,
                              0.05,
                              0.9),
                    nCPh = 0,
                    nGV = 0,
                    pmr = FALSE,
                    burnIn = 0.2,
                    iterations = 1000,
                    thinTo = 200,
                    progress = TRUE,
                    threads = 1,
                    maxParents = NULL,
                    interAm = NULL) {

  # Preprocessing checks -------------------------------------------------------

  # baycnts specific
  # check if it is a time series input
  isDbn = is.list(data)
  
  # Check that data is a matrix.
  if (!is.matrix(data) && !is.list(data)) {

    stop ('data is not a matrix or a list of matrices')

  }
  
  # Check for interAm when in time series mode
  if (isDbn && is.null(interAm)) {
    
    stop('need to define interAm when you have multi-dim data')
    
  }

  # Check that adjMatrix is a square matrix.
  if (dim(adjMatrix)[1] != dim(adjMatrix)[2]) {

    stop ('adjMatrix must be a square matrix')

  }

  # Check that adjMatrix has all zeros on the diagonal.
  if (sum(diag(adjMatrix)) != 0) {

    stop ('adjMatrix must have zeros on the diagonal')

  }
  
  # Check interAm has the same dim as adjMatrix
  if (isDbn) {
    if (
      !is.matrix(interAm) ||
      nrow(interAm) != nrow(adjMatrix) ||
      ncol(interAm) != ncol(adjMatrix)
    ) {
      
      stop('interAm must be a square matrix with the same dimesions as adjMatrix')
      
    }
    
  }
  
  # baycnts specific
  # build master objects for dbn mode or otherwise no change needed
  if (isDbn) {
    dataRef = data[[1L]] # for the checks
    Tsteps = length(data) # data is a list here
    K = ncol(dataRef) # number of nodes (one time step)
    nNodes = K * Tsteps
    nGV = nGV * K
    nCPh = nCPh * K
    
    data = build_master_data(data_list = data, nGV = nGV, nCPh = nCPh) # data is a matrix here
  } else {
    dataRef = data
    Tsteps = 1
    K = ncol(data)
    nNodes = K
  }

  # Check that adjMatrix has the same number of columns as data.
  if (dim(adjMatrix)[1] != dim(dataRef)[2]) {

    stop ('data must have the same number of columns as adjMatrix')

  }

  # Check the prior probability vector has three elements.
  if (length(prior) != 3) {

    stop ('prior must have three elements')

  }

  # Check the prior probability sums to 1.
  if (sum(prior) != 1) {

    stop ('prior must sum to 1')

  }

  # Check that burnIn, thinTo, and iterations agree.
  if ((1 - burnIn) * iterations < thinTo) {

    stop ('thinTo is greater than the number of iterations after the burnIn')

  }

  # Check that nGV > 0 if pmr is TRUE.
  if (pmr && nGV == 0) {

    stop('nGV must be at least 1 when pmr = TRUE')

  }

  # Check that the number of genetic variants is less than or equal to the
  # number of nodes.
  if (nGV > dim(dataRef)[2]) {

    stop ('nGV must be less than or equal to the number of columns in data')

  }

  # Check that the number of clinical phenotypes is less than or equal to the
  # number of nodes.
  if (nCPh > dim(dataRef)[2]) {

    stop ('nCPh must be less than or equal to the number of columns in data')

  }

  # Save the time when baycn starts.
  startTime <- Sys.time()

  # Set up ---------------------------------------------------------------------

  # Determine the number of nodes in the graph.
  nNodes <- ncol(adjMatrix)

  # Get the coordinates of the non zero elements in the adjacency matrix.
  coord <- coordinates(adjMatrix = adjMatrix)

  # Determine the number of edges in the network
  nEdges <- dim(coord)[2]

  # Determine the edge type for each edge in the graph (gv-ge or ge-ge, gv-gv)
  edgeType <- detType(coord = coord,
                      nEdges = nEdges,
                      nCPh = nCPh,
                      nGV = nGV,
                      nNodes = nNodes,
                      pmr = pmr)

  # Calculate the probability of choosing from 1 to nEdges.
  ztbProb <- dbinom(x = 1:nEdges,
                    size = nEdges,
                    prob = 1 / nEdges) / (1 - (1 - 1/nEdges)^nEdges)

  # A vector of integers that indicate which iterations will be kept.
  keep <- whichIt(burnIn = burnIn,
                  iterations = iterations,
                  thinTo = thinTo)

  # Display runtime message for finding potential cycles.
  if (progress) {

    cat('Identifying potential cycles: ')

    begin_time <- Sys.time()

  }

  # Display runtime message for finding potential cycles.
  if (progress) {

    end_time <- Sys.time()

    cat(as.double(round(end_time - begin_time, 3),
                  units = 'secs'), 'seconds', '\n')

  }

  # Initialize vectors, lists, and matrices ------------------------------------

  # A counter used to fill in the MarkovChain matrix.
  counter <- 0

  # Create a matrix to hold the accepted graph at each iteration of the MH
  # algorithm.
  MarkovChain <- matrix(nrow = thinTo,
                        ncol = nEdges)

  # Create a vector to hold the likelihood for the whole graph.
  likelihood <- vector(mode = 'numeric',
                       length = thinTo)

  # Create a vector to hold the decimal number for the accepted graph.
  graphDecimal <- vector(mode = 'integer',
                         length = thinTo)

  # Create a vector to hold the log likelihood for each node in the current
  # graph.
  currentLL <- vector(mode = 'numeric',
                      length = nNodes)

  # Create the likelihood environment ------------------------------------------

  # Create a new environment that has the log likelihood for all possible parent
  # combinations for each node. (LOG Likelihood Environment)
  # This will work no matter the dim of the data
  logle <- logLikEnv(data = data,
                     nCPh = nCPh,
                     nGV = nGV,
                     nNodes = nNodes)

  # Generate a starting graph --------------------------------------------------

  # Randomly generate a starting graph.
  currentES <- sample(x = c(0, 1, 2),
                      size = nEdges,
                      replace = TRUE,
                      prob = prior)
  #currentES <- c(1, 1, 1, 1, 1, 1)

  # Check if any gv nodes have ge node parents.
  if (pmr || nCPh >= 1) {

    # Reverse the direction of any ge -> gv, gv -> cph, or ge -> cph edges.
    currentES[edgeType == 1 & currentES == 1] <- 0

  }

  # Check if there are potential directed cycles in the graph and
  # remove them if they are.
  currentES <- cycleRmvr(currentES = currentES,
                         nNodes = nNodes,
                         coord = coord,
                         prior = prior,
                         edgeType = edgeType, 
                         nCPh = nCPh, 
                         pmr = pmr)

  # Turn the current edge states in to an adjacency matrix.
  currentAM <- toAdjMatrix(coordinates = coord,
                           graph = currentES,
                           nEdges = nEdges,
                           nNodes = nNodes)

  # Enforce maxParents limit on initial graph
  if (!is.null(maxParents)) {
    nParents <- colSums(currentAM)
    while (any(nParents > maxParents)) {
      # Find nodes with too many parents
      overloaded <- which(nParents > maxParents)

      # For each overloaded node, randomly remove edges until under limit
      for (node in overloaded) {
        parents <- which(currentAM[, node] == 1)
        nToRemove <- length(parents) - maxParents

        if (nToRemove > 0) {
          toRemove <- sample(parents, nToRemove)
          # Find edge indices for these parents and set to 2 (no edge)
          for (parent in toRemove) {
            edgeIdx <- which(coord[1, ] == parent & coord[2, ] == node)
            if (length(edgeIdx) > 0) {
              currentES[edgeIdx] <- 2
            }
          }
        }
      }

      # Rebuild adjacency matrix
      currentAM <- toAdjMatrix(coordinates = coord,
                               graph = currentES,
                               nEdges = nEdges,
                               nNodes = nNodes)
      nParents <- colSums(currentAM)
    }
  }
  
  if (isDbn) {
    tsCurrentAM = build_master_AM(intra_AM = currentAM,
                                  inter_AM = interAm,
                                  T_steps = Tsteps,
                                  nGV = nGV,
                                  nCPh = nCPh,
                                  K = K)
  } else {
    tsCurrentAM = currentAM
  }

  # Look up the log likelihood for each node.
  currentLL <- lull(data = data,
                    am = tsCurrentAM,
                    likelihood = currentLL,
                    llenv = logle,
                    nCPh = nCPh,
                    nGV = nGV,
                    nNodes = nNodes,
                    wNodes = 1:nNodes)

  # Copy the currentLL vector to the proposedLL vector. They need to be the same
  # before starting the MH algorithm because only the nodes that have different
  # parents in the proposed graph will have the likelihood changed in the
  # proposedLL vector.
  proposedLL <- currentLL

  # Start the Metropolis-Hastings algorithm ------------------------------------

  # Create progress bar.
  if (progress) {

    pb <- txtProgressBar(min = 0, max = iterations, style = 3)

  }

  # Run through the iterations of the Metropolis Hastings algorithm.
  for (e in 1:iterations) {

    # Change the current graph to a different graph at each iteration.
    proposedES <- mutate(edgeType = edgeType,
                         graph = currentES,
                         nCPh = nCPh,
                         nEdges = nEdges,
                         pmr = pmr,
                         prior = prior,
                         ztbProb = ztbProb)

    # If there are potential directed cycles in the graph remove them
    # Does not include: edgetype, nCPh, pmr, 
    proposedES <- cycleRmvr(currentES = proposedES,
                           nNodes = nNodes,
                           coord = coord,
                           prior = prior,
                           edgeType = edgeType, 
                           nCPh = nCPh, 
                           pmr = pmr)

    # Get the indices for the edge states that are different between the current
    # and proposed edge state vectors.
    difference <- divergent(coordinates = coord,
                            currentES = currentES,
                            proposedES = proposedES)

    proposedAM <- convertAM(adjMatrix = currentAM,
                            coordinates = coord,
                            edgeStates = proposedES,
                            wEdges = difference$dEdges)

    # Check if any node exceeds maxParents limit
    if (!is.null(maxParents)) {
      # Count parents for each node (column sums of adjacency matrix)
      nParents <- colSums(proposedAM)

      # If any node has too many parents, reject this proposal immediately
      if (any(nParents > maxParents)) {
        # Skip to next iteration, keeping currentES
        if (progress) {
          setTxtProgressBar(pb, e)
        }
        next
      }
    }
    
    # if T > 1, need to modify the adj matrix
    if (isDbn) {
      tsProposedAM = build_master_AM(intra_AM = proposedAM,
                                     inter_AM = interAm,
                                     T_steps = Tsteps,
                                     nGV = nGV,
                                     nCPh = nCPh,
                                     K = K)
      
      # Track the nodes that changed
      dNodes = as.integer(
        unlist(
          lapply(difference$dNodes, function(j) {
            vapply(seq_len(Tsteps), function(t) {
              master_node_idx(j, t, K, Tsteps, nGV, nCPh)
            }, integer(1L))
          })
        )
      )
      
    } else {
      tsProposedAM = proposedAM
      dNodes = difference$dNodes
    }

    # Look up the log likelihood for each node whose parents have changed.
    proposedLL <- lull(data = data,
                       am = tsProposedAM,
                       likelihood = currentLL,
                       llenv = logle,
                       nCPh = nCPh,
                       nGV = nGV,
                       nNodes = nNodes,
                       wNodes = dNodes)

    # Determine whether to accept the proposed graph proposedES or the current
    # graph currentES.
    decision <- caRatio(current = currentES,
                        currentLL = currentLL,
                        edgeType = edgeType,
                        nCPh = nCPh,
                        pmr = pmr,
                        proposed = proposedES,
                        proposedLL = proposedLL,
                        prior = prior,
                        wEdges = difference$dEdges)

    # If the proposed graph is accepted then change the currentES and currentLL
    # vectors to the proposedES and proposedLL
    if (decision == 'proposed') {

      currentES <- proposedES
      currentLL <- proposedLL
      currentAM <- proposedAM

    }

    # Store the accepted graph after burnin and thinning.
    if (e %in% keep) {

      # Update the counter if e is in the keep vector
      counter <- counter + 1

      # Store the accepted graph in the MarkovChain matrix
      MarkovChain[counter, ] <- currentES

      # Store the log likelihood of the accepted graph.
      likelihood[[counter]] <- sum(currentLL)

      # Calculate the decimal for the accepted graph.
      graphDecimal[[counter]] <- sum(currentES * 3^(0:(nEdges - 1)))

    }

    # Update progress bar.
    if (progress) {

      setTxtProgressBar(pb, e)

    }

  }

  # Close the progress bar.
  if (progress) {

    close(pb)

  }

  # Posterior probability matrix -----------------------------------------------

  # Create a matrix to hold the posterior probability of each edge.
  posteriorES <- as.data.frame(matrix(nrow = nEdges,
                                      ncol = 4))

  # zero - the edge points from the node with a smaller index to the
  # node with a larger index
  # one - the edge points from the node with a larger index to the
  # node with a smaller index
  # two - the edge is absent between the two nodes
  colnames(posteriorES) <- c('nodes', 'zero', 'one', 'two')

  rownames(posteriorES) <- paste('edge', 1:nEdges, sep = '')

  # Calculate the proportion of 0, 1, and 2 for each edge in the
  # network.
  for (e in 1:nEdges) {

    # calculate the posterior probability for state 0, 1, and 2
    edge_0 <- sum(MarkovChain[, e] == 0)
    edge_1 <- sum(MarkovChain[, e] == 1)
    edge_2 <- sum(MarkovChain[, e] == 2)

    posteriorES[e, 2:4] <- c(round(edge_0 / thinTo, 4),
                             round(edge_1 / thinTo, 4),
                             round(edge_2 / thinTo, 4))

    # Determine which nodes the edge is between
    posteriorES[e, 1] <- paste0(colnames(dataRef)[coord[1, e]],
                                '-',
                                colnames(dataRef)[coord[2, e]])

  }

  # Posterior probability adjacency matrix -------------------------------------

  # Probabilities for upper triangular matrix
  upper <- posteriorES[, 2]

  # Probabilities for lower triangular matrix
  lower <- posteriorES[, 3]

  # Initialize a pxp matrix to fill in later with the estimated probabilities
  # from bacyn.
  posteriorPM <- matrix(0, nrow = nNodes, ncol = nNodes)

  # loop through the row column coordinates of the true edge indicies to fill in
  # the posterior edge probabilities.
  for (v in 1:nEdges) {

    # Fill in the upper triangular matrix from the '0' column.
    posteriorPM[coord[1, v], coord[2, v]] <- upper[[v]]

    # Fill in hte lower triangular matrix from the '1' column.
    posteriorPM[coord[2, v], coord[1, v]] <- lower[[v]]

  }

  # name the rows and columns with the names from the data matrix.
  colnames(posteriorPM) <- colnames(dataRef)
  rownames(posteriorPM) <- colnames(dataRef)

  # Calculate the stepSize based on the burnIn, iterations, and thinTo
  stepSize <- ifelse(burnIn == 0,
                     ceiling(iterations / thinTo),
                     ceiling(((1 - burnIn) * iterations) / thinTo))

  # Save the time when baycn ends.
  endTime <- Sys.time()

  # Create a baycn object ------------------------------------------------------
  baycnObj <- new('baycn',
                  burnIn = burnIn * 100,
                  chain = MarkovChain,
                  decimal = graphDecimal,
                  iterations = iterations,
                  posteriorES = posteriorES,
                  posteriorPM = posteriorPM,
                  likelihood = likelihood,
                  stepSize = stepSize,
                  time = as.double((endTime - startTime),
                                   units = 'secs'))

  return (baycnObj)

}

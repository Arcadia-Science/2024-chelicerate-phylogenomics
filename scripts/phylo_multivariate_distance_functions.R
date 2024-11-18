require(phytools)
require(Rcpp)
require(RcppArmadillo)

# The following R-script contains a number of scripts that have been written
# to, as efficiently as possible, calculate multivariate distances between
# species or proteins, given a set of trait measurements and a corresponding
# species tree or gene family tree.

# Numerous functions have been written in Rcpp to take advantage of C's
# efficiency, particularly for the mahalanobis distance calculations among
# thousands of proteins and thousands of gene families.

# These functions get to be a bit complex, so I've attempted to walk through
# what they are, and how they work together below. The two core functions are:
# 1) "phylo_gls_transform":
#    - This function carries out the phylogenetic GLS transformation described
#      in Butler et al., 2000 (described in more detail below). In brief,
#      this function takes trait data and a corresponding phylogeny and returns
#      transformed data that represent the residual trait variation after
#      accounting for variance explained by phylogeny.
#    - NOTE - this is effectively a wrapper function around the Rcpp
#             function "phylo_correction"
# 2) "pairwise_phylo_mahalanobis":
#    - This takes the transformed data and the tree, and calculates the pairwise
#      mahalanobis distances. Distances are calculated using data corresponding
#      to terminal nodes (e.g. species) by default. Otherwise distances can be
#      calculated using data corresponding to all terminal AND internal nodes
#      by specifying include_internal = FALSE. Doing so is appropriate when
#      calculating distances among gene families using gene family evolutionary
#      event counts.
#    - NOTE - Internally this calls a number of internal functions (described
#             in more detail later):
#             - "compute_inverse_covariance"
#             - "chunk_pairs"
#             - "process_chunks"
# 3) "phylo_correction":
#    - This Rcpp function takes as input the phylogenetic covariance matrix and
#      traits, returning either a vector or matrix of phylogenetic GLS
#      transformed data if one or multiple traits are provided respectively
# 4) "compute_inverse_covariance":
#    - This Rcpp function calculates the inverse of the covariance matrix of
#      transformed data for use in calculation of the mahalanobis distances.
#      It is implemented in Rcpp for the increased efficiency of matrix algebra
# 5) "chunk_pairs":
#    - Helper function to facilitate the calculation of mahalanobis distances.
#    - Generates all pairwise comparisons (among indices), and then splits
#      these pairs into smaller groups (chunks) of a specified size.
#      Each chunk is a matrix of pairs, and the function returns a list of
#      these chunk matrices.
# 6) "process_chunks":
#    - Takes the list of "chunked" pairwise comparisons produced by
#      "chunk_pairs" and the inverse phenotypic covariance matrix produced by
#      "compute_inverse_covariance", and calculates the pairwise mahalanobis
#      distances in parallel using mclapply returning the full matrix.
#    - Internally, calls the Rcpp function "compute_chunk_distances" that
#      calculates the actual distances for each chunk.
# 7) "compute_chunk_distances":
#     - This function calculates the actual pairwise mahalanobis distances using
#       the matrix of data, inverse of the phenotypic covariance matrix, and
#       matrix of indices indicating the pair of observations we want to
#       calculate distances among (i.e. the chunk)

# Define a function that conducts the phylogenetic GLS transformation,
# returning the residual trait variation after accounting for the phylogenetic
# effect
phylo_gls_transform <-
  function(data, tree, include_internal = FALSE) {
    # Obtain the phylogenetic VCV for both species and internal nodes:
    phylo_vcv <- vcvPhylo(tree, anc.nodes = include_internal)

    if (include_internal == TRUE) {
      colnames(phylo_vcv) <- c(tree$tip.label, tree$node.label[-1])
      rownames(phylo_vcv) <- c(tree$tip.label, tree$node.label[-1])

      # Identify the root node:
      root_node_name <- tree$node.label[1]
      root_node <-
        which(c(tree$tip.label, tree$node.label) == root_node_name)
    }

    if (is.list(data)) {
      # Assume data is a list containing s_counts, sl_counts, d_counts, and t_counts
      # Remove the root node from the data:
      transf_data <- list()
      for (i in 1:length(data)) {
        # Remove the root node from the data:
        data[[i]] <-
          data[[i]][, -which(colnames(data[[i]]) == root_node_name)]
        # Ensure the data matches the order of the species tree
        data[[i]] <-
          data[[i]][, match(c(tree$tip.label, tree$node.label[-1]), colnames(data[[i]]))]
        dim_names <- list(rownames(data[[i]]), colnames(data[[i]]))
        transf_data[[i]] <-
          t(do.call(cbind, phylo_correction(t(data[[i]]), phylo_vcv)))
        dimnames(transf_data[[i]]) <- dim_names
      }
      transf_data <- do.call(cbind, transf_data)
    } else {
      if (include_internal == TRUE) {
        # Remove the root node from the data:
        data <- data[, -which(colnames(data) == root_node_name)]
        # Ensure the data matches the order of the species tree
        data <-
          data[, match(c(tree$tip.label, tree$node.label[-1]), colnames(data))]
      } else {
        # Ensure the data matches the order of the species tree
        data <- data[, match(tree$tip.label, colnames(data))]
      }
      dim_names <- list(rownames(data), colnames(data))

      # Now, obtain the phylogenetic GLS transformed data
      transf_data <-
        t(do.call(cbind, phylo_correction(t(data), phylo_vcv)))
      dimnames(transf_data) <- dim_names
    }
    # Return the phylogenetically transformed data
    return(transf_data)
  }

# Define the function to calculate all pairwise Mahalanobis distances
pairwise_phylo_mahalanobis <-
  function(transf_data, tree, include_internal = FALSE) {
    # Calculate the inverse of the covariance matrix for calculation of #
    # mahalanobis distances
    iV <- compute_inverse_covariance(transf_data)

    # Break generate a list of all possible unique pairwise comparisons,
    # dividing them into chunks of 10000
    pair_chunks <- chunk_pairs(nrow(transf_data), 10000)

    # Now, calculate distances for each chunk with Rcpp, parallelizing across
    # chunks and returning a symmetrical distance matrix
    distance_matrix <- process_chunks(pair_chunks, transf_data, iV)
    rownames(distance_matrix) <- rownames(transf_data)
    colnames(distance_matrix) <- rownames(transf_data)
    return(distance_matrix)
  }

# A function to, given a matrix of trait observations and a corresponding
# phylogenetic variance-covariance matrix, transform the data using a
# phylogenetic GLS transformation using the method described in
# Butler et al., 2000 - https://doi.org/10.1111/j.0014-3820.2000.tb00026.x
# Note, by default, this assumes a model of brownian motion wherein:
# "the elements of the G matrix (phylogenetic vcv) are simply the
# amount of time (Tbm) from the root of the phylogeny to the most recent
# common ancestor of the pair of taxa. The diagonal entries (species variances)
# are the amounts of time from the root of the tree to each species (the depth
# of the tree if a molecular clock is assumed)"
cppFunction(
  depends = "RcppArmadillo",
  code = '
List phylo_correction(NumericMatrix traits, NumericMatrix phylo_vcv) {
  // Get the phylogenetic variance-covariance matrix
  arma::mat Gmatrix = as<arma::mat>(phylo_vcv);

  // Compute the Cholesky decomposition of the inverse of the Phylogenetic VCV (Gmatrix)
  arma::mat correction_factor;
  try {
    correction_factor = arma::chol(arma::inv(Gmatrix));
  } catch (std::runtime_error& e) {
    stop("Error in Cholesky decomposition: ", e.what());
  }

  // Get the number of traits and samples
  int num_traits = traits.ncol();
  int num_samples = traits.nrow();

  // Initialize a list to store the corrected traits
  List U(num_traits);

  // Loop through each trait and perform the correction
  for (int t = 0; t < num_traits; t++) {

    // Get the current trait
    arma::vec current_trait(num_samples);
    for (int i = 0; i < num_samples; i++) {
      current_trait(i) = traits(i, t);
    }

    // Create a model matrix for the current trait
    arma::mat M = arma::ones(num_samples, 1);
    M.col(0) = current_trait;

    // Perform the correction by post-multiplying the model matrix
    // by the correction factor matrix
    arma::mat corrected_traits = correction_factor * M;

    // Store the corrected trait in the list
    U[t] = wrap(corrected_traits);
  }

  return U;
}
'
)

# A function to calculate the inverse covariance matrix of the phylogenetically
# transformed data for use in calculation of the mahalanobis distance
cppFunction(
  depends = "RcppArmadillo",
  code = "
  arma::mat compute_inverse_covariance(const arma::mat& transformed_data) {
    // Compute the corrected VCV
    arma::mat corrected_vcv = arma::cov(transformed_data) + arma::eye(transformed_data.n_cols, transformed_data.n_cols) * 1e-6;
    // Compute the inverse of the corrected VCV
    arma::mat iV = arma::inv(corrected_vcv);
    return iV;
  }
"
)

# A function to facilitate the parallelization of mahalanobis distances
# by breaking the pairwise sets of comparisons into a list of "chunks"
chunk_pairs <- function(n, chunk_size) {
  total_pairs <- combn(1:n, 2)
  num_chunks <- ceiling(ncol(total_pairs) / chunk_size)
  split_cols <-
    split(
      1:ncol(total_pairs),
      rep(1:num_chunks, each = chunk_size, len = ncol(total_pairs))
    )
  lapply(split_cols, function(cols) {
    total_pairs[, cols, drop = FALSE]
  })
}

process_chunks <- function(pair_chunks, data, iV) {
  n <- nrow(data)
  distance_matrix <- matrix(NA, n, n)

  process_chunk <- function(pair_chunk) {
    chunk_distances <- compute_chunk_distances(data, iV, pair_chunk)
    list(pair_chunk = pair_chunk, distances = chunk_distances)
  }

  chunk_results <-
    mclapply(pair_chunks, process_chunk, mc.cores = detectCores())

  for (result in chunk_results) {
    pair_chunk <- result$pair_chunk
    distances <- result$distances
    row_indices <- pair_chunk[1, ]
    col_indices <- pair_chunk[2, ]
    for (i in 1:length(row_indices)) {
      distance_matrix[row_indices[i], col_indices[i]] <- distances[i]
      distance_matrix[col_indices[i], row_indices[i]] <-
        distances[i]
    }
  }
  diag(distance_matrix) <- 0
  return(distance_matrix)
}

# This function calculates the actual mahalanobis distances between the
# set of observating within a matrix of paired indices (i.e. a chunk) given
# the phenotypic data, the phenotypic covariance matrix, and the matrix
# of indices among which we'd like to make the calculations.
cppFunction(
  depends = "RcppArmadillo",
  code = "
  arma::mat compute_chunk_distances(const arma::mat& data, const arma::mat& iV, const arma::umat& pair_indices) {
    int num_pairs = pair_indices.n_cols;
    arma::mat distances(num_pairs, 1);

    for(int p = 0; p < num_pairs; ++p) {
      int i = pair_indices(0, p) - 1;  // Adjust for 0-based indexing in C++
      int j = pair_indices(1, p) - 1;  // Adjust for 0-based indexing in C++

      // Compute the difference between the two observations
      arma::rowvec diff = data.row(i) - data.row(j);

      // Ensure the matrix operation is correct
      arma::mat temp = diff * iV * arma::trans(diff);  // This should result in a 1x1 matrix
      distances(p) = std::sqrt(temp(0,0));  // Extract the scalar value from the 1x1 matrix
    }

    return distances;
  }
"
)

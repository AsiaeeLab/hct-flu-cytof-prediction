# ==============================================================================
# Clustering Methods
# ==============================================================================
# This file contains functions for different clustering approaches and
# assigning test cells to training clusters.

library(stats)  # For hclust, kmeans, dist, cutree
library(igraph)  # For Leiden clustering

# Load FlowSOM if available
if (requireNamespace("FlowSOM", quietly = TRUE)) {
  library(FlowSOM)
} else {
  warning("FlowSOM package not installed. FlowSOM clustering will not be available.")
}

# Load ConsensusClusterPlus if available (needed for FlowSOM)
if (requireNamespace("ConsensusClusterPlus", quietly = TRUE)) {
  library(ConsensusClusterPlus)
}


#' Perform hierarchical clustering
#'
#' @param data Cell-level feature matrix
#' @param k Number of clusters
#' @param method Linkage method (e.g., "complete", "average", "ward.D2")
#' @return List with cluster assignments and model object
hierarchical_clustering <- function(data, k, method = "complete") {
  hc <- hclust(dist(data), method = method)
  clusters <- cutree(hc, k = k)

  return(list(
    clusters = clusters,
    model = hc,
    method = "hierarchical"
  ))
}


#' Perform k-means clustering
#'
#' @param data Cell-level feature matrix
#' @param k Number of clusters
#' @param nstart Number of random starts
#' @param iter_max Maximum iterations
#' @return List with cluster assignments and model object
kmeans_clustering <- function(data, k, nstart = 25, iter_max = 100) {
  set.seed(123)  # For reproducibility
  km <- kmeans(data, centers = k, nstart = nstart, iter.max = iter_max)

  return(list(
    clusters = km$cluster,
    model = km,
    method = "kmeans"
  ))
}


#' Calculate Euclidean distance and assign to nearest cluster
#'
#' Assigns each test cell to the nearest training cluster based on
#' Euclidean distance to cluster centroids.
#'
#' @param test_data Test cell feature matrix
#' @param centroids Matrix of cluster centroids (k × p)
#' @return Vector of cluster assignments for test cells
euclidean_assignment <- function(test_data, centroids) {
  apply(test_data, 1, function(row) {
    dists <- apply(centroids, 1, function(centroid) {
      sum((row - centroid)^2)
    })
    which.min(dists)
  })
}


#' Assign test cells to training clusters
#'
#' Calculates centroids from training data and assigns test cells to
#' nearest cluster. This prevents information leakage.
#'
#' @param k Number of clusters
#' @param train_data Training cell feature matrix
#' @param train_clusters Training cluster assignments
#' @param test_data Test cell feature matrix
#' @return Vector of cluster assignments for test cells
assign_to_clusters <- function(k, train_data, train_clusters, test_data) {
  # Calculate centroids from training data only
  centroids <- sapply(1:k, function(cluster) {
    colMeans(train_data[train_clusters == cluster, , drop = FALSE])
  })

  # Assign test cells to nearest centroid
  test_clusters <- euclidean_assignment(test_data, t(centroids))

  return(test_clusters)
}


#' Perform FlowSOM clustering
#'
#' Self-Organizing Map followed by hierarchical clustering.
#' Popular method from flow cytometry literature.
#'
#' @param data Cell-level feature matrix
#' @param k Number of meta-clusters
#' @param xdim X dimension of SOM grid
#' @param ydim Y dimension of SOM grid
#' @return List with cluster assignments and model object
flowsom_clustering <- function(data, k, xdim = NULL, ydim = NULL) {
  if (!requireNamespace("FlowSOM", quietly = TRUE)) {
    stop("FlowSOM package is required but not installed. Install with: BiocManager::install('FlowSOM')")
  }

  # Adaptive SOM grid size: make grid larger than k to allow meaningful meta-clustering
  # Use a grid that gives us ~2-3x more SOM nodes than final clusters
  # This allows the SOM to capture structure while k controls final granularity
  if (is.null(xdim) || is.null(ydim)) {
    target_nodes <- max(100, k * 3)  # At least 100 nodes, or 3x k
    grid_size <- ceiling(sqrt(target_nodes))
    xdim <- grid_size
    ydim <- grid_size
  }

  # FlowSOM expects a flowFrame or matrix
  # Build SOM
  set.seed(123)
  fsom <- FlowSOM::FlowSOM(data,
                           colsToUse = 1:ncol(data),
                           xdim = xdim,
                           ydim = ydim,
                           nClus = k,
                           seed = 123)

  # Get cluster assignments
  clusters <- FlowSOM::GetClusters(fsom)

  return(list(
    clusters = clusters,
    model = fsom,
    method = "flowsom",
    som_grid_size = c(xdim, ydim)
  ))
}


#' Perform Leiden clustering
#'
#' Graph-based community detection using the Leiden algorithm.
#' Improved version of Louvain clustering.
#'
#' @param data Cell-level feature matrix
#' @param k Number of clusters (approximate, via resolution parameter)
#' @param n_neighbors Number of neighbors for kNN graph
#' @param resolution Resolution parameter (higher = more clusters)
#' @return List with cluster assignments and model object
leiden_clustering <- function(data, k, n_neighbors = 30, resolution = NULL) {
  # Build kNN graph
  # For large datasets, prefer RANN for speed and compatibility.
  if (requireNamespace("RANN", quietly = TRUE)) {
    # RANN includes the query point as the nearest neighbor.
    knn_result <- RANN::nn2(data = data, query = data, k = n_neighbors + 1)
    nn_index <- knn_result$nn.idx[, -1, drop = FALSE]
    edges <- cbind(
      rep(1:nrow(data), each = n_neighbors),
      as.vector(t(nn_index))
    )
  } else if (requireNamespace("FNN", quietly = TRUE)) {
    knn_result <- FNN::get.knn(data, k = n_neighbors)
    edges <- cbind(
      rep(1:nrow(data), each = n_neighbors),
      as.vector(t(knn_result$nn.index))
    )
  } else {
    # Fallback: simple distance-based
    dists <- as.matrix(dist(data))
    edges <- c()
    for (i in 1:nrow(data)) {
      neighbors <- order(dists[i, ])[2:(n_neighbors + 1)]
      edges <- rbind(edges, cbind(i, neighbors))
    }
  }

  # Create igraph
  g <- graph_from_edgelist(edges, directed = FALSE)
  g <- simplify(g)  # Remove duplicates

  # If resolution not provided, estimate from desired k
  if (is.null(resolution)) {
    # Heuristic: resolution ≈ k / (n_nodes / 10)
    # This is approximate and may need tuning
    resolution <- k / (nrow(data) / 100)
    resolution <- max(0.1, min(resolution, 2.0))  # Bound it
  }

  # Perform Leiden clustering
  set.seed(123)
  leiden_result <- cluster_leiden(g,
                                  resolution = resolution,
                                  objective_function = "modularity")

  clusters <- membership(leiden_result)

  # Check if we got approximately k clusters
  n_clusters_found <- length(unique(clusters))
  if (abs(n_clusters_found - k) > k * 0.3) {
    warning(sprintf("Leiden found %d clusters, target was %d. Consider adjusting resolution parameter.",
                    n_clusters_found, k))
  }

  return(list(
    clusters = clusters,
    model = leiden_result,
    n_clusters_found = n_clusters_found,
    resolution = resolution,
    method = "leiden"
  ))
}


#' Perform clustering with specified method
#'
#' Unified interface for different clustering methods.
#'
#' @param data Cell-level feature matrix
#' @param k Number of clusters
#' @param method Clustering method: "hierarchical", "kmeans", "flowsom", "leiden"
#' @param ... Additional parameters passed to specific clustering function
#' @return List with cluster assignments and model object
perform_clustering <- function(data, k, method = "hierarchical", ...) {
  if (method == "hierarchical") {
    return(hierarchical_clustering(data, k, ...))
  } else if (method == "kmeans") {
    return(kmeans_clustering(data, k, ...))
  } else if (method == "flowsom") {
    return(flowsom_clustering(data, k, ...))
  } else if (method == "leiden") {
    return(leiden_clustering(data, k, ...))
  } else {
    stop("Unknown clustering method: ", method)
  }
}

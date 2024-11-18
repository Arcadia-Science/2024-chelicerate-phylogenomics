# Load required packages and source custom functions used herein.
require(igraph)
require(uwot)
require(dbscan)
require(ggrepel)

source("./scripts/phylo_multivariate_distance_functions.R")
source("./scripts/misc_phylo_profiling_functions.R")

# The following function calls a number of custom functions internally and runs
# the entire phylogenetic profiling pipeline from scratch.
# Starting from a dataframe of gene family evolutionary event counts (i.e.
# speciations, duplications, transfers, or losses) and a species tree, the
# pipeline carries out the following steps:
# 1) It conducts a phylogenetic GLS transformation on the count data,
#    effectively returning the residual variation after accounting for the
#    variation explained by shared evolutionary history (phylogeny) alone
# 2) It calculates pairwise distances among all gene families using either
#    mahalanobis distances (using the transformed count data) or euclidean
#    distances (using principal components obtained from the transformed data)
# 3) It then identifies, for each gene family the K-nearest neighbors. It does
#    so using a K of 15 for use in visualization in umap embedding space, as
#    well as the user-defined K to be used to construct a reduced network later
#    used for network-based cluster inference (community detection)
# 4) Using this KNN network, it then runs one or several community detection
#    algorithms, clustering gene families based on their phylogenetic profiles
# 5) Gene families and the corresponding gene family clusters are then
#    visualized in a UMAP embedding space using both a static plot (for the
#    best-performing clustering algorithm based on modularity) or an
#    interactive plot that includes visualization of the results of all
#    clustering algorithms
run_phylo_profiling <-
  function(data = NULL,
           tree = NULL,
           cluster_method = c(
             "none",
             "all",
             "leiden",
             "infomap",
             "walktrap-2",
             "walktrap-4",
             "walktrap-6"
           ),
           k = 5,
           distance_metric = "mahalanobis",
           use_anc_nodes = TRUE,
           profile_title = NULL,
           annotations = NULL,
           seed = 1234) {
    # Set the seed for reproducibility:
    set.seed(seed)

    # Set the names/methods of different clustering methods to be accessed
    # later on
    cluster_names <- c(
      "Infomap Cluster",
      "Walktrap-2 Cluster",
      "Walktrap-4 Cluster",
      "Walktrap-6 Cluster",
      "Leiden Cluster"
    )
    cluster_methods <- c(
      "infomap", "walktrap2",
      "walktrap4", "walktrap6", "leiden"
    )

    # Conduct the phylogenetic GLS transformation, returning the residual trait
    # variation after accounting for and removing the phylogenetic effect
    transf_data <-
      phylo_gls_transform(data, tree, include_internal = use_anc_nodes)
    data <- transf_data
    pcs <- prcomp(transf_data)$x

    if (distance_metric == "mahalanobis") {
      # Compute Mahalanobis distances in parallel
      print("Calculating mahalanobis distances using transformed data...")
      distance_matrix <-
        pairwise_phylo_mahalanobis(transf_data, tree, include_internal = use_anc_nodes)
      umap_dat <- transf_data
    } else {
      print(
        "Calculating euclidean distances using principal components obtained from transformed data..."
      )
      distance_matrix <- as.matrix(dist(pcs))
      umap_dat <- pcs
    }

    # If we choose to cluster gene families based on their phylogenetic profiles
    # (the primary utility of this pipeline outside of simply calculating
    # distances among gene families), we then will construct KNN networks,
    # run user-specified network-based community detection algorithms, and
    # then visualize these gene family clusters in a UMAP embedding space.
    if (cluster_method != "none") {
      print(paste0("Determining the ", k, " nearest-neighbors..."))
      # Note we will construct two adjacency matrices: one for
      # clustering, using the user-specified value of k, and one using,
      # a k of 15 for use with UMAP
      # Calculate K nearest neighbors for each node:
      # Here we use k, the number of nearest neighbors defined in the
      # function call
      nearest_neighbors_umap <-
        dbscan::kNN(as.dist(distance_matrix),
          k = 15,
          search = "linear"
        )
      sparse_dists_umap <- convert_to_dgCMatrix(
        cbind(
          1:nrow(nearest_neighbors_umap$id),
          nearest_neighbors_umap$id
        ),
        cbind(
          rep(
            1e-10, nrow(nearest_neighbors_umap$dist)
          ),
          nearest_neighbors_umap$dist
        )
      )
      knn_adj_matrix <-
        convert_to_dgCMatrix(
          nearest_neighbors_umap$id[, 1:k],
          exp(-nearest_neighbors_umap$dist[, 1:k])
        )

      # Remove isolated (singleton) components, as this will just pose problems
      # down the line
      n_components <-
        igraph::components(igraph::graph_from_adjacency_matrix(sparse_dists_umap))
      singletons <-
        which(n_components$membership %in% which(n_components$csize == 1))
      sparse_dists_umap <-
        sparse_dists_umap[-singletons, -singletons]
      knn_adj_matrix <- knn_adj_matrix[-singletons, -singletons]
      distance_matrix <- distance_matrix[-singletons, -singletons]
      annotations <- annotations[-singletons, ]

      # Construct the igraph network from the adjacency matrix we contructed
      # using the user-specified value of K
      g <-
        igraph::graph_from_adjacency_matrix(knn_adj_matrix, weighted = T, mode = "undirected")

      # Now, if clustering was requested while calling the function, go ahead
      # and run the algorithm(s)
      if (cluster_method == "all") {
        # run all algorithms and compare performance based on modularity
        print("Performing infomap clustering...")
        infomap_parts <- igraph::cluster_infomap(g, nb.trials = 20)
        print("Performing walktrap clustering with a step size of 2...")
        walk2_parts <- igraph::cluster_walktrap(g, steps = 2)
        print("Performing walktrap clustering with a step size of 4...")
        walk4_parts <- igraph::cluster_walktrap(g, steps = 4)
        print("Performing walktrap clustering with a step size of 6...")
        walk6_parts <- igraph::cluster_walktrap(g, steps = 6)
        print("Performing leiden clustering...") # Test multiple different resolution parameters
        leiden_parts <- list(
          igraph::cluster_leiden(
            g,
            objective_function = "modularity",
            n_iterations = 50,
            resolution_parameter = 0.25
          ),
          igraph::cluster_leiden(
            g,
            objective_function = "modularity",
            n_iterations = 50,
            resolution_parameter = 0.5
          ),
          igraph::cluster_leiden(
            g,
            objective_function = "modularity",
            n_iterations = 50,
            resolution_parameter = 0.75
          ),
          igraph::cluster_leiden(
            g,
            objective_function = "modularity",
            n_iterations = 50,
            resolution_parameter = 1
          ),
          igraph::cluster_leiden(
            g,
            objective_function = "modularity",
            n_iterations = 50,
            resolution_parameter = 1.25
          ),
          igraph::cluster_leiden(
            g,
            objective_function = "modularity",
            n_iterations = 50,
            resolution_parameter = 1.5
          ),
          igraph::cluster_leiden(
            g,
            objective_function = "modularity",
            n_iterations = 50,
            resolution_parameter = 1.75
          )
        )
        leiden_mods <- c(
          igraph::modularity(g, leiden_parts[[1]]$membership),
          igraph::modularity(g, leiden_parts[[2]]$membership),
          igraph::modularity(g, leiden_parts[[3]]$membership),
          igraph::modularity(g, leiden_parts[[4]]$membership),
          igraph::modularity(g, leiden_parts[[5]]$membership),
          igraph::modularity(g, leiden_parts[[6]]$membership),
          igraph::modularity(g, leiden_parts[[7]]$membership)
        )
        leiden_parts <-
          leiden_parts[[which(leiden_mods == max(leiden_mods))]]

        # Plot modularity for each algorithm, so we can better understand their
        # relative performance
        cluster_modularity <-
          data.frame(
            clustering_method = cluster_methods,
            modularity = c(
              igraph::modularity(g, infomap_parts$membership),
              igraph::modularity(g, walk2_parts$membership),
              igraph::modularity(g, walk4_parts$membership),
              igraph::modularity(g, walk6_parts$membership),
              igraph::modularity(g, leiden_parts$membership)
            )
          )
        cluster_modularity_plt <-
          ggplot(data = cluster_modularity, aes(x = clustering_method, y = modularity)) +
          xlab("Community Detection Algorithm") +
          ylab("Modularity") +
          geom_point(size = 3) +
          geom_line(aes(group = 1)) +
          theme_classic(base_size = 16)

        # Produce a dataframe of gene family cluster memberships for each algorithm
        cluster_ids <- list(
          "infomap_cluster_id" = infomap_parts$membership,
          "walktrap2_cluster_id" = walk2_parts$membership,
          "walktrap4_cluster_id" = walk4_parts$membership,
          "walktrap6_cluster_id" = walk6_parts$membership,
          "leiden_cluster_id" = leiden_parts$membership
        )

        # And order these by performance based on modularity
        cluster_orders <-
          order(cluster_modularity$modularity, decreasing = T)
        cluster_methods <- cluster_methods[cluster_orders]
        cluster_names <- cluster_names[cluster_orders]
        cluster_ids <- cluster_ids[cluster_orders]
        best_cluster_name <- cluster_names[1]

        # and construct a list to contain the full set of network objects
        final_networks <- list(
          infomap = infomap_parts,
          walktrap2 = walk2_parts,
          walktrap4 = walk4_parts,
          walktrap6 = walk6_parts,
          leiden = leiden_parts
        )
        # If only one algorithm was requested, run those....
      } else if (cluster_method == "infomap") {
        print("Performing infomap clustering...")
        final_networks <- igraph::cluster_infomap(g)
        cluster_ids <- final_networks$membership
        best_cluster_name <- "Infomap Cluster"
      } else if (cluster_method == "walktrap-2") {
        print("Performing walktrap clustering with a step size of 2...")
        final_networks <- igraph::cluster_walktrap(g, steps = 2)
        cluster_ids <- final_networks$membership
        best_cluster_name <- "Walktrap (2) Cluster"
      } else if (cluster_method == "walktrap-4") {
        print("Performing walktrap clustering with a step size of 4...")
        final_networks <- igraph::cluster_walktrap(g, steps = 4)
        cluster_ids <- final_networks$membership
        best_cluster_name <- "Walktrap (4) Cluster"
      } else if (cluster_method == "walktrap-6") {
        print("Performing walktrap clustering with a step size of 6...")
        final_networks <- igraph::cluster_walktrap(g, steps = 6)
        cluster_ids <- final_networks$membership
        best_cluster_name <- "Walktrap (6) Cluster"
      } else {
        print("Performing leiden clustering...")
        final_networks <- list(
          igraph::cluster_leiden(
            g,
            objective_function = "modularity",
            n_iterations = 50,
            resolution_parameter = 0.25
          ),
          igraph::cluster_leiden(
            g,
            objective_function = "modularity",
            n_iterations = 50,
            resolution_parameter = 0.5
          ),
          igraph::cluster_leiden(
            g,
            objective_function = "modularity",
            n_iterations = 50,
            resolution_parameter = 0.75
          ),
          igraph::cluster_leiden(
            g,
            objective_function = "modularity",
            n_iterations = 50,
            resolution_parameter = 1
          ),
          igraph::cluster_leiden(
            g,
            objective_function = "modularity",
            n_iterations = 50,
            resolution_parameter = 1.25
          ),
          igraph::cluster_leiden(
            g,
            objective_function = "modularity",
            n_iterations = 50,
            resolution_parameter = 1.5
          ),
          igraph::cluster_leiden(
            g,
            objective_function = "modularity",
            n_iterations = 50,
            resolution_parameter = 1.75
          )
        )
        leiden_mods <- c(
          igraph::modularity(g, final_networks[[1]]$membership),
          igraph::modularity(g, final_networks[[2]]$membership),
          igraph::modularity(g, final_networks[[3]]$membership),
          igraph::modularity(g, final_networks[[4]]$membership),
          igraph::modularity(g, final_networks[[5]]$membership),
          igraph::modularity(g, final_networks[[6]]$membership),
          igraph::modularity(g, final_networks[[7]]$membership)
        )
        final_networks <-
          final_networks[[which(leiden_mods == max(leiden_mods))]]
        cluster_ids <- final_networks$membership
        best_cluster_name <- "Leiden Cluster"
      }

      print("Obtaining UMAP projections...")
      if (cluster_method == "all") {
        # identify singleton clusters so we can exclude them:
        keep_clusts <-
          lapply(
            seq_along(cluster_ids),
            function(i) {
              table(cluster_ids[[i]])[which(table(cluster_ids[[i]]) > 1)]
            }
          )

        # determine plotting colors (for each clustering algorithm) and prepare for plotting:
        cols <- list()
        unique_clusters <- list()
        color_mappings <- list()
        for (i in 1:length(cluster_ids)) {
          cols[[i]] <-
            arcadia_color_discrete(accent_v4, length(keep_clusts[[i]]))

          names(cols[[i]]) <-
            names(keep_clusts[[i]][order(keep_clusts[[i]], decreasing = T)])
          # Get unique cluster IDs
          unique_clusters[[i]] <- names(keep_clusts[[i]])
          # Map unique cluster IDs to colors
          color_mappings[[i]] <- cols[[i]]
        }
        best_cols <- cols[[1]]
        # Get unique cluster IDs
        unique_best_clusters <- unique_clusters[[1]]
        # Map unique cluster IDs to colors
        best_color_mapping <- color_mappings[[1]]
      } else {
        # identify singleton clusters so we can exclude them:
        keep_clusts <-
          table(cluster_ids)[which(table(cluster_ids) > 1)]

        # determine plotting colors:
        best_cols <-
          arcadia_color_discrete(accent_v4, length(keep_clusts))
        names(best_cols) <-
          names(keep_clusts[order(keep_clusts, decreasing = T)])
        # Get unique cluster IDs
        unique_best_clusters <- names(keep_clusts)
        # Map unique cluster IDs to colors
        color_mappings <- list()
        color_mappings[[1]] <- best_cols
      }

      # Based on exploratory use of different approaches, use uwot's
      # t-Distributed UMAP with a slightly higher value for local connectivity
      # (here, using a local neighborhood size of 5)
      umap_results <-
        uwot::tumap(
          as.dist(distance_matrix),
          n_threads = detectCores(),
          local_connectivity = 5,
          verbose = T
        )
      umap_results_axes <- data.frame(umap_results[, 1:2])

      # Combine cluster IDs with the UMAP results
      if (cluster_method == "all") {
        umap_results_axes <-
          cbind(umap_results_axes, do.call(cbind, cluster_ids))
        umap_results_axes[, -c(1:2)] <-
          apply(umap_results_axes[, -c(1:2)], 2, as.factor)
      } else {
        umap_results_axes <- cbind(umap_results_axes, cluster_ids)
        colnames(umap_results_axes)[3] <-
          paste0(cluster_method, "_cluster_id")
        umap_results_axes[, 3] <- as.factor(umap_results_axes[, 3])
      }

      print("Plotting: best-performing algorithm...")
      # Temporarily drop singleton clusters for plotting
      umap_results_axes_final <- umap_results_axes
      keep_clusts <-
        names(table(umap_results_axes[, 3]))[which(table(umap_results_axes[, 3]) > 1)]
      umap_results_axes <-
        umap_results_axes[which(umap_results_axes[, 3] %in% keep_clusts), ]

      # Plot in 2D using the best-performing algorithm
      umap_results_axes$best_cluster_ids <-
        as.factor(as.numeric(umap_results_axes[, 3]))
      umap_plt <-
        ggplot(
          umap_results_axes,
          aes(x = X1, y = X2, fill = best_cluster_ids)
        ) +
        scale_fill_manual(values = best_cols) +
        geom_point(
          pch = 21,
          stroke = 0.25,
          alpha = 0.85,
          size = 1.5,
          color = "DarkSlateGrey"
        ) +
        guides(fill = guide_legend(title = element_text(best_cluster_name))) +
        theme_bw()

      # And visualize interactively with plotly, labeling by cluster and
      # including provided annotations.
      print("Plotting: best-performing algorithm...")
      # First, combine the annotations with the umap layout/cluster ids
      umap_results_axes_final$protein_name <-
        gsub("\\s*\\([^()]*\\)$", "", annotations[, 2][[1]])

      # Return the umap_result_axes back to its original state
      # we'll do this each time we add a new plot to the dropdown
      umap_results_axes <- umap_results_axes_final

      # Initialize the plotly object
      umap_pltly <- plot_ly()
      nclusts <- c()
      # Loop through each clustering method
      for (method_idx in seq_along(cluster_methods)) {
        method <- cluster_methods[method_idx]
        color_map <- color_mappings[[method_idx]]

        # Only plot non-singleton clusters:
        clust_freqs <-
          table(umap_results_axes[[paste0(method, "_cluster_id")]])
        keep_clusts <- names(clust_freqs[which(clust_freqs > 5)])
        nclusts[method_idx] <- length(keep_clusts)
        umap_results_axes <-
          umap_results_axes[which(umap_results_axes[[paste0(method, "_cluster_id")]] %in% keep_clusts), ]

        # Convert cluster IDs to a factor for consistent coloring
        cluster_factor <-
          as.factor(umap_results_axes[[paste0(method, "_cluster_id")]])
        cluster_colors <- color_map[as.character(cluster_factor)]

        # Create hover text for each point
        hover_text <-
          paste(
            "Gene Family ID: ",
            rownames(umap_results_axes),
            "<br>",
            paste0(
              str_to_title(method), " Cluster ",
              umap_results_axes[[paste0(method, "_cluster_id")]]
            ),
            "<br>",
            paste0("Prot. Name: ", umap_results_axes$protein_name)
          )

        # Add a trace for the current method
        umap_pltly <- umap_pltly %>%
          add_trace(
            data = umap_results_axes,
            x = ~X1,
            y = ~X2,
            type = "scatter",
            mode = "markers",
            marker = list(
              size = 4,
              line = list(width = 0.25, color = "DarkSlateGrey"),
              opacity = 1,
              color = cluster_colors
            ),
            text = hover_text,
            hoverinfo = "text",
            name = method,
            showlegend = FALSE,
            # Hide legend to avoid overcrowding; we're focusing on the dropdown
            visible = (method_idx == 1) # Only the first method's trace is visible initially
          )
        umap_results_axes <- umap_results_axes_final
      }

      # Create dropdown menu items to toggle visibility
      visibility_lists <-
        lapply(seq_along(cluster_methods), function(idx) {
          visibility_vector <- rep(FALSE, length(cluster_methods))
          visibility_vector[idx] <- TRUE

          list(
            method = "restyle",
            args = list("visible", visibility_vector),
            label = paste0(
              gsub("Cluster", "Clusters", cluster_names[idx]),
              ": N = ",
              nclusts[idx]
            )
          )
        })

      # Finalize the plot layout
      umap_pltly <- umap_pltly %>%
        layout(
          title = paste0("UMAP: ", profile_title, " Profiles"),
          xaxis = list(title = "UMAP Axis 1"),
          yaxis = list(title = "UMAP Axis 2"),
          updatemenus = list(
            list(
              y = 1.075,
              xref = "paper",
              yref = "paper",
              yanchor = "top",
              xanchor = "left",
              active = 0,
              buttons = visibility_lists
            )
          )
        )

      print("Pipeline finished!")
      return(
        list(
          graph = g,
          phylo_corrected_data = data,
          umap_clusters = umap_results_axes,
          static_umap = umap_plt,
          interactive_umap = umap_pltly,
          cluster_alg_modularity = cluster_modularity_plt
        )
      )
    } else {
      print("Plotting...")
      # No clusters were inferred, so simply obtain the t-UMAP embeddings and
      # plot statically.
      umap_results <-
        tumap(
          as.dist(distance_matrix),
          n_threads = detectCores(),
          local_connectivity = 5,
          verbose = T
        )
      umap_results_axes <- data.frame(umap_results[, 1:2])

      umap_plt <-
        ggplot(
          umap_results_axes,
          aes(x = X1, y = X2)
        ) +
        geom_point(alpha = 0.5, size = 0.5) +
        theme_bw()

      print("Pipeline finished!")
      return(
        list(
          graph = g,
          phylo_corrected_data = data,
          umap_layout = umap_results_axes,
          static_umap = umap_plt
        )
      )
    }
  }

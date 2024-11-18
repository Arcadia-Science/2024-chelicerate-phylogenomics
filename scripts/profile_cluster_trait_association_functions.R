require(phyr)
require(logistf)

# Define a function to calculate event means for these species for each gene
# family, and assess how strongly the mean event counts within each gene family
# cluster correlates with a binary trait.
cluster_assoc_test <-
  function(count_dat,
           cluster_res,
           suppressors,
           nonsuppressors,
           method = c("rank", "ols", "logistic"),
           min_clust_size = 5) {
    if (length(method) == 2) {
      method <- "rank"
    }
    # Reduce the count data down to only include those for which we have cluster
    # IDs
    cluster_ids <- cluster_res$umap_clusters
    count_dat <-
      count_dat[which(rownames(count_dat) %in% rownames(cluster_ids)), ]
    
    # Calculate means
    suppressor_mean <-
      apply(count_dat[, which(colnames(count_dat) %in%
                                suppressors)], 1, mean)
    nonsuppressor_mean <-
      apply(count_dat[, which(colnames(count_dat) %in%
                                nonsuppressors)], 1, mean)
    clusters <- cluster_ids$leiden_cluster_id
    
    ps <- c()
    coef <- c()
    clust_freqs <- table(clusters)
    cls <-
      sort(as.numeric(names(clust_freqs)[which(clust_freqs >= min_clust_size)]))
    idx <- 1
    if (method == "ols") {
      # Conduct a simple OLS to test for a difference in event counts between
      # suppressing and non-suppressing species.
      for (i in cls) {
        cl_ids <- which(clusters == i)
        dat <-
          data.frame(
            "count" = c(suppressor_mean[cl_ids], nonsuppressor_mean[cl_ids]),
            "detection_suppression" = c(
              rep(1, length(cl_ids)),
              rep(0, length(cl_ids))
            )
          )
        
        fit <- lm(count ~ detection_suppression, data = dat)
        ps[idx] <- summary(fit)$coefficients[8]
        coef[idx] <- summary(fit)$coefficients[2]
        idx <- idx + 1
      }
    } else if (method == "logistic") {
      # Conduct a logistic regression to identify which gene families have
      # event counts distributions that predict whether a species is
      # detection-suppressing or not.
      for (i in cls) {
        cl_ids <- which(clusters == i)
        dat <-
          data.frame(
            "count" = c(suppressor_mean[cl_ids], nonsuppressor_mean[cl_ids]),
            "detection_suppression" = c(
              rep(1, length(cl_ids)),
              rep(0, length(cl_ids))
            )
          )
        # Use Firth's logistic regression to reduce the bias in the ML
        # estimates that arise due to separation by adding a penalty to the
        # likelihood, preventing infinite estimates.
        fit <-
          suppressWarnings(logistf(
            detection_suppression ~ count,
            data = dat,
            family = "binomial"
          ))
        log <- capture.output(fit <- summary(fit))
        ps[idx] <- fit$prob[[2]]
        coef[idx] <- fit$coefficients[[2]]
        idx <- idx + 1
      }
    } else {
      # Conduct a simple wilcox rank sum test to see if counts in one population
      # differs from that in the other
      for (i in cls) {
        cl_ids <- which(clusters == i)
        test_res <-
          wilcox.test(
            suppressor_mean[cl_ids],
            nonsuppressor_mean[cl_ids],
            exact = FALSE,
            conf.int = TRUE
          )
        ps[i] <- test_res$p.value
        coef[i] <- test_res$estimate[[1]]
        idx <- idx + 1
      }
    }
    return(list(
      pvals = ps,
      coefficients = coef,
      clusters = cls
    ))
  }

# Define a function to take results from the above function, and plot the
# results, returning a table of summarized results, as well as an intuitive
# figure that indicates which gene families are the strongest predictors,
# ordering them by their respective coefficients.
plot_cluster_traitcorr <-
  function(pvals,
           coefficients,
           profile_type,
           cluster_data,
           cluster_ids,
           data_type = c("genefam_clusters", "gene_families"),
           method = c("rank", "ols", "logistic")) {
    if (length(method) == 2) {
      method <- "rank"
    }
    if (method %in% c("rank", "ols")) {
      ylab_string <- "Mean-shift by\nDetection-suppression"
    } else {
      ylab_string <- "Log-Odds: Detection Suppression"
    }
    if (length(data_type) == 2) {
      data_type <- "genefam_clusters"
      xlab_string <-
        paste(str_to_title(profile_type), "Profile Cluster")
    } else {
      xlab_string <- paste("Gene Family:", str_to_title(profile_type))
    }
    # Remove any gene families for which we could not conduct tests:
    to_remove <- which(is.na(coefficients))
    if (length(to_remove) > 0) {
      pvals <- pvals[-to_remove]
      coefficients <- coefficients[-to_remove]
    }
    res <-
      data.frame(
        x = 1:length(pvals),
        profile_type = str_to_lower(profile_type),
        cluster_data = str_to_lower(cluster_data),
        cluster_id = cluster_ids[order(coefficients)],
        coefficient = coefficients[order(coefficients)],
        pval = pvals[order(coefficients)],
        qval = p.adjust(pvals[order(coefficients)], method = "fdr")
      )
    extremes10_pos <-
      quantile(res$coefficient[which(res$pval <= 0.05 &
                                       res$coefficient >= 0)], c(0.9))
    extremes_pos <-
      quantile(res$coefficient[which(res$pval <= 0.05 &
                                       res$coefficient >= 0)], c(0.95))
    extremes10_neg <-
      quantile(res$coefficient[which(res$pval <= 0.05 &
                                       res$coefficient <= 0)], c(0.1))
    extremes_neg <-
      quantile(res$coefficient[which(res$pval <= 0.05 &
                                       res$coefficient <= 0)], c(0.05))
    extremes10 <- c(extremes10_neg, extremes10_pos)
    extremes <- c(extremes_neg, extremes_pos)
    pval_breaks <- c(0, 0.01, 0.05)
    
    # Add a column indicating the significance level
    res$signif_level <- "No"
    res[which(res$pval <= 0.05), "signif_level"] <- "Yes"
    res[which(
      res$pval <= 0.05 &
        res$coefficient <= extremes10[1] |
        res$pval <= 0.05 &
        res$coefficient >= extremes10[2]
    ), "signif_level"] <- "Yes & Top 10% Coef."
    res[which(
      res$pval <= 0.05 &
        res$coefficient <= extremes[1] |
        res$pval <= 0.05 &
        res$coefficient >= extremes[2]
    ), "signif_level"] <- "Yes & Top 5% Coef."
    res$signif_level <-
      factor(res$signif_level,
             levels = c("No", "Yes", "Yes & Top 10% Coef.", "Yes & Top 5% Coef.")
      )
    res$signif_fdr <- "No"
    res$signif_fdr[which(res$qval <= 0.05)] <- "Yes"
    res$signif_fdr <-
      factor(res$signif_fdr, levels = c("No", "Yes"))
    
    # Label the clusters with the ten most extreme positive and negative
    # coeffients
    sig_clusts <-
      res[res$pval <= 0.05, ]
    label_neg <-
      which(sig_clusts$coefficient <= 0)[order(sig_clusts$coefficient[sig_clusts$coefficient <= 0], decreasing = FALSE)][1:10]
    label_pos <-
      which(sig_clusts$coefficient >= 0)[order(sig_clusts$coefficient[sig_clusts$coefficient >= 0], decreasing = TRUE)][1:10]
    to_label <- na.omit(c(label_neg, label_pos))
    if (length(to_label) > 0) {
      sig_clusts <- sig_clusts[to_label, ]
    }
    
    # Now plot, highlighting those
    plt <- suppressWarnings(
      ggplot(
        res,
        aes(
          x = x,
          y = coefficient,
          size = signif_level,
          alpha = signif_level,
          color = signif_level,
          shape = signif_fdr
        )
      ) +
        geom_point() +
        geom_text_repel(
          data = sig_clusts,
          aes(label = cluster_id, size = NULL),
          size = 3,
          min.segment.length = 0.05,
          box.padding = 1,
          nudge_x = 5,
          show.legend = FALSE,
          max.overlaps = 30
        ) +
        scale_alpha_discrete(
          range = c(
            "No" = 0.75,
            "Yes" = 1,
            "Yes & Top 10% Coef." = 1,
            "Yes & Top 5% Coef." = 1
          ),
          drop = FALSE
        ) +
        scale_size_discrete(
          range = c(
            "No" = 1.5,
            "Yes" = 2,
            "Yes & Top 10% Coef." = 2.5,
            "Yes & Top 5% Coef." = 3
          ),
          drop = FALSE
        ) +
        scale_shape_manual(
          values = c("No" = 16, "Yes" = 15),
          drop = FALSE
        ) +
        scale_color_manual(
          values = c(
            "No" = "#292928",
            "Yes" = "#5088C5",
            "Yes & Top 10% Coef." = "#F28360",
            "Yes & Top 5% Coef." = "#C85152"
          ),
          drop = FALSE
        ) +
        geom_hline(yintercept = 0, linetype = "solid") +
        theme_classic(base_size = 14) +
        guides(
          size = FALSE,
          alpha = FALSE,
          color = guide_legend(
            title = expression(P <= 0.05),
            nrow = 2,
            title.position = "top"
          ),
          shape = guide_legend(
            title = expression(FDR <= 0.05),
            nrow = 2,
            title.position = "top"
          )
        ) +
        xlab(xlab_string) +
        ylab(ylab_string) +
        theme(legend.position = "top")
    )
    
    # Finally, remove the first column which was just used for plotting
    res <- res[, -1]
    if (data_type == "gene_families") {
      colnames(res)[3] <- "gene_family_id"
    }
    return(list(
      correlation_summary = res,
      correlation_plot = plt
    ))
  }

# A function to assess, for each cluster that strongly predicts a binary trait,
# which constitutent gene families are the strongest
# predictors
test_fam_corrs_per_clust <-
  function(count_dat,
           clusters,
           phylo_profiling_res,
           spptree) {
    counter <- 1
    clstr_plts <- list()
    for (clstr in clusters$cluster_id) {
      print(paste0(
        "Processing cluster ",
        clstr,
        ": ",
        counter,
        " of ",
        length(clusters$cluster_id)
      ))
      clust_ogs <-
        rownames(phylo_profiling_res$umap_clusters)[which(phylo_profiling_res$umap_clusters$leiden_cluster_id == clstr)]
      clstr_dat <-
        count_dat[which(rownames(count_dat) %in% clust_ogs), ]
      
      # Reformat and add the trait:
      clstr_dat <- as.data.frame(t(clstr_dat))
      # If there are any columns where the predictor is invariant (i.e. all 0), remove
      invar <- which(colSums(clstr_dat) == 0)
      if (length(invar) > 0) {
        clstr_dat <- clstr_dat[, -invar]
      }
      clstr_dat$detection_suppression <- 0
      clstr_dat$detection_suppression[which(rownames(clstr_dat) %in% suppressors)] <-
        1
      
      # Now, fit phylogenetic logistic regressions for each gene family, testing whether
      # they are significant predictors
      # Prepare the data:
      clstr_dat <-
        clstr_dat[match(spptree$tip.label, rownames(clstr_dat)), ]
      
      # Fit the models
      clust_res <-
        pbapply::pblapply(1:c(ncol(clstr_dat) - 1), function(i) {
          newdat <- clstr_dat[, c(i, ncol(clstr_dat))]
          colnames(newdat)[1] <- "gf"
          fit <-
            summary(
              phylolm::phyloglm(
                detection_suppression ~ gf,
                phy = spptree,
                data = newdat,
                method = "logistic_IG10",
                btol = 10
              )
            )
          coefs <- fit$coefficients[[2]]
          coef_ps <- fit$coefficients[[8]]
          return(c(coefs, coef_ps))
        }, cl = detectCores())
      clust_res <- do.call(rbind, clust_res)
      
      clstr_plts[[counter]] <-
        plot_cluster_traitcorr(
          pvals = clust_res[, 2],
          coefficients = clust_res[, 1],
          profile_type = "speciation",
          cluster_data = "combined",
          cluster_ids = colnames(clstr_dat[, -ncol(clstr_dat)]),
          data_type = "gene_families",
          method = "logistic"
        )
      clstr_plts[[counter]]$correlation_plot <-
        clstr_plts[[counter]]$correlation_plot +
        ggtitle(paste0("Cluster ", clstr)) +
        theme(
          plot.title = element_text(hjust = 0, vjust = -7),
          legend.title = element_text(hjust = 0),
          plot.title.position = "plot",
          plot.margin = margin(0, 1, 1, 1, "cm")
        )

      clstr_plts[[counter]]$cluster_data <- clstr_dat
      names(clstr_plts)[counter] <- paste0("Cluster_", clstr)
      counter <- counter + 1
    }
    return(clstr_plts)
  }

# And a function to save the results that are returned by the above function
save_per_fam_associations <-
  function(per_fam_results = NULL,
           output_dir = "./",
           event_type = c("speciation", "transfer", "loss")) {
    if (is.null(per_fam_results)) {
      print("ERROR: No result provided! Stopping...")
    }
    if (length(event_type) > 1 & !is.null(per_fam_results)) {
      event_type <- "speciation"
    }
    
    # Save these results out to file. Create a sub-directory for each event type,
    # since we will be saving a number of plots and tables for each.
    if (event_type == "speciation") {
      output_dir <- paste0(output_dir, "speciation_associations/")
      dir.create(output_dir,
                 recursive = TRUE,
                 showWarnings = FALSE
      )
    } else if (event_type == "transfer") {
      output_dir <- paste0(output_dir, "transfer_associations/")
      dir.create(output_dir,
                 recursive = TRUE,
                 showWarnings = FALSE
      )
    } else {
      output_dir <- paste0(output_dir, "loss_associations/")
      dir.create(output_dir,
                 recursive = TRUE,
                 showWarnings = FALSE
      )
    }

    # Loop through significant clusters and save outputs.
    for (clstr in 1:length(ls(per_fam_results))) {
      clstr_name <- str_to_lower(names(per_fam_results)[clstr])
      # Save the per-family correlation summaries
      write.table(
        per_fam_results[[clstr]]$correlation_summary,
        file = paste0(
          output_dir,
          clstr_name,
          "_per_family_correlation_results.tsv"
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE
      )
      # Save the data (species event counts per-cluster)
      write.table(
        per_fam_results[[clstr]]$cluster_data,
        file = paste0(
          output_dir,
          clstr_name,
          "_per_species_speciation_counts.tsv"
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE
      )
      # And save the plot visualizing the per-family correlation test results
      ggsave(
        per_fam_results[[clstr]]$correlation_plot,
        file = paste0(
          output_dir,
          clstr_name,
          "_per_family_correlation_results.pdf"
        ),
        height = 6,
        width = 7
      )
    }
  }

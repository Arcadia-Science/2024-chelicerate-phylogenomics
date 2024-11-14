# This contains a series of functions used to summarize outputs produced from
# running the NovelTree workflow
require(tidyverse) # Version 2.0.0
require(ggtree) # Version 3.6.2
require(ggnewscale) # Version 0.4.8
require(parallel) # Version 4.2.1
require(reshape2) # Version 1.4.4
require(cowplot) # Version 1.1.1
require(phytools) # Version 1.5.1
require(ape) # Version 5.7.1
library(ggforce) # Version 0.4.1
library(plotly) # Version 4.10.1
library(cogeqc) # Version 1.2.1

####################################
# Function to read busco results
read_busco_summary <-
  function(busco_summary = NULL) {
    classes <-
      c("Complete_SC",
        "Complete_duplicate",
        "Fragmented",
        "Missing")
    busco_res <-
      read_table(busco_summary,
                 progress = F,
                 show_col_types = F)[,-c(3, 8)]
    busco_res[1, 1] <- gsub(".fasta", "", busco_res[1, 1])
    colnames(busco_res) <-
      c(
        "Species",
        "Lineage",
        "Complete_SC",
        "Complete_duplicate",
        "Fragmented",
        "Missing"
      )
    busco_res <-
      reshape2::melt(busco_res, id = c("Species", "Lineage"))
    colnames(busco_res) <-
      c("File", "Lineage", "Class", "Frequency")
    busco_res <-
      busco_res[, c("Class", "Frequency", "Lineage", "File")]
    busco_res$Class <- factor(busco_res$Class, levels = classes)
    return(busco_res)
  }

# Function to calculate the bootstrapped support values from
# a distribution of trees, given a user-specified fixed tree.
# This will just be used to recalculate BS support values for
# the re-rooted asteroid tree.
calc_bs_support <- function(bs_trees = NULL,
                            species_tree = NULL) {
  # Compute bootstrap support
  support <- prop.clades(species_tree, bs_trees)
  
  # Add support to the tree
  species_tree$node.label <-
    round(support, 2) # Multiply by 100 to convert to percentages
  
  return(species_tree)
}

# Function to plot the inferred species trees
plot_spp_tree <-
  function(tree = NULL,
           tree_fpath = NULL,
           samplesheet = NULL,
           samplesheet_fpath = NULL,
           nodelab_name = "Bootstrap",
           grp_name = "Group",
           tip_grp_cols = NULL,
           outgroup = NULL,
           cladogram = F) {
    if (is.null(tree)) {
      tree <- read.tree(tree_fpath)
    }
    if (is.null(tree) & is.null(tree_fpath)) {
      print("Error! Provide tree object or filepath to tree")
    }
    if (is.null(samplesheet)) {
      samplesheet <- read.delim(samplesheet_fpath, sep = ",")
    }
    if (is.null(samplesheet) & is.null(samplesheet_fpath)) {
      print("Error! Provide samplesheet object or filepath to samplesheet")
    }
    
    # If plotting the asteroid tree, we need to do some light reformatting of
    # names that were induced through the use of disco to decompose multi-copy
    # gene family trees (replacement of underscores with hyphens in species names).
    tree$tip.label[order(tree$tip.label)] <-
      samplesheet$species[order(samplesheet$species)]
    
    # If an outgroup is provided, root using those taxa.
    if (!is.null(outgroup)) {
      mrca <- getMRCA(tree, outgroup)
      tree <-
        reroot(tree, node = mrca, position = 0.5)
    }
    
    # Now, handle the topological support values
    suppressWarnings(tree$node.label <- as.numeric(tree$node.label))
    if (nodelab_name == "Bootstrap") {
      # If we rerooted the tree, recalculate the support values
      # from the bootstrapped distribution of trees
      if (!is.null(outgroup)) {
        # read in the bootstrapped trees
        workflow_outdir <- gsub("([^/]*/[^/]*)$", "", tree_fpath)
        bs_trees <-
          read.newick(paste0(workflow_outdir, "asteroid/asteroid.bsTrees.newick"))
        tree <- calc_bs_support(bs_trees, tree)
      }
      support_tibble <- tibble(node = 1:Nnode(tree) + Ntip(tree),
                               support = tree$node.label)
      support_cols <- 6
      support_lims <- c(25, 100)
      support_breaks <- c(25, 50, 75, 100)
    } else {
      support_tibble <- tibble(node = 1:Nnode(tree) + Ntip(tree),
                               support = as.numeric(round(tree$node.label, 2)))
      support_cols <- 3
      support_lims <- c(-1, 1)
      support_breaks <- c(-1,-0.5, 0.0, 0.5, 1)
    }
    taxa_tibble <- tibble(tip = samplesheet$species,
                          grp = samplesheet$taxonomy)
    
    
    # Determine if we're plotting as a cladogram or not (e.g. forcing
    # the tree to be ultrametric, ignoring branch lengths)
    if (cladogram == T) {
      bl <- "none"
    } else {
      bl <- "branch.length"
    }
    
    tree_plt <- suppressWarnings(
      ggtree(tree, branch.length = bl) %<+% support_tibble +
        geom_text(
          aes(label = support),
          hjust = -0.5,
          size = 3,
          na.rm = TRUE
        ) +
        geom_nodepoint(
          aes(size = support, fill = support),
          shape = 21,
          color = "black",
          na.rm = TRUE
        ) +
        scale_size_continuous(
          range = c(0.5, 3),
          limits = support_lims,
          breaks = support_breaks
        ) +
        scale_fill_distiller(
          direction = 1,
          palette = support_cols,
          limits = support_lims,
          breaks = support_breaks
        )
    )
    
    tree_plt <- suppressMessages(
      suppressWarnings(
        tree_plt %<+% taxa_tibble +
          guides(
            fill = guide_legend(title = nodelab_name, nrow = 2),
            size = guide_legend(title = nodelab_name, nrow = 2)
          ) +
          new_scale_fill() +
          geom_tiplab(
            aes(label = label),
            size = 10,
            hjust = -.5,
            as_ylab = T
          ) +
          geom_tippoint(
            aes(fill = grp),
            pch = 22,
            size = 3,
            color = "black"
          ) +
          scale_fill_manual(values = tsar_cols) +
          theme(
            legend.position = "top",
            legend.justification = "left"
          ) +
          guides(fill = guide_legend(title = grp_name, nrow = 2))
      )
    )
    
    return(tree_plt)
  }

###########################################################
# OrthoFinder summary functions
###########################################################
# Function to read in and get per-species gene copy number
# for each orthogroup
get_per_spp_og_counts <-
  function(results_dir = NULL,
           out_dir = NULL) {
    orthogroup_dir <-
      list.files(
        path = paste0(results_dir, "orthofinder/complete_dataset/"),
        pattern = "Results_Inflation",
        full.names = T
      )
    
    # read in per-species gene copy number per gene family
    og_counts <-
      read.table(
        paste0(orthogroup_dir, "/Orthogroups/Orthogroups.GeneCount.tsv"),
        header = T,
        check.names = F
      )
    colnames(og_counts) <- gsub("\\..*", "", colnames(og_counts))
    
    # And calculate the number of species in each gene family
    og_counts$NumSpecies <-
      rowSums(og_counts[,-c(1, ncol(og_counts))] > 0)
    
    # Write out to file if an output directory is provided
    if (!is.null(out_dir)) {
      # Create the directory if it doesn't yet exist
      dir.create(out_dir, showWarnings = F, recursive = T)
      # And write out to file
      write.table(
        og_counts,
        file = paste0(out_dir, "og_copy_num_per_spp.tsv"),
        sep = "\t",
        quote = F,
        row.names = F,
        col.names = T
      )
    }
    return(og_counts)
  }

# Function to generally summarize (producing both plots and tables)
# the results of orthofinders "comparative genomics statistics"
get_ofinder_summaries <-
  function(results_dir = NULL,
           show_plots = F,
           out_dir = "summarized-results/orthofinder/") {
    # Store relevant fpaths to variables
    orthogroup_dir <-
      list.files(
        path = paste0(results_dir, "orthofinder/complete_dataset/"),
        pattern = "Results_Inflation",
        full.names = T
      )
    phylohog_dir <-
      paste0(results_dir, "orthofinder/complete_dataset/Results_HOGs/")
    
    # read in hierarchical orthogroup pairwise species overlaps
    hog_spp_overlaps <-
      read.table(
        paste0(
          phylohog_dir,
          "/Comparative_Genomics_Statistics/OrthologuesStats_Totals.tsv"
        ),
        check.names = F
      )
    
    # read in orthogroup (gene family) pairwise species overlaps
    og_spp_overlaps <-
      read.table(
        paste0(
          orthogroup_dir,
          "/Comparative_Genomics_Statistics/Orthogroups_SpeciesOverlaps.tsv"
        ),
        check.names = F
      )
    
    # Read in overall OG stats: these will be reformatted in
    # a more useful way and plotted
    overall_stats_fpath <-
      list.files(
        orthogroup_dir,
        pattern = "Statistics_Overall.tsv",
        recursive = T,
        full.names = T
      )
    og_gene_freqs <-
      read_tsv(
        overall_stats_fpath,
        skip = 26,
        n_max = 20,
        show_col_types = FALSE,
        col_types = c("cnnnn"),
        col_names = c(
          "avg_genes_per_spp_in_og",
          "number_of_ogs",
          "percentage_of_ogs",
          "number_of_genes",
          "percentage_of_genes"
        )
      )
    
    # Exclude rows for which all entries = 0
    og_gene_freqs <-
      og_gene_freqs[which(rowSums(og_gene_freqs[,-1]) > 0),]
    og_gene_freqs[, 1] <-
      factor(og_gene_freqs$avg_genes_per_spp_in_og,
             levels = og_gene_freqs$avg_genes_per_spp_in_og)
    
    # Now begin to plot
    og_gene_hists <- list()
    og_gene_hists[[1]] <- suppressWarnings(
      ggplot(data = og_gene_freqs[-1,],
             aes(x = avg_genes_per_spp_in_og,
                 y = number_of_ogs)) +
        geom_histogram(
          stat = "identity",
          color = "black",
          fill = "grey"
        ) +
        theme_classic(base_size = 12) +
        theme(axis.text.x.bottom =
                element_text(
                  angle = 45,
                  vjust = 1,
                  hjust = 1
                )) +
        xlab("Avg. # genes-per-spp in GF") +
        ylab("# Families")
    )
    
    og_gene_hists[[2]] <- suppressWarnings(
      ggplot(
        data = og_gene_freqs[-1,],
        aes(x = avg_genes_per_spp_in_og,
            y = percentage_of_ogs)
      ) +
        geom_histogram(
          stat = "identity",
          color = "black",
          fill = "grey"
        ) +
        theme_classic(base_size = 12) +
        theme(axis.text.x.bottom =
                element_text(
                  angle = 45,
                  vjust = 1,
                  hjust = 1
                )) +
        xlab("Avg. # genes-per-spp in GF") +
        ylab("% Families")
    )
    
    og_gene_hists[[3]] <- suppressWarnings(
      ggplot(data = og_gene_freqs[-1,],
             aes(x = avg_genes_per_spp_in_og,
                 y = number_of_genes)) +
        geom_histogram(
          stat = "identity",
          color = "black",
          fill = "grey"
        ) +
        theme_classic(base_size = 12) +
        theme(axis.text.x.bottom =
                element_text(
                  angle = 45,
                  vjust = 1,
                  hjust = 1
                )) +
        xlab("Avg. # genes-per-spp in GF") +
        ylab("# Genes")
    )
    
    og_gene_hists[[4]] <- suppressWarnings(
      ggplot(
        data = og_gene_freqs[-1,],
        aes(x = avg_genes_per_spp_in_og,
            y = percentage_of_genes)
      ) +
        geom_histogram(
          stat = "identity",
          color = "black",
          fill = "grey"
        ) +
        theme_classic(base_size = 12) +
        theme(axis.text.x.bottom =
                element_text(
                  angle = 45,
                  vjust = 1,
                  hjust = 1
                )) +
        xlab("Avg. # genes-per-spp in GF") +
        ylab("% Genes")
    )
    
    # Combine into one
    og_gene_hists_plts <-
      plot_grid(og_gene_hists[[1]],
                og_gene_hists[[2]],
                og_gene_hists[[3]],
                og_gene_hists[[4]],
                ncol = 2)
    
    # print the plots if that is requested
    if (show_plots == T) {
      print(og_gene_hists_plts)
    }
    
    ###################################################
    # Now plot the distribution of orthogroups across #
    # species                                         #
    ###################################################
    og_spp_counts <-
      read_tsv(
        overall_stats_fpath,
        skip = 48,
        show_col_types = FALSE,
        col_types = c("nn"),
        col_names = c("number_of_species_in_og",
                      "number_of_ogs")
      )
    og_spp_hist <- suppressWarnings(
      ggplot(data = og_spp_counts,
             aes(x = number_of_species_in_og,
                 y = number_of_ogs)) +
        geom_histogram(
          stat = "identity",
          color = "black",
          fill = "grey",
          size = 0.5
        ) +
        theme_classic(base_size = 14) +
        xlab("# Species in Gene Family") +
        ylab("# Gene Families")
    )
    
    # And plot if requested
    if (show_plots == T) {
      print(og_spp_hist)
    }
    
    # Write out to file if an output directory is provided
    if (!is.null(out_dir)) {
      # Create the directory if it doesn't yet exist
      dir.create(out_dir, showWarnings = F, recursive = T)
      
      # Write out to file, both the histograms of gene count per og/genes per
      # species per orthogroup
      ggsave(
        og_gene_hists_plts,
        filename = paste0(out_dir, "genes_per_spp_in_ogs_histograms.pdf"),
        height = 9,
        width = 13
      )
      write.table(
        og_gene_freqs,
        file = paste0(out_dir, "genes_per_spp_in_ogs_histograms.tsv"),
        sep = "\t",
        quote = F,
        col.names = T,
        row.names = F
      )
      
      # Same for histograms of eh # of species per orthogroup.
      ggsave(
        og_spp_hist,
        filename = paste0(out_dir, "num_species_in_og_hist.pdf"),
        height = 6,
        width = 6
      )
      write.table(
        og_spp_counts,
        file = paste0(out_dir, "num_species_in_og_hist.tsv"),
        row.names = F,
        col.names = T,
        quote = F
      )
    }
  }

# Function to summarize gene family content, both with respect to number of
# species included in each gene family, and the per-species mean copy number
summ_genefam_composition <-
  function(results_dir = NULL,
           show_plots = F,
           out_dir = "summarized-results/orthofinder/") {
    # Read in the gene copy number and species count for each gene family
    res <-
      read.table(
        paste0(results_dir, "filtered_orthogroups/all_ogs_counts.csv"),
        sep = ",",
        header = T
      )
    
    # Plot the mean species gene copy number against the number of species
    # in each gene family
    copynum_v_numspp <-
      ggplot(data = res, aes(y = mean_copy_num, x = num_spp)) +
      geom_point(position = position_jitter(width = 0.2),
                 alpha = 0.25,
                 size = 0.75) +
      theme_classic(base_size = 14) +
      scale_y_log10() +
      annotation_logticks(sides = "l") +
      xlab("# Species in gene family") +
      ylab("Mean per-species\ngene copy #")
    
    # And plot the histogram of species counts across all gene families
    fam_sppcount_dist <-
      ggplot(data = res, aes(x = num_spp)) +
      geom_histogram(fill = "lightgrey", color = "black") +
      theme_classic(base_size = 14) +
      scale_y_log10() +
      annotation_logticks(sides = "l") +
      xlab("# Species in gene family") +
      ylab("# Gene families")
    
    # Combine them
    genefam_composition_plt <-
      suppressMessages(plot_grid(copynum_v_numspp, fam_sppcount_dist, ncol = 2))
    
    # print the plots if that is requested
    if (show_plots == T) {
      print(genefam_composition_plt)
    }
    
    # And save
    ggsave(
      genefam_composition_plt,
      filename = paste0(out_dir, "genefamily_compostion.png"),
      height = 6,
      width = 12,
      dpi = 600
    )
    ggsave(
      genefam_composition_plt,
      filename = paste0(out_dir, "genefamily_compostion.pdf"),
      height = 6,
      width = 12
    )
  }

# Summarize and plot orthogroup statistics per species
get_per_spp_ofinder_stats <-
  function(tree = NULL,
           tree_fpath = NULL,
           results_dir = NULL,
           show_plots = T,
           samplesheet = NULL,
           samplesheet_fpath = NULL,
           grp_name = "Group",
           tip_grp_cols = NULL,
           count_cols = arcadia_cividis,
           prop_cols = arcadia_viridis,
           outgroup = NULL,
           out_dir = "summarized-results/orthofinder/") {
    if (is.null(tree)) {
      tree <- read.tree(tree_fpath)
    }
    if (is.null(tree) & is.null(tree_fpath)) {
      print("Error! Provide tree object or filepath to tree")
    }
    if (is.null(samplesheet)) {
      samplesheet <- read.delim(samplesheet_fpath, sep = ",")
    }
    if (is.null(samplesheet) & is.null(samplesheet_fpath)) {
      print("Error! Provide samplesheet object or filepath to samplesheet")
    }
    
    ###########################################
    # Orthogroup statistics per species       #
    ###########################################
    species <- tree$tip.label[order(tree$tip.label)]
    
    og_stats_perspp_fpath <-
      list.files(
        path = paste0(results_dir, "/orthofinder/complete_dataset/"),
        pattern = "Statistics_PerSpecies.tsv",
        recursive = T,
        full.names = T
      )
    
    # First read in the overall stats per species and write out as a table
    og_stats_perspp <-
      read_tsv(
        og_stats_perspp_fpath,
        skip = 1,
        n_max = 10,
        col_types = c("f", rep("n", length(species))),
        col_names = c("Statistic", species)
      )
    
    # Transpose and then plot stats as a heatmap alongside the the species tree
    og_stats_perspp <- og_stats_perspp %>%
      gather(key = species, value = value, 2:ncol(og_stats_perspp)) %>%
      spread(key = names(og_stats_perspp)[1], value = "value")
    # Check the above for use of "spread_"
    og_stats_perspp <-
      data.frame(og_stats_perspp[order(match(og_stats_perspp$species,
                                             tree$tip.label)),])
    rownames(og_stats_perspp) <- tree$tip.label
    
    # Pull out and rename relevant stats
    plot_dat_nums <- og_stats_perspp[c(2, 3, 4, 7, 9, 10)]
    plot_dat_props <- og_stats_perspp[c(5, 6, 8, 11)]
    colnames(plot_dat_nums) <-
      c(
        "# genes",
        "# genes in OGs",
        "# unassigned genes",
        "# OGs containing spp.",
        "# spp. specific OGs",
        "# genes in spp. specific OGs"
      )
    colnames(plot_dat_props) <-
      c(
        "% genes in OGs",
        "% unassigned genes",
        "% OGs containing spp.",
        "% genes in spp. specific OGs"
      )
    
    # Now prepare to plot
    # Get metadata for plotting
    taxa_tibble <- tibble(tip = samplesheet$species,
                          grp = samplesheet$taxonomy)
    
    # Root if using the asteroid tree and a specified outgroup
    if (!is.null(outgroup)) {
      # find the most recent common ancestor of the outgroup
      mrca_out <- getMRCA(tree, tip = outgroup)
      # find the parent of this MRCA node
      parent_node <- which(tree$edge[, 2] == mrca_out)
      # the parent edge number is the second column value at the corresponding row
      parent_edge_number <- tree$edge[parent_node, 2]
      # root the tree at the parent node
      tree <-
        root(tree,
             edgelabel = T,
             node = parent_edge_number,
             position = 0.25)
      
      # Make the basic tree plot
      base_p <- suppressMessages(
        ggtree(tree, size = 0.75, branch.length = "none") %<+% taxa_tibble +
          geom_tiplab(size = 7, as_ylab = T)
      )
    } else {
      # Make the basic tree plot
      base_p <- suppressMessages(ggtree(tree, size = 0.75) %<+% taxa_tibble +
                                   geom_tiplab(size = 7, as_ylab = T))
    }
    
    # Now one with the tip labels
    tree_p <- suppressMessages(
      base_p +
        geom_tippoint(
          aes(fill = grp),
          pch = 22,
          size = 3,
          color = "black"
        ) +
        scale_fill_manual(values = tsar_cols) +
        guides(
          fill = guide_legend(
            title = grp_name,
            nrow = 2,
            title.position = "top",
            title.hjust = 0
          )
        ) +
        theme(
          legend.position = "top",
          legend.justification = "left"
        ) +
        new_scale_fill()
    )
    
    # Combine with heatmap of per-species orthogroup count stats
    count_p <- suppressMessages(
      gheatmap(
        base_p,
        plot_dat_nums,
        offset = 0.0,
        font.size = 2.5,
        width = 0.6,
        hjust = 1,
        colnames_angle = 25,
        colnames_offset_y = 0.25
      ) +
        scale_fill_gradientn(
          colors = count_cols$color_dict,
          values = count_cols$values,
          trans = "log10"
        ) +
        theme(legend.position = "none")
    )
    
    # Get legends...
    # For the tree & tip labels
    leg1 <-
      get_legend(
        tree_p +
          theme(
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 10),
            legend.key.width = unit(1, "cm")
          )
      )
    
    # For the count stats
    leg2 <- get_legend(
      count_p +
        labs(fill = "Count") +
        guides(fill = guide_colorbar(
          title.position = "top",
          title.hjust = 0
        )) +
        theme(
          legend.position = "top",
          legend.justification = "left",
          legend.key.width = unit(1, "cm"),
          plot.margin = margin(10, 10, 15, 10),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 10)
        )
    )
    
    # And finally, proportional stats
    leg3 <- suppressMessages(get_legend(
      gheatmap(
        base_p,
        plot_dat_props,
        offset = 0.5,
        font.size = 2.5,
        width = 0.6,
        hjust = 1,
        colnames_angle = 25,
        colnames_offset_y = 0.25
      ) +
        scale_fill_gradientn(
          colors = prop_cols$color_dict,
          values = prop_cols$values,
          trans = "log10",
          labels = c(0, 5, 10, 25, 50, 75, 100),
          breaks = c(0, 5, 10, 25, 50, 75, 100)
        ) +
        theme(
          legend.position = "top",
          legend.justification = "left",
          legend.key.width = unit(1, "cm"),
          plot.margin = margin(10, 10, 15, 10),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 10)
        ) +
        guides(fill = guide_colorbar(
          title.position = "top",
          title.hjust = 0
        )) +
        labs(fill = "Percent")
    ))
    
    og_spp_stats_p <- suppressMessages(
      gheatmap(
        tree_p,
        plot_dat_nums,
        offset = 0.0,
        font.size = 2.5,
        width = 0.6,
        hjust = 1,
        colnames_angle = 25,
        colnames_offset_y = 0.25
      ) +
        scale_fill_gradientn(
          colors = count_cols$color_dict,
          values = count_cols$values,
          trans = "log10"
        ) +
        theme(legend.position = "none") +
        new_scale_fill()
    )
    
    og_spp_stats_p <- suppressMessages(
      gheatmap(
        og_spp_stats_p,
        plot_dat_props,
        offset = 3,
        font.size = 2.5,
        width = 0.4,
        hjust = 1,
        colnames_angle = 25,
        colnames_offset_y = 0.25
      ) +
        scale_fill_gradientn(
          colors = prop_cols$color_dict,
          values = prop_cols$values,
          trans = "log10",
          labels = c(0, 5, 10, 25, 50, 75, 100),
          breaks = c(0, 5, 10, 25, 50, 75, 100)
        ) +
        theme(
          legend.position = "none",
          plot.margin = margin(-10, 10, 35, 10)
        ) +
        coord_cartesian(clip = "off")
    )
    
    legs <-
      plot_grid(leg1, leg2, leg3, ncol = 3)
    
    og_spp_stats_p_final <-
      plot_grid(legs,
                og_spp_stats_p,
                ncol = 1,
                rel_heights = c(0.15, 1))
    
    # Plot out if requested
    if (show_plots == T) {
      print(og_spp_stats_p_final)
    }
    
    # Write out to file if an output directory is provided
    if (!is.null(out_dir)) {
      # Create the directory if it doesn't yet exist
      dir.create(out_dir, showWarnings = F, recursive = T)
      # And write out to file c( both the tabular summary and the figure)
      write.table(
        og_stats_perspp,
        file = paste0(out_dir, "per_spp_og_stats.tsv"),
        row.names = F,
        col.names = T,
        quote = F
      )
      
      # That was involved.... Save this to a pdf
      ggsave(
        og_spp_stats_p_final,
        filename = paste0(out_dir, "species_tree_gene_og_stat_heatmaps.pdf"),
        height = 7,
        width = 8
      )
      ggsave(
        og_spp_stats_p_final,
        filename = paste0(out_dir, "species_tree_gene_og_stat_heatmaps.png"),
        height = 7,
        width = 8,
        dpi = 600
      )
    }
    
    # Return both the dataframe and plot of per-species stats (in the inevitable
    # case that plotting size needs to tweaked for the given dataset)
    return(
      list(
        og_stats_per_spp_table = og_stats_perspp,
        og_stats_per_spp_plt = og_spp_stats_p_final
      )
    )
  }

###########################################################
# SpeciesRax/GeneRax summary functions
###########################################################
# A helper function to get per-gene family evolutionary events
# inferred from generax. Run internally to the "summarize_generax"
# function below
get_og_event_counts <-
  function(i,
           per_spp_og_counts = per_spp_og_counts,
           per_og_events = per_og_events) {
    # Populate an empty dataframe
    # Create a dataframe to store the per-OG event counts
    per_og_event_counts <-
      data.frame(
        gene_family = NA,
        speciation = NA,
        speciation_loss = NA,
        duplication = NA,
        transfer = NA,
        transfer_loss = NA,
        loss = NA,
        number_gene_copies = NA,
        number_species = NA
      )
    
    # Read in the the orthogroup-wide (across spp) event counts
    tmp <- read.table(per_og_events[i], sep = ":", check.names = F)
    gf <- gsub("_eventCounts.txt", "", per_og_events[i]) %>%
      gsub(".*reconciliations/", "", .)
    gf_col <- data.frame(gene_family = gf)
    
    # Get the per-species gene-count for this gene family
    species <-
      colnames(per_spp_og_counts)[-c(1, c((ncol(
        per_spp_og_counts
      ) - 1):ncol(per_spp_og_counts)))]
    
    counts <-
      per_spp_og_counts[which(per_spp_og_counts$Orthogroup == gf),]
    og_spps <- species[which(species %in% colnames(counts))]
    
    # Now fill
    per_og_event_counts[1,] <- c(gf, tmp$V2, counts$NumSpecies)
    
    # Return these event results as output
    return(per_og_event_counts)
  }

# Fill in the matrix of recipient-donor transfer events
get_og_event_rates <-
  function(i, per_og_rates_fpaths = per_og_rates_fpaths) {
    # Create a dataframe to store family-level rates
    per_og_event_rates <-
      data.frame(
        gene_family = NA,
        duplication = NA,
        loss = NA,
        transfer = NA
      )
    
    # Identify which gene family we"re dealing with
    gf <-
      gsub(".*OG", "OG", per_og_rates_fpaths[i]) %>%
      gsub("/stats.txt", "", .)
    
    # Read in and store those rates
    per_og_event_rates[1,] <-
      c(gf, scan(per_og_rates_fpaths[i],
                 what = "character", quiet = T)[6:8])
    per_og_event_rates[1, 2:4] <-
      as.numeric(per_og_event_rates[1, 2:4])
    
    # And return as output
    return(per_og_event_rates)
  }

get_spp_event_rates <-
  function(i, per_spp_rates_fpaths = per_spp_rates_fpaths) {
    rates <- read.table(per_spp_rates_fpaths[i], header = F, sep = " ")
    
    # Identify which gene family we're dealing with
    gf <-
      gsub(".*OG", "OG", per_spp_rates_fpaths[i]) %>%
      gsub("/per_species_rates.txt", "", .)
    
    # Create a dataframe to store species-level rates
    per_og_event_rates <-
      data.frame(
        gene_family = rep(gf, nrow(rates)),
        species = rates[, 1],
        duplication = rates[, 2],
        loss = rates[, 3],
        transfer = rates[, 4]
      )
    
    # And return as output
    return(per_og_event_rates)
  }

get_og_events_per_spp <-
  function(i,
           per_spp_og_counts = per_spp_og_counts,
           per_spp_og_events = per_spp_og_events) {
    # Get the names of species
    species <-
      colnames(per_spp_og_counts)[-c(1, c((ncol(
        per_spp_og_counts
      ) - 1):ncol(per_spp_og_counts)))]
    
    # Read in table of per-species event counts for this gene family
    tmp <-
      read.table(per_spp_og_events[i],
                 row.names = 1,
                 check.names = F)
    tmp <-
      data.frame(t(tmp[which(rownames(tmp) %in% species),]), check.names = F)
    
    # Identify which gene family we"re dealing with
    gf <-
      gsub("_speciesEventCounts.txt", "", per_spp_og_events[i]) %>%
      gsub(".*reconciliations/", "", .)
    gf_col <- data.frame(gene_family = gf)
    
    # Create one empty table for each parameter for the per-species counts per OG
    per_spp_og_speciation <-
      data.frame(matrix(ncol = length(species) + 1, nrow = 0))
    colnames(per_spp_og_speciation) <- c("gene_family", species)
    per_spp_og_duplication <-
      data.frame(matrix(ncol = length(species) + 1, nrow = 0))
    colnames(per_spp_og_duplication) <- c("gene_family", species)
    per_spp_og_transfer <-
      data.frame(matrix(ncol = length(species) + 1, nrow = 0))
    colnames(per_spp_og_transfer) <- c("gene_family", species)
    per_spp_og_loss <-
      data.frame(matrix(ncol = length(species) + 1, nrow = 0))
    colnames(per_spp_og_loss) <- c("gene_family", species)
    
    # And populate counts of duplication, transfer and loss
    specs <- cbind(gf_col, tmp[1,])
    dups <- cbind(gf_col, tmp[2,])
    loss <- cbind(gf_col, tmp[3,])
    transf <- cbind(gf_col, tmp[4,])
    
    # Now populate, allowing for species to not be present
    per_spp_og_speciation <-
      plyr::rbind.fill(per_spp_og_speciation, specs)
    per_spp_og_duplication <-
      plyr::rbind.fill(per_spp_og_duplication, dups)
    per_spp_og_loss <-
      plyr::rbind.fill(per_spp_og_loss, loss)
    per_spp_og_transfer <-
      plyr::rbind.fill(per_spp_og_transfer, transf)
    
    return(
      list(
        speciations = per_spp_og_speciation,
        duplications = per_spp_og_duplication,
        transfers = per_spp_og_transfer,
        losses = per_spp_og_loss
      )
    )
  }

# Now some quick helper functions to pull out speciations, duplications,
# transfers, and losses
get_speciations <-
  function(i, per_spp_events = per_spp_events) {
    speciactions <- per_spp_events[[i]]$speciations
    return(speciactions)
  }
get_duplications <-
  function(i, per_spp_events = per_spp_events) {
    duplications <- per_spp_events[[i]]$duplications
    return(duplications)
  }
get_losses <-
  function(i, per_spp_events = per_spp_events) {
    losses <- per_spp_events[[i]]$losses
    return(losses)
  }
get_transfers <-
  function(i, per_spp_events = per_spp_events) {
    transfers <- per_spp_events[[i]]$transfers
    return(transfers)
  }

# A function to get detailed info about transfer events (recipients, donors)
# for each gene family.
get_tranfer_donor_recips <-
  function(i,
           per_spp_og_counts = NULL,
           spp_tranf_rates_fpaths = NULL) {
    species <-
      colnames(per_spp_og_counts)[-c(1, (ncol(per_spp_og_counts) - 1):ncol(per_spp_og_counts))]
    
    per_spp_og_transfer_donor <-
      data.frame(matrix(ncol = length(species) + 1, nrow = 0))
    colnames(per_spp_og_transfer_donor) <- c("gene_family", species)
    per_spp_og_transfer_recip <-
      data.frame(matrix(ncol = length(species) + 1, nrow = 0))
    colnames(per_spp_og_transfer_recip) <- c("gene_family", species)
    
    transf_count_mat <-
      matrix(
        nrow = length(species),
        ncol = length(species),
        dimnames = list(species, species),
        data = 0
      )
    # Get the name of the gene family we're working with.
    gf <-
      gsub(".*OG", "OG", spp_tranf_rates_fpaths[i]) %>% gsub("_transfers.txt", "", .)
    # And the species that are in this gene family
    og_spps <- per_spp_og_counts[which(per_spp_og_counts$Orthogroup == gf),-c(1, (ncol(per_spp_og_counts) - 1):ncol(per_spp_og_counts))]
    og_spps <- colnames(og_spps[which(og_spps[1,] > 0),])
    # And then the table of recipient-donor events
    # But only if transfer events were inferred
    if (file.size(spp_tranf_rates_fpaths[i]) != 0L) {
      tmp <- read.table(spp_tranf_rates_fpaths[i], check.names = F)
      donor_spps <- og_spps[which(og_spps %in% tmp$V1)]
      donors <- summary(as.factor(tmp$V1))
      donors <-
        data.frame(as.list(donors[which(names(donors) %in% species)]), check.names = F)
      recip_spps <- og_spps[which(og_spps %in% tmp$V2)]
      recips <- summary(as.factor(tmp$V2))
      recips <-
        data.frame(as.list(recips[which(names(recips) %in% species)]), check.names = F)
      
      # And only if the transfer events occurred between tips
      # Make sure species included in the gene family have integer
      # counts, all other species NA
      if (nrow(donors) > 0) {
        donors$gene_family <- gf
        non_donors <- og_spps[-which(og_spps %in% donor_spps)]
        # If there are species who did not donate gene copies via transfer:
        if (length(non_donors) > 0) {
          non_donors <-
            data.frame(matrix(
              ncol = length(non_donors),
              dimnames = list(NULL, non_donors),
              data = 0
            ),
            check.names = F)
          donors <- cbind(donors, non_donors)
        }
        per_spp_og_transfer_donor <-
          plyr::rbind.fill(per_spp_og_transfer_donor, donors)
      } else {
        non_donors <-
          data.frame(matrix(
            ncol = length(og_spps),
            dimnames = list(NULL, og_spps),
            data = 0
          ),
          check.names = F)
        per_spp_og_transfer_donor <-
          plyr::rbind.fill(per_spp_og_transfer_donor, non_donors)
      }
      if (nrow(recips) > 0) {
        recips$gene_family <- gf
        non_recips <- og_spps[-which(og_spps %in% recip_spps)]
        # If there are species who did not receive gene copies via transfer:
        if (length(non_recips) > 0) {
          non_recips <-
            data.frame(matrix(
              ncol = length(non_recips),
              dimnames = list(NULL, non_recips),
              data = 0
            ),
            check.names = F)
          recips <- cbind(recips, non_recips)
        }
        per_spp_og_transfer_recip <-
          plyr::rbind.fill(per_spp_og_transfer_recip, recips)
      } else {
        non_recips <-
          data.frame(matrix(
            ncol = length(og_spps) + 1,
            dimnames = list(NULL, c("gene_family", og_spps)),
            data = c(gf, rep(0, length(og_spps)))
          ),
          check.names = F)
        per_spp_og_transfer_recip <-
          plyr::rbind.fill(per_spp_og_transfer_recip, non_recips)
      }
      # Now fill in the donor-recipient matrix
      tmp <- tmp[which(tmp$V1 %in% species & tmp$V2 %in% species),]
      for (x in 1:nrow(tmp)) {
        donor <- which(species == tmp$V1[x])
        recip <- which(species == tmp$V2[x])
        
        # y-axis: recipient, x-axis: donor
        transf_count_mat[donor, recip] <-
          transf_count_mat[donor, recip] + 1
      }
    } else {
      per_spp_og_transfer_donor <-
        data.frame(matrix(
          ncol = length(og_spps) + 1,
          dimnames = list(NULL, c("gene_family", og_spps)),
          data = c(gf, rep(0, length(og_spps)))
        ),
        check.names = F)
      per_spp_og_transfer_recip <-
        data.frame(matrix(
          ncol = length(og_spps) + 1,
          dimnames = list(NULL, c("gene_family", og_spps)),
          data = c(gf, rep(0, length(og_spps)))
        ),
        check.names = F)
    }
    return(
      list(
        "summed_matrix" = transf_count_mat,
        "gf_transfer_donors" = per_spp_og_transfer_donor,
        "gf_transfer_recips" = per_spp_og_transfer_recip
      )
    )
  }

# Read in and produce summaries for all generax results
# This included a multitude of steps
summarize_generax_per_family <-
  function(results_dir = "./",
           per_spp_og_counts = per_spp_og_counts,
           out_dir = NULL,
           nparallel = detectCores() - 1) {
    generax_dir <- paste0(results_dir, "generax/per_family_rates")
    
    # Get lists of files for each sort of analysis/output from generax that we"ll be summarizing
    per_og_events <-
      list.files(
        full.names = T,
        path = generax_dir,
        pattern = "eventCounts.txt",
        recursive = T
      )
    per_spp_og_events <-
      list.files(
        full.names = T,
        path = generax_dir,
        pattern = "speciesEventCounts.txt",
        recursive = T
      )
    per_og_rates_fpaths <-
      list.files(
        full.names = T,
        path = generax_dir,
        pattern = "stats.txt",
        recursive = T
      )
    per_og_rates_fpaths <-
      per_og_rates_fpaths[which(grepl("results/", per_og_rates_fpaths))]
    spp_tranf_rates_fpaths <-
      list.files(
        full.names = T,
        path = generax_dir,
        pattern = "transfers",
        recursive = T
      )
    per_spp_coverage <-
      list.files(
        full.names = T,
        path = generax_dir,
        pattern = "Coverage.txt",
        recursive = T
      )
    
    # Get the counts of each event type (duplications, transfers, losses, etc)
    # per-og, across all species
    message("Extracting event counts per-orthogroup, across all species.")
    per_og_event_res <-
      do.call(
        rbind,
        mclapply(
          X = 1:length(per_og_events),
          get_og_event_counts,
          per_spp_og_counts = per_spp_og_counts,
          per_og_events = per_og_events,
          mc.cores = nparallel
        )
      )
    
    # And the same, but for rates
    message("Extracting event rates per-orthogroup, across all species.")
    per_og_event_rates <-
      do.call(
        rbind,
        mclapply(
          X = 1:length(per_og_events),
          get_og_event_rates,
          per_og_rates_fpaths = per_og_rates_fpaths,
          mc.cores = nparallel
        )
      )
    
    # Get the counts of events per species, per orthogroup
    # Begin by first summarizing these event counts per orthogroup
    message("Extracting event counts per-species, per-orthogroup.")
    per_spp_events <-
      mclapply(
        1:length(per_spp_og_events),
        get_og_events_per_spp,
        per_spp_og_counts = per_spp_og_counts,
        per_spp_og_events = per_spp_og_events,
        mc.cores = nparallel
      )
    
    # And then pull out each event type individually
    message("Now, pulling out each event type individually.")
    per_spp_og_speciation <-
      do.call(
        rbind,
        mclapply(
          1:length(per_spp_events),
          get_speciations,
          per_spp_events = per_spp_events,
          mc.cores = nparallel
        )
      )
    per_spp_og_duplication <-
      do.call(
        rbind,
        mclapply(
          1:length(per_spp_events),
          get_duplications,
          per_spp_events = per_spp_events,
          mc.cores = nparallel
        )
      )
    per_spp_og_loss <-
      do.call(
        rbind,
        mclapply(
          1:length(per_spp_events),
          get_losses,
          per_spp_events = per_spp_events,
          mc.cores = nparallel
        )
      )
    
    # Clean up the large interim list
    rm(per_spp_events)
    
    # Now, focusing on transfers - get a summed matrix of transfers among species,
    # with donors along the x-axis, and recipients along the y.
    # y-axis: recipient, x-axis: donor
    message("Summarizing gene transfer recipient events.")
    transf_res <-
      transpose(
        mclapply(
          1:length(spp_tranf_rates_fpaths),
          get_tranfer_donor_recips,
          per_spp_og_counts = per_spp_og_counts,
          spp_tranf_rates_fpaths = spp_tranf_rates_fpaths,
          mc.cores = nparallel
        )
      )
    message("Summarizing gene transfer events into a matrix of donor-recipient species pairs.")
    transf_count_mat <- Reduce("+", transf_res$summed_matrix)
    message("Pulling out the count of transfer-donor events for each species per gene family")
    transf_donors <- do.call("rbind", transf_res$gf_transfer_donors)
    message("Pulling out the count of transfer-recipient events for each species per gene family")
    transf_recips <- do.call("rbind", transf_res$gf_transfer_recips)
    
    # Generate a list containing alll required outputs for plotting
    results <- list(
      rates_per_og = per_og_event_rates,
      events_per_og = per_og_event_res,
      lgt_count_mat = transf_count_mat,
      speciations_per_spp = per_spp_og_speciation,
      duplications_per_spp = per_spp_og_duplication,
      losses_per_spp = per_spp_og_loss,
      transfer_donor_counts = transf_donors,
      transfer_recip_counts = transf_recips
    )
    
    # Write each output to their own stand-alone summary table
    if (!is.null(out_dir)) {
      # Create the directory if it doesn't yet exist
      dir.create(out_dir, showWarnings = F, recursive = T)
      # And write out to file
      write.table(
        per_og_event_rates,
        file = paste0(out_dir, "event_rates_per_gene_family.tsv"),
        sep = "\t",
        quote = F,
        row.names = F,
        col.names = T
      )
      write.table(
        per_og_event_res,
        file = paste0(out_dir, "event_counts_per_gene_family.tsv"),
        sep = "\t",
        quote = F,
        row.names = F,
        col.names = T
      )
      write.table(
        transf_count_mat,
        file = paste0(out_dir, "hgt_summed_counts_recip_donor.tsv"),
        sep = "\t",
        quote = F,
        row.names = T,
        col.names = NA
      )
      write.table(
        per_spp_og_speciation,
        file = paste0(
          out_dir,
          "speciation_count_per_species_per_gene_family.tsv"
        ),
        sep = "\t",
        quote = F,
        row.names = F,
        col.names = T
      )
      write.table(
        per_spp_og_duplication,
        file = paste0(
          out_dir,
          "duplication_count_per_species_per_gene_family.tsv"
        ),
        sep = "\t",
        quote = F,
        row.names = F,
        col.names = T
      )
      write.table(
        per_spp_og_loss,
        file = paste0(out_dir, "loss_count_per_per_species_gene_family.tsv"),
        sep = "\t",
        quote = F,
        row.names = F,
        col.names = T
      )
      write.table(
        transf_donors,
        file = paste0(
          out_dir,
          "transfer_donor_count_per_species_per_gene_family.tsv"
        ),
        sep = "\t",
        quote = F,
        row.names = F,
        col.names = T
      )
      write.table(
        transf_recips,
        file = paste0(
          out_dir,
          "transfer_recipient_count_per_species_per_gene_family.tsv"
        ),
        sep = "\t",
        quote = F,
        row.names = F,
        col.names = T
      )
    }
    return(results)
  }

summarize_generax_per_species <-
  function(results_dir = "./",
           per_spp_og_counts = per_spp_og_counts,
           out_dir = NULL,
           nparallel = detectCores() - 1) {
    generax_dir <- paste0(results_dir, "generax/per_species_rates")
    
    # Get lists of files for each sort of analysis/output from generax that we"ll be summarizing
    per_og_events <-
      list.files(
        full.names = T,
        path = generax_dir,
        pattern = "eventCounts.txt",
        recursive = T
      )
    per_spp_og_events <-
      list.files(
        full.names = T,
        path = generax_dir,
        pattern = "speciesEventCounts.txt",
        recursive = T
      )
    per_spp_rates_fpaths <-
      list.files(
        full.names = T,
        path = generax_dir,
        pattern = "per_species_rates.txt",
        recursive = T
      )
    spp_tranfers_fpaths <-
      list.files(
        full.names = T,
        path = generax_dir,
        pattern = "_transfers.txt",
        recursive = T
      )
    per_spp_coverage <-
      list.files(
        full.names = T,
        path = generax_dir,
        pattern = "Coverage.txt",
        recursive = T
      )
    
    # Get the counts of each event type (duplications, transfers, losses, etc)
    # per-og, across all species
    message("Extracting event counts for each species per gene family.")
    per_og_event_res <-
      do.call(
        rbind,
        mclapply(
          X = 1:length(per_og_events),
          get_og_event_counts,
          per_spp_og_counts = per_spp_og_counts,
          per_og_events = per_og_events,
          mc.cores = nparallel
        )
      )
    
    # And the same, but for rates
    message("Extracting event rates per-species, for each gene family.")
    per_spp_event_rates <-
      do.call(
        rbind,
        mclapply(
          X = 1:length(per_spp_rates_fpaths),
          get_spp_event_rates,
          per_spp_rates_fpaths = per_spp_rates_fpaths,
          mc.cores = nparallel
        )
      )
    
    # Get the counts of events per species, per orthogroup
    # Begin by first summarizing these event counts per orthogroup
    message("Extracting event counts per-species, per-orthogroup.")
    per_spp_events <-
      mclapply(
        1:length(per_spp_og_events),
        get_og_events_per_spp,
        per_spp_og_counts = per_spp_og_counts,
        per_spp_og_events = per_spp_og_events,
        mc.cores = nparallel
      )
    
    # And then pull out each event type individually
    message("Now, pulling out each event type individually.")
    per_spp_og_speciation <-
      do.call(
        rbind,
        mclapply(
          1:length(per_spp_events),
          get_speciations,
          per_spp_events = per_spp_events,
          mc.cores = nparallel
        )
      )
    per_spp_og_duplication <-
      do.call(
        rbind,
        mclapply(
          1:length(per_spp_events),
          get_duplications,
          per_spp_events = per_spp_events,
          mc.cores = nparallel
        )
      )
    per_spp_og_loss <-
      do.call(
        rbind,
        mclapply(
          1:length(per_spp_events),
          get_losses,
          per_spp_events = per_spp_events,
          mc.cores = nparallel
        )
      )
    
    # Clean up the large interim list
    rm(per_spp_events)
    
    # Now, focusing on transfers - get a summed matrix of transfers among species,
    # with donors along the x-axis, and recipients along the y.
    # y-axis: recipient, x-axis: donor
    message("Summarizing gene transfer recipient events.")
    transf_res <-
      transpose(
        mclapply(
          1:length(spp_tranfers_fpaths),
          get_tranfer_donor_recips,
          per_spp_og_counts = per_spp_og_counts,
          spp_tranf_rates_fpaths = spp_tranfers_fpaths,
          mc.cores = nparallel
        )
      )
    message("Summarizing gene transfer events into a matrix of donor-recipient species pairs.")
    transf_count_mat <- Reduce("+", transf_res$summed_matrix)
    message("Pulling out the count of transfer-donor events for each species per gene family")
    transf_donors <- do.call("rbind", transf_res$gf_transfer_donors)
    message("Pulling out the count of transfer-recipient events for each species per gene family")
    transf_recips <- do.call("rbind", transf_res$gf_transfer_recips)
    
    # Generate a list containing alll required outputs for plotting
    results <- list(
      spp_rates_per_og = per_spp_event_rates,
      events_per_og = per_og_event_res,
      lgt_count_mat = transf_count_mat,
      speciations_per_spp = per_spp_og_speciation,
      duplications_per_spp = per_spp_og_duplication,
      losses_per_spp = per_spp_og_loss,
      transfer_donor_counts = transf_donors,
      transfer_recip_counts = transf_recips
    )
    
    # Write each output to their own stand-alone summary table
    if (!is.null(out_dir)) {
      # Create the directory if it doesn't yet exist
      dir.create(out_dir, showWarnings = F, recursive = T)
      # And write out to file
      write.table(
        per_spp_event_rates,
        file = paste0(out_dir, "species_event_rates_per_gene_family.tsv"),
        sep = "\t",
        quote = F,
        row.names = F,
        col.names = T
      )
      write.table(
        per_og_event_res,
        file = paste0(out_dir, "species_event_counts_per_gene_family.tsv"),
        sep = "\t",
        quote = F,
        row.names = F,
        col.names = T
      )
      write.table(
        transf_count_mat,
        file = paste0(out_dir, "hgt_summed_counts_recip_donor.tsv"),
        sep = "\t",
        quote = F,
        row.names = T,
        col.names = NA
      )
      write.table(
        per_spp_og_speciation,
        file = paste0(
          out_dir,
          "speciation_count_per_species_per_gene_family.tsv"
        ),
        sep = "\t",
        quote = F,
        row.names = F,
        col.names = T
      )
      write.table(
        per_spp_og_duplication,
        file = paste0(
          out_dir,
          "duplication_count_per_species_per_gene_family.tsv"
        ),
        sep = "\t",
        quote = F,
        row.names = F,
        col.names = T
      )
      write.table(
        per_spp_og_loss,
        file = paste0(out_dir, "loss_count_per_per_species_gene_family.tsv"),
        sep = "\t",
        quote = F,
        row.names = F,
        col.names = T
      )
      write.table(
        transf_donors,
        file = paste0(
          out_dir,
          "transfer_donor_count_per_species_per_gene_family.tsv"
        ),
        sep = "\t",
        quote = F,
        row.names = F,
        col.names = T
      )
      write.table(
        transf_recips,
        file = paste0(
          out_dir,
          "transfer_recipient_count_per_species_per_gene_family.tsv"
        ),
        sep = "\t",
        quote = F,
        row.names = F,
        col.names = T
      )
    }
    return(results)
  }

# This function will return the node number given the node name and a phylogenetic tree
update_node_names <- function(species_rates = NULL,
                              tree = NULL) {
  # Get the names of internal nodes from speciesrax
  nodelabs <- species_rates$node[-c(1:length(tree$tip.label))]
  # Get the original tip labels
  tiplabs <- tree$tip.label
  # Update the tip labels to only use hyphens instead of underscores
  tree$tip.label <- gsub("_", "-", tree$tip.label)
  # Now update the node labels accordingly
  for (i in 1:length(tiplabs)) {
    nodelabs <- gsub(tiplabs[i], tree$tip.label[i], nodelabs)
  }
  
  # Now identify the node number that corresponds to each of GeneRax's named nodes
  for (node_name in nodelabs) {
    # Now split the string on underscores to get the taxa names
    taxa <- unlist(strsplit(node_name, "_"))[2:3]
    
    # Find the node number (MRCA) corresponding to the taxa in the original node name
    nodelabs[which(nodelabs == node_name)] <-
      getMRCA(tree, match(taxa, tree$tip.label))
  }
  # Now update the rates table with these node ids
  species_rates$node[-c(1:length(tree$tip.label))] <- nodelabs
  rownames(species_rates) <- species_rates$node
  
  # And do the same for species (terminal nodes).
  # These are ordered intuitively.
  species_rates$node[1:length(tree$tip.label)] <-
    1:length(tree$tip.label)
  species_rates$node <- as.numeric(species_rates$node)
  
  return(species_rates)
}

# A function to plot GeneRax DTL rates per species
plot_rates_per_spp <-
  function(tree = NULL,
           species_rates = NULL,
           align_tips = T,
           tip_offset = 0.1,
           rate_cols = list(arcadia_dahlias, arcadia_pansies, arcadia_poppies)) {
    # Combine the DTL rates with the species tree for plotting.
    dtl_rates_tree <- full_join(tree, species_rates, by = "node")
    
    # define a function to get the plotting breakpoints:
    get_colbreaks <- function(x) {
      c(min(x), mean(range(x)), max(x))
    }
    
    # And plot:
    dup_rate_p <-
      ggtree(
        dtl_rates_tree,
        aes(color = duplication),
        ladderize = TRUE,
        continuous = "color",
        size = 1.5
      ) +
      scale_color_gradientn(
        colors = rate_cols[[1]]$color_dict[1:3],
        values = c(0, 0.3, 1.0),
        name = "Duplication\nRate",
        labels = function(values)
          sprintf("%.2f", values)
      ) +
      geom_tiplab(
        offset = tip_offset,
        align = align_tips,
        size = 3,
        color = "black"
      ) +
      theme(
        legend.position = "top",
        plot.margin = margin(5, 140, 10, 5),
        legend.text = element_text(size = 10),
        legend.key.width = unit(0.9, "cm"),
        legend.title = element_text(size = 12)
      ) +
      coord_cartesian(clip = "off")
    
    transf_rate_p <-
      ggtree(
        dtl_rates_tree,
        aes(color = transfer),
        ladderize = TRUE,
        continuous = "color",
        size = 1.5
      ) +
      scale_color_gradientn(
        colors = rate_cols[[2]]$color_dict[1:4],
        name = "Transfer\nRate",
        labels = function(values)
          sprintf("%.2f", values)
      ) +
      geom_tiplab(
        offset = tip_offset,
        align = align_tips,
        size = 3,
        color = "black"
      ) +
      theme(
        legend.position = "top",
        plot.margin = margin(5, 140, 10, 5),
        legend.text = element_text(size = 10),
        legend.key.width = unit(0.9, "cm"),
        legend.title = element_text(size = 12)
      ) +
      coord_cartesian(clip = "off")
    
    loss_rate_p <-
      ggtree(
        dtl_rates_tree,
        aes(color = loss),
        ladderize = TRUE,
        continuous = "color",
        size = 1.5
      ) +
      scale_color_gradientn(
        colors = rev(rate_cols[[3]]$color_dict[4:8]),
        name = "Loss\nRate",
        labels = function(values)
          sprintf("%.2f", values)
      ) +
      geom_tiplab(
        offset = tip_offset,
        align = align_tips,
        size = 3,
        color = "black"
      ) +
      theme(
        legend.position = "top",
        plot.margin = margin(5, 140, 10, 5),
        legend.text = element_text(size = 10),
        legend.key.width = unit(0.9, "cm"),
        legend.title = element_text(size = 12)
      ) +
      coord_cartesian(clip = "off")
    
    dt_ratio_p <-
      ggtree(
        dtl_rates_tree,
        aes(color = dt_ratio),
        ladderize = TRUE,
        continuous = "color",
        size = 1.5
      ) +
      scale_color_gradientn(
        colors = rate_cols[[3]]$color_dict[1:4],
        name = "Duplication /\nTransfer Rate",
        labels = function(values)
          sprintf("%.2f", values)
      ) +
      geom_tiplab(
        offset = tip_offset,
        align = align_tips,
        size = 3,
        color = "black"
      ) +
      theme(
        legend.position = "top",
        plot.margin = margin(5, 140, 10, 5),
        legend.text = element_text(size = 10),
        legend.key.width = unit(0.9, "cm"),
        legend.title = element_text(size = 12)
      ) +
      coord_cartesian(clip = "off")
    
    rate_plts <- list(
      "rate_tree_data" = dtl_rates_tree,
      "duplication" = dup_rate_p,
      "transfer" = transf_rate_p,
      "loss" = loss_rate_p,
      "dt_ratio" = dt_ratio_p
    )
    return(rate_plts)
  }

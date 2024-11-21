# clustering gene families based on per-species/node event counts
source("./scripts/noveltree-summary-functions.R")
source("./scripts/core_phylo_profiling_pipeline.R")
source("./scripts/profile_cluster_trait_association_functions.R")

# Pull in NovelTree outputs from the Chelicerata dataset
result_directory <-
  "chelicerata-v1-10062023/" # available on Zenodo
generax_dir <-
  paste0(result_directory, "generax/per_species_rates/")

# Read in the species tree, retaining the node labels. This can be found in the
# result directories of any gene family.
spptree <-
  read.tree(paste0(
    generax_dir,
    "OG0000017/species_trees/inferred_species_tree.newick"
  ))
spptree <- ladderize(spptree)

# Read in the GeneRax event counts. Use the
per_spp_og_events <-
  list.files(
    full.names = TRUE,
    path = generax_dir,
    pattern = "speciesEventCounts.txt",
    recursive = TRUE
  )

og_counts <- get_per_spp_og_counts(results_dir = result_directory)

# Define a function to read in the event counts per-node (terminal and internal)
get_og_events_per_node <-
  function(i, per_spp_og_events = per_spp_og_events, spptree) {
    # Read in table of per-species event counts for this gene family
    tmp <- read.table(per_spp_og_events[i], check.names = FALSE)
    colnames(tmp) <- c("node", "s", "sl", "d", "t")

    return(tmp)
  }

# Use this function to read in the actual event count data for each gene family
gf_node_events <- list()
for (i in 1:length(per_spp_og_events)) {
  gf_node_events[[i]] <-
    get_og_events_per_node(i, per_spp_og_events, spptree)
}

# Get the name of orthogroups
ogs <-
  gsub(".*/", "", per_spp_og_events) |>
  gsub(pattern = "_.*", replacement = "")
names(gf_node_events) <- ogs

# Create empty dataframes to hold each count variable, one row
# per gene family
s_counts <-
  matrix(
    NA,
    nrow = length(gf_node_events),
    ncol = nrow(gf_node_events[[1]]),
    dimnames = list(
      names(gf_node_events),
      c(gf_node_events[[1]]$node)
    )
  )
sl_counts <- s_counts
d_counts <- s_counts
t_counts <- s_counts

# Populate these empty count dataframes with the observed counts of each
# respective event type
for (i in 1:length(gf_node_events)) {
  t_dat <- t(gf_node_events[[i]][, -1])
  t_dat <- t_dat
  s_counts[i, ] <- log10(t_dat[1, ] + 1)
  sl_counts[i, ] <- log10(t_dat[2, ] + 1)
  d_counts[i, ] <- log10(t_dat[3, ] + 1)
  t_counts[i, ] <- log10(t_dat[4, ] + 1)
}

# Create a list of these four event count datasets - rather than providing
# a single dataframe of all counts, we provide a list so that we can
# individually conduct Cholensky transformations on each dataset to obtain
# a "phylogenetically corrected" VCV to use in mahalanobis distance calculation
all_counts <- cbind(s_counts, sl_counts, d_counts, t_counts)

# Read in the gene family annotations so we can include
# them in plots
prot_annots <-
  read_tsv("./metadata/final_genefam_protein_name_protein_name_summary.tsv")
prot_annots <-
  prot_annots[which(prot_annots$Orthogroup %in% names(gf_node_events)), ]

prot_gos <-
  read_tsv("./metadata/final_genefam_go_term_go_id_summary.tsv")
prot_gos <-
  prot_gos[which(prot_gos$Orthogroup %in% names(gf_node_events)), ]
colnames(prot_gos)[2:3] <- c("best_gos", "second_best_gos")

# Now, go ahead and conduct the phylogenetic profiling, clustering gene families
# according to different types of gene family evolutionary events

# Counts of speciation events
s_dists_clusts <-
  run_phylo_profiling(
    data = s_counts,
    tree = spptree,
    cluster_method = "all",
    profile_title = "Speciation",
    annotations = prot_annots,
    distance_metric = "mahalanobis",
    k = 10
  )
# Counts of duplication events
d_dists_clusts <-
  run_phylo_profiling(
    d_counts,
    spptree,
    cluster_method = "all",
    profile_title = "Duplication",
    annotations = prot_annots,
    distance_metric = "mahalanobis",
    k = 10
  )
# Counts of loss events
sl_dists_clusts <-
  run_phylo_profiling(
    sl_counts,
    spptree,
    cluster_method = "all",
    profile_title = "Loss",
    annotations = prot_annots,
    distance_metric = "mahalanobis",
    k = 10
  )
# Counts of horizontal transfer events
t_dists_clusts <-
  run_phylo_profiling(
    t_counts,
    spptree,
    cluster_method = "all",
    profile_title = "Transfer",
    annotations = prot_annots,
    distance_metric = "mahalanobis",
    k = 10
  )

# And using the combination of all event types
all_dists_clusts <-
  run_phylo_profiling(
    list(
      s_counts = s_counts,
      sl_counts = sl_counts,
      t_counts = t_counts,
      d_counts = d_counts
    ),
    spptree,
    cluster_method = "all",
    profile_title = "Combined",
    annotations = prot_annots,
    distance_metric = "mahalanobis",
    k = 10
  )

# Now, save the constituent results to their respective subdirectories. Again,
# create these directories if you haven't already done so.
dir.create("chelicerate-results/plots/",
  showWarnings = FALSE,
  recursive = TRUE
)
dir.create(
  "chelicerate-results/phylo-transformed-profiles/",
  showWarnings = FALSE,
  recursive = TRUE
)
dir.create(
  "chelicerate-results/umap-layout-cluster-ids/",
  showWarnings = FALSE,
  recursive = TRUE
)

# First the PDFs of the UMAPs depicting the clustering of gene families based
# on their phylogenetic profiles using the best-performing clustering algorithm
# according to modularity score
ggsave(
  s_dists_clusts$static_umap,
  file = "chelicerate-results/plots/speciation_profile_umap_best_clusters.pdf",
  height = 10,
  width = 10
)
ggsave(
  d_dists_clusts$static_umap,
  file = "chelicerate-results/plots/duplication_profile_umap_best_clusters.pdf",
  height = 10,
  width = 10
)
ggsave(
  sl_dists_clusts$static_umap,
  file = "chelicerate-results/plots/loss_profile_umap_best_clusters.pdf",
  height = 10,
  width = 10
)
ggsave(
  t_dists_clusts$static_umap,
  file = "chelicerate-results/plots/transfer_profile_umap_best_clusters.pdf",
  height = 10,
  width = 10
)
ggsave(
  all_dists_clusts$static_umap,
  file = "chelicerate-results/plots/combined_profile_umap_best_clusters.pdf",
  height = 10,
  width = 10
)

# And now save the HTML widgets of the plotly UMAP plots depicting the results
# of each clustering algorithm. Points are labeled according to their cluster
# IDs and the protein names from genefinder
htmlwidgets::saveWidget(
  as_widget(s_dists_clusts$interactive_umap),
  "chelicerate-results/plots/speciation_profile_umap_clustering.html"
)
htmlwidgets::saveWidget(
  as_widget(d_dists_clusts$interactive_umap),
  "chelicerate-results/plots/duplication_profile_umap_clustering.html"
)
htmlwidgets::saveWidget(
  as_widget(sl_dists_clusts$interactive_umap),
  "chelicerate-results/plots/loss_profile_umap_clustering.html"
)
htmlwidgets::saveWidget(
  as_widget(t_dists_clusts$interactive_umap),
  "chelicerate-results/plots/transfer_profile_umap_clustering.html"
)
htmlwidgets::saveWidget(
  as_widget(all_dists_clusts$interactive_umap),
  "chelicerate-results/plots/combined_profile_umap_clustering.html"
)

# Save out the phylogenetic GLS transformed profiles
# Don't worry about saving out the transformed data from the combined analysis -
# these are already contained in the other outputs
write.table(
  s_dists_clusts$phylo_corrected_data,
  file = "chelicerate-results/phylo-transformed-profiles/speciation_phylo_gls_transformed_profiles.tsv",
  quote = FALSE,
  sep = ",",
  col.names = TRUE,
  row.names = TRUE
)
write.table(
  d_dists_clusts$phylo_corrected_data,
  file = "chelicerate-results/phylo-transformed-profiles/duplication_phylo_gls_transformed_profiles.tsv",
  quote = FALSE,
  sep = ",",
  col.names = TRUE,
  row.names = TRUE
)
write.table(
  sl_dists_clusts$phylo_corrected_data,
  file = "chelicerate-results/phylo-transformed-profiles/loss_phylo_gls_transformed_profiles.tsv",
  quote = FALSE,
  sep = ",",
  col.names = TRUE,
  row.names = TRUE
)
write.table(
  t_dists_clusts$phylo_corrected_data,
  file = "chelicerate-results/phylo-transformed-profiles/transfer_phylo_gls_transformed_profiles.tsv",
  quote = FALSE,
  sep = ",",
  col.names = TRUE,
  row.names = TRUE
)

# And lastly, save out to file the umap layout that includes cluster IDs as well
# as protein annotations. Before we do so, go ahead and combine the gene
# ontology annotations with these.
s_dists_clusts$umap_clusters <-
  cbind(s_dists_clusts$umap_clusters, prot_gos[which(prot_gos$Orthogroup %in% rownames(s_dists_clusts$umap_clusters)), 2:3])
d_dists_clusts$umap_clusters <-
  cbind(d_dists_clusts$umap_clusters, prot_gos[which(prot_gos$Orthogroup %in% rownames(d_dists_clusts$umap_clusters)), 2:3])
sl_dists_clusts$umap_clusters <-
  cbind(sl_dists_clusts$umap_clusters, prot_gos[which(prot_gos$Orthogroup %in% rownames(sl_dists_clusts$umap_clusters)), 2:3])
t_dists_clusts$umap_clusters <-
  cbind(t_dists_clusts$umap_clusters, prot_gos[which(prot_gos$Orthogroup %in% rownames(t_dists_clusts$umap_clusters)), 2:3])
all_dists_clusts$umap_clusters <-
  cbind(all_dists_clusts$umap_clusters, prot_gos[which(prot_gos$Orthogroup %in% rownames(all_dists_clusts$umap_clusters)), 2:3])
write.table(
  s_dists_clusts$umap_clusters,
  file = "chelicerate-results/umap-layout-cluster-ids/speciation_profile_clusters_annotations.tsv",
  quote = FALSE,
  sep = ",",
  col.names = TRUE,
  row.names = TRUE
)
write.table(
  d_dists_clusts$umap_clusters,
  file = "chelicerate-results/umap-layout-cluster-ids/duplication_profile_clusters_annotations.tsv",
  quote = FALSE,
  sep = ",",
  col.names = TRUE,
  row.names = TRUE
)
write.table(
  sl_dists_clusts$umap_clusters,
  file = "chelicerate-results/umap-layout-cluster-ids/loss_profile_clusters_annotations.tsv",
  quote = FALSE,
  sep = ",",
  col.names = TRUE,
  row.names = TRUE
)
write.table(
  t_dists_clusts$umap_clusters,
  file = "chelicerate-results/umap-layout-cluster-ids/transfer_profile_clusters_annotations.tsv",
  quote = FALSE,
  sep = ",",
  col.names = TRUE,
  row.names = TRUE
)
write.table(
  all_dists_clusts$umap_clusters,
  file = "chelicerate-results/umap-layout-cluster-ids/combined_profile_clusters_annotations.tsv",
  quote = FALSE,
  sep = ",",
  col.names = TRUE,
  row.names = TRUE
)


################################################################################
# Alright, so now that we've generated our clusters from phylogenetic profiles,
# let's go another step further and see if we can determine which clusters are
# most strongly associated with the trait of interest.
# To do so, we will take the following approach. Looping through each gene
# family, we will:
# 1) Conduct phylogenetic GLS transformations of gene family evolutionary event
#    counts for each species and obtain per-species means
# 2) Pair these data with the detection-suppression trait data
# 3) Summarize across states, obtaining mean event counts for species in state
#    1 or 0 for detection suppression.
# 4) Then, for each cluster, conduct non-parametric tests to determine whether
#    detection-suppressing species have significantly different (or greater)
#    event counts than non-detection suppressing species

# Start by obtaining the transformed data for tips. Keep only the data for
# terminal branches:
terminal <-
  which(colnames(s_counts) %in% spptree$tip.label)
s_dat <- s_counts[, terminal]
d_dat <- d_counts[, terminal]
t_dat <- t_counts[, terminal]
l_dat <- sl_counts[, terminal]

s_dat <-
  as.data.frame(t(do.call(cbind, phylo_correction(
    t(s_dat), vcvPhylo(spptree, anc.nodes = FALSE)
  ))))
d_dat <-
  as.data.frame(t(do.call(cbind, phylo_correction(
    t(d_dat), vcvPhylo(spptree, anc.nodes = FALSE)
  ))))
t_dat <-
  as.data.frame(t(do.call(cbind, phylo_correction(
    t(t_dat), vcvPhylo(spptree, anc.nodes = FALSE)
  ))))
l_dat <-
  as.data.frame(t(do.call(cbind, phylo_correction(
    t(l_dat), vcvPhylo(spptree, anc.nodes = FALSE)
  ))))

dimnames(s_dat) <-
  list(rownames(s_counts[, terminal]), colnames(s_counts[, terminal]))
dimnames(d_dat) <-
  list(rownames(d_counts[, terminal]), colnames(d_counts[, terminal]))
dimnames(t_dat) <-
  list(rownames(t_counts[, terminal]), colnames(t_counts[, terminal]))
dimnames(l_dat) <-
  list(rownames(sl_counts[, terminal]), colnames(sl_counts[, terminal]))

# Read in trait data:
metadat <- read_tsv("./metadata/chelicerate-species-metadata.txt")

# Update the host detection suppression trait so NA's and Unknown's are No
metadat$host_detection_suppression[which(is.na(metadat$host_detection_suppression))] <-
  "No"
metadat$host_detection_suppression[which(metadat$host_detection_suppression == "Unknown")] <-
  "No"
metadat$host_detection_suppression <-
  factor(metadat$host_detection_suppression, levels = c("No", "Yes"))

# Determine which species can and cannot suppress detection
suppressors <-
  metadat$id[which(metadat$host_detection_suppression == "Yes")]
nonsuppressors <-
  metadat$id[which(metadat$host_detection_suppression == "No")]

### Remove clusters with small samples sizes (<25 OGs), which led to noisy parameter estimates, possibly due to a greater sensitivity to model initialization conditions and low statistical power
s_cluster_counts <- table(s_dists_clusts$umap_clusters$leiden_cluster_id)
t_cluster_counts <- table(t_dists_clusts$umap_clusters$leiden_cluster_id)
sl_cluster_counts <- table(sl_dists_clusts$umap_clusters$leiden_cluster_id)
all_cluster_counts <- table(all_dists_clusts$umap_clusters$leiden_cluster_id)

s_valid_clusters <- names(s_cluster_counts[s_cluster_counts >= 25])
t_valid_clusters <- names(t_cluster_counts[t_cluster_counts >= 25])
sl_valid_clusters <- names(sl_cluster_counts[sl_cluster_counts >= 25])
all_valid_clusters <- names(all_cluster_counts[all_cluster_counts >= 25])

s_dists_clusts$umap_clusters <- s_dists_clusts$umap_clusters[
  s_dists_clusts$umap_clusters$leiden_cluster_id %in% s_valid_clusters, 
]
t_dists_clusts$umap_clusters <- t_dists_clusts$umap_clusters[
  t_dists_clusts$umap_clusters$leiden_cluster_id %in% t_valid_clusters, 
]
sl_dists_clusts$umap_clusters <- sl_dists_clusts$umap_clusters[
  sl_dists_clusts$umap_clusters$leiden_cluster_id %in% sl_valid_clusters, 
]
all_dists_clusts$umap_clusters <- all_dists_clusts$umap_clusters[
  all_dists_clusts$umap_clusters$leiden_cluster_id %in% all_valid_clusters, 
]

# Conduct the tests, for each gene family and each event type of whether event
# counts differ between detection-suppressing species
# NOTE: We cannot use duplications here, as they are not inferred to have
# occurred in terminal branches
s_res <- cluster_assoc_test(
  count_dat = s_dat,
  cluster_res = s_dists_clusts,
  suppressors,
  nonsuppressors,
  method = "logistic"
)
t_res <- cluster_assoc_test(t_dat, t_dists_clusts,
  suppressors, nonsuppressors,
  method = "logistic"
)
l_res <- cluster_assoc_test(l_dat, sl_dists_clusts,
  suppressors, nonsuppressors,
  method = "logistic"
)

# Now do the same, but using the clusters inferred from the full set of event types
s_comb_res <- cluster_assoc_test(s_dat, all_dists_clusts,
  suppressors, nonsuppressors,
  method = "logistic"
)
t_comb_res <- cluster_assoc_test(t_dat, all_dists_clusts,
  suppressors, nonsuppressors,
  method = "logistic"
)
l_comb_res <- cluster_assoc_test(l_dat, all_dists_clusts,
  suppressors, nonsuppressors,
  method = "logistic"
)

# Plot the results for each event type, using the clusters inferred from the
# same corresponding event types
s_corr_res <- plot_cluster_traitcorr(
  pvals = s_res$pvals,
  coefficients = s_res$coefficients,
  profile_type = "speciation",
  cluster_data = "speciation",
  cluster_ids = s_res$clusters,
  method = "logistic"
)
t_corr_res <-
  plot_cluster_traitcorr(t_res$pvals,
    t_res$coefficient,
    "transfer",
    "transfer",
    t_res$clusters,
    method = "logistic"
  )
l_corr_res <-
  plot_cluster_traitcorr(l_res$pvals,
    l_res$coefficient,
    "loss",
    "loss",
    l_res$clusters,
    method = "logistic"
  )

# And do the same for each event type, but using the clusters inferred from
# the full set of event data
s_comb_corr_res <-
  plot_cluster_traitcorr(
    s_comb_res$pvals,
    s_comb_res$coefficient,
    "speciation",
    "combined",
    s_comb_res$clusters,
    method = "logistic"
  )
t_comb_corr_res <-
  plot_cluster_traitcorr(
    t_comb_res$pvals,
    t_comb_res$coefficient,
    "transfer",
    "combined",
    t_comb_res$clusters,
    method = "logistic"
  )
l_comb_corr_res <-
  plot_cluster_traitcorr(
    l_comb_res$pvals,
    l_comb_res$coefficient,
    "loss",
    "combined",
    l_comb_res$clusters,
    method = "logistic"
  )

# For now (and for the sake of simplicity), let's just focus on the clusters
# inferred from all profiles
# Combine the results from the above analyses, creating a table that includes
# our ranking of gene families based on their mean-shifts
final_clust_associations <-
  rbind(
    s_comb_corr_res$correlation_summary,
    t_comb_corr_res$correlation_summary,
    l_comb_corr_res$correlation_summary
  )

# Now, go ahead and write this out to file, in addition to the four summary
# plots, and one that combined the four into separate panels
dir.create(
  "./chelicerate-results/host_detection_suppression_association_results/",
  showWarnings = FALSE
)
write.table(
  final_clust_associations,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  file = "chelicerate-results/host_detection_suppression_association_results/combined_profile_cluster_host_detection_suppression_associations.tsv"
)

comb_cls_plt <-
  plot_grid(
    s_comb_corr_res$correlation_plot,
    t_comb_corr_res$correlation_plot,
    l_comb_corr_res$correlation_plot,
    ncol = 2,
    nrow = 2,
    labels = c("a.", "b.", "c.")
  )

ggsave(comb_cls_plt,
  height = 12,
  width = 14,
  file = "./chelicerate-results/host_detection_suppression_association_results/combined_profile_cluster_host_detection_suppression_associations.pdf"
)
ggsave(
  s_comb_corr_res$correlation_plot,
  height = 6,
  width = 7,
  file = "./chelicerate-results/host_detection_suppression_association_results/speciation_profile_cluster_host_detection_suppression_associations.pdf"
)
ggsave(
  t_comb_corr_res$correlation_plot,
  height = 6,
  width = 7,
  file = "./chelicerate-results/host_detection_suppression_association_results/transfer_profile_cluster_host_detection_suppression_associations.pdf"
)
ggsave(
  l_comb_corr_res$correlation_plot,
  height = 6,
  width = 7,
  file = "./chelicerate-results/host_detection_suppression_association_results/loss_profile_cluster_host_detection_suppression_associations.pdf"
)

# Now take the most strongly predictive clusters and assess which of the
# constituent gene families are the strongest predictors.
# Note, we still are treating these as univariate analyses.
# Since we will be conducting analyses per-gene family, use the
# untransformed count data per-species - we will conduct phylogenetic logistic
# regressions using the speciation, loss, and transfer counts. As before,
# we do not look at duplications within species, as these are not observed in
# terminal branches.

# Create empty dataframes to hold each count variable, one row
# per gene family
terminal <-
  which(gf_node_events[[1]]$node %in% spptree$tip.label)
s_counts <-
  matrix(
    NA,
    nrow = length(gf_node_events),
    ncol = nrow(gf_node_events[[1]][terminal, ]),
    dimnames = list(
      names(gf_node_events),
      gf_node_events[[1]]$node[terminal]
    )
  )
sl_counts <- s_counts
d_counts <- s_counts
t_counts <- s_counts

# Populate these empty count dataframes with the observed counts of each
# respective event type
for (i in 1:length(gf_node_events)) {
  t_dat <- t(gf_node_events[[i]][, -1])
  s_counts[i, ] <- t_dat[1, terminal]
  sl_counts[i, ] <- t_dat[2, terminal]
  d_counts[i, ] <- t_dat[3, terminal]
  t_counts[i, ] <- t_dat[4, terminal]
}

# Identify which gene families we should do a post-hoc assessment of per-family
# association with host detection suppression
post_hoc_clusts <-
  final_clust_associations[which(final_clust_associations$signif_level %in%
                                   c("Yes & Top 10% Coef.", "Yes & Top 5% Coef.")), ]

# subset by event type
s_final_clusts <-
  post_hoc_clusts[which(post_hoc_clusts$profile_type == "speciation"), ]
t_final_clusts <-
  post_hoc_clusts[which(post_hoc_clusts$profile_type == "transfer"), ]
l_final_clusts <-
  post_hoc_clusts[which(post_hoc_clusts$profile_type == "loss"), ]

# Start with speciation
s_gf_clust_res <-
  test_fam_corrs_per_clust(
    count_dat = s_counts,
    clusters = s_final_clusts,
    phylo_profiling_res = all_dists_clusts,
    spptree = spptree
  )
# transfer counts
t_gf_clust_res <-
  test_fam_corrs_per_clust(
    count_dat = t_counts,
    clusters = t_final_clusts,
    phylo_profiling_res = all_dists_clusts,
    spptree = spptree
  )
# and loss counts
l_gf_clust_res <-
  test_fam_corrs_per_clust(
    count_dat = sl_counts,
    clusters = l_final_clusts,
    phylo_profiling_res = all_dists_clusts,
    spptree = spptree
  )

# Now save everything out to file
save_per_fam_associations(
  per_fam_results = s_gf_clust_res,
  output_dir = "chelicerate-results/host_detection_suppression_association_results/combined_profile_clusters/",
  event_type = "speciation"
)
save_per_fam_associations(
  per_fam_results = t_gf_clust_res,
  output_dir = "chelicerate-results/host_detection_suppression_association_results/combined_profile_clusters/",
  event_type = "transfer"
)
save_per_fam_associations(
  per_fam_results = l_gf_clust_res,
  output_dir = "chelicerate-results/host_detection_suppression_association_results/combined_profile_clusters/",
  event_type = "loss"
)

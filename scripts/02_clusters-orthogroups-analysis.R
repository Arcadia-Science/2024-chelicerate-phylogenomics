library(tidyverse)

######################################
# Organize files from trait association analysis, NovelTree run, and annotations
######################################

# Load in orthogroups TSV
# this file links, for each orthogroup, the locus tag of protein belonging in that
# orthogroup. This file is within the NovelTree results downloaded from Zenodo (chelicerata-v1-10062023)
# at this location: chelicerata-v1-10062023/orthofinder/complete_dataset/Results_Inflation_1.5/Orthogroups/Orthogroups.tsv
# Adjust your path below as needed

orthogroups <-
  read_tsv(
    "chelicerata-v1-10062023/orthofinder/complete_dataset/Results_Inflation_1.5/Orthogroups/Orthogroups.tsv"
  )

# pivot longer and separate rows in locus tag by "," so each row is a single
# entry for a locus tag to then combine with annotation information
orthogroups_long <- orthogroups %>%
  pivot_longer(!Orthogroup, names_to = "species", values_to = "locus_tag") %>%
  drop_na()

orthogroups_separated <- orthogroups_long %>%
  separate_rows(locus_tag, sep = ",\\s*") %>%
  mutate(orthogroup = Orthogroup) %>%
  select(orthogroup, species, locus_tag)

# Annotation files are too big for github.
# Update the path below to wherever you have stored these files.
# Downloaded from Zenodo (annotated/)

annotation_directory <- "./annotated"
annotation_files <-
  list.files(annotation_directory,
    pattern = "*.tsv",
    full.names = TRUE
  )

all_annotations <- map_df(annotation_files, function(file_name) {
  read_tsv(file_name) %>%
    mutate(source_file = basename(file_name))
})

annotations_df <- all_annotations %>%
  mutate(species_name = gsub("_combined.tsv", "", source_file)) %>%
  mutate(locus_tag = gene_name) %>%
  select(
    species_name,
    locus_tag,
    Length,
    egg_Description,
    KO,
    KO_definition,
    deepsig_feature,
    deepsig_start,
    deepsig_end
  )

# orthogroups with annotations
orthogroup_annotations <-
  left_join(orthogroups_separated, annotations_df) %>%
  select(-species_name)

# Cluster results
# Read in all cluster results for each noveltree model -
# speciation, transfer, loss, in subdirectories.
# These files come from the host_detection_suppression_association directory,
# where significant clusters of gene families were identified based on different
# models
profile_clusters_dir <-
  "./chelicerate-results/host_detection_suppression_association_results/combined_profile_clusters"

profile_files <-
  list.files(
    profile_clusters_dir,
    "*_per_family_correlation_results.tsv",
    full.names = TRUE,
    recursive = TRUE
  )

all_profiles <- map_df(profile_files, function(file_name) {
  read_tsv(file_name) %>%
    mutate(filename = paste(basename(dirname(file_name)),
      basename(file_name),
      sep = "_"
    ))
})

# Model comes from the filename, bug where "speciation" is carried through for
# profile_type for all files but is incorrect
all_profiles_df <- all_profiles %>%
  mutate(cluster = gsub(
    "_per_family_correlation_results.tsv",
    "",
    filename
  )) %>%
  mutate(model = cluster) %>%
  mutate(model = gsub("_.+$", "", model)) %>%
  mutate(cluster = gsub("^.*?_.*?_", "", cluster)) %>%
  mutate(orthogroup = gene_family_id) %>%
  select(
    cluster,
    model,
    profile_type,
    orthogroup,
    signif_level,
    signif_fdr,
    coefficient,
    pval
  )

# Select specific clusters that are the top 10% and only those under the
# speciation model
file_path <-
  "chelicerate-results/host_detection_suppression_association_results/combined_profile_cluster_host_detection_suppression_associations.tsv"

final_clust_associations <- read.table(
  file_path,
  sep = "\t",
  header = TRUE,
  quote = "",
  stringsAsFactors = FALSE
)

# Filter based on signif_level and positive coefficient
cluster_list <- final_clust_associations %>%
  filter(signif_level %in% c("Yes & Top 10% Coef.", "Yes & Top 5% Coef.") &
           coefficient > 0)

speciation_select_clusters <- all_profiles_df %>%
  filter(cluster %in% paste0("cluster_", cluster_list$cluster_id)) %>%
  filter(model == "speciation")

speciation_select_clusters %>%
  select(orthogroup) %>%
  unique() %>%
  count() # 832 orthogroups from the top 10% clusters that are positively
# predictive for host detection suppression under the speciation model

# Filter for significant clusters, which meets pval < 0.05
signf_speciation_select_clusters <- speciation_select_clusters %>%
  filter(signif_level != "No")

signf_speciation_select_clusters %>%
  select(orthogroup) %>%
  unique() %>%
  count() # 86 orthogroups that are significantly, positively associated with
# host detection suppression

# combine clusters/orthogroups with orthogroups/locus tags
signf_clusters_orthogroups_annotations <-
  left_join(signf_speciation_select_clusters,
    orthogroup_annotations,
    relationship = "many-to-many"
  ) # 3461 proteins

# counts of genes across species
signf_clusters_orthogroups_counts <-
  signf_clusters_orthogroups_annotations %>%
  group_by(species, cluster, orthogroup) %>%
  count() %>%
  pivot_wider(names_from = "species", values_from = "n") %>%
  mutate(across(where(is.numeric), ~ coalesce(., 0)))

# filter for Amblyomma locus tags to check differential expression profiles
amblyomma_signf_clusters_orthogroups <-
  signf_clusters_orthogroups_annotations %>%
  filter(species == "Amblyomma-americanum")

# Write out main files.
# Only writing out files for the filtered down sets of orthogroups/proteins that
# are within the top 10% clusters, are positively and significantly associated
# with host detection suppression under the speciation model.

# Annotation table
output_path <- "chelicerate-results/clusters-orthogroups-analysis"

if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

write.table(
  signf_clusters_orthogroups_annotations,
  "chelicerate-results/clusters-orthogroups-analysis/top-positive-significant-clusters-orthogroups-annotations.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Counts table
write.table(
  signf_clusters_orthogroups_counts,
  "chelicerate-results/clusters-orthogroups-analysis/top-positive-significant-clusters-orthogroups-counts.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Amblyomma cluster genes, also filtered for positive coefficient and speciation
# model.
write.table(
  amblyomma_signf_clusters_orthogroups,
  "chelicerate-results/clusters-orthogroups-analysis/top-positive-significant-amblyomma-clusters-orthogroups.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

######################################
# Filter down clusters based on %s of secreted proteins, those that are
# expressed in Amblyomma salivary transcriptome
######################################

# Breakdown of % of proteins that are secreted in each orthogroup.
# Filter to get orthogroups where at least 50% of the proteins in that group are
# predicted to have signal peptides by deepsig

majority_secreted_orthogroups <-
  signf_clusters_orthogroups_annotations %>%
  group_by(orthogroup) %>%
  mutate(total_count = n()) %>%
  ungroup() %>%
  group_by(orthogroup, deepsig_feature, total_count) %>%
  summarise(n = n()) %>%
  mutate(percentage = (n / total_count) * 100) %>%
  filter(deepsig_feature == "Signal peptide") %>%
  arrange(desc(percentage)) %>%
  filter(percentage > 50)

majority_secreted_orthogroups_list <-
  majority_secreted_orthogroups %>%
  pull(orthogroup) # 15 orthogroups after secretion filter

# Expressed amblyomma proteins in orthogroups
# Per orthogroup, get number of Amblyomma genes in that orthogroup
# filter by expression = 1, percentage over total genes in orthogroup to get
# number that is expressed or "present" in the Isoseq salivary transcriptome

# namemap of Isoseq reads to transcripts
isoseq_namemap <-
  read.csv("metadata/SRR19070014_orthofuser_namemap.csv")
# presence/absence of expression for each transcript. On Zenodo as tx2gene.tsv.
transcript_mapping <-
  read.table(
    "chelicerate-results/clusters-orthogroups-analysis/tx2gene.tsv",
    col.names = c("transcript", "presence")
  )

# Filtered transcript mapping to only those with present genes
transcript_mapping_present <- transcript_mapping %>%
  filter(presence != "-1") %>%
  mutate(presence = gsub(".TU", "", presence)) %>%
  mutate(gene = gsub("_", "-", presence)) %>%
  select(transcript, gene)

transcripts_mapped_list <- transcript_mapping_present %>%
  select(gene) %>%
  unique() %>%
  pull(gene)

amblyomma_clusters_expression <-
  amblyomma_signf_clusters_orthogroups %>%
  mutate(gene = gsub("Amblyomma-americanum_", "", locus_tag)) %>%
  mutate(gene = gsub(".model", "", gene)) %>%
  mutate(expression = if_else(gene %in% transcripts_mapped_list, 1, 0))

expressed_signal_peptides <- amblyomma_clusters_expression %>%
  filter(expression == 1) %>%
  filter(deepsig_feature == "Signal peptide")

majority_amblyomma_orthogroups_expressed <-
  amblyomma_clusters_expression %>%
  group_by(orthogroup) %>%
  mutate(total_count = n()) %>%
  ungroup() %>%
  group_by(orthogroup, expression, total_count) %>%
  summarise(n = n()) %>%
  filter(expression == 1) %>%
  mutate(percent_expressed = (n / total_count) * 100) %>%
  select(orthogroup, expression, n, total_count, percent_expressed) %>%
  arrange(desc(percent_expressed)) %>%
  filter(percent_expressed > 25)

# filter the expressed genes based on secretion filter above
filtered_orthogroups_by_secretion_expression <-
  majority_amblyomma_orthogroups_expressed %>%
  filter(orthogroup %in% majority_secreted_orthogroups_list) %>%
  left_join(majority_secreted_orthogroups, by = "orthogroup") %>%
  mutate(total_amblyomma_genes = total_count.x) %>%
  mutate(total_amblyomma_genes_expressed = n.x) %>%
  mutate(total_proteins_in_orthogroup = total_count.y) %>%
  mutate(total_secreted_proteins = n.y) %>%
  mutate(percent_secreted = percentage) %>%
  ungroup() %>%
  select(
    orthogroup,
    total_amblyomma_genes,
    total_amblyomma_genes_expressed,
    percent_expressed,
    total_proteins_in_orthogroup,
    total_secreted_proteins,
    percent_secreted
  )

filtered_orthogroups_list <-
  filtered_orthogroups_by_secretion_expression %>%
  pull(orthogroup) # 13 orthogroups left

# filter the annotations and counts tables to this filtered set of orthogroups
filtered_annotations <- signf_clusters_orthogroups_annotations %>%
  filter(orthogroup %in% filtered_orthogroups_list)

filtered_counts <- signf_clusters_orthogroups_counts %>%
  filter(orthogroup %in% filtered_orthogroups_list)

# Further filter for groups that have at least 10 total counts across 6 or more
# tick species.
# get list of ticks from the chelicerate metadata used for phylo profiling

chelicerate_metadata <-
  read_tsv("metadata/chelicerate-species-metadata.txt")

tick_metadata <- chelicerate_metadata %>%
  filter(str_detect(common_name, "(?i)tick"))

tick_list <- tick_metadata %>%
  pull(id)

# filter the clusters table to only include tick species
tick_columns <- c("orthogroup", tick_list)
tick_columns <- tick_columns %>% intersect(names(filtered_counts))

# filter for orthogroups with non-zero counts for at least 6 tick species.
filtered_tick_orthogroups_6species <- filtered_counts %>%
  ungroup() %>%
  select(all_of(tick_columns)) %>%
  mutate(across(-orthogroup, ~ if_else(. > 0, 1, 0))) %>% # present vs absent
  mutate(tick_count = rowSums(select(., -orthogroup))) %>%
  filter(tick_count >= 6) %>%
  pull(orthogroup)

# filter for orthogroups with non-zero counts for at least 6 tick species where
# there are at least 10 counts across tick species.
filtered_counts_filtered_6species10counts <- filtered_counts %>%
  filter(orthogroup %in% filtered_tick_orthogroups_6species) %>%
  select(all_of(tick_columns)) %>%
  rowwise() %>%
  mutate(sum_tick_columns = sum(c_across(everything()), na.rm = TRUE)) %>%
  filter(sum_tick_columns >= 10) %>%
  select(-sum_tick_columns)

species6counts10orthogroups_list <-
  unique(filtered_counts_filtered_6species10counts$orthogroup)

# Filtered dataframes based on 50% of proteins in orthogroup are secreted, have
# to have an Amblyomma representative and 25% expressed in salivary
# transcriptome, and orthogroup has non-zero counts for at least 6 tick species
# where there are at least 10 counts across tick species.

final_filtered_counts <- filtered_counts %>%
  filter(orthogroup %in% species6counts10orthogroups_list)

final_filtered_counts_tickonly <- filtered_counts %>%
  filter(orthogroup %in% species6counts10orthogroups_list) %>%
  select(orthogroup, all_of(tick_columns))

final_filtered_annotations <- filtered_annotations %>%
  filter(orthogroup %in% species6counts10orthogroups_list)

# write out files
write.table(
  final_filtered_annotations,
  "chelicerate-results/clusters-orthogroups-analysis/filtered-orthogroups-annotations-for-tm-prediction.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  amblyomma_clusters_expression,
  "chelicerate-results/clusters-orthogroups-analysis/amblyomma_clusters_expression.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  final_filtered_counts,
  "chelicerate-results/clusters-orthogroups-analysis/filtered-orthogroups-counts-for-tm-prediction.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# plots and analysis
final_filtered_annotations %>%
  ggplot(aes(x = Length)) +
  geom_histogram()

final_filtered_annotations %>%
  filter(deepsig_feature == "Signal peptide") %>%
  ggplot(aes(x = Length)) +
  geom_histogram()

locus_tags_list <- final_filtered_annotations %>%
  filter(deepsig_feature == "Signal peptide") %>%
  pull(locus_tag)

write.table(
  locus_tags_list,
  "metadata/filtered-host-detection-association-protein-locus-tags-for-tm-prediction.txt",
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  col.names = FALSE
)

######################################
# The next 3 steps were done outside of this R script:
# Set of concatenated proteins from all chelicerate species pulled from those
# that were directly input to the NovelTree run
# 1. pulled out protein seqs with this command/script:
# python3 ./scripts/grab_proteins_from_locus_tags.py \
#   2024-06-24-all-chelicerate-noveltree-proteins.fasta \
#   ./metadata/filtered-host-detection-association-protein-locus-tags-for-tm-prediction.txt
#   ./chelicerate-results/deeptmhmm-results/filtered-proteins-for-tm-prediction.fasta
# The FASTA file 2024-06-24-all-chelicerate-noveltree-proteins.fasta is
# available on Zenodo.
# 2. Transmembrane domain predictions from deepTMHMM tool with
# chelicerate-results/deeptmhmm-results/filtered-proteins-for-tm-prediction.fasta
# on this webserver https://dtu.biolib.com/DeepTMHMM
# 3. Downloaded results from webserver and used this command to make table:
# awk -F'\\| ' '/^>/{gsub(/^>/, "", $1); print $1 "\t" $2}' \
# chelicerate-results/deeptmhmm-results/predicted_topologies.3line > \
# chelicerate-results/deeptmhmm-results/tmpredictions.txt
# 4. We only had proteins left that had just a SP prediction, no TM or TM+SP,
# but the code below would remove them if needed.
######################################

# select proteins in orthogroups with SP predictions, toss out SP+TM or just TM
transmembrane_predictions <-
  read_tsv("chelicerate-results/deeptmhmm-results/tmpredictions.txt",
    col_names = FALSE
  ) %>%
  mutate(locus_tag = X1) %>%
  mutate(tm_prediction = X2) %>%
  select(locus_tag, tm_prediction)

filtered_annotations_tm_predictions <-
  final_filtered_annotations %>%
  filter(deepsig_feature == "Signal peptide") %>%
  left_join(transmembrane_predictions)

orthogroup_tmpreds <- filtered_annotations_tm_predictions %>%
  filter(species %in% tick_list) %>%
  group_by(orthogroup, tm_prediction) %>%
  count()

# Identify orthogroups with "SP+TM" or "TM" in tm_prediction
orthogroups_to_exclude <- orthogroup_tmpreds %>%
  filter(tm_prediction %in% c("SP+TM", "TM")) %>%
  pull(orthogroup)

# Filter the original dataframe to exclude these orthogroups
# Also filter to remove rows where signif_level=No
filtered_orthogroup_tmpreds <-
  filtered_annotations_tm_predictions %>%
  filter(!orthogroup %in% orthogroups_to_exclude) %>%
  filter(signif_level != "No")

# counts of proteins in only tick species in the 10 remaining orthogroups
filtered_tick_SP_proteins_counts <- filtered_orthogroup_tmpreds %>%
  group_by(species, orthogroup) %>%
  count() %>%
  pivot_wider(names_from = "species", values_from = "n") %>%
  mutate(across(where(is.numeric), ~ coalesce(., 0)))

filtered_annotations_tm_predictions %>%
  group_by(orthogroup) %>%
  select(orthogroup, locus_tag) %>%
  count() # 10 orthogroups with 275 proteins to select candidates from

# Merge expression data into filtered_orthogroup_tmpreds for use in candidate
# selection.
# Create a vector to store the expression values
expression_values <- rep(NA, nrow(filtered_orthogroup_tmpreds))

# Find the indices where locus_tag is present in amblyomma_clusters_expression
match_indices <-
  match(
    filtered_orthogroup_tmpreds$locus_tag,
    amblyomma_clusters_expression$locus_tag
  )

# Update the expression values at corresponding indices
expression_values[!is.na(match_indices)] <-
  amblyomma_clusters_expression$expression[match_indices[!is.na(match_indices)]]

# Append the expression column to filtered_orthogroup_tmpreds
filtered_orthogroup_tmpreds$expression <- expression_values

all_locus_tags <- filtered_orthogroup_tmpreds %>%
  pull(locus_tag)

# write out final filtered annotation table, counts table, locus tags

# final annotation table
write.table(
  filtered_orthogroup_tmpreds,
  "chelicerate-results/clusters-orthogroups-analysis/final-filtered-SP-orthogroups-annotations.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# final counts table
write.table(
  filtered_tick_SP_proteins_counts,
  "chelicerate-results/clusters-orthogroups-analysis/final-filtered-SP-orthogroups-counts.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# final locus tag list
write.table(
  all_locus_tags,
  "metadata/final-filtered-host-detection-association-SP-protein-locus-tags.txt",
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  col.names = FALSE
)

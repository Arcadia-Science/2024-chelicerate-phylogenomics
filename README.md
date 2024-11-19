# 2024-chelicerate-phylogenomics

[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/projects/miniconda/en/latest/)

## Purpose

We used phylogenomics to investigate patterns of gene family evolution across ticks and other chelicerates, which include a diverse array of parasites. We used phylogenetic profiling and trait-association tests to predict gene families that may enable parasitic species to feed on hosts undetected for prolonged periods (>1 day). This repository contains the scripts used to generate phylogenetic profiles for each gene family using the outputs of [NovelTree](https://doi.org/10.57844/arcadia-9602-3351), cluster families by profile similarity, and identify clusters that predict suppression of host detection.

## Installation and Setup

This repository uses conda to manage software environments and installations. You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/projects/miniconda/en/latest/). After installing conda and [mamba](https://mamba.readthedocs.io/en/latest/), run the following command to create the pipeline run environment.

```{bash}
mamba env create -n chelicerate --file envs/dev.yml
conda activate chelicerate
```
The `phyr` package is not available through conda but is required by this repo. To install the package, run the following from within your activated environment to open R:

```{bash}
R
```
Once you are in an R session, run the following:

```{bash}
# Install remotes if not already installed
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}

# Set repositories
options(repos = c(
  phyr = 'https://daijiang.r-universe.dev',
  CRAN = 'https://cloud.r-project.org'
))

# Install specific version of phyr
remotes::install_version("phyr", version = "1.1.2", repos = "https://daijiang.r-universe.dev")
```

Quit the R session by running `quit()` and then enter `n` to not save the workspace.

This repository contains R scripts which can be run in Rstudio, which can be installed following the instructions [here](https://posit.co/download/rstudio-desktop/). Once Rstudio is installed, run the following from command line to open Rstudio from your activated conda environment:

```{bash}
open -a Rstudio
```

You can now move on to running the necessary [scripts](scripts) as described below. Make sure your working directory is set to this repo when in Rstudio.

## Data

Some of the data needed to run these scripts is too large for GitHub and has been deposited on [Zenodo](10.5281/zenodo.14113178). This includes:
- Outputs of the NovelTree run that are used as input into [phylo_profiling_genefam_evol_counts.R](scripts/phylo_profiling_genefam_evol_counts.R): chelicerata-v1-10062023.zip
- Chelicerate gene annotations used to do orthogroup filtering in [clusters_orthogroups-analysis.R](scripts/clusters_orthogroups-analysis.R): annotated.zip 
- Presence/absence of expression for each transcript in *A. americanum* salivary transcriptome used in [clusters_orthogroups-analysis.R](scripts/clusters_orthogroups-analysis.R): tx2gene.tsv
- Chelicerate protein sequences used as input in [clusters_orthogroups-analysis.R](scripts/clusters_orthogroups-analysis.R) and needed for [deepTMHMM webserver](https://dtu.biolib.com/DeepTMHMM) prediction: 2024-06-24-all-chelicerate-noveltree-proteins.fasta

## Overview

### The [`clusters_orthogroups-analysis.R`](scripts/clusters_orthogroups-analysis.R) script will perform the following:

1. Conduct phylogenetic profiling to identify groups (clusters) of gene families that have very similar patterns of gene duplication, transfer, and loss (gene family evolutionary events)
    - These results are contained within [`chelicerate-results/umap-layout-cluster-ids/*_profile_clusters_annotations.tsv`](chelicerate-results/umap-layout-cluster-ids/)
        - If you ever want to visualize these profile clusters, there are static plots (pdfs) of the UMAP projections of each gene family, and for each event type, with each gene family/point colored according to their (leiden) cluster ID (e.g. [`speciation_profile_umap_best_clusters.pdf`](chelicerate-results/plots/speciation_profile_umap_best_clusters.pdf))
        - Additionally, there are interactive HTMLs where you can hover over to see the name of the gene family and corresponding annotation, as well as see the cluster assignments from the other clustering algorithms used (e.g. [`speciation_profile_umap_clustering.html`](chelicerate-results/plots/speciation_profile_umap_clustering.html)
2. Within those groups/profile clusters, ask whether the mean counts of each event type (excluding duplication, as these events are occur deeper in the tree) is associated with suppression of host detection (i.e. a logistic regression, asking whether suppression of host detection, a binary trait, is predicted by the mean count of each event type).
    - These results are contained within [`chelicerate-results/host_detection_suppression_association_results/combined_profile_cluster_host_detection_suppression_associations.tsv`](chelicerate-results/host_detection_suppression_association_results/combined_profile_cluster_host_detection_suppression_associations.tsv)
    - The pdfs (e.g. [`loss_profile_cluster_detection_suppression_associations.pdf`](chelicerate-results/host_detection_suppression_association_results/loss_profile_cluster_host_detection_suppression_associations.pdf) are plots of these results, where the Y-axis is the fitted coefficient of the logistic regression (here, log-odds, since we’re doing logistic regression), and each point is a profile cluster (i.e. cluster of gene families)
    - These are ordered by their coefficient, and the most significantly associated clusters are labeled according to their cluster ID (for that specific event type)
3. Then, for those profile clusters that have significant associations between suppression of host detection and counts of each event type, conduct a similar association test, this time within each profile cluster
    - Here, we ask (for each gene family) if the counts of each event type observed in each species predict whether they suppress detection by the host (again, a logistic regression).
    - These results are contained within [`chelicerate-results/host_detection_suppression_association_results/combined_profile_clusters`](chelicerate-results/host_detection_suppression_association_results/combined_profile_clusters/)
        - The subdirectories `loss_associations`, `speciation_associations`, and `transfer_associations` contain the association test results using each respective event type, and the tables/plots are interpreted the same as for the profile cluster association test results
    - The key difference between this step and the association tests done in `1.` is that now the unit being tested is the specific gene family - the annotations for which can be obtained from [`umap-layout-cluster-ids/*_profile_clusters_annotations.tsv`](chelicerate-results/umap-layout-cluster-ids/)
    - One key note: the per-family functional annotations in these `*annotations.tsv` files do not differ, but the different cluster IDs will.
        - Because the statistical tests were conducted on each individual event type, you’d actually need to look at the corresponding table to pull out those IDs (e.g. [`speciation_profile_cluster_host_detection_suppression_annotations.tsv`](chelicerate-results/umap-layout-cluster-ids/speciation_profile_clusters_annotations.tsv))
    - Regardless, for now it's best to focus on the `leiden_cluster_id`, since these are the clusters the statistical tests were conducted on.

#### A few quick notes on interpretation:
1. "Speciations" in the context of gene family evolution captures multiple features of gene family evolution and is closely ties to gene copy number - it thus is a correlate of the other event types - gene duplication, transfer, and loss.
    - For instance, increasing gene duplication increases the observed number of speciation events, whereas elevated rates of gene loss will decrease the number of speciation events.
    - For these reasons, it is very information rich and variable among species, making it particularly suitable for these statistical analyses of trait association.
    - Thus, if you were to choose one event type to focus on, it's recommended you focus on speciations.
2. It's crucially important that you interpret estimated coefficients in the context of the underlying biology for each event type.
    - For instance, a positive association between gene losses and host detection suppression means that these gene families likely _prevent_ or _limit_ the capacity of species (in some way) to suppress detection by the host.
        - A reasonable hypothesis here would be that it was selectively beneficial for those genes to be lost in species that suppress host detection.
        - Instead, the gene families that exhibit a negative association between gene losses and suppression of host detection would be hypothesized to be important for suppression of host detection, as it may be selectively advantageous for them to be retained in these species.
    - This is in contrast to the count of speciations, which again is a strong correlate of gene copy number.
        - Here, a positive association between speciation count and suppression of host detection would suggest that having a greater number of or diversity in these genes contributes to the capacity of species to suppress host detection.

### Analysis of Clusters of Orthogroups Positively Associated with Suppression of Host Detection with [`scripts/clusters-orthogroups-analysis.R`](scripts/clusters-orthogroups-analysis.R)

We next sought to narrow down our pool of orthogroups to families of putative salivary effector proteins that could be direct mediators of the suppression of host detection. Gathering the different files together and filtering for these orthogroups is performed in the script `scripts/clusters-orthogroups-analysis.R`. First the identified clusters of orthogroups are combined with orthgroup information for what proteins are in that orthogroup and the annotations of those proteins. From combining this information, tables were created of both annotation descriptions for each protein and counts of orthogroups for each species. Then orthogroups were filtered based on:
1. Selecting top 10% correlated clusters of orthogroups under the speciation model and are positively and significantly associated with suppression of host detection. Then within these clusters, removing any orthogroup that was individually negatively associated with suppression of host detection, and also any orthogroup with an initial p value > 0.05.

    Orthogroups/proteins at this level of filtering included 76 orthogroups and 3414 proteins, which are described in [`chelicerate-results/clusters-orthogroups-analysis/top-positive-significant-clusters-orthogroups-annotations.tsv`](chelicerate-results/clusters-orthogroups-analysis/top-positive-significant-clusters-orthogroups-annotations.tsv).

    It is important to note that after FDR testing no individual orthogroup was individually significantly associated with suppression of host detection. From the orthogroups that remained, we performed additional filtering as follows below.

2. Filtering for orthogroups where at least 50% of the proteins in that orthogroup were predicted to contain signal peptides by DeepSig.
3. Requiring that the orthogroup contain an *Amblyomma* representative, and at least 25% of the *Amblyomma* genes in that orthogroup were counted as "expressed" in the IsoSeq salivary transcriptome (presence/absence, where presence counts as "expressed"). After this filter we had 12 orthogroups remaining.
4. Orthogroup must have at least 10 total counts across 6 or more tick species. After this filter we had 10 orthogroups remaining.
5. Removed orthogroups with proteins predicted to have transmembrane domain predictions with the [deepTMHMM webserver](https://dtu.biolib.com/DeepTMHMM). Protein sequences were extracted from the full set of chelicerate proteins with `scripts/grab_proteins_from_locus_tags.py`, creating [this file](chelicerate-results/deeptmhmm-results/filtered-proteins-for-tm-prediction.fasta) that was used for prediction.
We ran the concatenated set of proteins in the webserver since the command-line program has license restrictions for commercial entities. We kept proteins designated as only signal peptide (SP), meaning we removed proteins that had predicted transmembranes (TM) along with those designated as transmembrane and signal peptide (SP + TM). From the remaining 10 orthogroups, all proteins in all orthogroups only had deepTHMM predictions with `SP`.

This left us with 10 orthogroups and a total of 275 proteins, which are described in [`chelicerate-results/clusters-orthogroups-analysis/final-filtered-SP-orthogroups-annotations.tsv`](chelicerate-results/clusters-orthogroups-analysis/final-filtered-SP-orthogroups-annotations.tsv).

## Contributing

See how we recognize [feedback and contributions to our code](https://github.com/Arcadia-Science/arcadia-software-handbook/blob/main/guides-and-standards/guide-credit-for-contributions.md).

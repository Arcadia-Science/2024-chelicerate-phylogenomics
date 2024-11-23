# Load required libraries or install them if not present
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Set options to avoid installation prompts
options(repos = c(CRAN = "https://cran.rstudio.com"))
options(Ncpus = parallel::detectCores())

# Bioconductor package versions are pulled from the Bioconductor version.
# The versions of ggtree and cogeqc that we need were both available in
# Bioconductor version 3.19.
bioc_packages <- c(
  "ggtree" = "3.19", # for ggtree 3.12.0
  "cogeqc" = "3.19" # for cogeqc 1.8.0
)

cran_packages <- c(
  # required packages w/o the correct version or with version conflicts on Conda
  "uwot" = "0.2.2",
  "Matrix" = "1.7-0",
  "phyr" = "1.1.0",
  "phytools" = "2.1-1",
  "RcppArmadillo" = "0.12.8.3.0",
  "pbapply" = "1.7-2",
  "plyr" = "1.8.9",
  "phylolm" = "2.6.2",
  # installed by conda but double check the versions weren't changed by bioconductor
  "dbscan" = "1.1.12",
  "ape" = "5.8",
  "cowplot" = "1.1.3",
  "dbscan" = "1.1-12",
  "ggforce" = "0.4.2",
  "ggnewscale" = "0.4.10",
  "ggrepel" = "0.9.5",
  "gridExtra"="2.3",
  "igraph"="2.0.3",
  "jpeg"="0.1-10",
  "latticeExtra"="0.6-30",
  "logistf"="1.26.0",
  "mvtnorm"="1.3-2",
  "plotly"="4.10.4",
  "Rcpp"="1.0.12",
  "reshape2"="1.4.4",
  "tidyverse"="2.0.0"
)

install_cran_packages <- function(packages) {
  for (pkg in names(packages)) {
    required_version <- packages[pkg]
    # Packages with a dash in their version are switched to only periods once
    # installed. Switch the version for this function so things won't be
    # unnecessarily re-installed.
    required_version <- gsub("\\-", "\\.", required_version)
    # Check if package is installed
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(sprintf("Package '%s' is not installed. Installing version %s.", pkg, required_version))
      remotes::install_version(pkg, version = required_version, upgrade = "never")
    } else {
      # Check the version of the installed package
      current_version <- as.character(packageVersion(pkg))
      if (current_version != required_version) {
        message(sprintf("Package '%s' is installed but at version %s. Required version is %s. Updating...", pkg, current_version, required_version))
        remotes::install_version(pkg, version = required_version, upgrade = "never")
      } else {
        message(sprintf("Package '%s' is already installed with the correct version %s.", pkg, current_version))
      }
    }
  }
}

# Function to install Bioconductor packages with specific versions
install_bioc_packages <- function(packages) {
  for (pkg in names(packages)) {
    BiocManager::install(pkg, version = packages[pkg], ask = FALSE, update = FALSE)
  }
}

# Install packages
# Bioconductor packages need to be installed first as BioManager will twiddle
# with the versions of other packages
BiocManager::install(version = "3.19", ask=FALSE)
install_bioc_packages(bioc_packages)
install_cran_packages(cran_packages)

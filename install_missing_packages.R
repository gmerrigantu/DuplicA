# Install missing R packages for DuplicA

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install Bioconductor packages
BiocManager::install("clusterProfiler")
BiocManager::install("Biostrings")
BiocManager::install("biomaRt")

cat("Packages installed successfully\n")

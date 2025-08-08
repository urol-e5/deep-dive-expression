# Install required packages for miRNA analysis

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install CRAN packages
cran_packages <- c("dplyr", "tidyr", "ggplot2", "stringr", "reshape2", 
                   "corrplot", "patchwork", "viridis", "RColorBrewer", "readr",
                   "igraph", "ggraph", "tidygraph", "network", "sna", "intergraph")

for (pkg in cran_packages) {
  if (!require(pkg, character.only = TRUE)) {
    tryCatch({
      install.packages(pkg, dependencies = TRUE)
      cat("Installed:", pkg, "\n")
    }, error = function(e) {
      cat("Failed to install:", pkg, "-", e$message, "\n")
    })
  } else {
    cat("Package already installed:", pkg, "\n")
  }
}

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (requireNamespace("BiocManager", quietly = TRUE)) {
  bioc_packages <- c("edgeR", "limma")
  for (pkg in bioc_packages) {
    if (!require(pkg, character.only = TRUE)) {
      tryCatch({
        BiocManager::install(pkg, ask = FALSE)
        cat("Installed Bioconductor package:", pkg, "\n")
      }, error = function(e) {
        cat("Failed to install Bioconductor package:", pkg, "-", e$message, "\n")
      })
    } else {
      cat("Bioconductor package already installed:", pkg, "\n")
    }
  }
}

cat("Package installation completed!\n")

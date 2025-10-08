#!/usr/bin/env Rscript

# REDAC Setup Script
# This script installs all required packages for the REDAC Shiny application
# Authors: Francesco Russo, IEOMI-CNR, Napoli; Giovanni Maria De Filippis, University of Naples Federico II

cat("=============================================================\n")
cat("         REDAC Package Installation Script\n")
cat("=============================================================\n")
cat("Installing required packages for REDAC application...\n\n")

# Function to check if a package is installed
is_package_installed <- function(pkg) {
  return(pkg %in% installed.packages()[, "Package"])
}

# Function to install CRAN packages if not already installed
install_cran_package <- function(pkg) {
  if (!is_package_installed(pkg)) {
    cat(paste("Installing CRAN package:", pkg, "\n"))
    tryCatch({
      install.packages(pkg, dependencies = TRUE)
      cat(paste("✓ Successfully installed:", pkg, "\n"))
    }, error = function(e) {
      cat(paste("✗ Failed to install", pkg, ":", e$message, "\n"))
    })
  } else {
    cat(paste("✓ Package already installed:", pkg, "\n"))
  }
}

# Function to install Bioconductor packages if not already installed
install_bioc_package <- function(pkg) {
  if (!is_package_installed(pkg)) {
    cat(paste("Installing Bioconductor package:", pkg, "\n"))
    tryCatch({
      BiocManager::install(pkg, dependencies = TRUE, update = FALSE)
      cat(paste("✓ Successfully installed:", pkg, "\n"))
    }, error = function(e) {
      cat(paste("✗ Failed to install", pkg, ":", e$message, "\n"))
    })
  } else {
    cat(paste("✓ Package already installed:", pkg, "\n"))
  }
}

# First, install BiocManager if not available
if (!is_package_installed("BiocManager")) {
  cat("Installing BiocManager...\n")
  install.packages("BiocManager")
} else {
  cat("✓ BiocManager already installed\n")
}

# Load BiocManager
library(BiocManager)

# Set repositories
options(repos = BiocManager::repositories())

cat("\n--- Installing CRAN packages ---\n")

# CRAN packages required by REDAC
cran_packages <- c(
  # Shiny and web interface
  "shiny",
  "shinydashboard", 
  "shinythemes",
  "DT",
  
  # Data manipulation and processing
  "dplyr",
  "tidyr",
  "tibble",
  "stringr",
  "reshape2",
  "janitor",
  
  # Plotting and visualization
  "ggplot2",
  "plotly",
  "gplots",
  "vioplot",
  "ggdendro",
  "ggraph",
  "ggrepel",
  "plot3D",
  "factoextra",
  
  # Statistical analysis
  "mgcv",
  "lattice",
  
  # Network analysis
  "igraph",
  
  # Multivariate analysis
  "FactoMineR",
  
  # Web requests and JSON
  "httr",
  "jsonlite",
  
  # Document generation and reporting
  "markdown",
  "knitr",
  "tinytex",
  
  # File I/O
  "readxl",
  
  # Financial data (for quantmod if needed)
  "quantmod"
)

# Install CRAN packages
for (pkg in cran_packages) {
  install_cran_package(pkg)
}

cat("\n--- Installing Bioconductor packages ---\n")

# Bioconductor packages required by REDAC
bioc_packages <- c(
  # Core Bioconductor infrastructure
  "BiocManager",
  "SummarizedExperiment",
  "AnnotationDbi",
  
  # RNA-seq analysis
  "DESeq2",
  "edgeR",
  "preprocessCore",
  
  # Functional enrichment analysis
  "clusterProfiler",
  "enrichplot",
  "DOSE",
  "ReactomePA",
  "pathview",
  "org.Hs.eg.db",
  
  # Additional annotation packages
  "GO.db",
  "KEGG.db",

  # Other
  "rmdHelpers"
)

# Install Bioconductor packages
for (pkg in bioc_packages) {
  install_bioc_package(pkg)
}

cat("\n--- Verifying installation ---\n")

# Test loading critical packages
critical_packages <- c(
  "shiny", "DT", "ggplot2", "plotly", "dplyr", 
  "BiocManager", "DESeq2", "edgeR", "clusterProfiler", 
  "org.Hs.eg.db", "httr", "jsonlite"
)

failed_packages <- c()

for (pkg in critical_packages) {
  tryCatch({
    library(pkg, character.only = TRUE)
    cat(paste("✓", pkg, "loaded successfully\n"))
  }, error = function(e) {
    cat(paste("✗ Failed to load", pkg, ":", e$message, "\n"))
    failed_packages <- c(failed_packages, pkg)
  })
}

cat("\n=============================================================\n")

if (length(failed_packages) == 0) {
  cat("🎉 SUCCESS: All packages installed and loaded successfully!\n")
  cat("REDAC is ready to use.\n")
} else {
  cat("⚠️  WARNING: Some packages failed to install or load:\n")
  for (pkg in failed_packages) {
    cat(paste("   -", pkg, "\n"))
  }
  cat("\nPlease try installing these packages manually:\n")
  cat("For CRAN packages: install.packages('package_name')\n")
  cat("For Bioconductor packages: BiocManager::install('package_name')\n")
}

cat("=============================================================\n")

# Create .Renviron file template if it doesn't exist
if (!file.exists(".Renviron")) {
  cat("\nCreating .Renviron template file...\n")
  writeLines(c(
    "# REDAC Environment Variables",
    "# Add your API key for the AI models below:",
    "# API_KEY=your_api_key_here"
  ), ".Renviron")
  cat("✓ .Renviron template created. Please add your API_KEY.\n")
} else {
  cat("✓ .Renviron file already exists.\n")
}

cat("\n📋 SETUP COMPLETE!\n")
cat("To run REDAC, execute: shiny::runApp()\n")
cat("Don't forget to set your API_KEY in the .Renviron file.\n\n")
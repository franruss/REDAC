# REDAC: RNA-seq Expression Data Analysis Chatbot

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R](https://img.shields.io/badge/R-4.0%2B-blue.svg)](https://www.r-project.org/)
[![Shiny](https://img.shields.io/badge/Shiny-Web%20App-brightgreen.svg)](https://shiny.rstudio.com/)

REDAC is a web-based R Shiny application that offers an interactive, queryable platform designed to simplify RNA-seq expression data exploration and analysis. It provides a straightforward approach to perform differential RNA-seq analysis rapidly, easily, and transparently through natural language queries powered by Large Language Models (LLMs).

## Features

- **Interactive Differential Expression Analysis** using edgeR with natural language queries
- **Dual LLM Integration** with Gemma and LLaMA for code generation and biological interpretation
- **Comprehensive Visualizations** including PCA plots, volcano plots, MA plots, heatmaps, and network visualizations
- **KEGG Pathway Enrichment Analysis** with automated biological interpretation
- **Automated Report Generation** in both HTML (interactive) and Word (static) formats
- **Expert and Novice-Friendly** interface with alternative R code generation for developers

## Online Access

🌐 **Try REDAC Online**: [https://frusso.shinyapps.io/REDAC](https://frusso.shinyapps.io/REDAC)

## Installation

### Prerequisites

- R (≥ 4.0.0)
- RStudio (recommended)
- Internet connection for package installation

### Local Installation

1. **Clone the repository:**
```bash
git clone https://github.com/franruss/REDAC.git
cd REDAC
```

2. **Install required R packages:**
```r
# Install required CRAN packages
install.packages(c("shiny", "shinydashboard", "shinythemes", "ggplot2", 
                   "plotly", "gplots", "lattice", "mgcv", "vioplot", 
                   "httr", "jsonlite", "stringr", "knitr", "tidyr", 
                   "dplyr", "igraph", "FactoMineR", "factoextra", 
                   "plot3D", "ggraph", "ggrepel", "markdown"))

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("SummarizedExperiment", "edgeR", "preprocessCore", 
                       "AnnotationDbi", "clusterProfiler", "org.Hs.eg.db", 
                       "enrichplot", "DOSE"))
```

3. **Run the application:**
```r
# Load required libraries
source("setup.R")

# Launch the app
shiny::runApp()
```

## Usage

REDAC consists of three main modules:

### 1. Perform a Complete Analysis
- Upload your RNA-seq count matrix (TSV format: genes as rows, samples as columns)
- Specify your comparison using natural language (e.g., "execute an rnaseq analysis between treated 3,4 and wt 1,2 samples, de regulated")
- Get automated quality assessment plots (PCA, dendrograms, density plots)
- Receive differential expression results with interactive volcano and MA plots
- Access AI-powered analysis explanations and alternative R code

### 2. Enrichment Analysis and Result Interpretation
- Upload edgeR differential expression results
- Perform KEGG pathway enrichment analysis
- Get dual AI interpretations from Gemma and LLaMA models
- Visualize results with dot plots and network diagrams

### 3. Plot Generation
- Generate individual plots from count matrices or DE results
- Available plots: boxplot, violin plot, heatmap, PCA, network visualization, and more
- Interactive exploration with zoom and gene identification features

## Input Data Format

**Count Matrix** (for complete analysis):
- Tab-separated values (TSV) file
- Genes as rows, samples as columns
- Raw count values (non-negative integers)
- Gene identifiers as row names

**DE Results** (for enrichment analysis):
- Output from edgeR or similar tools
- Required columns: logFC/log2FoldChange, P.Value/pvalue, adj.P.Val/padj

## Example Use Case

Using publicly available Gefitinib resistance data:
1. Upload RNA-seq count data from HCC827 cell lines
2. Compare "HCC827-GR-Pulse vs HCC827-parental, UP-genes"
3. Identify 1,291 up-regulated genes including resistance markers (ALDH1A1, ABCB1, CD44)
4. Discover enriched pathways: TNF signaling, MAPK signaling, NF-kappa B signaling
5. Get expert-level biological interpretation from dual AI models

**Case Study Data**: [Zenodo Record](https://zenodo.org/records/11057181)

## Authors

- **Francesco Russo** - Consiglio Nazionale delle Ricerche (CNR), Naples, Italy
- **Giovanni Maria De Filippis** - University of Naples Federico II, Naples, Italy

[Author Profile](https://scholar.google.it/citations?user=CrlFjHsAAAAJ&hl=it)

<!-- ## Citation

If you use REDAC in your research, please cite:

> De Filippis GM, Sahu P, Ambrosio P, et al. REDAC: RNA-seq Expression Data Analysis Chatbot. *Bioinformatics Advances* (2025). DOI: [To be published] -->

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Funding

This work has been funded by the **MUR (PRIN2022PNRR-2022MZJR9X, P2022YA3LL)**.

## Contact

📧 **francesco.russo@cnr.it**

## Acknowledgments

- Ministero dell'Università e della Ricerca (MUR), Italy
- Department of Electrical Engineering and Information Technology, University of Naples Federico II
- IRCCS Istituto Romagnolo per lo studio dei Tumori "Dino Amadori" (IRST)
- Cell Biology and Biotherapy Unit, Istituto Nazionale Tumori-IRCCS-Fondazione G. Pascale

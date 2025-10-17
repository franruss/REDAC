
# Install and load the edgeR package
#install.packages("edgeR")
library(edgeR)

# Configurazioni dataset-specifiche
dataset_configs <- list(
    file = "./data/Raw-data-RNAseq-for-submission-with-gene-name-3-1-24.csv",
    file = "./data/GSE164073_Eye_count_matrix.tsv.txt",
    file = "./data/GSE146458_raw_counts_GRCh38.p13_NCBI.txt",
    file = "./data/chol_tcga_gdc_read_counts.txt",
    file = "./data/coad_tcga_gdc_read_counts.txt",
    file = "./data/acc_tcga_gdc_read_counts.txt",
    file = "./data/dlbclnos_tcga_gdc_read_counts.txt",
    file = "./data/ucs_tcga_gdc_read_counts.txt",
    file = "./data/plmeso_tcga_gdc_seq_read_counts.txt",
    file = "./data/GSE279712_human.merged.counts.txt"
)
length(dataset_configs)

# Codice generato da LLAMA e racchiuso in una funzione
function_generated_by_LLM <- function(df=df){
  group <- factor(c(rep("wt", 2), rep("treated", 2)))

  # Create a DGEList object
  y <- DGEList(counts = df[, 1:4], group = group)

  # Filter out lowly expressed genes
  keep <- rowSums(cpm(y) > 1) >= 2
  y <- y[keep, ]
  number_of_genes_after_filtering = dim(y)[1]

  # Normalize the data
  y <- calcNormFactors(y, method = "TMM")

  # Estimate the dispersion
  design <- model.matrix(~0 + group, data = y$samples)
  y <- estimateDisp(y, design)

  # Fit the model
  fit <- glmFit(y, design)

  # Perform the differential expression analysis
  lrt <- glmLRT(fit, coef = "groupwt")

  # Get the results
  results <- as.data.frame(lrt$table)

  # Filter for de-regulated genes
  significant_genes <- results[results$PValue < 0.05, ]
  return(list(significant_genes=significant_genes, number_of_genes_after_filtering=number_of_genes_after_filtering))
}

REDAC_de_function <- function(df=df){
  which_results = "de"
  # Load your count matrix
  counts = df
  
  # Define the sample groups
  group <- factor(c(rep("wt", 2), rep("treated", 2)))
  dge <- DGEList(counts = df[, 1:4], group = group)

  # Filter low counts,
  keep <- rowSums(counts) >= dim(counts)[2]
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  number_of_genes_after_filtering = dim(dge)[1]
  # Normalize and estimate dispersion,
  dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge)

  # Differential expression,
  et <- exactTest(dge, pair = c('wt', 'treated'))
  res <- topTags(et, n = nrow(dge))$table
  res <- res[order(res$FDR), ]
  res$log2FoldChange <- round(res$logFC, 4)
  res$logFC <- NULL
  DEres <- na.omit(DEres)
  significant_genes = DEres[DEres$FDR < 0.05, ]
  return(list(significant_genes=significant_genes, number_of_genes_after_filtering=number_of_genes_after_filtering))
}

read_and_clean_colnames <- function(file_path, sep = sep, header = header, rownames = rownames ){
  # Read the first line (header)
  header_line <- readLines(file_path, n = 1)
  
  # Split column names
  raw_colnames <- strsplit(header_line, split = sep, fixed = TRUE)[[1]]
  
  # Clean column names: replace spaces and invalid characters with underscores
  cleaned_colnames <- make.names(raw_colnames)
  cleaned_colnames <- gsub("\\.+", "_", cleaned_colnames)  # collapse multiple dots to underscores
  cleaned_colnames <- gsub("_$", "", cleaned_colnames)     # remove trailing underscores if any
  cleaned_colnames <- cleaned_colnames[-1]    # remove trailing underscores if any
  
  
  # Read the full data (skip the header already processed)
  df <- read.table(file_path, sep = sep, header = header, row.names = rownames, skip = 1)
  
  # Assign cleaned column names
  colnames(df) <- cleaned_colnames
  
  return(df)
}

# Primo Ciclo di Testing della funzione generata da LLAMA: function_generated_by_LLM

for (i in 1:length(dataset_configs)){
  df <- read_and_clean_colnames(dataset_configs[i]$file, sep = "\t", header = TRUE, rownames = 1)
  print("executing")
  print(dataset_configs[i]$file)
  print("total genes found in the input file")
  print(dim(df)[1])
  res = function_generated_by_LLM(df) 
  print("number_of_genes_after_filtering")
  print(res$number_of_genes_after_filtering)
  print("significant genes found")
  print(dim(res$significant_genes)[1])
  if(res$number_of_genes_after_filtering==(dim(res$significant_genes)[1])){
    print("FAILED: Something wrong in the execution! Number of genes after filtering and number of significant genes are the same!")
  }
  if(res$number_of_genes_after_filtering>(dim(res$significant_genes)[1])){
    print("SUCCESS: number of genes after filtering greather than number of significant genes!")
  }
  print("_________________________________")
}

# Primo Ciclo di Testing della funzione scritta da me in REDAC: REDAC_de_function

for (i in 1:length(dataset_configs)){
  df <- read_and_clean_colnames(dataset_configs[i]$file, sep = "\t", header = TRUE, rownames = 1)
  print("executing")
  print(dataset_configs[i]$file)
  print("total genes found in the input file")
  print(dim(df)[1])
  res = REDAC_de_function(df) 
  print("number_of_genes_after_filtering")
  print(res$number_of_genes_after_filtering)
  print("significant genes found")
  print(dim(res$significant_genes)[1])
  if(res$number_of_genes_after_filtering==(dim(res$significant_genes)[1])){
    print("FAILED: Something wrong in the execution! Number of genes after filtering and number of significant genes are the same!")
  }
  if(res$number_of_genes_after_filtering>(dim(res$significant_genes)[1])){
    print("SUCCESS: number of genes after filtering greather than number of significant genes!")
  }
  print("_________________________________")
}




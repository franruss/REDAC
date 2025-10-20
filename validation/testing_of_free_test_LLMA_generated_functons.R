
# Install and load the edgeR package
#install.packages("edgeR")
library(edgeR)
library(DESeq2)

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
function_generated_by_LLM_1 <- function(df=df){
  group <- factor(c(rep("wt", 2), rep("treated", 2)))

  # Create a DGEList object
  y <- DGEList(counts = df[, 1:4], group = group)

  # Filter out lowly expressed genes
  keep <- rowSums(cpm(y) > 1) >= 2
  y <- y[keep, ]
  #number_of_genes_after_filtering = dim(y)[1]

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
  return(list(significant_genes=significant_genes, number_of_genes_after_filtering=dim(y)[1]))
}

### Here’s an example of how you can perform an RNA-seq analysis using the edgeR package in R. In this example, we’ll 
# compare the treated samples (HCC827_GR_High_A and HCC827_GR_High_B) to the wild-type samples (HCC827_Parental_A and HCC827_Parental_B).
# Codice generato da LLAMA e racchiuso in una funzione
function_generated_by_LLM_2 <- function(df=df){
# Create a group vector
  # Create a DGEList object
  y <- DGEList(counts = df, group = c(rep("wt", 2), rep("treated", 2), rep("low", 2), rep("pulse", 2), rep("wt2", 2), rep("treated2", 2), rep("low2", 2), rep("pulse2", 2)))
  
  # Filter out lowly expressed genes
  keep <- rowSums(cpm(y) > 1) >= 2
  y <- y[keep, ]
  
  # Normalize the data
  y <- calcNormFactors(y, method = "TMM")
  
  # Design matrix
  design <- model.matrix(~0 + group, data = y$samples)
  
  # Fit the model
  y <- estimateDisp(y, design)
  
  # Perform the differential expression analysis
  fit <- glmFit(y, design)
  
  # Contrast of interest
  contrast <- makeContrasts(groupwt - grouptreated, levels = design)
  
  # Perform the differential expression analysis
  lrt <- glmLRT(fit, contrast = contrast)
  
  # Get the results
  results <- lrt$table
  significant_genes <- results[results$PValue < 0.05, ]
  return(list(significant_genes=significant_genes, number_of_genes_after_filtering= dim(y)[1]))
}

##prompt: Example: perform an rnaseq analysis between treated 3,4 and wt 1,2 samples, de regulated and use exactTest of edgeR package

# Codice generato da LLAMA e racchiuso in una funzione
function_generated_by_LLM_3 <- function(df=df){
group <- factor(c(rep("wt", 2), rep("treated", 2)))

# Create a DGEList object
y <- DGEList(counts = df[, 1:4], group = group)

# Filter out lowly expressed genes
keep <- rowSums(cpm(y) > 1) >= 2
y <- y[keep, ]

# Normalize the data
y <- calcNormFactors(y, method = "TMM")

# Perform the exact test
design <- model.matrix(~0 + group, data = y$samples)
colnames(design) <- levels(y$samples$group)
y <- estimateDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = "treated")

# Get the results
results <- lrt$table

# Sort the results by p-value
results <- results[order(results$PValue), ]
significant_genes <- results[results$PValue < 0.05, ]
return(list(significant_genes=significant_genes, number_of_genes_after_filtering= dim(y)[1]))
}

# Codice generato da LLAMA e racchiuso in una funzione
function_generated_by_LLM_4 <- function(df=df){
# Define the sample groups
group <- factor(c(rep("wt", 2), rep("treated", 2)))

# Create a DGEList object
y <- DGEList(counts = df[, 1:4], group = group)

# Filter out lowly expressed genes
keep <- rowSums(cpm(y) > 1) >= 2
y <- y[keep, ]

# Normalize the data
y <- calcNormFactors(y, method = "TMM")

# Perform the exact test
design <- model.matrix(~0 + group, data = y$samples)
colnames(design) <- levels(y$samples$group)
y <- estimateDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = "treated")

# Get the results
results <- lrt$table

# Sort the results by p-value
results <- results[order(results$PValue), ]
significant_genes <- results[results$PValue < 0.05, ]
return(list(significant_genes=significant_genes, number_of_genes_after_filtering= dim(y)[1]))
}

# Codice generato da LLAMA e racchiuso in una funzione
function_generated_by_LLM_5 <- function(df=df){
# Define the sample groups
group <- factor(c(rep("wt", 2), rep("treated", 2)))

# Create a DGEList object
y <- DGEList(counts = df[, 1:4], group = group)

# Filter out lowly expressed genes
keep <- rowSums(cpm(y) > 1) >= 2
y <- y[keep, ]

# Normalize the data
y <- calcNormFactors(y, method = "TMM")

# Perform the exact test
design <- model.matrix(~0 + group, data = y$samples)
colnames(design) <- levels(y$samples$group)
y <- estimateDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = "treated")

# Get the results
results <- lrt$table

# Sort the results by p-value
results <- results[order(results$PValue), ]
significant_genes <- results[results$PValue < 0.05, ]
return(list(significant_genes=significant_genes, number_of_genes_after_filtering=number_of_genes_after_filtering))
}

# Codice generato da LLAMA e racchiuso in una funzione
function_generated_by_LLM_6 <- function(df=df){
# Define the sample groups
group <- factor(c("wt", "wt", "treated", "treated"))

# Create a DGEList object
y <- DGEList(counts = df[, 1:4], group = group)

# Filter out lowly expressed genes
keep <- rowSums(cpm(y) > 1) >= 2
y <- y[keep, ]

# Normalize the data
y <- calcNormFactors(y, method = "TMM")

# Perform the exact test
design <- model.matrix(~0 + group, data = y$samples)
colnames(design) <- levels(y$samples$group)
y <- estimateDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = "treated")

# Get the results
results <- lrt$table

# Sort the results by p-value
results <- results[order(results$PValue), ]
significant_genes <- results[results$PValue < 0.05, ]
return(list(significant_genes=significant_genes, number_of_genes_after_filtering=number_of_genes_after_filtering))
}

## prompt: Example: perform an rnaseq analysis by using exactTest of edgeR package between treated 3,4 and wt 1,2 samples, de regulated

# Codice generato da LLAMA e racchiuso in una funzione
function_generated_by_LLM_7 <- function(df=df){
# We'll assume that the first two columns are the wt samples and the next two columns are the treated samples
group <- factor(c(rep("wt", 2), rep("treated", 2), rep("wt", 2), rep("treated", 2), rep("wt", 2), rep("treated", 2), rep("wt", 2), rep("treated", 2)))
y <- DGEList(counts = df[, -1], group = group)

# Normalize the data
y <- calcNormFactors(y, method = "TMM")

# Estimate the dispersion
design <- model.matrix(~0 + group, data = y$samples)
y <- estimateDisp(y, design)

# Perform the exact test
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = "groupwt")

# Get the results
results <- as.data.frame(lrt$table)

# Filter the results to get the differentially expressed genes
de_results <- results[results$PValue < 0.05, ]

group <- factor(c(rep("wt", 2), rep("treated", 2), rep("wt", 2), rep("treated", 2), rep("wt", 2), rep("treated", 2), rep("wt", 2), rep("treated", 2)))
y <- DGEList(counts = df[, -1], group = group)

# Normalize the data
y <- calcNormFactors(y, method = "TMM")

# Estimate the dispersion
design <- model.matrix(~0 + group, data = y$samples)
y <- estimateDisp(y, design)

# Perform the exact test
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = "groupwt")

# Get the results
results <- as.data.frame(lrt$table)

# Filter the results to get the differentially expressed genes
de_results <- results[results$PValue < 0.05, ]

return(list(significant_genes=de_results, number_of_genes_after_filtering=number_of_genes_after_filtering))
}

# Codice generato da LLAMA e racchiuso in una funzione
function_generated_by_LLM_8 <- function(df=df){
# Define the group variable
group <- factor(c(rep("wt", 2), rep("treated", 2), rep("wt", 2), rep("treated", 2), rep("wt", 2), rep("treated", 2), rep("wt", 2), rep("treated", 2)))
group <- factor(c("wt", "wt", "treated", "treated", "wt", "wt", "treated", "treated", "wt", "wt", "treated", "treated", "wt", "wt", "treated", "treated"))

# Create a DGEList object
y <- DGEList(counts = df[, -1], group = group)

# Normalize the data
y <- calcNormFactors(y, method = "TMM")

# Estimate the dispersion
design <- model.matrix(~0 + group, data = y$samples)
y <- estimateDisp(y, design)

# Perform the exact test
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = "groupwt")

# Get the results
results <- as.data.frame(lrt$table)

# Filter the results to get the differentially expressed genes
de_results <- results[results$PValue < 0.05, ]

return(list(significant_genes=de_results, number_of_genes_after_filtering=number_of_genes_after_filtering))
}

# Codice generato da LLAMA e racchiuso in una funzione
function_generated_by_LLM_9 <- function(df=df){
# Define the group variable
group <- factor(c("wt", "wt", "treated", "treated", "wt", "wt", "treated", "treated", "wt", "wt", "treated", "treated", "wt", "wt", "treated", "treated"))
names(group) <- colnames(df[, -1])

# Create a DGEList object
y <- DGEList(counts = df[, -1], group = group)

# Normalize the data
y <- calcNormFactors(y, method = "TMM")

# Estimate the dispersion
design <- model.matrix(~0 + group, data = y$samples)
y <- estimateDisp(y, design)

# Perform the exact test
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = "groupwt")

# Get the results
results <- as.data.frame(lrt$table)

# Filter the results to get the differentially expressed genes
de_results <- results[results$PValue < 0.05, ]

return(list(significant_genes=de_results, number_of_genes_after_filtering=number_of_genes_after_filtering))
}

### prompt: Here’s an example of how you can perform an RNA-seq analysis using the edgeR package in R. Since you also mentioned using DESeq2, I’ll provide an example for that as well.

function_generated_by_LLM_10 <- function(df=df){
# Create a DGEList object for edgeR
group <- factor(c("wt", "wt", "treated", "treated"))
df=df[,1:4]
# Now, let's use DESeq2
# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = as.matrix(df), colData = data.frame(condition = c("wt", "wt", "treated", "treated")), design = ~ condition)

# Perform the differential expression analysis
dds <- DESeq(dds)

# Get the results
res <- results(dds)
# Filter for de-regulated genes
res <- na.omit(res)
significant_genes <- res[res$padj<0.05, ]
return(list(significant_genes=significant_genes, number_of_genes_after_filtering=dim(df)[1]))
}

list_of_functions_generated_by_LLM = list(function_generated_by_LLM_1,function_generated_by_LLM_2,function_generated_by_LLM_3,function_generated_by_LLM_4,function_generated_by_LLM_5,
                                       function_generated_by_LLM_6,function_generated_by_LLM_7,function_generated_by_LLM_8,function_generated_by_LLM_9,function_generated_by_LLM_10)

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
  DEres <- res
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

no_execution_errors = 0
non_reliable_results = 0
success = 0
execution_errors = 0
for (i in 1:length(dataset_configs)){
   for (j in 1:length(list_of_functions_generated_by_LLM)){
     print(paste("test on data: ",i," on function: ",j),sep="")
     error_message=NULL
     tryCatch({
      df <- read_and_clean_colnames(dataset_configs[i]$file, sep = "\t", header = TRUE, rownames = 1)
      print("executing")
      print(dataset_configs[i]$file)
      print("total genes found in the input file")
      print(dim(df)[1])
      res = list_of_functions_generated_by_LLM[[j]](df) 
      print("number_of_genes_after_filtering")
      print(res$number_of_genes_after_filtering)
      print("significant genes found")
      print(dim(res$significant_genes)[1])
     }, error = function(e) {
       error_message <<- e$message
     })
     if(is.null(error_message)){
         print("Execution Completed!")
         no_execution_errors = no_execution_errors + 1
         if(res$number_of_genes_after_filtering==(dim(res$significant_genes)[1])){
             print("FAILED: Something wrong in the results! Number of genes after filtering and number of significant genes are the same!")
             non_reliable_results = non_reliable_results + 1
         }
         if(res$number_of_genes_after_filtering>(dim(res$significant_genes)[1])){
             print("SUCCESS: number of genes after filtering greather than number of significant genes!")
             success = success + 1
         }
     }else{
        print("An execution error found!")
        execution_errors = execution_errors + 1
     }
     print("_________________________________")
   
   }
}
no_execution_errors 
non_reliable_results 
success 
execution_errors 






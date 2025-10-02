# Function definitions
distToSim  <- function(x) {1-2*x^2}
simToDist <- function(x) {sqrt(1/2*(1-(x)))}
simFun <- function(x, y) {corCosine(x, y)} # Cosine similarity between two vectors or conformable matrices.
# For centered data, cosine similarity equals Pearson correlation.
# Definition: cos(theta) = (X · Y) / ( ||X|| * ||Y|| )

ntp_mine <- function (emat, templates, nPerm = 1000, distance = "cosine", 
                      nCores = 1, seed = NULL, verbose = getOption("verbose"), 
                      doPlot = FALSE) 
{
  if (class(emat)[1] == "ExpressionSet") 
    emat <- suppressPackageStartupMessages(Biobase::exprs(emat))
  if (is.data.frame(emat) | is.vector(emat)) 
    emat <- as.matrix(emat)
  if (is.null(rownames(emat))) 
    stop("missing emat rownames, check input")
  if (is.null(templates$class) | is.null(templates$probe)) {
    stop("missing columns in templates, check input")
  }
  if (is.character(templates$class)) 
    templates$class <- factor(templates$class)
  if (is.factor(templates$probe) | is.integer(templates$probe)) {
    warning("templates$probe coerced to character", call. = FALSE)
    templates$probe <- as.character(templates$probe)
  }
  keepP <- stats::complete.cases(emat)
  if (sum(!keepP) > 0) {emat <- emat[keepP, , drop = FALSE]}
  keepT <- templates$probe %in% rownames(emat)
  if (sum(!keepT) > 0) {
    if (isTRUE(verbose)){
      message(paste0(sum(!keepT), "/", length(keepT), 
                     " templates features not in emat, discarded"))
    }
    templates <- templates[keepT, ]
  }
  if (min(table(templates$class)) < 2) {
    message("<2 matched features/class")
    stop("check templates$probe is matchable against rownames(emat)", 
         call. = FALSE)
  }
  if (min(table(templates$class)) < 5) {
    warning("<5 matched features/class - unstable predictions", 
            call. = FALSE)
  }
  N <- ncol(emat)
  K <- nlevels(templates$class)
  S <- nrow(templates)
  P <- nrow(emat)
  class.names <- levels(templates$class)
  templates$class <- as.numeric(templates$class) # Convert class factor to numeric (needed to build tmat design matrix)
  emat.mean <- round(mean(emat), 2)
  if (abs(emat.mean) > 1) {
    isnorm <- " <- check feature centering!"
    emat.sd <- round(stats::sd(emat), 2)
    warning(paste0("emat mean=", emat.mean, "; sd=", 
                   emat.sd, isnorm), call. = FALSE)
  }
  feat.class <- paste(range(table(templates$class)), collapse = "-")
  mm <- match(templates$probe, rownames(emat), nomatch = 0) # Indices of predictor (template) probes in expression matrix
  if (!all(rownames(emat)[mm] == templates$probe)) {stop("error matching probes, check rownames(emat) and templates$probe")}
  pReplace <- length(templates$probe) > length(unique(templates$probe))
  tmat <- matrix(rep(templates$class, K), ncol = K)
  for (k in (1:K)){tmat[, k] <- as.numeric(tmat[, k] == k)}
  tmat # design matrix
  if (K == 2){tmat[tmat == 0] <- -1}
  res=NULL
  #ntpFUN <- function(n) { # n corre su tutti i campioni
  for (n in (1:N)){
  n.sim <- as.vector(simFun( emat[mm, n, drop = FALSE] ,tmat)) # Correlate sample with design matrix columns (class templates)
  # Rationale: identifies which class template (up genes close to 1, down genes close to 0 or -1) best matches the sample.
  # plot(emat[mm, n, drop = FALSE] , tmat[,1])  # (Optional exploratory plot)
  n.sim.perm.max <- apply(simFun(matrix(emat[, n][sample.int(P, S * nPerm, replace = TRUE)], ncol = nPerm), tmat), 1, max) 
  # Permutation: build random expression matrices (sampling with replacement) and compute class similarities; take per-class maxima.
  n.ntp <- which.max(n.sim) # Selected class with the highest observed similarity score
    n.sim.ranks <- rank(-c(n.sim[n.ntp], (n.sim.perm.max)))
    n.pval <- n.sim.ranks[1]/length(n.sim.ranks)
    res[[n]] <- (c(n.ntp, simToDist(n.sim), n.pval))
  }
  #res <- lapply(seq_len(N), ntpFUN)
  res <- data.frame(do.call(rbind, res))
  colnames(res) <- c("prediction", paste0("d.", class.names), "p.value")
  res$prediction <- factor(class.names[res$prediction], levels = class.names)
  rownames(res) <- colnames(emat)
  res$p.value[res$p.value < 1/nPerm] <- 1/nPerm
  res$FDR <- stats::p.adjust(res$p.value, "fdr")
  return(res)
}

#  treated_samples = c(1,2)
#  wt_samples = c(3,4)
#  counts  <- read.table("C:/Users/f.russo.ENGIBBC/Desktop/REDAC_submission/expression_file.txt", header = TRUE, sep = "\t", row.names = 1)
#  treated_samples = colnames(counts)[as.numeric(treated_samples)]
#  wt_samples =      colnames(counts)[as.numeric(wt_samples)]
#  which_results = "de"
setClassUnion("ExpData", c("matrix", "SummarizedExperiment"))

de_analysis_old <- function (counts,treated_samples,wt_samples,which_results) {
  DEres=NULL
  condition = rep("treated_samples",length(treated_samples))
  treated_samples_condition = cbind( treated_samples, condition)
  condition = rep("wt_samples",length(wt_samples))
  wt_samples_condition = cbind(wt_samples, condition)
  # Build the samples metadata table containing comparison info
  samples = rbind(wt_samples_condition, treated_samples_condition)
  colnames(samples) = c("patients","condition")
  rownames(samples) = samples[,1]
  dataset = counts 
  dim(dataset)
  dataset =   dataset[, which(colnames(dataset) %in% samples[,1]) ] #selezioni solo colonne di interesse dal file dataset
  head(dataset)
  ### CONTROLLA che deve essere TRUE
  if(dim(dataset)[2] == dim(samples)[1]){
    samples = as.data.frame(samples)
  dataset3 = dataset[,samples[,1]] # Reorder dataset columns to match sample ordering
    dataset3[1:3,]
    if(all(colnames(dataset3)==rownames(samples))){ # Sanity check: column order must match metadata rownames
      ###################################
      #### Ensure ALL values above are TRUE
      ##################################
      name = paste("comparison_",unique(samples[,2])[1],"_vs_",unique(samples[,2])[2],sep="")
      print(name)
      samples$condition <- factor(samples$condition, levels = unique(samples$condition))

  # IMPORTANT: ensure matching dimensions between metadata and count matrix
      nrow(samples) == ncol(dataset3)

      dataset3 = as.matrix(dataset3)
      dim(dataset3)
      colnames(dataset3)
      print(samples, row.names=FALSE)
      dds <- DESeqDataSetFromMatrix(countData=dataset3, colData=samples, design=~condition, tidy = FALSE)
      # Tidy = TRUE. For matrix input: whether the first column of countData is the rownames for the count matrix

  # Remove genes with total counts < number of samples (low-information features)
      keep <- rowSums(counts(dds)) >= dim(counts(dds))[2]
      dds <- dds[keep,]

  # Run DESeq2 (may take some time)
      dds <- DESeq(dds)
      #counts are divided by those:
      sizeFactors(dds)
      normalized_counts = counts(dds, normalized=TRUE)
      dim(normalized_counts)

  # Potentially time-consuming step (minutes on larger datasets)
      res = results(dds, tidy=TRUE)
      head(res)

      res <- res[order(res$padj),]
      head(res)
      res$baseMean = round(res$baseMean,1)
      res$log2FoldChange = round(res$log2FoldChange,6)

      #write.table(res[,c(1,3,7)],file = paste("global_results_on_filtered_",name,".txt",sep=""), quote = FALSE,row.names = FALSE,col.names = TRUE,sep = "\t")
      DEres = res
      if(which_results=="down"){DEres=DEres[DEres$log2FoldChange<0,]}else if(which_results=="up"){DEres=DEres[DEres$log2FoldChange>0,]}
  DEres = na.omit(DEres) # Remove rows with NA
      #print(DEres)
    }else{print("check the columns!")}
  }else{print("check the columns in the question you wrote!")}
  DEres
}

de_analysis_old_2 <- function (counts,treated_samples,wt_samples,which_results) {
  DEres <- NULL
  
  # Build sample metadata
  condition <- rep("treated_samples", length(treated_samples))
  treated_samples_condition <- cbind(treated_samples, condition)
  
  condition <- rep("wt_samples", length(wt_samples))
  wt_samples_condition <- cbind(wt_samples, condition)
  
  samples <- rbind(wt_samples_condition, treated_samples_condition)
  colnames(samples) <- c("patients", "condition")
  rownames(samples) <- samples[, 1]
  
  # Subset and order count matrix
  dataset <- counts[, which(colnames(counts) %in% samples[, 1])]
  
  if (ncol(dataset) == nrow(samples)) {
    samples <- as.data.frame(samples)
    dataset <- dataset[, samples[, 1]]  # reorder columns to match sample metadata
    
    if (all(colnames(dataset) == rownames(samples))) {
      name <- paste("comparison_", unique(samples[, 2])[1], "_vs_", unique(samples[, 2])[2], sep = "")
      print(name)
      
      samples$condition <- factor(samples$condition, levels = unique(samples$condition))
      
      # Create DGEList
      group <- samples$condition
      dge <- DGEList(counts = dataset, group = group)
      
      # Filter low counts
      dim(dge$counts)
      keep <- rowSums(counts) >= dim(counts)[2]
      dge <- dge[keep, , keep.lib.sizes=FALSE]
      dim(dge$counts)
      
      # Normalize and estimate dispersion
      dge <- calcNormFactors(dge)
      dge <- estimateDisp(dge)
      
      # Perform differential expression test
      et <- exactTest(dge, pair = rev(levels(group)))  # reversed for treated vs wt
      res <- topTags(et, n = nrow(dge))$table
      
      # Clean up results
      res <- res[order(res$FDR), ]
      res$log2FoldChange <- round(res$logFC, 6)
      res$logFC <- NULL
      
      # Filter based on user input
      if (which_results == "down") {
        DEres <- res[res$log2FoldChange < 0, ]
      } else if (which_results == "up") {
        DEres <- res[res$log2FoldChange > 0, ]
      } else {
        DEres <- res
      }
      
      DEres <- na.omit(DEres)
    } else {
      print("Check the column order of count data vs sample metadata.")
    }
  } else {
    print("Mismatch between number of samples in counts and metadata.")
  }
  
  return(DEres)
}


de_analysis <- function (counts,treated_samples,wt_samples,which_results) {
  DEres <- NULL
# Build metadata
  treated_samples = as.numeric(gsub("\\D", "", treated_samples))
  wt_samples      = as.numeric(gsub("\\D", "", wt_samples))
treated_condition <- data.frame(patients = treated_samples, condition = 'treated')
print(treated_condition)
wt_condition <- data.frame(patients = wt_samples, condition = 'wt')
samples <- rbind(wt_condition, treated_condition)
print("samples")
print(samples)
rownames(samples) <- samples$patients

# Subset and order count matrix
dataset <- counts[, samples$patients]
samples$condition <- factor(samples$condition, levels = c('wt', 'treated'))
print("dataset")
print(head(dataset))
# Create DGEList object
dge <- DGEList(counts = dataset, group = samples$condition)

# Filter low-count genes
keep <- rowSums(counts) >= dim(counts)[2]
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Normalize and estimate dispersion
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge)

# Differential expression testing
et <- exactTest(dge, pair = c('wt', 'treated'))
res <- topTags(et, n = nrow(dge))$table
res <- res[order(res$FDR), ]
res$log2FoldChange <- round(res$logFC, 4)
res$logFC <- NULL

# Filter DE results according to which_results parameter
if (which_results == 'down') {
  DEres <- res[res$log2FoldChange < 0, ]
} else if (which_results == 'up') {
  DEres <- res[res$log2FoldChange > 0, ]
} else {
  DEres <- res
}
DEres <- na.omit(DEres)
return(DEres)
}

read_and_clean_colnames2 <- function(file_path, sep = sep, header = header, rownames = rownames ) {
  # Read the Excel sheet
  df <- readxl::read_excel(file_path, sheet = 1)
  # Clean column names (e.g., replaces spaces with underscores)
  df <- janitor::clean_names(df)
  df = as.data.frame(df)
  gene_names= df[,1]
  df=as.matrix(as.matrix.data.frame(df[,-1]))
  row.names(df) = gene_names
  return(df)
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

calculate_correlation <- function (counts,treated_samples,wt_samples,which_results) {
# Compute Pearson correlation for all genes versus patient overall survival (OS) and return indices above a threshold.
patientlist <- read.table(file="coadread_data_bcr_clinical_data_patient.txt",header=TRUE,sep="\t") # Clinical metadata
coadread_expression <- read.table(file="coadread_data_RNA_Seq_v2_expression_median.txt",header=TRUE,sep="\t") # Expression data
coad_expr = (coadread_expression[,c(-1,-2)])
# Optional: quantile normalization step was considered but is commented out.
myfqua = coad_expr
colnames(myfqua) = colnames(coadread_expression)[3:length((coadread_expression))]
rownames(myfqua) = coadread_expression[,1]
patients = patientlist$PATIENT_ID
patients = strsplit(as.character(patients),"-")
patientnames=NULL
# Build patient identifiers using dots instead of hyphens (TCGA.xx.yy.01)
for (p in 1:length(patients)){
  patientnames[p] = paste("TCGA",patients[[p]][2],patients[[p]][3],"01",sep=".")
}
patientlist2 = cbind(patientnames,patientlist)
patientlist3=NULL
list_absent_patients=NULL
k=0
for (i in 1:length(colnames(myfqua))){
  temp = subset(patientlist2,patientlist2$patientnames==colnames(myfqua)[i]) # Check if patient exists in clinical file
  if (dim(temp)[1]>0 && as.character(temp$OS_MONTHS)!="[Not Available]") {
    patientlist3 = rbind(patientlist3,temp)
  }else{
    k=k+1
    print("patient not found!")
    list_absent_patients[k] = colnames(myfqua)[i]
  }
}
list_absent_patients
# Remove patients lacking clinical data from expression matrix
patients_to_be_eliminated = which(colnames(myfqua) %in% list_absent_patients) 
gene_expression2 = myfqua[,-patients_to_be_eliminated]
dim(patientlist3)
dim(gene_expression2)
cor_res = vector(mode="numeric",length=length(rownames(gene_expression2)))
cor_test_res = vector(mode="numeric",length=length(rownames(gene_expression2)))
method_used="pearson"
genes = rownames(myfqua)
for(j in 1:length(rownames(gene_expression2))){
  print(j)
  gene_expression_line = subset(gene_expression2,as.character(rownames(gene_expression2))==as.character(genes[j]) )
  cor_res[j]      = round(cor(as.numeric(gene_expression_line),as.numeric(patientlist3$OS_MONTHS), method=method_used),3)
  cor_test_res[j] = round(cor.test(as.numeric(gene_expression_line),as.numeric(patientlist3$OS_MONTHS),method = method_used)[[3]],7)
}
which(cor_res>0.2)
}

##### END definitions
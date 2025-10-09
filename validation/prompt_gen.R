# Configurazioni dataset-specifiche
dataset_configs <- list(
  "RNAseq_submission" = list(
    file = "./data/Raw-data-RNAseq-for-submission-with-gene-name-3-1-24.csv",
    treated_cols = "1, 2",
    wt_cols = "3, 4"
  ),

  "GSE164073_Eye" = list(
    file = "./data/GSE164073_Eye_count_matrix.tsv.txt",
    treated_cols = "MW4_cornea_CoV2_1, MW5_cornea_CoV2_2",
    wt_cols = "MW1_cornea_mock_1, MW2_cornea_mock_2"
  ),  #https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164073

  "GSE146458_raw" = list(
    file = "./data/GSE146458_raw_counts_GRCh38.p13_NCBI.txt",
    treated_cols = "GSM4386209, GSM4386210",
    wt_cols = "GSM4386207, GSM4386208"
  ) #https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146458
)

# Funzione per generare prompts basati sulle colonne
generate_prompts <- function(treated_cols, wt_cols) {
  treated_underscore <- gsub(", ", "_", treated_cols)
  wt_underscore <- gsub(", ", "_", wt_cols)

  prompts <- c(
    sprintf("perform rnaseq analysis between treated %s and wt %s samples, up regulated", treated_cols, wt_cols),
    sprintf("compare treated %s vs control %s down-regulated genes", treated_cols, wt_cols),
    sprintf("analyze treated columns %s against wild type columns %s, only genes down regulated", treated_cols, wt_cols),
    sprintf("differential expression analysis: treated samples %s vs wild type samples %s, all differentially expressed genes", treated_cols, wt_cols),
    sprintf("compare columns %s (treated) with columns %s (control), show upregulated", treated_cols, wt_cols),
    sprintf("DE analysis treated_%s vs control_%s up", treated_underscore, wt_underscore),
    sprintf("rnaseq comparison: treated samples %s versus wt samples %s, downregulated genes", treated_cols, wt_cols),
    sprintf("perform differential gene expression analysis between treatment group (samples %s) and control group (samples %s), focus on upregulated genes", treated_cols, wt_cols)
  )

  return(prompts)
}

# Generare prompts per ogni dataset
dataset_prompts <- lapply(dataset_configs, function(config) {
  list(
    file = config$file,
    treated_cols = config$treated_cols,
    wt_cols = config$wt_cols,
    prompts = generate_prompts(config$treated_cols, config$wt_cols)
  )
})

# Accesso ai prompts generati
# dataset_prompts$RNAseq_submission$prompts[1]
# dataset_prompts$GSE164073_Eye$prompts
# dataset_prompts$GSE146458_raw$prompts[3]

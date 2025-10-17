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
  ), #https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146458
  "chol_tcga" = list(
    file = "./data/chol_tcga_gdc_read_counts.txt",
    treated_cols = "TCGA-3X-AAV9-01A",	"TCGA-3X-AAVA-01A",	
    wt_cols = "TCGA-3X-AAVB-01A",	"TCGA-3X-AAVC-01A"
  ), #https://https://www.cbioportal.org/
  "coad_tcga" = list(
    file = "./data/coad_tcga_gdc_read_counts.txt",
    treated_cols = "TCGA-3L-AA1B-01A", "TCGA-4N-A93T-01A",	
    wt_cols = "TCGA-OR-A5J3-01A", "TCGA-OR-A5J5-01A"
  ), #https://https://www.cbioportal.org/
  "acc_tcga" = list(
    file = "./data/acc_tcga_gdc_read_counts.txt",
    treated_cols = "TCGA-OR-A5J1-01A", "TCGA-OR-A5J2-01A",	
    wt_cols = "TCGA-4T-AA8H-01A", "TCGA-5M-AAT4-01A"
  ), #https://https://www.cbioportal.org/
  "dlbclnos_tcga" = list(
    file = "./data/dlbclnos_tcga_gdc_read_counts.txt",
    treated_cols = "TCGA-FA-8693-01A", "TCGA-FA-A4BB-01A",	
    wt_cols = "TCGA-FA-A4XK-01A", "TCGA-FA-A6HN-01A"
  ), #https://https://www.cbioportal.org/
  "ucs_tcga" = list(
    file = "./data/ucs_tcga_gdc_read_counts.txt",
    treated_cols = "TCGA-N5-A4R8-01A", "TCGA-N5-A4RA-01A",	
    wt_cols = "TCGA-N5-A4RD-01A", "TCGA-N5-A4RF-01A"
  ), #https://https://www.cbioportal.org/
  "plmeso_tcga" = list(
    file = "./data/plmeso_tcga_gdc_seq_read_counts.txt",
    treated_cols = "TCGA-3H-AB3K-01A", "TCGA-3H-AB3L-01A",	
    wt_cols = "	TCGA-3H-AB3M-01A", "TCGA-3H-AB3O-01A"
  ), #https://https://www.cbioportal.org/
  "GSE279712" = list(
    file = "./data/GSE279712_human.merged.counts.txt",
    treated_cols = "104307-001-001R", "104307-001-072",	
    wt_cols = "104307-001-040", "104307-001-081"
  ) #https://https://www.cbioportal.org/
)
length(dataset_configs)
# Funzione per generare prompts basati sulle colonne
generate_prompts <- function(treated_cols, wt_cols) {
  treated_underscore <- gsub(", ", "_", treated_cols)
  wt_underscore <- gsub(", ", "_", wt_cols)

  prompts <- c(
    sprintf("perform rnaseq analysis between treated %s and wt %s samples, up regulated", treated_cols, wt_cols),
    sprintf("compare treated %s vs control %s down regulated genes", treated_cols, wt_cols),
    sprintf("analyze treated columns %s against wild type columns %s, only genes down regulated", treated_cols, wt_cols),
    sprintf("differential expression analysis: treated samples %s vs wild type samples %s, all de regulated genes", treated_cols, wt_cols),
    sprintf("compare columns %s (treated) with columns %s (control), show upregulated", treated_cols, wt_cols),
    sprintf("DE analysis treated_%s vs control_%s up regulated", treated_underscore, wt_underscore),
    sprintf("rnaseq comparison: treated samples %s versus wt samples %s, down regulated genes", treated_cols, wt_cols),
    sprintf("perform differential gene expression analysis between treatment group (samples %s) and control group (samples %s), focus on upregulated genes", treated_cols, wt_cols),
    sprintf("make the comparison of the treated patients %s against wt ones %s, only down genes", treated_cols, wt_cols),
    sprintf("calculate significant differential genes between treatment samples (samples %s) and control ones (samples %s), expecially upregulated genes", treated_cols, wt_cols),
    sprintf("do the analysis of the %s treated and of the %s control. I want the up expressed genes only", treated_cols, wt_cols),
    sprintf(" %s treated  %s wt, de regulated", treated_cols, wt_cols),
    sprintf("generate only de between %s treated  %s wild type", treated_cols, wt_cols),
    sprintf("comparison between %s and %s, that are treated and wild types respectively, just up regulated genes", treated_cols, wt_cols),
    sprintf("I want to know de regulated between treated samples %s vs wild type samples %s", treated_cols, wt_cols)
    
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

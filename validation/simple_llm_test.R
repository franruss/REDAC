# Script di Test Semplificato per REDAC
# Replica della funzione `output$chat_output2` per testing automatizzato
# Basato sul codice in validation/validation.r

# ==============================================================================
# SETUP
# ==============================================================================
library(markdown)
library(shiny)
library(shinydashboard)
library(gplots)
library(ggplot2)
library(shinythemes)
library(mgcv)
library(lattice)
library(SummarizedExperiment)
library(edgeR)
library(vioplot)
library(preprocessCore)
library(AnnotationDbi)
library(clusterProfiler)
library(org.Hs.eg.db)
library(httr)
library(jsonlite)
library(BiocManager)
library(enrichplot)
library(DOSE)
library(stringr)
library(knitr)
library(tidyr)
library(dplyr, warn.conflicts = FALSE)
library(igraph)
library(FactoMineR)
library(factoextra)
library(plotly)
library(plot3D)
library(ggraph)
library(ggrepel)
library(DT)
library(quantmod)
library(tibble)
library(reshape2)
library(ggdendro) 
library(tinytex)
library(readxl)
library(janitor)
library(rmdHelpers)


# Carica ambiente
source("./definitions.R")
source("./.Renviron")
readRenviron("./.Renviron")

# Librerie
library(httr)
library(jsonlite)
library(markdown)

# Configurazione API
api_key <- Sys.getenv("API_KEY")
provider_url <- "https://api.together.xyz/v1/chat/completions"

# ==============================================================================
# FUNZIONE PRINCIPALE (REPLICA DA validation.r)
# ==============================================================================

test_llm_response <- function(dataset_file, prompt_text, model_name = "meta-llama/Llama-3.3-70B-Instruct-Turbo-Free") {
  
  cat("========================================\n")
  cat("Testing LLM Response\n")
  cat("Dataset:", dataset_file, "\n")
  cat("Prompt:", prompt_text, "\n") 
  cat("Model:", model_name, "\n")
  cat("========================================\n")
  
  # Verifica file dataset
  if (!file.exists(dataset_file)) {
    stop("Dataset file not found: ", dataset_file)
  }
  
  # Verifica API key
  if (api_key == "" || is.na(api_key)) {
    stop("NO API key found! Check .Renviron file")
  }
  
  # Leggi dati
  my_data <- read_and_clean_colnames(dataset_file, sep = "\t", header = TRUE, rownames = 1)
  
  # Prepara testo dati (prime 3 righe)
  data_text <- paste(capture.output(head(my_data, 3)), collapse = "\n")
  
  # Prepara richiesta API (esattamente come nel codice originale)
  body <- list(
    model = model_name,
    messages = list(
      list(role = "system", content = "Write an R code that uses edgeR package"),
      list(role = "user", content = paste(prompt_text, " on a dataset that have: ", 
                                         dim(my_data)[1], " rows, and ", dim(my_data)[2], 
                                         " columns. Here are the first 3 rows of my dataframe: \n", 
                                         data_text, " use edgeR package for this request: ", prompt_text))
    ),
    max_tokens = 1500
  )
  
  # Fai richiesta
  response <- POST(
    provider_url,
    add_headers(
      Authorization = paste("Bearer", api_key),
      `Content-Type` = "application/json"
    ),
    body = toJSON(body, auto_unbox = TRUE),
    encode = "json"
  )
  
  # Processa risposta
  if (status_code(response) != 200) {
    stop("API request failed with status: ", status_code(response))
  }
  
  content <- httr::content(response, as = "text")
  parsed <- fromJSON(content)
  
  # Estrai contenuto
  llm_response <- parsed$choices$message$content
  
  # Converti da markdown a HTML (come nel codice originale)
  html_output <- markdown::markdownToHTML(text = llm_response)
  
  # Analizza la risposta
  analysis <- analyze_llm_response(llm_response, my_data)
  
  # Stampa risultati
  cat("\n--- LLM RESPONSE ---\n")
  cat(llm_response)
  cat("\n\n--- ANALYSIS ---\n")
  cat("Contains R code:", analysis$has_r_code, "\n")
  cat("Contains edgeR:", analysis$mentions_edger, "\n")
  cat("Code blocks found:", analysis$code_blocks_count, "\n")
  cat("Estimated quality:", analysis$quality_score, "/10\n")
  cat("Issues found:", paste(analysis$issues, collapse = ", "), "\n")
  
  # Ritorna tutto
  return(list(
    dataset_info = list(
      file = dataset_file,
      dimensions = dim(my_data),
      columns = colnames(my_data)
    ),
    prompt = prompt_text,
    model = model_name,
    raw_response = llm_response,
    html_response = html_output,
    analysis = analysis,
    timestamp = Sys.time()
  ))
}

# ==============================================================================
# FUNZIONE DI ANALISI DELLA RISPOSTA
# ==============================================================================

analyze_llm_response <- function(response_text, dataset) {
  
  # Inizializza analisi
  analysis <- list(
    has_r_code = FALSE,
    mentions_edger = FALSE,
    code_blocks_count = 0,
    quality_score = 0,
    issues = c(),
    code_snippets = c()
  )
  
  # Verifica presenza di codice R
  if (grepl("```r|```R|library\\(|install\\.packages|<-|function\\(", response_text)) {
    analysis$has_r_code <- TRUE
    analysis$quality_score <- analysis$quality_score + 3
  }
  
  # Verifica menzione di edgeR
  if (grepl("edgeR|EdgeR|exactTest|DGEList", response_text, ignore.case = TRUE)) {
    analysis$mentions_edger <- TRUE
    analysis$quality_score <- analysis$quality_score + 3
  }
  
  # Conta blocchi di codice
  code_blocks <- regmatches(response_text, gregexpr("```[^`]*```", response_text))[[1]]
  analysis$code_blocks_count <- length(code_blocks)
  analysis$code_snippets <- code_blocks
  
  if (analysis$code_blocks_count > 0) {
    analysis$quality_score <- analysis$quality_score + 2
  }
  
  # Verifica problemi comuni
  if (grepl("I cannot|I can't|I'm sorry|I don't", response_text)) {
    analysis$issues <- c(analysis$issues, "Refusal to help")
    analysis$quality_score <- analysis$quality_score - 3
  }
  
  if (grepl("hallucination|made up|fictional", response_text, ignore.case = TRUE)) {
    analysis$issues <- c(analysis$issues, "Potential hallucination")
    analysis$quality_score <- analysis$quality_score - 2
  }
  
  if (!grepl("differential|expression|DE", response_text, ignore.case = TRUE)) {
    analysis$issues <- c(analysis$issues, "Missing key concepts")
    analysis$quality_score <- analysis$quality_score - 1
  }
  
  # Verifica sintassi R di base
  if (analysis$has_r_code) {
    if (!grepl("library\\(|require\\(", response_text)) {
      analysis$issues <- c(analysis$issues, "Missing library calls")
    }
    
    if (grepl("col[0-9]+|column[0-9]+", response_text) && 
        !any(grepl(paste0("col", 1:ncol(dataset), collapse = "|"), response_text))) {
      analysis$issues <- c(analysis$issues, "Invalid column references")
    }
  }
  
  # Assicura che il punteggio sia tra 0 e 10
  analysis$quality_score <- max(0, min(10, analysis$quality_score))
  
  if (length(analysis$issues) == 0) {
    analysis$issues <- "None detected"
  }
  
  return(analysis)
}

# ==============================================================================
# FUNZIONE PER TEST MULTIPLI
# ==============================================================================

run_multiple_tests <- function(dataset_files, prompts, model_name = "meta-llama/Llama-3.3-70B-Instruct-Turbo-Free") {
  
  results <- list()
  test_count <- 0
  
  for (dataset in dataset_files) {
    for (prompt in prompts) {
      test_count <- test_count + 1
      
      cat("\n", "="*60, "\n")
      cat("TEST", test_count, "of", length(dataset_files) * length(prompts), "\n")
      
      tryCatch({
        result <- test_llm_response(dataset, prompt, model_name)
        results[[paste0("test_", test_count)]] <- result
      }, error = function(e) {
        cat("ERROR in test", test_count, ":", e$message, "\n")
        results[[paste0("test_", test_count)]] <- list(
          error = e$message,
          dataset = dataset,
          prompt = prompt,
          timestamp = Sys.time()
        )
      })
      
      # Pausa tra richieste
      Sys.sleep(2)
    }
  }
  
  # Genera summary
  generate_test_summary(results)
  
  return(results)
}

# ==============================================================================
# FUNZIONE SUMMARY
# ==============================================================================

generate_test_summary <- function(results) {
  
  cat("\n", "="*80, "\n")
  cat("TEST SUMMARY\n")
  cat("="*80, "\n")
  
  total_tests <- length(results)
  successful_tests <- sum(sapply(results, function(x) is.null(x$error)))
  
  cat("Total tests:", total_tests, "\n")
  cat("Successful:", successful_tests, "\n")
  cat("Failed:", total_tests - successful_tests, "\n")
  
  if (successful_tests > 0) {
    # Analizza qualità delle risposte
    quality_scores <- sapply(results, function(x) {
      if (is.null(x$error) && !is.null(x$analysis$quality_score)) {
        return(x$analysis$quality_score)
      } else {
        return(NA)
      }
    })
    
    quality_scores <- quality_scores[!is.na(quality_scores)]
    
    if (length(quality_scores) > 0) {
      cat("Average quality score:", round(mean(quality_scores), 2), "/10\n")
      cat("Best score:", max(quality_scores), "\n")
      cat("Worst score:", min(quality_scores), "\n")
    }
    
    # Analizza problemi comuni
    all_issues <- unlist(sapply(results, function(x) {
      if (is.null(x$error) && !is.null(x$analysis$issues)) {
        return(x$analysis$issues)
      } else {
        return(NULL)
      }
    }))
    
    if (length(all_issues) > 0 && !all(all_issues == "None detected")) {
      cat("\nCommon issues:\n")
      issue_table <- table(all_issues[all_issues != "None detected"])
      for (i in 1:length(issue_table)) {
        cat("-", names(issue_table)[i], ":", issue_table[i], "times\n")
      }
    }
  }
  
  cat("="*80, "\n")
}

# ==============================================================================
# CONFIGURAZIONE TEST PREDEFINITI
# ==============================================================================

# Dataset di test (modifica i percorsi secondo la tua struttura)
test_datasets <- c(
  "../data/GSE164073_Eye_count_matrix.tsv.txt"
)

# Prompts di test (semplici per iniziare)
test_prompts <- c(
  "Plot a heatmap",
  "Perform differential expression analysis", 
  "Create a volcano plot",
  "Generate a PCA plot",
  "Make a boxplot of the data",
  "Show me the data distribution"
)

# ==============================================================================
# ESEMPIO DI UTILIZZO
# ==============================================================================

cat("Simple LLM Testing Script Loaded!\n")
cat("Available functions:\n")
cat("- test_llm_response(dataset_file, prompt_text, model_name)\n")
cat("- run_multiple_tests(dataset_files, prompts, model_name)\n\n")

cat("Example usage:\n")
cat('result <- test_llm_response("../data/GSE164073_Eye_count_matrix.tsv.txt", "Plot a heatmap")\n')
cat('multi_results <- run_multiple_tests(test_datasets, test_prompts)\n')

# Uncomment to run a quick test:
quick_test <- test_llm_response("./data/GSE164073_Eye_count_matrix.tsv.txt", "Plot a heatmap")

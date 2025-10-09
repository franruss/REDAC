# Script di Test per Valutazione Prompts REDAC
# Trasformazione della funzione Shiny in script standalone per testare 
# diversi prompts su diversi datasets senza interfaccia grafica

# ==============================================================================
# SETUP E LIBRERIE
# ==============================================================================
library(markdown)
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

# ==============================================================================
# Carica le definizioni e l'ambiente
source("./definitions.R")
source("./.Renviron")
readRenviron("./.Renviron")

# Librerie necessarie
library(BiocManager)
options(repos = BiocManager::repositories())
library(httr)
library(jsonlite)
library(DT)
library(edgeR)
library(DESeq2)
library(markdown)

# Configurazione API
provider <- "together"
provider_url <- "https://api.together.xyz/v1/chat/completions"
api_key <- Sys.getenv("API_KEY")

# Modelli disponibili
MODELS_ <- c("Gemma-3n" = "google/gemma-3n-E4B-it", 
             "Llama-3.3-70B" = "meta-llama/Llama-3.3-70B-Instruct-Turbo-Free")

url <- provider_url

# ==============================================================================
# FUNZIONI DI SUPPORTO
# ==============================================================================

# Funzione per testare un singolo prompt su un dataset
test_prompt_on_dataset <- function(prompt, dataset_path, model_name, test_id = NULL) {
  
  cat("=============================================================\n")
  cat("Test ID:", ifelse(is.null(test_id), "N/A", test_id), "\n")
  cat("Prompt:", prompt, "\n")
  cat("Dataset:", dataset_path, "\n")
  cat("Model:", model_name, "\n")
  cat("Timestamp:", Sys.time(), "\n")
  cat("=============================================================\n")
  
  # Inizializza risultato del test
  test_result <- list(
    test_id = test_id,
    prompt = prompt,
    dataset_path = dataset_path,
    model_name = model_name,
    timestamp = Sys.time(),
    success = FALSE,
    error_type = NULL,
    error_message = NULL,
    execution_time = NULL,
    result_summary = NULL,
    parsed_json = NULL
  )
  
  start_time <- Sys.time()
  
  tryCatch({
    # Verifica che esista il dataset
    if (!file.exists(dataset_path)) {
      stop("Dataset file does not exist: ", dataset_path)
    }
    
    # Leggi i dati
    my_data <- read_and_clean_colnames(dataset_path, sep = "\t", header = TRUE, rownames = 1)
    
    # Prepara i dati per il prompt (prime 3 righe)
    data_text <- paste(capture.output(head(my_data, 3)), collapse = "\n")
    
    # Verifica API key
    if (api_key == "" || is.na(api_key)) {
      stop("NO API key found! Check .Renviron file")
    }
    
    # Prepara la richiesta API
    body <- list(
      model = model_name,
      messages = list(
        list(role = "system", content = "Respond only in JSON format. The JSON must have this structure: 
             {\"functiontoberun\":[\"analysis\"],
             \"treated\":[[\"col1\"],[\"col2\"]],\"wt\":[[\"col3\"],[\"col4\"]],
             \"regulated\":[\"de\",\"down\",\"up\"]}. Do not add explanations or other text."),
        list(role = "user", content = prompt)
      ),
      max_tokens = 500
    )
    
    # Fai la richiesta API
    response <- POST(
      url,
      add_headers(
        Authorization = paste("Bearer", api_key),
        `Content-Type` = "application/json"
      ),
      body = toJSON(body, auto_unbox = TRUE),
      encode = "json"
    )
    
    # Verifica risposta HTTP
    if (status_code(response) != 200) {
      stop("API request failed with status: ", status_code(response))
    }
    
    content <- httr::content(response, as = "text")
    parsed <- fromJSON(content)
    
    # Gestisci diversi formati di risposta
    json_content <- parsed$choices$message$content
    
    # Se la risposta è in formato markdown con ```json, estrai il JSON
    if (grepl("```json", json_content)) {
      json_lines <- strsplit(json_content, "\n")[[1]]
      json_line_idx <- which(json_lines != "" & !grepl("```", json_lines))
      if (length(json_line_idx) > 0) {
        json_content <- json_lines[json_line_idx[1]]
      }
    }
    
    # Prova a parsare il JSON
    parsed_json <- fromJSON(json_content)
    test_result$parsed_json <- parsed_json
    
    # Verifica struttura JSON
    validation_result <- validate_json_structure(parsed_json, my_data)
    
    if (validation_result$valid) {
      # Se JSON è valido, prova ad eseguire l'analisi
      if (!is.null(parsed_json$functiontoberun) && 
          length(parsed_json$treated) > 0 && 
          length(parsed_json$wt) > 0) {
        
        # Esegui l'analisi differenziale
        analysis_result <- execute_differential_analysis(
          my_data, 
          parsed_json$treated, 
          parsed_json$wt, 
          parsed_json$regulated[1]
        )
        
        test_result$success <- TRUE
        test_result$result_summary <- paste(
          "Analysis completed successfully.",
          "Genes analyzed:", nrow(analysis_result),
          "Significant genes (FDR < 0.05):", sum(analysis_result$FDR < 0.05, na.rm = TRUE)
        )
        
      } else {
        test_result$error_type <- "Invalid JSON structure"
        test_result$error_message <- "Missing required fields in JSON response"
      }
    } else {
      test_result$error_type <- validation_result$error_type
      test_result$error_message <- validation_result$error_message
    }
    
  }, error = function(e) {
    test_result$error_type <<- classify_error(e$message)
    test_result$error_message <<- e$message
  })
  
  test_result$execution_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  # Stampa risultato
  cat("Result:", ifelse(test_result$success, "SUCCESS", "FAILURE"), "\n")
  if (!test_result$success) {
    cat("Error Type:", test_result$error_type, "\n")
    cat("Error Message:", test_result$error_message, "\n")
  } else {
    cat("Summary:", test_result$result_summary, "\n")
  }
  cat("Execution Time:", round(test_result$execution_time, 2), "seconds\n")
  cat("\n")
  
  return(test_result)
}

# Funzione per validare la struttura JSON
validate_json_structure <- function(parsed_json, dataset) {
  
  result <- list(valid = FALSE, error_type = NULL, error_message = NULL)
  
  if (is.null(parsed_json) || !is.list(parsed_json)) {
    result$error_type <- "Parsing/formatting errors"
    result$error_message <- "Failed to parse JSON response or not a valid list"
    return(result)
  }
  
  # Verifica campi richiesti
  required_fields <- c("functiontoberun", "treated", "wt", "regulated")
  missing_fields <- setdiff(required_fields, names(parsed_json))
  
  if (length(missing_fields) > 0) {
    result$error_type <- "Invalid/unrecognized function"
    result$error_message <- paste("Missing required fields:", paste(missing_fields, collapse = ", "))
    return(result)
  }
  
  # Verifica che le colonne specificate esistano nel dataset
  all_cols <- c(unlist(parsed_json$treated), unlist(parsed_json$wt))
  all_cols <- as.numeric(gsub("\\D", "", all_cols))  # Estrai numeri
  
  max_col <- ncol(dataset)
  invalid_cols <- all_cols[all_cols > max_col | all_cols < 1]
  
  if (length(invalid_cols) > 0) {
    result$error_type <- "Nonexistent sample reference"
    result$error_message <- paste("Referenced column(s) do not exist in dataset:", 
                                 paste(invalid_cols, collapse = ", "),
                                 ". Dataset has", max_col, "columns")
    return(result)
  }
  
  result$valid <- TRUE
  return(result)
}

# Funzione per eseguire l'analisi differenziale
execute_differential_analysis <- function(counts, treated_samples, wt_samples, regulated) {
  
  # Estrai numeri dalle colonne
  treated_samples <- as.numeric(gsub("\\D", "", unlist(treated_samples)))
  wt_samples <- as.numeric(gsub("\\D", "", unlist(wt_samples)))
  
  # Usa la funzione de_analysis dal definitions.R
  result <- de_analysis(counts, treated_samples, wt_samples, regulated)
  
  return(result)
}

# Funzione per classificare il tipo di errore
classify_error <- function(error_message) {
  
  if (grepl("parse|JSON|format", error_message, ignore.case = TRUE)) {
    return("Parsing/formatting errors")
  }
  
  if (grepl("function|method|recognized", error_message, ignore.case = TRUE)) {
    return("Invalid/unrecognized function")
  }
  
  if (grepl("column|sample|exist", error_message, ignore.case = TRUE)) {
    return("Nonexistent sample reference")
  }
  
  if (grepl("API|key|request", error_message, ignore.case = TRUE)) {
    return("API/Connection error")
  }
  
  return("Execution error in R")
}

# ==============================================================================
# FUNZIONE PRINCIPALE PER TESTARE MULTIPLE COMBINAZIONI
# ==============================================================================

run_comprehensive_test <- function(prompts_list, datasets_list, models_list = MODELS_, output_file = NULL) {
  
  cat("Starting comprehensive prompt testing...\n")
  cat("Prompts:", length(prompts_list), "\n")
  cat("Datasets:", length(datasets_list), "\n") 
  cat("Models:", length(models_list), "\n")
  cat("Total tests:", length(prompts_list) * length(datasets_list) * length(models_list), "\n\n")
  
  all_results <- list()
  test_counter <- 0
  
  for (model_name in models_list) {
    for (dataset_path in datasets_list) {
      for (prompt in prompts_list) {
        
        test_counter <- test_counter + 1
        test_id <- paste0("test_", sprintf("%04d", test_counter))
        
        result <- test_prompt_on_dataset(
          prompt = prompt,
          dataset_path = dataset_path, 
          model_name = model_name,
          test_id = test_id
        )
        
        all_results[[test_id]] <- result
        
        # Pausa tra richieste per evitare rate limiting
        Sys.sleep(1)
      }
    }
  }
  
  # Salva risultati se specificato
  if (!is.null(output_file)) {
    saveRDS(all_results, output_file)
    # save also as txt
    # write.table(do.call(rbind, lapply(all_results, as.data.frame)), 
    #                 file = sub("\\.rds$", ".txt", output_file), 
    #                 sep = "\t", row.names = FALSE, quote = FALSE)
    cat("Results saved to:", output_file, "\n")
  }
  
  # Genera report summary
  generate_summary_report(all_results)
  
  return(all_results)
}

# Funzione per generare report riassuntivo
generate_summary_report <- function(results) {
  
  cat("\n==============================================================================\n")
  cat("SUMMARY REPORT\n")
  cat("==============================================================================\n")
  
  total_tests <- length(results)
  successful_tests <- sum(sapply(results, function(x) x$success))
  failed_tests <- total_tests - successful_tests
  
  cat("Total tests:", total_tests, "\n")
  cat("Successful tests:", successful_tests, "(", round(successful_tests/total_tests*100, 1), "%)\n")
  cat("Failed tests:", failed_tests, "(", round(failed_tests/total_tests*100, 1), "%)\n\n")
  
  # Analisi errori per tipo
  error_types <- sapply(results, function(x) if(!x$success) x$error_type else NA)
  error_types <- error_types[!is.na(error_types)]
  
  if (length(error_types) > 0) {
    cat("Error breakdown:\n")
    error_table <- table(error_types)
    for (i in 1:length(error_table)) {
      cat("-", names(error_table)[i], ":", error_table[i], 
          "(", round(error_table[i]/failed_tests*100, 1), "% of failures)\n")
    }
    cat("\n")
  }
  
  # Analisi per modello
  models_success <- tapply(sapply(results, function(x) x$success), 
                          sapply(results, function(x) x$model_name), 
                          function(x) sum(x)/length(x)*100)
  
  cat("Success rate by model:\n")
  for (model in names(models_success)) {
    cat("-", model, ":", round(models_success[model], 1), "%\n")
  }
  
  cat("\n==============================================================================\n")
}

# ==============================================================================
# CONFIGURAZIONE TEST PREDEFINITI
# ==============================================================================

# Carica la generazione dei prompts
source("validation/prompt_gen.R")

# Lista di prompts e datasets da testare
test_prompts <- c(
  dataset_prompts$RNAseq_submission$prompts,
  dataset_prompts$GSE164073_Eye$prompts,
  dataset_prompts$GSE146458_raw$prompts
)

test_datasets <- c(
  dataset_configs$RNAseq_submission$file,
  dataset_configs$GSE164073_Eye$file,
  dataset_configs$GSE146458_raw$file
)

# ==============================================================================
# ESECUZIONE TEST
# ==============================================================================

# Funzione per eseguire test semplificato (solo un modello)
run_simple_test <- function() {
  
  # Test su un singolo dataset e modello
  dataset_path <- "./data/GSE164073_Eye_count_matrix.tsv.txt"
  model_name <- MODELS_["Llama-3.3-70B"]
  
  cat("Running simple test with", length(test_prompts), "prompts...\n\n")
  
  results <- list()
  for (i in 1:length(test_prompts)) {
    result <- test_prompt_on_dataset(
      prompt = test_prompts[i],
      dataset_path = dataset_path,
      model_name = model_name,
      test_id = paste0("simple_test_", i)
    )
    results[[i]] <- result
    Sys.sleep(1) # Pausa tra richieste
  }
  
  generate_summary_report(results)
  return(results)
}

# ==============================================================================
# ESEMPIO DI UTILIZZO
# ==============================================================================

# Per testare singolo prompt:
# result <- test_prompt_on_dataset(
#   prompt = "perform rnaseq analysis between treated 1,2 and wt 3,4 samples, up regulated",
#   dataset_path = "./data/GSE164073_Eye_count_matrix.tsv.txt", 
#   model_name = MODELS_["Llama-3.3-70B"]
# )

# Per test semplificato:
# simple_results <- run_simple_test()

# Per test completo (attenzione: molto lungo!):
# comprehensive_results <- run_comprehensive_test(
#   prompts_list = test_prompts,
#   datasets_list = test_datasets,
#   models_list = MODELS_,
#   output_file = "test_results.rds"
# )

# Per test completo (Specifico):
comprehensive_results <- run_comprehensive_test(
  prompts_list = dataset_prompts$RNAseq_submission$prompts,
  datasets_list = test_datasets[1],
  models_list = MODELS_["Llama-3.3-70B"],
  output_file = "validation/test_results.rds"
)
# generate_summary_report(comprehensive_results)

cat("Script loaded successfully!\n")
cat("Available functions:\n")
cat("- test_prompt_on_dataset(): Test single prompt\n") 
cat("- run_simple_test(): Quick test with multiple prompts\n")
cat("- run_comprehensive_test(): Full test suite\n\n")
cat("Example usage:\n")
cat("simple_results <- run_simple_test()\n")


# read test results.rds
# test_results <- readRDS("validation/test_results.rds")

test_prompts[1]

# Script di Test per REDAC - Validazione Prompts

Questi script trasformano le funzioni dell'app Shiny REDAC in strumenti standalone per testare automaticamente diversi prompts su diversi datasets senza interfaccia grafica.

## File Creati

### 1. `test_prompts_script.R` - Script Completo di Test
Script principale che replica la funzione di analisi differenziale dall'app Shiny e la adatta per test automatizzati.

**Caratteristiche:**
- Test di prompt multipli su dataset multipli
- Supporto per modelli LLM diversi
- Validazione della struttura JSON delle risposte
- Classificazione automatica degli errori
- Report riassuntivi dettagliati
- Esecuzione dell'analisi differenziale quando il parsing è corretto

**Tipi di errore rilevati:**
- Parsing/formatting errors (errori di parsing JSON)
- Invalid/unrecognized function (funzioni non riconosciute)
- Execution error in R (errori durante l'esecuzione)
- Nonexistent sample reference (riferimenti a campioni inesistenti)
- API/Connection errors (errori di connessione API)

### 2. `simple_llm_test.R` - Script Semplificato
Script più semplice che replica la funzione `output$chat_output2` dall'app per testare la generazione di codice R.

**Caratteristiche:**
- Test singoli o multipli
- Analisi qualitativa delle risposte LLM
- Rilevamento di codice R e menzioni di edgeR
- Punteggio di qualità automatico
- Identificazione di problemi comuni

## Configurazione

### 1. File .Renviron
Assicurati di avere un file `.Renviron` nella directory principale con:
```
API_KEY=your_together_api_key_here
```

### 2. Struttura Directory
Assicurati che la struttura delle directory sia:
```
REDAC/
├── validation/
│   ├── test_prompts_script.R
│   ├── simple_llm_test.R
│   └── validation.r (originale)
├── definitions.R
├── .Renviron
└── data/
    ├── GSE164073_Eye_count_matrix.tsv.txt
    ├── GSE146458_raw_counts_GRCh38.p13_NCBI.txt
    └── Raw-data-RNAseq-for-submission-with-gene-name-3-1-24.csv
```

## Utilizzo

### Script Completo (`test_prompts_script.R`)

```r
# Carica lo script
source("test_prompts_script.R")

# Test singolo
result <- test_prompt_on_dataset(
  prompt = "perform rnaseq analysis between treated 1,2 and wt 3,4 samples, up regulated",
  dataset_path = "../data/GSE164073_Eye_count_matrix.tsv.txt", 
  model_name = MODELS_["Llama-3.3-70B"]
)

# Test rapido con più prompts
simple_results <- run_simple_test()

# Test completo (attenzione: può durare molto!)
comprehensive_results <- run_comprehensive_test(
  prompts_list = test_prompts,
  datasets_list = test_datasets,
  models_list = MODELS_,
  output_file = "test_results.rds"
)
```

### Script Semplificato (`simple_llm_test.R`)

```r
# Carica lo script
source("simple_llm_test.R")

# Test singolo
result <- test_llm_response(
  dataset_file = "../data/GSE164073_Eye_count_matrix.tsv.txt",
  prompt_text = "Plot a heatmap"
)

# Test multipli
multi_results <- run_multiple_tests(
  dataset_files = test_datasets, 
  prompts = test_prompts
)
```

## Prompts di Test Predefiniti

Lo script include diversi prompts di test che variano per:
- Formato del comando
- Specificazione dei campioni (numeri vs nomi)
- Tipo di regolazione richiesta (up/down/de)
- Dettaglio della descrizione

Esempi:
```r
test_prompts <- c(
  "perform rnaseq analysis between treated 1,2 and wt 3,4 samples, up regulated",
  "compare treated 1,2,3 vs control 4,5,6, down-regulated genes", 
  "analyze first two treated against second two wild type, only genes down regulated",
  "differential expression analysis: treated samples 1,2 vs wild type samples 3,4, all differentially expressed genes"
)
```

## Output e Report

### Report Riassuntivo
Lo script genera automaticamente report che includono:
- Numero totale di test
- Tasso di successo/fallimento
- Breakdown degli errori per tipo
- Prestazioni per modello
- Tempi di esecuzione

### Esempio di Output:
```
==============================================================================
SUMMARY REPORT
==============================================================================
Total tests: 24
Successful tests: 18 (75.0%)
Failed tests: 6 (25.0%)

Error breakdown:
- Parsing/formatting errors: 3 (50.0% of failures)
- Nonexistent sample reference: 2 (33.3% of failures)
- Execution error in R: 1 (16.7% of failures)

Success rate by model:
- Llama-3.3-70B: 80.0%
- Gemma-3n: 70.0%
```

## Personalizzazione

### Aggiungere Nuovi Dataset
```r
test_datasets <- c(
  test_datasets,
  "path/to/your/new/dataset.txt"
)
```

### Aggiungere Nuovi Prompts
```r
custom_prompts <- c(
  "your custom prompt here",
  "another custom prompt"
)
```

### Modificare Modelli
```r
custom_models <- c(
  "model1" = "provider/model1-name",
  "model2" = "provider/model2-name"
)
```

## Troubleshooting

### Errori Comuni
1. **"NO API key found!"** - Verifica il file .Renviron
2. **"Dataset file does not exist"** - Controlla i percorsi dei file
3. **"API request failed"** - Verifica connessione internet e validità API key
4. **Rate limiting** - Aumenta la pausa tra richieste (Sys.sleep)

### Debug
- I test stampano informazioni dettagliate durante l'esecuzione
- Usa `tryCatch` per gestire errori specifici
- Salva i risultati con `saveRDS()` per analisi posteriori

## Estensioni Future

Possibili miglioramenti:
- Integrazione con diversi provider LLM
- Test di performance e velocità
- Analisi semantica più approfondita delle risposte
- Test di robustezza con dataset malformati
- Integrazione con framework di testing automatico
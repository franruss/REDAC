# REDAC Validation Scripts

Standalone testing tools for REDAC Shiny app prompts validation without GUI interface.

## Scripts

### `test_prompts_script.R` - Comprehensive Testing
Main script for automated differential expression analysis testing.

**Features:**
- Multi-prompt, multi-dataset testing
- Multiple LLM model support
- JSON response validation
- Automatic error classification
- Detailed summary reports

### `simple_llm_test.R` - Simplified Testing
Basic script for R code generation testing.

**Features:**
- Single/multiple test runs
- Response quality scoring
- R code detection

## Setup

### 1. Environment File
Create `.Renviron` in the main directory:
```
API_KEY=your_together_api_key_here
```

### 2. Directory Structure
```
REDAC/
├── validation/
│   ├── test_prompts_script.R
│   └── simple_llm_test.R
├── definitions.R
├── .Renviron
└── data/
    └── [dataset files]
```

## Usage

### Comprehensive Testing
```r
source("test_prompts_script.R")

# Single test
result <- test_prompt_on_dataset(
  prompt = "perform rnaseq analysis between treated 1,2 and wt 3,4 samples",
  dataset_path = "../data/GSE164073_Eye_count_matrix.tsv.txt"
)

# Multiple tests
results <- run_simple_test()
```

### Simple Testing
```r
source("simple_llm_test.R")

result <- test_llm_response(
  dataset_file = "../data/GSE164073_Eye_count_matrix.tsv.txt",
  prompt_text = "Plot a heatmap"
)
```

## Test Results

Scripts automatically generate reports including:
- Success/failure rates
- Error type breakdown
- Model performance comparison
- Execution times

## Common Issues

1. **"NO API key found!"** - Check `.Renviron` file
2. **"Dataset file does not exist"** - Verify file paths
3. **"API request failed"** - Check internet connection and API key validity
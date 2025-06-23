source("definitions.R")
library(BiocManager)
options(repos = BiocManager::repositories())

shinyServer(function(input, output, session) {
  
  # Function to check RNA-Seq count table
  check_rnaseq_counts <- function(data) {
    
    # Check if the file has at least two columns (for multiple samples)
    if (ncol(data) < 2) {
      stop("The file does not contain multiple samples. RNA-Seq count data should have multiple sample columns.")
    }
    
    # Check if row names resemble gene IDs (basic check: non-numeric row names)
    if (all(grepl("^[0-9]+$", rownames(data)))) {
      stop("Row names appear to be numeric. Gene identifiers should be non-numeric (e.g., Ensembl IDs or gene symbols).")
    }
    
    # Check if all values are numeric
    if (!all(sapply(data, is.numeric))) {
      stop("Some values are non-numeric. RNA-Seq count data should only contain numeric values.")
    }
    
    # Check if values are non-negative integers
    if (!all(data >= 0)) {
      stop("RNA-Seq count data should not contain negative values.")
    }
    
    # Additional check: If values are ALL integers
    if (mean(data %% 1 == 0) != 1) {
      warning("There are non-integers in the file that you provided. RNA-Seq count data must contain integers.")
    }
    
    cat("The file appears to be an RNA-Seq count table.\n")
    return(TRUE)
  }
  
  # Function to check if the file is suitable for a plot on results
  check_de_result_data <- function(data) {
    
    # Required columns
    required_cols <- c("logFC","log2FoldChange","FC","log2FC","FoldChange","Fold","Change","P.Value", "pvalue", "pval", "adj.P.Val", "padj","Padj","Pval","Pvalue")
    
    # Check if any column matches logFC
    logfc_col <- grep("logFC|log2FoldChange|FC|log2FC|FoldChange|Fold Change", colnames(data), ignore.case = TRUE, value = TRUE)
    if (length(logfc_col) == 0) {
      stop("No log fold-change (logFC) column found. This is required to generate the plot.")
    }
    
    # Check if any column matches p-value
    pval_col <- grep("pval|p.value|pvalue|P.Value|pval|adj.P.Val|padj|Padj|Pval|Pvalue", colnames(data), ignore.case = TRUE, value = TRUE)
    if (length(pval_col) == 0) {
      stop("No p-value column found. This is required for a volcano plot.")
    }
    
    # Check if logFC values are numeric
    if (!is.numeric(data[[logfc_col[1]]])) {
      stop("logFC values are not numeric. Ensure the file contains valid numerical log fold-change values.")
    }
    
    # Check if p-value values are numeric and between 0 and 1
    if (!is.numeric(data[[pval_col[1]]]) || any(data[[pval_col[1]]] < 0 | data[[pval_col[1]]] > 1, na.rm = TRUE)) {
      stop("P-value column contains non-numeric values or values outside the range [0,1].")
    }
    
    cat("The file contains logFC and p-value columns and can be used to generate the plot.\n")
    return(TRUE)
  }
  
  # Function to check if the file is suitable for a plot on emrichment results
  check_enrichment_data <- function(data) {
    
    # Required columns
    required_cols <- c("Description",	"p.adjust",	"qvalue")
    
    # Check if any column matches Description
    logfc_col <- grep("Description", colnames(data), ignore.case = TRUE, value = TRUE)
    if (length(logfc_col) == 0) {
      stop("No Description column found. This is required to generate the plot.")
    }
    
    # Check if any column matches p.adjust
    pval_col <- grep("p.adjust", colnames(data), ignore.case = TRUE, value = TRUE)
    if (length(pval_col) == 0) {
      stop("No p.adjust column found. This is required for a volcano plot.")
    }
    
    cat("The file contains Description and p.adjust columns and can be used to generate the plot.\n")
    return(TRUE)
  }
  
 # answer_fun_plot <- eventReactive(input$run, {
 output$answer_fun_plot <- renderPlot({
   #req(input$text1, input$file5)
   prompt = input$text1
   inDataFile = input$file5
   #print(input$run)
      if (is.null(inDataFile) || (prompt=="")){return(NULL)}else{
        my_data <- read_and_clean_colnames(inDataFile$datapath, sep = "\t", header = TRUE, rownames = 1 )
        #   my_data <- read.table("expression_file.txt", header = TRUE, sep = "\t", row.names = 1)
        #   prompt = "could you make an violin (if you are able to) on columns 1,2,3?"
        api_key <- "..."  # Replace with your API key
        url <- "https://api.together.xyz/v1/chat/completions"
        
        body <- list(
          model = "google/gemma-2-27b-it",
          #model = "meta-llama/Llama-3.3-70B-Instruct-Turbo-Free",
          messages = list(
            list(role = "system", content = "Respond only in JSON format. The JSON must have this structure: 
               {\"functiontoberun\":[\"box\",\"violin\",\"heat\",\"corr heat\",\"volcano\",\"pca\",\"cluster\",
               \"dendro\",\"dens\",\"integration\",\"component\",
               \"MAplot\",\"dot\",\"net\",\"3D pca\"],
               \"columns\":[[\"col1\"],[\"col2\"],[\"col3\"]]}. Do not add explanations or other text."),
            list(role = "user", content = prompt)
          ),
          max_tokens = 500
        )
        response <- POST(
          url,
          add_headers(
            Authorization = paste("Bearer", api_key),
            `Content-Type` = "application/json"
          ),
          body = toJSON(body, auto_unbox = TRUE),
          encode = "json"
        )
        content <- httr::content(response, as = "text")
        #print(content)
        parsed <- fromJSON(content)
        #parsed <- tryCatch(fromJSON(content), error = function(e) NULL)
        #print(parsed$choices$message$content)
        # the output is in this format:
        # parsed$choices$message$content = "```json\n{\"function to be run\": [\"violin\"], \"columns to use\": [[\"1\", \"2\", \"3\"]]}\n```"
        # but it should be in this format:
        # parsed$choices$message$content = '{\"function to be run\": [\"violin\"], \"columns to use\": [[\"col1\"], [\"col2\"], [\"col3\"]]}'
        content = strsplit(parsed$choices$message$content,"\n" )[[1]][2] ## questa linea solo se usi Gemma LLM
        #print(content)
        #[1] "{\"function to be run\": [\"violin\"], \"columns\": [[\"1\"], [\"2\"], [\"3\"]]}"
        parsed <- fromJSON(content)
        print(parsed)
        #$`function to be run`
        #[1] "violin"
        if (is.null(parsed) || !is.list(parsed)) {
          ("Error: Failed to parse JSON response")
        }
        
        if (!is.null(parsed$functiontoberun) && !is.null(content)) {
          # Supponiamo che il dataframe sia chiamato df
          counts <- my_data  # Sostituiscilo con il dataframe effettivo
          # Eseguire le funzioni sui dati specificati
          for (i in seq_along(parsed$functiontoberun)) {
            func_name <- parsed$functiontoberun[i]
            cols <- unlist(parsed$columnstouse)
            #scegli qui le funzioni da usare
            if (func_name == "box") {
              if (check_rnaseq_counts(counts)){
                boxplot(log10(counts+1), main=paste("Boxplot"), 
                        col=c('red', 'blue', 'green','brown','purple','darkgreen','pink','orange','gold','darkblue','cyan','darkred'  ) , las=2  )
              }
            } else if (func_name == "violin") {
             # if (check_rnaseq_counts(counts)){
                vioplot(log10(counts+1), main=paste("Violin plot"),
                        col=c('red', 'blue', 'green','brown','purple','darkgreen','pink','orange','gold','darkblue','cyan','darkred' )  , las=2,  cexRow=0.4 )
            #  }
            } else if (func_name == "heat") {
              #if (check_rnaseq_counts(counts)){
              print("heat")
                if(dim(counts)[1]<50){
                  heatmap.2(log(as.matrix(counts+1)),main=paste("Heatmap of the most variable 50 genes",sep = ""),
                            scale='col',trace='none', col=colorRampPalette(c("red",'white',"blue"))(500), Colv = NA,Rowv = NA,
                            cexRow  = 0.4,cexCol  = 0.9, margins=c(8,6),dendrogram = "both",cex.main = 0.5)
                  }else{
                  vars <- sort(apply(counts, 1, var, na.rm = TRUE), decreasing = TRUE) #highest gene-wise variances
                  data <- counts[names(vars)[1:50], ]
                  heatmap.2(log(as.matrix(data+1)),main=paste("Heatmap of the most variable 50 genes",sep = ""),
                          scale='col',trace='none', col=colorRampPalette(c("red",'white',"blue"))(500), Colv = NA,Rowv = NA,
                          cexRow  = 0.4,cexCol  = 0.9, margins=c(8,6),dendrogram = "both",cex.main = 0.5)
                }
            #  }
            }else if (func_name == "cluster" || func_name == "dendro") {
             # if (check_rnaseq_counts(counts)){
                print("Clustering!")
                x = counts
                x = log10(x+1)
                #x = normalize.quantiles(as.matrix(x),copy=TRUE)
                # metodi per clusterizzare: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"
                colnames(x) = colnames(counts)
                rownames(x) = rownames(counts)
                clustering_used_method = "ward.D"
                distance_used_method = "euclidean"
                hclustfunc <- function(x)  hclust(x, method = clustering_used_method) # definisco la funzione hclustfunc da hclust
                distfunc <-   function(x)    dist(x, method = distance_used_method)      # definisco la funzione distfunc da dist
                fit <- hclustfunc(distfunc((t(x)))) #le innesto una dentro l'atra, chiamando prima distfunc
                plot(fit)
             # }
            }else if (func_name == "pca") {
             # if (check_rnaseq_counts(counts)){
                print("PCA")
                x = counts
                x_t <- t(x)
                # Perform PCA
                pca_result <- prcomp(x_t)
                # Extract PCA scores
                pca_data <- as.data.frame(pca_result$x)
                # Create a dataframe with sample names and first two principal components
                pca_data$Sample <- rownames(pca_data)
                # Plot PCA results (PC1 vs. PC2)
                plot(ggplot(pca_data, aes(x = PC1, y = PC2, label = Sample)) +
                       geom_point(color = "blue", size = 3) + 
                       geom_text(vjust = 1.5, hjust = 1.5, size = 3) + 
                       theme_minimal() +
                       labs(title = "PCA", x = "PC1", y = "PC2"))
             # }
            }else if (func_name == "volcano") {
              if (check_de_result_data(counts)){
                #Volcano Plot
                res_df = counts
                #plot(res$log2FoldChange,-log2(res$padj), pch=20, cex=0.2 ,ylab="-LOG2(P-VALUE ADJUSTED)", xlab="LOG2(FOLD CHANGE)",main= paste("Volcano Plot",sep=""))
                #abline(h=-log2(0.05),col='red')
                res_df$significant <- res_df$padj < 0.05
                ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
                  geom_point(alpha = 0.6) +
                  theme_minimal() +
                  labs(title = "Volcano Plot", x = "log2 Fold Change", y = "-log10 p-value")
              }
            }else if (func_name == "component") {
              #if (check_rnaseq_counts(counts)){
                screeplot(princomp(counts),main="Principal Component Histogram", col="blue", las=2)
              #}
            }else if (func_name == "dens") {
              if (check_rnaseq_counts(counts)){
                df = list(data.frame())
                x=counts
                for ( i in 1:length(colnames(x)) ) {
                  pp=x[,i]
                  colnames(x)[i]
                  col=rep(colnames(x)[i],length(pp))
                  df[[i]]=data.frame(counts=pp,idsample <- col)
                  if (i>1) {df[[i]] = rbind(df[[i-1]],df[[i]])}
                }
                df=data.frame(df[[length(colnames(counts))]])
                samples = df$idsample....col
                print(qplot(log(counts+1), main="Densities", data = df, geom = "density",colour=samples))
              }
            }else if (func_name == "MAplot") {
              if (check_de_result_data(counts)){
                res = counts
                plot(log2(res$baseMean + 1) , res$log2FoldChange ,col = "black", main="MA plot", xlab='log2(res$baseMean + 1)', ylab='log2FoldChange',pch=19,cex=0.3) 
                DE_genes = subset(res, padj<0.05)
                points(log2(DE_genes$baseMean + 1) , DE_genes$log2FoldChange, pch=19, col='red', cex=0.5)
                abline(h=0,col='red')
                
              }
            }else if (func_name == "dot") {
              if (check_de_result_data(counts)){
                #counts <- read.table("results_expression_file.txt", header = TRUE, sep = "\t", row.names = 1)
                print("dot performing")
                gene_symbol_list = rownames(subset(counts,as.numeric(counts$FDR)<0.05))
                length(gene_symbol_list)
                # Convert gene symbols to Entrez IDs
                entrez_ids <- mapIds(org.Hs.eg.db, keys = as.character(gene_symbol_list), column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
                entrez_ids <- na.omit(entrez_ids)
                de <- as.character(entrez_ids)
                x <- enrichKEGG(gene = de,  organism = 'hsa', pAdjustMethod="BH", qvalueCutoff = 0.001)
                plot(dotplot(x))
              }
            }else if (func_name == "net") {
              if (check_de_result_data(counts)){
                #counts <- read.table("results_expression_file.txt", header = TRUE, sep = "\t", row.names = 1)
                gene_symbol_list = rownames(subset(counts,as.numeric(counts$FDR)<0.05))
                length(gene_symbol_list)
                # Convert gene symbols to Entrez IDs
                entrez_ids <- mapIds(org.Hs.eg.db, keys = as.character(gene_symbol_list), column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
                entrez_ids <- na.omit(entrez_ids)
                de <- as.character(entrez_ids)
                x <- enrichKEGG(gene = de,  organism = 'hsa', pAdjustMethod="BH",qvalueCutoff = 0.001)
                edox <- setReadable(x, "org.Hs.eg.db", "ENTREZID")
                plot(cnetplot(edox, layout = igraph::layout_nicely,
                              showCategory = 5,
                              color_category = "black",
                              size_category = 0.5,
                              color_item = "red",
                              size_item = 0.5,
                              color_edge = "blue",
                              size_edge = 0.5,
                              node_label = "all",
                              foldChange = NULL,
                              hilight = "none",
                              hilight_alpha = 0.5))
              }
            }else if (func_name == "corr heat") {
              if (check_rnaseq_counts(counts)){
                # Calculate the correlation matrix
                cor_matrix <- cor(data)
                # Plot the correlation heatmap
                heatmap.2(cor_matrix, main = "Correlation Heatmap")
              }
            }else if (func_name == "3D pca") {
              if (check_rnaseq_counts(counts)){
                # Perform PCA (transpose if samples are columns!)
                pca_result <- PCA(t(counts), graph = FALSE)  # Note the transpose 't()' if samples are columns
                
                # Extract the PCA scores (rows are now samples after transposing)
                pca_data <- data.frame(
                  Sample = colnames(counts),
                  PC1 = pca_result$ind$coord[, 1],
                  PC2 = pca_result$ind$coord[, 2],
                  PC3 = pca_result$ind$coord[, 3]
                )
                
                # Create interactive 3D plot with sample labels as hover text
                p <- plot_ly(
                  data = pca_data,
                  x = ~PC1, y = ~PC2, z = ~PC3,
                  type = 'scatter3d',
                  mode = 'markers+text',
                  text = ~Sample,
                  hoverinfo = 'text',
                  marker = list(size = 5, color = 'blue')
                )
                
                p <- p %>% layout(
                  title = "3D PCA Plot",
                  scene = list(
                    xaxis = list(title = "PC1"),
                    yaxis = list(title = "PC2"),
                    zaxis = list(title = "PC3")
                  )
                )
                
                p
              }
            }else if (func_name == "integration") {
              print("integration")
              #gene_expr = read.table(file = gene_expr$datapath, header = TRUE,sep = "\t",row.names = 1)
              #ppi_data =  read.table(file = ppi_data$datapath, header = TRUE,sep = "\t", stringsAsFactors = FALSE)
              #counts <- read.table("C:/Users/f.russo.ENGIBBC/Desktop/REDAC_submission/", header = TRUE, sep = "\t", row.names = 1)
              #ppi_data  <- read.table("protein_protein_interaction_data.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
              ppi_data  <- read.table(file ="9606.protein.physical.links.detailed.v12.0_with_gene_symbols_.txt", header = TRUE,sep = "\t")
              protein_info <- read.table(file ="9606.protein.info.v12.0.txt", header = TRUE,sep = "\t")
              colnames(ppi_data) <- c("protein_code_2", "protein_code_1", "experimental", "Protein1", "Protein2")
              # colnames(ppi_data) <- c("Protein1", "Protein2")
              #Adesso elimino i conteggi nulli o quelli inferiori a 10 su tutti i campioni
              #keep <- rowSums(counts >= 10
              #dds <- dds[keep,]
              gene_expr = counts
              interaction_threshold = 970
              ppi_data = subset(ppi_data, ppi_data$experimental > interaction_threshold)
              print(dim(ppi_data))
              # Ensure genes in PPI match those in the gene expression dataset
              valid_genes <- rownames(gene_expr)
              ppi_filtered <- ppi_data %>% filter(Protein1 %in% valid_genes & Protein2 %in% valid_genes)
              
              # Function to compute Pearson correlation between two interacting proteins
              compute_correlation2 <- function(protein1, protein2, gene_expr) {
                expr1 <- gene_expr[protein1, ]
                expr2 <- gene_expr[protein2, ]
                
                if (sd(expr1) == 0 || sd(expr2) == 0) {
                  return(NA)  # Avoid zero variance issue
                }
                
                cor_value <- cor(as.numeric(expr1), as.numeric(expr2), method = "pearson", use = "complete.obs")
                cor_value 
              }
              
              ppi_filtered <- ppi_filtered[,3:5]
              ppi_filtered <- ppi_filtered[,c(2,3,1)]
              
              # Compute correlation for each interaction
              ppi_filtered <- ppi_filtered %>%
                rowwise() %>%
                mutate(Correlation = compute_correlation2(Protein1,Protein2, gene_expr) ) %>%
                drop_na()  # Remove NA values due to zero variance
              
              # Create a weighted graph
              ppi_graph <- graph_from_data_frame(ppi_filtered, directed = FALSE)
              
              # Plot the network using ggraph
              print("ggraph")
              p <- ggraph(ppi_graph, layout = "fr") +  
                geom_edge_link(aes(edge_alpha = abs(Correlation), 
                                   edge_width = abs(Correlation), 
                                   color = Correlation)) +
                geom_node_point(size = 1, color = "blue") +
                geom_text_repel(aes(x = x, y = y, label = name), size = 1, max.overlaps = 1000, segment.size = 0.2) +  
                scale_edge_color_gradient2(low = "blue", mid = "gray", high = "red", midpoint = 0) +  
                theme_minimal() +
                ggtitle("Gene Expression Weighted PPI Network") +
                theme(legend.position = "bottom")
              
              plot(p)
            }
          }   
          #return(parsed$choices$message$content) ##  Gemma output
        } else {(paste("Error: Unexpected response format -", content))}
      }
    })
  
  get_response_table <- function(prompt,inDataFile) {
   
    if (is.null(inDataFile) || (prompt=="")){return(NULL)}else{
      my_data <- read_and_clean_colnames(inDataFile$datapath, sep = "\t", header = TRUE, rownames = 1 )
      #   my_data <- read.table("C:/Users/f.russo.ENGIBBC/Desktop/REDAC_submission/expression_file.txt", header = TRUE, sep = "\t", row.names = 1)
      #   prompt = "Can you perform a de analysis of treated 1,2,5 against wt 3,7,8 ?"
      api_key <- "..."  # Replace with your API key
      url <- "https://api.together.xyz/v1/chat/completions"
      
      body <- list(
        model = "google/gemma-2-27b-it",
        #model = "meta-llama/Llama-3.3-70B-Instruct-Turbo-Free",
        messages = list(
          list(role = "system", content = "Respond only in JSON format. The JSON must have this structure: 
               {\"functiontoberun\":[\"analysis\"],
                \"treated\":[[\"col1\"],[\"col2\"]],\"wt\":[[\"col3\"],[\"col4\"]],
                \"regulated\":[\"down\",\"up\",\"de\"]}. Do not add explanations or other text."),
          list(role = "user", content = prompt)
        ),
        max_tokens = 500
      )
      response <- POST(
        url,
        add_headers(
          Authorization = paste("Bearer", api_key),
          `Content-Type` = "application/json"
        ),
        body = toJSON(body, auto_unbox = TRUE),
        encode = "json"
      )
      content <- httr::content(response, as = "text")
      #print(content)
      parsed <- fromJSON(content)
      #parsed <- tryCatch(fromJSON(content), error = function(e) NULL)
      #print(parsed$choices$message$content)
      # the output is in this format:
      # parsed$choices$message$content = "```json\n{\"function to be run\": [\"violin\"], \"columns to use\": [[\"1\", \"2\", \"3\"]]}\n```"
      # but it should be in this format:
      # parsed$choices$message$content = '{\"function to be run\": [\"violin\"], \"columns to use\": [[\"col1\"], [\"col2\"], [\"col3\"]]}'
      content = strsplit(parsed$choices$message$content,"\n" )[[1]][2] ## questa linea solo se usi Gemma LLM
      #print(content)
      #[1] "{\"function to be run\": [\"violin\"], \"columns\": [[\"1\"], [\"2\"], [\"3\"]]}"
      parsed <- fromJSON(content)
      print(parsed)
      #$`function to be run`
      #[1] "violin"
      if (is.null(parsed) || !is.list(parsed)) {
        ("Error: Failed to parse JSON response")
      }
      if (!is.null(parsed$functiontoberun) && length(parsed$treated) > 0  && length(parsed$wt) > 0 && !is.null(content)) {
       withProgress(message = 'Running, Please Wait!', detail = 'This may take a while...', value = 0, {     
        # Supponiamo che il dataframe sia chiamato df
        counts <- my_data  # Sostituiscilo con il dataframe effettivo
        
        # Eseguire le funzioni sui dati specificati
        for (i in seq_along(parsed$functiontoberun)) {
          func_name       <- parsed$functiontoberun[i]
          treated_samples <- unlist(str_replace_all(parsed$treated, c("col" = "")))
          wt_samples      <- unlist(str_replace_all(parsed$wt, c("col" = "")))
          regulated       <- parsed$regulated[i]
          #scegli qui le funzioni da usare
          if (func_name == "analysis") {
            print("rnaseqanalysis")
            final_result  = de_analysis(counts,colnames(counts)[as.numeric(treated_samples)],colnames(counts)[as.numeric(wt_samples)],regulated)
          } else if (func_name == "enrichment") {
            final_result = gene_pathway_enrichment(as.character(rownames(counts)))
          }
        }
        n <- 100
        for (i in 1:n) {
          # Increment the progress bar, and update the detail text.
          incProgress(1/n, detail = paste("\n  ",i, "% completed!"))
          # Pause for 0.1 seconds to simulate a long computation.
          Sys.sleep(0.0005)
        }
       })     
        #return(parsed$choices$message$content) ##  LLM output
      } else {
        (paste("Error: Unexpected response format -", content))
      }
    }
    list(result = final_result, treated_samples = treated_samples, wt_samples = wt_samples, which_results = regulated)
  }
  
  output$download <- downloadHandler(
    filename = function(){paste("results_of_the_analysis_of_",input$file2$name,sep="")},
    content = function(fname){write.table(get_response_table(input$text2, input$file2)$result, fname, col.names = TRUE,row.names = FALSE,quote=FALSE, sep = "\t")}
  )


  
  output$chat_output <- renderUI({
    inDataFile = input$file6
    prompt = input$text6
    if ((is.null(inDataFile)) || (prompt=="")){return(NULL)}else{
      print("answer")
      my_data <- read_and_clean_colnames(inDataFile$datapath, sep = "\t", header = TRUE, rownames = 1 )
      #   prompt = "Plot a heatmap"
      data_text <- paste(capture.output(head(my_data, 3)), collapse = "\n")  # Only show top 20 rows
      api_key <- "..."  # Replace with your API key
      url <- "https://api.together.xyz/v1/chat/completions"
      body <- list(
        #model = "google/gemma-2-27b-it",
        model = "meta-llama/Llama-3.3-70B-Instruct-Turbo-Free",
        messages = list(
          list(role = "system", content = "use a Chain of Thought and Act as an expert data analyst in R programming language"),
          list(role = "user", content = paste(prompt," on a dataset that have: ", dim(my_data)[1]," rows, and ",dim(my_data)[2]," columns. 
                                              Here are the first 3 rows:\n", data_text)
          )
        ),
        max_tokens = 1500
      )
      response <- POST(
        url,
        add_headers(
          Authorization = paste("Bearer", api_key),
          `Content-Type` = "application/json"
        ),
        body = toJSON(body, auto_unbox = TRUE),
        encode = "json"
      )
      content <- httr::content(response, as = "text")
      parsed <- fromJSON(content)
      HTML(markdown::markdownToHTML(parsed$choices$message[2]$content))
    }
  })
  
  output$chat_output_short_advice <- renderUI({
    inDataFile = input$file2
    prompt = input$text2
    if ((is.null(inDataFile)) || (prompt=="")){return(NULL)}else{
      my_data <- read_and_clean_colnames(inDataFile$datapath, sep = "\t", header = TRUE, rownames = 1 )
      #   prompt = "Can you perform a de analysis of treated 1,2,5 against wt 3,7,8 ?"
      data_text <- paste(capture.output(head(my_data, 3)), collapse = "\n")  # Only show top 20 rows
      api_key <- "..."  # Replace with your API key
      url <- "https://api.together.xyz/v1/chat/completions"
      body <- list(
        model = "google/gemma-2-27b-it",
        #model = "meta-llama/Llama-3.3-70B-Instruct-Turbo-Free",
        messages = list(
          list(role = "system", content = "Act as an expert data analyst suggesting a general way to perform the request. Do not write code."),
          list(role = "user", content = paste(prompt," on a dataset that have ", dim(my_data)[1]," rows, and ",dim(my_data)[2]," columns. 
                                              Here are the first 3 rows:\n", data_text)
          )
        ),
        max_tokens = 1700
      )
      response <- POST(
        url,
        add_headers(
          Authorization = paste("Bearer", api_key),
          `Content-Type` = "application/json"
        ),
        body = toJSON(body, auto_unbox = TRUE),
        encode = "json"
      )
      content <- httr::content(response, as = "text")
      parsed <- fromJSON(content)
      HTML(markdown::markdownToHTML(parsed$choices$message[2]$content))
    }
  })
  
  output$chat_output_interpretationGemma <- renderUI({
      my_data <- enrichment_function()$results
      print(my_data[,1:5])
      my_data= my_data$Description
      data_text = my_data
      data_text = paste((data_text), collapse = ", ")
      api_key <- "..."  # Replace with your API key
      url <- "https://api.together.xyz/v1/chat/completions"
      body <- list(
        model = "google/gemma-2-27b-it",
        messages = list(
          list(role = "system", content = "Interpret data by giving possible explanation of the biological processes involved"),
          list(role = "user", content = paste("Give a unified summary about these biological processes together with exaplanations: ", data_text)
          )
        ),
        max_tokens = 2500
      )
      response <- POST(
        url,
        add_headers(
          Authorization = paste("Bearer", api_key),
          `Content-Type` = "application/json"
        ),
        body = toJSON(body, auto_unbox = TRUE),
        encode = "json"
      )
      content <- httr::content(response, as = "text")
      parsed <- fromJSON(content)
      HTML(markdown::markdownToHTML(parsed$choices$message[2]$content))
  })
  
  output$chat_output_interpretationLlama <- renderUI({
      my_data <- enrichment_function()$results
      print(my_data[,1:5])
      my_data= my_data$Description
      data_text = my_data
      data_text = paste((data_text), collapse = ", ")
      api_key <- "..."  # Replace with your API key
      url <- "https://api.together.xyz/v1/chat/completions"
      body <- list(
        model = "meta-llama/Llama-3.3-70B-Instruct-Turbo-Free",
        messages = list(
          list(role = "system", content = "Interpret data by giving possible explanation of the biological processes involved"),
          list(role = "user", content = paste("Give a unified summary about these biological processes together with exaplanations: ", data_text)
          )
        ),
        max_tokens = 2500
      )
      response <- POST(
        url,
        add_headers(
          Authorization = paste("Bearer", api_key),
          `Content-Type` = "application/json"
        ),
        body = toJSON(body, auto_unbox = TRUE),
        encode = "json"
      )
      content <- httr::content(response, as = "text")
      parsed <- fromJSON(content)
      HTML(markdown::markdownToHTML(parsed$choices$message[2]$content))
  })
  
  output$show_input_fun <- DT::renderDT({  
    inFile <- input$file3
    if (is.null(inFile)){return(NULL)}else{
      # by clicking the 'Run' button the code below is executed!
      counts = read_and_clean_colnames(inFile$datapath, header=T, sep="\t", rownames=1)
      counts
    }
  })
  
  output$show_input_fun_plot <- DT::renderDT({  
    inFile <- input$file5
    if (is.null(inFile)){return(NULL)}else{
      # by clicking the 'Run' button the code below is executed!
      counts = read_and_clean_colnames(inFile$datapath, header=T, sep="\t", rownames=1)
      counts
    }
  })
  
  # Function for pathway enrichment analysis
  gene_pathway_enrichment <- function(gene_list) {
    
    if (all(grepl("^[A-Za-z0-9.-]+$", gene_list[1:10]))) { #we have gene symbols
    
      # Convert gene symbols to Entrez IDs
      entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_list, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
      entrez_ids <- na.omit(entrez_ids)
      #print(entrez_ids)
    
    if (length(entrez_ids) == 0) {
      stop("No valid Entrez IDs found for the given gene symbols.")
    }
    
    # Perform KEGG pathway enrichment analysis
    kegg_enrichment <- enrichKEGG(gene = entrez_ids,  organism = 'hsa', pAdjustMethod="BH",qvalueCutoff = 0.001)
    #kegg_enrichment <- enrichGO(gene = entrez_ids, OrgDb= org.Hs.eg.db, keyType= "ENTREZID", ont= "BP", pAdjustMethod = "BH", qvalueCutoff  = 0.00001)
    #kegg_enrichment <- enrichPathway(gene = entrez_ids, organism = "human", pAdjustMethod = "BH", qvalueCutoff  = 0.0001, readable = TRUE)
    
    if (is.null(kegg_enrichment)) {
      return("No significant pathways found.")
    }
    
    }else{
      # Perform KEGG pathway enrichment analysis
      kegg_enrichment <- enrichKEGG(gene = entrez_ids,  organism = 'hsa', pAdjustMethod="BH",qvalueCutoff = 0.001)
    }
      # Convert results to a data frame
      enrichmentDE <- as.data.frame(kegg_enrichment)
      enrichmentDE[,c(4,2,5,11,12,13,14)]
      #enrichmentDE
  }
  
  enrichment_function <- eventReactive(input$run3, {
    inFile3 = input$file3
    result_de <- read_and_clean_colnames(inFile3$datapath, header = TRUE, sep = "\t", rownames = 1)
    gene_list <- rownames(subset(result_de, as.numeric(result_de$FDR) < 0.05)) #dagli solo quelli significativi
    print(head(gene_list))
    final_result = gene_pathway_enrichment(as.character(gene_list))
    print(head(final_result))
    write.table(final_result, paste("enrichment_on_",input$file3$name ,sep=""), row.names = TRUE, sep="\t",quote=FALSE, col.names = TRUE)
    return(list(results = final_result, counts = result_de))
  })
  
  output$download_results3 <- downloadHandler(
    filename = function() { paste("enrichment_on_",input$file3$name ,sep="") },
    content = function(file) {
      file.copy(paste("enrichment_on_",input$file3$name ,sep=""), file)
    }
  )
  
  output$resultsTableEnrich <- renderDT({
    datatable(enrichment_function()$results)
  })
  
  output$inputTableEnrich <- renderDT({
    #datatable(enrichment_function()$counts)
    inDataFile = input$file3
    if (is.null(inDataFile)){
      counts = matrix(data=" ",nrow=2,ncol=2)
      colnames(counts) = c(" "," ")
      datatable(counts)
    }else{
      counts <- read.table(inDataFile$datapath, header = TRUE, sep = "\t", row.names = 1)
      datatable(counts)
    }
  })
  
  output$generate_dotplot <- renderPlotly({
   inFile3 = input$file3
   #counts <- read.table("C:/Users/f.russo.ENGIBBC/Desktop/REDAC_submission/results_expression_file.txt", header = TRUE, sep = "\t", row.names = 1)
   counts <- read_and_clean_colnames(inFile3$datapath, header = TRUE, sep = "\t", rownames = 1)
   print(head(counts))
   gene_symbol_list <- rownames(subset(counts, as.numeric(counts$FDR) < 0.05))
   entrez_ids <- mapIds(
    org.Hs.eg.db,
    keys = as.character(gene_symbol_list),
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
   )
   entrez_ids <- na.omit(entrez_ids)
   de <- as.character(entrez_ids)
   print(head(de))
   # KEGG enrichment
   #x = gene_pathway_enrichment(as.character(de))
   x <- enrichKEGG(gene = de,organism = 'hsa',pAdjustMethod = "BH",qvalueCutoff = 0.001)
  
  # Convert to data frame
  x_df <- as.data.frame(x)
  print(x_df[,c(1,2,3,4,5,6,7,8,9,10,11)])
  # Keep only top N results (optional)
  x_df <- x_df %>% arrange(qvalue) #%>% head(20)
  
  # Create Plotly dot plot
  fig = plot_ly(
    data = x_df,
    x = ~GeneRatio,
    y = ~reorder(Description, -qvalue),
    type = 'scatter',
    mode = 'markers',
    marker = list(
      size = ~Count * 0.7,
      color = ~qvalue,
      colorscale = 'Magma',
      colorbar = list(title = "qvalue"),
      reversescale = TRUE,
      showscale = TRUE
    ),
    text = ~paste(
      "Pathway:", Description,
      "<br>Gene Ratio:", GeneRatio,
      "<br>Count:", Count,
      "<br>qvalue:", signif(qvalue, 4)
    ),
    hoverinfo = "text"
  ) %>%
    layout(
      title = "KEGG Pathway Enrichment",
      xaxis = list(title = "Gene Ratio"),
      yaxis = list(title = ""),
      margin = list(l = 150)  # for long pathway names
    )
  fig
  })
  
  
  
  dea_results2 <- eventReactive(input$run2, {
    req(input$file2)
    req(input$text2)
    #print(input$file2)
    #print(input$text2)
    res = get_response_table(input$text2, input$file2)$result
    res_df <- as.data.frame(res) #%>% rownames_to_column("Gene")
    # Save locally
    write.table(res_df, paste("result_of_",input$text2,"_on_",input$file2$name,sep=""), row.names = TRUE, sep="\t",quote=FALSE, col.names = TRUE)
    counts <- read_and_clean_colnames(input$file2$datapath, header=T,rownames = 1, sep="\t")
    return(list(results = res_df, counts = counts, text=input$text2))
  })
  
  output$download_results2 <- downloadHandler(
    filename = function() { paste("result_of_",input$text2,"_on_",input$file2$name,sep="") },
    content = function(file) {
      file.copy(paste("result_of_",input$text2,"_on_",input$file2$name,sep=""), file)
    }
  )
  
  output$resultsTable2 <- renderDT({
    #datatable(dea_results2()$results)
    resulttable = dea_results2()$results
    resulttable = as.data.frame(resulttable)
    row.names(resulttable) <- sapply(row.names(resulttable), function(x) 
      toString(tags$a(href=paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",x,"&keywords=",x), x)))
    resulttable = datatable(as.data.frame(resulttable), escape=FALSE)
    resulttable
  })
  
  output$show_input_fun2 <- renderDT({
      inDataFile = input$file2
      if (is.null(inDataFile)){
        counts = matrix(data=" ",nrow=2,ncol=2)
        colnames(counts) = c(" "," ")
        datatable(counts)
      }else{
        counts <- read_and_clean_colnames(inDataFile$datapath, sep = "\t", header = TRUE, rownames = 1 )
        #counts <- read.table(inDataFile$datapath, header = TRUE, sep = "\t", row.names = 1)
        #datatable(dea_results2()$counts)
        datatable(counts)
      }
  })
  
  output$chat_output2 <- renderUI({
   # inDataFile = input$file2
  #  prompt = input$text2
   # if ((is.null(inDataFile)) || (prompt=="")){return(NULL)}else{
      print("answer")
    my_data <- dea_results2()$counts
    prompt <- dea_results2()$text
      #my_data <- read_and_clean_colnames(inDataFile$datapath, sep = "\t", header = TRUE, rownames = 1 )
      #   prompt = "Plot a heatmap"
      data_text <- paste(capture.output(colnames(my_data)), collapse = "\n")  # Only show top 20 rows
      api_key <- "..."  # Replace with your API key
      url <- "https://api.together.xyz/v1/chat/completions"
      body <- list(
        #model = "google/gemma-2-27b-it",
        model = "meta-llama/Llama-3.3-70B-Instruct-Turbo-Free",
        messages = list(
          list(role = "system", content = "use a Chain of Thought and Act as an expert data analyst in R programming language that uses DESEQ2 package"),
          list(role = "user", content = paste(prompt," on a dataset that have: ", dim(my_data)[1]," rows, and ",dim(my_data)[2]," columns. 
                                              Here are the sample names:\n", data_text, " use DESEQ2 package for this request",prompt)
          )
        ),
        max_tokens = 1500
      )
      response <- POST(
        url,
        add_headers(
          Authorization = paste("Bearer", api_key),
          `Content-Type` = "application/json"
        ),
        body = toJSON(body, auto_unbox = TRUE),
        encode = "json"
      )
      content <- httr::content(response, as = "text")
      parsed <- fromJSON(content)
      HTML(markdown::markdownToHTML(parsed$choices$message[2]$content))
    #}
  })
  output$chat_output_short_advice2 <- renderUI({
   # inDataFile = input$file2
  #  prompt = input$text2
    #if ((is.null(inDataFile)) || (prompt=="")){return(NULL)}else{
    #my_data <- read_and_clean_colnames(inDataFile$datapath, sep = "\t", header = TRUE, rownames = 1 )
    my_data <- dea_results2()$counts
    prompt <- dea_results2()$text
      #   prompt = "Can you perform a de analysis of treated 1,2,5 against wt 3,7,8 ?"
      data_text <- paste(capture.output(head(my_data, 3)), collapse = "\n")  # Only show top 20 rows
      api_key <- "..."  # Replace with your API key
      url <- "https://api.together.xyz/v1/chat/completions"
      body <- list(
        model = "google/gemma-2-27b-it",
        #model = "meta-llama/Llama-3.3-70B-Instruct-Turbo-Free",
        messages = list(
          list(role = "system", content = "Act as an expert rnaseq data analyst suggesting a general way to perform the request (suggesting the use of DESEQ2, EdgeR, limma, for de genes identification) 
                    and the generation of several plots such as PCA, densities, heatmap, dendrogram, violin, et cetera to explore input data and MAplot, volcano plot, et cetera to explore results. Do not write code."),
          list(role = "user", content = paste(prompt," on a dataset that have ", dim(my_data)[1]," rows, and ",dim(my_data)[2]," columns. 
                                              Here are the first 3 rows:\n", data_text)
          )
        ),
        max_tokens = 1700
      )
      response <- POST(
        url,
        add_headers(
          Authorization = paste("Bearer", api_key),
          `Content-Type` = "application/json"
        ),
        body = toJSON(body, auto_unbox = TRUE),
        encode = "json"
      )
      content <- httr::content(response, as = "text")
      parsed <- fromJSON(content)
      HTML(markdown::markdownToHTML(parsed$choices$message[2]$content))
    #}
  })
  output$volcanoPlot2 <- renderPlotly({
    res_df <- dea_results2()$results
    res_df$significant <- res_df$FDR < 0.05
    res_df$gene <- rownames(res_df)
    volcano_plot <- plot_ly(
      data = res_df,
      x = ~log2FoldChange,
      y = ~-log10(FDR),
      type = 'scatter',
      mode = 'markers',
      color = ~significant,
      colors = c("gray", "red"),  # or customize your color scheme
      text = ~paste(
        "Gene: ", gene,
        "<br>log2FC: ", round(log2FoldChange, 2),
        "<br>-log10(p): ", round(-log10(FDR), 2)
      ),
      hoverinfo = 'text',
      marker = list(size = 6, opacity = 0.6)
    ) %>%
      layout(
        title = list(text = "Volcano Plot", font = list(size = 16)),
        xaxis = list(title = "log2 Fold Change"),
        yaxis = list(title = "-log10 FDR"),
        legend = list(title = list(text = "Significance"))
      )
    
    volcano_plot
  })
  
  output$foldchangePlot2 <- renderPlotly({
    res_df <- dea_results2()$results
    res_df$significant <- res_df$FDR < 0.05
    #plot(log2(res$baseMean + 1) , res$log2FoldChange ,col = "black", main="DESeq2 Fold Change Plot", xlab='Mean of Normalized Counts', ylab='log2FoldChange',pch=19,cex=0.3) 
    res_df$gene <- rownames(res_df)
    print(head(res_df))

    # Step 3: Create interactive MA plot
    ma_plot <- plot_ly(
      data = res_df,
      x = ~logCPM,
      y = ~log2FoldChange,
      type = 'scatter',
      mode = 'markers',
      color = ~significant,
      colors = c("gray", "red"),  # Customize colors if needed
      text = ~paste(
        "Gene: ", gene,
        "<br>logCPM: ", round(logCPM, 2),
        "<br>log2FC: ", round(log2FoldChange, 2)
      ),
      hoverinfo = 'text',
      marker = list(size = 6, opacity = 0.6)
    ) %>%
      layout(
        title = list(text =  "Fold Change vs logCPM Plot", font = list(size = 16)),
        xaxis = list(title = "logCPM"),
        yaxis = list(title = "log2 Fold Change"),
        legend = list(title = list(text = "Significance"))
      )
    
    ma_plot
  })
  
  output$pcaPlot2 <- renderPlotly({
    #datatable(dea_results2()$counts)
    #x <- dea_results2()$counts
    inDataFile = input$file2
    if (is.null(inDataFile)){
      # Create an empty plot
      fig <- plot_ly() %>%
        layout(
          xaxis = list(showticklabels = FALSE, zeroline = FALSE),
          yaxis = list(showticklabels = FALSE, zeroline = FALSE),
          annotations = list(
            text = "",
            xref = "", yref = "",
            showarrow = FALSE,
            font = list(size = 20),
            x = 0.5, y = 0.5
          )
        )
      
      fig
    }else{
      #x <- dea_results2()$counts
    x =  read_and_clean_colnames(inDataFile$datapath, sep = "\t", header = TRUE, rownames = 1 )
    
    # Step 1: Remove genes (rows) with zero variance
    x <- x[apply(x, 1, var) != 0, ]
    
    # Step 2: Transpose and run PCA
    x_t <- t(x)
    pca_result <- prcomp(x_t, scale. = TRUE)
    
    # Step 3: Extract PCA scores
    pca_data <- as.data.frame(pca_result$x)
    pca_data$Sample <- rownames(pca_data)
    
    # Step 4: Get variance explained for labels
    variance_explained <- summary(pca_result)$importance[2, 1:2] * 100
    
    # Step 5: Plotly PCA plot
    fig <- plot_ly(
      data = pca_data,
      x = ~PC1,
      y = ~PC2,
      type = 'scatter',
      mode = 'markers+text',
      text = ~Sample,
      hoverinfo = 'text',
      textposition = 'top right',
      marker = list(size = 6, color = 'blue')
    ) %>%
      layout(
        title = "PCA: PC1 vs PC2",
        xaxis = list(title = paste0("PC1 (", round(variance_explained[1], 1), "%)")),
        yaxis = list(title = paste0("PC2 (", round(variance_explained[2], 1), "%)")),
        showlegend = FALSE
      )
    fig
    }
  })
  
  output$integration_network <- renderPlotly({
    datatable(dea_results2()$counts)
    # counts <- read.table("C:/Users/f.russo.ENGIBBC/Desktop/REDAC_submission/expression_file.txt", header = TRUE, sep = "\t", row.names = 1)
    # ppi_data  <- read.table("C:/Users/f.russo.ENGIBBC/Desktop/REDAC_submission/protein_protein_interaction_data.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    ppi_data  <- read.table(file ="9606.protein.physical.links.detailed.v12.0_with_gene_symbols_.txt", header = TRUE,sep = "\t")
    protein_info <- read.table(file ="9606.protein.info.v12.0.txt", header = TRUE,sep = "\t")
    colnames(ppi_data) <- c("protein_code_2", "protein_code_1", "experimental", "Protein1", "Protein2")
    # colnames(ppi_data) <- c("Protein1", "Protein2")
    #Adesso elimino i conteggi nulli o quelli inferiori a 10 su tutti i campioni
    #keep <- rowSums(counts > dim(counts)[2])
    #counts <- counts[keep,]
    gene_expr = counts
    interaction_threshold = 980
    ppi_data = subset(ppi_data, ppi_data$experimental > interaction_threshold)
    print(dim(ppi_data))
    # Ensure genes in PPI match those in the gene expression dataset
    valid_genes <- rownames(gene_expr)
    ppi_filtered <- ppi_data %>% filter(Protein1 %in% valid_genes & Protein2 %in% valid_genes)
    
    # Function to compute Pearson correlation between two interacting proteins
    compute_correlation2 <- function(protein1, protein2, gene_expr) {
      expr1 <- gene_expr[protein1, ]
      expr2 <- gene_expr[protein2, ]
      
      if (sd(expr1) == 0 || sd(expr2) == 0) {
        return(NA)  # Avoid zero variance issue
      }
      
      cor_value <- cor(as.numeric(expr1), as.numeric(expr2), method = "pearson", use = "complete.obs")
      cor_value 
    }
    
    ppi_filtered <- ppi_filtered[,3:5]
    ppi_filtered <- ppi_filtered[,c(2,3,1)]
    
    # Compute correlation for each interaction
    ppi_filtered <- ppi_filtered %>%
      rowwise() %>%
      mutate(Correlation = compute_correlation2(Protein1,Protein2, gene_expr) ) %>%
      drop_na()  # Remove NA values due to zero variance
    
    # Create a weighted graph
    ppi_graph <- graph_from_data_frame(ppi_filtered, directed = FALSE)
    
    # Generate layout
    layout <- layout_with_fr(ppi_graph)
    V(ppi_graph)$x <- layout[,1]
    V(ppi_graph)$y <- layout[,2]
    
    # Get edge and node data
    edges <- get.data.frame(ppi_graph)
    nodes <- data.frame(id = V(ppi_graph)$name,
                        x = V(ppi_graph)$x,
                        y = V(ppi_graph)$y)
    
    # Start plot
    fig <- plot_ly()
    
    # Add edge traces
    for (i in 1:nrow(edges)) {
      edge <- edges[i, ]
      v0 <- which(V(ppi_graph)$name == edge$from)
      v1 <- which(V(ppi_graph)$name == edge$to)
      
      fig <- fig %>% add_trace(
        type = "scatter",
        mode = "lines",  
        x = c(V(ppi_graph)$x[v0], V(ppi_graph)$x[v1]),
        y = c(V(ppi_graph)$y[v0], V(ppi_graph)$y[v1]),
        line = list(
          width = abs(edge$Correlation) * 2,
          color = ifelse(edge$Correlation > 0, "red", "blue"),
          opacity = 0.5
        ),
        hoverinfo = "text",
        text = paste0(edge$from, " - ", edge$to, "<br>Correlation: ", round(edge$Correlation, 2)),
        showlegend = TRUE
      )
    }
    
    # Add node trace
    fig <- fig %>% add_trace(
      type = "scatter",
      mode = "markers+text",
      x = nodes$x,
      y = nodes$y,
      text = nodes$id,
      textposition = "top center",
      marker = list(size = 5, color = 'skyblue'),
      hoverinfo = "text",
      showlegend = TRUE
    )
    
    # Final layout
    fig <- fig %>% layout(
      title = 'Gene Expression Weighted PPI Network',
      xaxis = list(title = "", showgrid = FALSE, zeroline = FALSE),
      yaxis = list(title = "", showgrid = FALSE, zeroline = FALSE)
    )
    
    fig
  })
  
  output$pcaCompPlot2 <- renderPlotly({
    #datatable(dea_results2()$counts)
    #counts <- dea_results2()$counts
    inDataFile = input$file2
    if (is.null(inDataFile)){
      # Create an empty plot
      fig <- plot_ly() %>%
        layout(
          xaxis = list(showticklabels = FALSE, zeroline = FALSE),
          yaxis = list(showticklabels = FALSE, zeroline = FALSE),
          annotations = list(
            text = "",
            xref = "", yref = "",
            showarrow = FALSE,
            font = list(size = 20),
            x = 0.5, y = 0.5
          )
        )
      
      fig
    }else{
    #x <- dea_results2()$counts
    counts =  read_and_clean_colnames(inDataFile$datapath, sep = "\t", header = TRUE, rownames = 1 )
    
    # Step 1: Remove genes (rows) with zero variance
    counts <-  counts[apply( counts, 1, var) != 0, ]
    pca_result <- prcomp(t(counts), scale. = TRUE)  # Transpose if samples are columns
    
    # Calculate variance explained
    var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
    
    # Create Plotly screeplot
    fig <- plot_ly(
      x = paste0("PC", 1:length(var_explained)),
      y = var_explained,
      type = "bar",
      marker = list(color = 'blue'),
      text = ~paste0(round(var_explained, 2), "%"),
      hoverinfo = 'text'
    ) %>%
      layout(
        title = "Principal Component Histogram",
        xaxis = list(title = "Principal Components"),
        yaxis = list(title = "Variance Explained (%)")
      )
    fig
    }
  })
  
  output$violinPlot2 <- renderPlotly({
    #datatable(dea_results2()$counts)
    #counts <- dea_results2()$counts
    inDataFile = input$file2
    if (is.null(inDataFile)){
      # Create an empty plot
      fig <- plot_ly() %>%
        layout(
          xaxis = list(showticklabels = FALSE, zeroline = FALSE),
          yaxis = list(showticklabels = FALSE, zeroline = FALSE),
          annotations = list(
            text = "",
            xref = "", yref = "",
            showarrow = FALSE,
            font = list(size = 20),
            x = 0.5, y = 0.5
          )
        )
      
      fig
    }else{
      #x <- dea_results2()$counts
      counts =  read_and_clean_colnames(inDataFile$datapath, sep = "\t", header = TRUE, rownames = 1 )
    # Step 1: Log-transform counts
    log_counts <- log10(counts + 1)
    # Step 2: Add gene names as a column before melting
    log_counts$Gene <- rownames(log_counts)
    # Step 3: Melt into long format
    log_counts_long <- melt(log_counts, id.vars = "Gene", variable.name = "Sample", value.name = "Expression")
      violin_plot <- plot_ly(      
        data = log_counts_long,
        x = ~Sample,
        y = ~Expression,
        type = 'violin',
        split = ~Sample,
        box = list(
          visible = T
        ),
        meanline = list(
          visible = T
        )
      )%>%
      layout(
        title = "Violin Plot of log10(counts + 1)",
        xaxis = list(title = "Sample", tickangle = 90),
        yaxis = list(title = "log10(Expression + 1)")
      )
    violin_plot
    }
  })
  
  output$densityPlot2 <- renderPlotly({
    #datatable(dea_results2()$counts)
    #counts <- dea_results2()$counts
    inDataFile = input$file2
    if (is.null(inDataFile)){
      # Create an empty plot
      fig <- plot_ly() %>%
        layout(
          xaxis = list(showticklabels = FALSE, zeroline = FALSE),
          yaxis = list(showticklabels = FALSE, zeroline = FALSE),
          annotations = list(
            text = "",
            xref = "", yref = "",
            showarrow = FALSE,
            font = list(size = 20),
            x = 0.5, y = 0.5
          )
        )
      
      fig
    }else{
      #x <- dea_results2()$counts
      counts =  read_and_clean_colnames(inDataFile$datapath, sep = "\t", header = TRUE, rownames = 1 )
    print("density")
    # Prepare the dataframe
    df_list <- lapply(seq_len(ncol(counts)), function(i) {
      data.frame(
        counts = counts[, i],
        idsample = colnames(counts)[i]
      )
    })
    
    df <- bind_rows(df_list)
    
    # Log-transform counts
    df$log_counts <- log(df$counts + 1)
    
    # Initialize plotly object
    p <- plot_ly()
    
    # Add a density curve per sample
    for (sample_name in unique(df$idsample)) {
      sample_data <- df %>% filter(idsample == sample_name)
      d <- density(sample_data$log_counts)
      
      p <- p %>%
        add_trace(
          x = d$x,
          y = d$y,
          type = "scatter",
          mode = "lines",
          name = sample_name
        )
    }
    
    # Add layout
    p <- p %>%
      layout(
        title = "Log Gene Count Densities by Sample",
        xaxis = list(title = "log(count + 1)"),
        yaxis = list(title = "Density"),
        legend = list(title = list(text = "Sample"))
      )
    p
    }
  })
  
  output$dendroPlot2 <- renderPlotly({
    #datatable(dea_results2()$counts)
    inDataFile = input$file2
    if (is.null(inDataFile)){
      # Create an empty plot
      fig <- plot_ly() %>%
        layout(
          xaxis = list(showticklabels = FALSE, zeroline = FALSE),
          yaxis = list(showticklabels = FALSE, zeroline = FALSE),
          annotations = list(
            text = "",
            xref = "", yref = "",
            showarrow = FALSE,
            font = list(size = 20),
            x = 0.5, y = 0.5
          )
        )
      
      fig
    }else{
    #x <- dea_results2()$counts
    x =  read_and_clean_colnames(inDataFile$datapath, sep = "\t", header = TRUE, rownames = 1 )
    x = log10(x+1)
    counts <- x
    print("Clustering!")
    #x <- normalize.quantiles(as.matrix(x), copy = TRUE)
    colnames(x) <- colnames(counts)
    rownames(x) <- rownames(counts)
    
    clustering_used_method <- "ward.D"
    distance_used_method <- "euclidean"
    
    hclustfunc <- function(x) hclust(x, method = clustering_used_method)
    distfunc <- function(x) dist(x, method = distance_used_method)
    
    fit <- hclustfunc(distfunc(t(x)))
    
    # Convert hclust to dendrogram and then to dendrogram segments
    dend <- as.dendrogram(fit)
    p <- ggdendrogram(dend , rotate = FALSE, size = 2)
    fig=ggplotly(p)
    fig
    }
  })
  
  output$heatmapPlot2 <- renderPlotly({
    #datatable(dea_results2()$counts)
    #counts <- dea_results2()$counts
    inDataFile = input$file2
    if (is.null(inDataFile)){
      # Create an empty plot
      fig <- plot_ly() %>%
        layout(
          xaxis = list(showticklabels = FALSE, zeroline = FALSE),
          yaxis = list(showticklabels = FALSE, zeroline = FALSE),
          annotations = list(
            text = "",
            xref = "", yref = "",
            showarrow = FALSE,
            font = list(size = 20),
            x = 0.5, y = 0.5
          )
        )
      
      fig
    }else{
      #x <- dea_results2()$counts
      counts =  read_and_clean_colnames(inDataFile$datapath, sep = "\t", header = TRUE, rownames = 1 )
      if(dim(counts)[1]<1000){
        # Step 1: Compute log-transformed expression matrix
        tog_counts <- log(as.matrix(counts + 1))
        # Step 3: Scale by row (z-score normalization per gene)
        scale_rows <- function(x) {
          t(scale(t(x)))
        }
        scaled_counts <- scale_rows(top_counts)
        
        # Step 4: Create Plotly heatmap
        heatmap_plot <- plot_ly(
          y = colnames(scaled_counts),
          x = rownames(scaled_counts),
          z = t(scaled_counts),
          type = "heatmap",
          colors = colorRampPalette(c("red", "white", "blue"))(500),
          showscale = TRUE
        ) %>%
          layout(
            title = list(text = "Heatmap of the Most Variable 1000 Genes", font = list(size = 12)),
            xaxis = list(title = "", tickangle = 90 ,tickfont = list(size = 10)),
            yaxis = list(title = "", tickfont = list(size = 10)),
            margin = list(l = 100, r = 50, b = 100, t = 50)
          )
        
        fig = heatmap_plot
      }else{
        # Step 1: Compute log-transformed expression matrix
        log_counts <- log(as.matrix(counts + 1))
        
        # Step 2: Select top 50 most variable genes
        gene_variances <- apply(log_counts, 1, var)
        top_genes <- names(sort(gene_variances, decreasing = TRUE))[1:1000]
        top_counts <- log_counts[top_genes, ]
        
        # Step 3: Scale by row (z-score normalization per gene)
        scale_rows <- function(x) {
          t(scale(t(x)))
        }
        scaled_counts <- scale_rows(top_counts)
        
        # Step 4: Create Plotly heatmap
        heatmap_plot <- plot_ly(
          y = colnames(scaled_counts),
          x = rownames(scaled_counts),
          z = t(scaled_counts),
          type = "heatmap",
          colors = colorRampPalette(c("red", "white", "blue"))(500),
          showscale = TRUE
        ) %>%
          layout(
            title = list(text = "Heatmap of the Most Variable 1000 Genes", font = list(size = 12)),
            xaxis = list(title = "", tickangle = 90 ,tickfont = list(size = 10)),
            yaxis = list(title = "", tickfont = list(size = 10)),
            margin = list(l = 100, r = 50, b = 100, t = 50)
          )
        
        fig = heatmap_plot
      } 
     fig
    }
  })
  
  output$pca3DPlot2 <- renderPlotly({
   # datatable(dea_results2()$counts)
   # counts <- dea_results2()$counts
    inDataFile = input$file2
    if (is.null(inDataFile)){
      # Create an empty plot
      p <- plot_ly() %>%
        layout(
          xaxis = list(showticklabels = FALSE, zeroline = FALSE),
          yaxis = list(showticklabels = FALSE, zeroline = FALSE),
          annotations = list(
            text = "",
            xref = "", yref = "",
            showarrow = FALSE,
            font = list(size = 20),
            x = 0.5, y = 0.5
          )
        )
      
      p
    }else{
      #x <- dea_results2()$counts
      counts =  read_and_clean_colnames(inDataFile$datapath, sep = "\t", header = TRUE, rownames = 1 )
    # Perform PCA (transpose if samples are columns!)
    pca_result <- PCA(t(counts), graph = FALSE)  # Note the transpose 't()' if samples are columns
    # Extract the PCA scores (rows are now samples after transposing)
    pca_data <- data.frame(
      Sample = colnames(counts),
      PC1 = pca_result$ind$coord[, 1],
      PC2 = pca_result$ind$coord[, 2],
      PC3 = pca_result$ind$coord[, 3]
    )
    # Create interactive 3D plot with sample labels as hover text
    p <- plot_ly(
      data = pca_data,
      x = ~PC1, y = ~PC2, z = ~PC3,
      type = 'scatter3d',
      mode = 'markers+text',
      text = ~Sample,
      hoverinfo = 'text',
      marker = list(size = 5, color = 'blue')
    )
    p <- p %>% layout(
      title = "3D PCA Plot",
      scene = list(
        xaxis = list(title = "PC1"),
        yaxis = list(title = "PC2"),
        zaxis = list(title = "PC3")
      )
    )
    p
    }
  })
  
  
  
  
  
  
  get_response_table3 <- function(prompt,inDataFile) {
    my_data <- read_and_clean_colnames(inDataFile$datapath, sep = "\t", header = TRUE, rownames = 1 )
    my_data
  }

  dea_results3 <- eventReactive(input$run6, {
    req(input$file6)
    counts <- read_and_clean_colnames(input$file6$datapath, header=T,rownames = 1)
    return(list(counts = counts))
  })
  
  output$show_input_fun3 <- renderDT({
    datatable(dea_results3()$counts)
  })
  
  output$generate_plot3 <- renderPlotly({
      req(input$text6)
      prompt <-input$text6
      inDataFile = input$file6
      #counts <- datatable(dea_results3()$counts)
      if (is.null(inDataFile)){
        # Create an empty plot
        fig <- plot_ly() %>%
          layout(
            xaxis = list(showticklabels = FALSE, zeroline = FALSE),
            yaxis = list(showticklabels = FALSE, zeroline = FALSE),
            annotations = list(
              text = "",
              xref = "", yref = "",
              showarrow = FALSE,
              font = list(size = 20),
              x = 0.5, y = 0.5
            )
          )
        
        fig
      }else{
        my_data <- read_and_clean_colnames(inDataFile$datapath, sep = "\t", header = TRUE, rownames = 1 )
        counts = my_data
      print(head(counts))
      api_key <- "..."  # Replace with your API key
      url <- "https://api.together.xyz/v1/chat/completions"
      body <- list(
        model = "google/gemma-2-27b-it",
        #model = "meta-llama/Llama-3.3-70B-Instruct-Turbo-Free",
        messages = list(
          list(role = "system", content = "Respond only in JSON format. The JSON must have this structure: 
               {\"functiontoberun\":[\"box\",\"surface\",\"density\",\"violin\",\"heat\",\"corrheat\",\"volcano\",\"pca\",\"cluster\",\"dendro\",\"integration\",\"pcacomp\",\"MAplot\",\"dot\",\"KEGG\",\"net\",\"3Dpca\"],
               \"columns\":[[\"col1\"],[\"col2\"],[\"col3\"]]}. Do not add explanations or other text."),
          list(role = "user", content = prompt)
        ),
        max_tokens = 500
      )
      response <- POST(
        url,
        add_headers(
          Authorization = paste("Bearer", api_key),
          `Content-Type` = "application/json"
        ),
        body = toJSON(body, auto_unbox = TRUE),
        encode = "json"
      )
      content <- httr::content(response, as = "text")
      #print(content)
      parsed <- fromJSON(content)
      #parsed <- tryCatch(fromJSON(content), error = function(e) NULL)
      #print(parsed$choices$message$content)
      # the output is in this format:
      # parsed$choices$message$content = "```json\n{\"function to be run\": [\"violin\"], \"columns to use\": [[\"1\", \"2\", \"3\"]]}\n```"
      # but it should be in this format:
      # parsed$choices$message$content = '{\"function to be run\": [\"violin\"], \"columns to use\": [[\"col1\"], [\"col2\"], [\"col3\"]]}'
      content = strsplit(parsed$choices$message$content,"\n" )[[1]][2] ## questa linea solo se usi Gemma LLM
      #print(content)
      #[1] "{\"function to be run\": [\"violin\"], \"columns\": [[\"1\"], [\"2\"], [\"3\"]]}"
      parsed <- fromJSON(content)
      print(parsed)
      #$`function to be run`
      #[1] "violin"
      if (is.null(parsed) || !is.list(parsed)) {
        ("Error: Failed to parse JSON response")
      }
      
      withProgress(message = 'Running, Please Wait!', detail = 'This may take a while...', value = 0, {     
        
      if (!is.null(parsed$functiontoberun) && !is.null(content)) {
        # Eseguire le funzioni sui dati specificati
          func_name <- parsed$functiontoberun[1]
          cols <- unlist(parsed$columnstouse)
          #scegli qui le funzioni da usare
          if      (func_name == "pca")      {
          # Step 1: Remove genes (rows) with zero variance
          x <- counts
          x <- x[apply(x, 1, var) != 0, ]
          
          # Step 2: Transpose and run PCA
          x_t <- t(x)
          pca_result <- prcomp(x_t, scale. = TRUE)
          
          # Step 3: Extract PCA scores
          pca_data <- as.data.frame(pca_result$x)
          pca_data$Sample <- rownames(pca_data)
          
          # Step 4: Get variance explained for labels
          variance_explained <- summary(pca_result)$importance[2, 1:2] * 100
          
          # Step 5: Plotly PCA plot
          library(plotly)
          
          pca_plot <- plot_ly(
            data = pca_data,
            x = ~PC1,
            y = ~PC2,
            type = 'scatter',
            mode = 'markers+text',
            text = ~Sample,
            hoverinfo = 'text',
            textposition = 'top right',
            marker = list(size = 6, color = 'blue')
          ) %>%
            layout(
              title = "PCA: PC1 vs PC2",
              xaxis = list(title = paste0("PC1 (", round(variance_explained[1], 1), "%)")),
              yaxis = list(title = paste0("PC2 (", round(variance_explained[2], 1), "%)")),
              showlegend = FALSE
            )
          print("plot")
          fig=pca_plot
          } 
          else if (func_name == "3Dpca")   {
            # Perform PCA (transpose if samples are columns!)
            # Perform PCA (transpose if samples are columns!)
            pca_result <- PCA(t(counts), graph = FALSE)  # Note the transpose 't()' if samples are columns
            
            # Extract the PCA scores (rows are now samples after transposing)
            pca_data <- data.frame(
              Sample = colnames(counts),
              PC1 = pca_result$ind$coord[, 1],
              PC2 = pca_result$ind$coord[, 2],
              PC3 = pca_result$ind$coord[, 3]
            )
            
            # Create interactive 3D plot with sample labels as hover text
            p <- plot_ly(
              data = pca_data,
              x = ~PC1, y = ~PC2, z = ~PC3,
              type = 'scatter3d',
              mode = 'markers+text',
              text = ~Sample,
              hoverinfo = 'text',
              marker = list(size = 5, color = 'blue')
            )
            
            p <- p %>% layout(
              title = "3D PCA Plot",
              scene = list(
                xaxis = list(title = "PC1"),
                yaxis = list(title = "PC2"),
                zaxis = list(title = "PC3")
              )
            )
            print("PCA 3dim")
            fig = p
          }
          else if (func_name == "heat")     {
            if(dim(counts)[1]<1000){
              # Step 1: Compute log-transformed expression matrix
              log_counts <- log(as.matrix(counts + 1))
              
              # Step 2: Select top 50 most variable genes
              #gene_variances <- apply(log_counts, 1, var)
              #top_genes <- names(sort(gene_variances, decreasing = TRUE))[1:50]
              top_counts <- log_counts
              
              # Step 3: Scale by row (z-score normalization per gene)
              scale_rows <- function(x) {
                t(scale(t(x)))
              }
              scaled_counts <- scale_rows(top_counts)
              
              # Step 4: Create Plotly heatmap
              heatmap_plot <- plot_ly(
                y = colnames(scaled_counts),
                x = rownames(scaled_counts),
                z = t(scaled_counts),
                type = "heatmap",
                colors = colorRampPalette(c("red", "white", "blue"))(500),
                showscale = TRUE
              ) %>%
                layout(
                  title = list(text = "Heatmap of the Most Variable 1000 Genes", font = list(size = 12)),
                  xaxis = list(title = "", tickangle = 90 ,tickfont = list(size = 10)),
                  yaxis = list(title = "", tickfont = list(size = 10)),
                  margin = list(l = 100, r = 50, b = 100, t = 50)
                )
              
              fig = heatmap_plot
            }else{
              # Step 1: Compute log-transformed expression matrix
              log_counts <- log(as.matrix(counts + 1))
              
              # Step 2: Select top 50 most variable genes
              gene_variances <- apply(log_counts, 1, var)
              top_genes <- names(sort(gene_variances, decreasing = TRUE))[1:1000]
              top_counts <- log_counts[top_genes, ]
              
              # Step 3: Scale by row (z-score normalization per gene)
              scale_rows <- function(x) {
                t(scale(t(x)))
              }
              scaled_counts <- scale_rows(top_counts)
              
              # Step 4: Create Plotly heatmap
              heatmap_plot <- plot_ly(
                y = colnames(scaled_counts),
                x = rownames(scaled_counts),
                z = t(scaled_counts),
                type = "heatmap",
                colors = colorRampPalette(c("red", "white", "blue"))(500),
                showscale = TRUE
              ) %>%
                layout(
                  title = list(text = "Heatmap of the Most Variable 1000 Genes", font = list(size = 12)),
                  xaxis = list(title = "", tickangle = 90 ,tickfont = list(size = 10)),
                  yaxis = list(title = "", tickfont = list(size = 10)),
                  margin = list(l = 100, r = 50, b = 100, t = 50)
                )
              
            fig = heatmap_plot
            } 
          }
          else if (func_name == "box")      {
            # Step 1: Log-transform counts
            log_counts <- log10(counts + 1)
            # Step 2: Add gene names as a column before melting
            log_counts$Gene <- rownames(log_counts)
            # Step 3: Melt into long format
            log_counts_long <- melt(log_counts, id.vars = "Gene", variable.name = "Sample", value.name = "Expression")
            box_plot <- plot_ly(      
              data = log_counts_long,
              x = ~Sample,
              y = ~Expression,
              type = 'box',
              split = ~Sample
            )%>%
              layout(
                title = "BoxPlot of log10(counts + 1)",
                xaxis = list(title = "Sample", tickangle = 90),
                yaxis = list(title = "log10(Expression + 1)")
              )
            fig = box_plot
          } 
          else if (func_name == "surface")  {
            if ( dim(counts)[1]>100 ){
            log_counts <- log10(counts + 1)
            gene_variances <- apply(log_counts, 1, var)
            top_genes <- names(sort(gene_variances, decreasing = TRUE))[1:100]
            gene_expression <- as.matrix(log_counts[top_genes, ])
            
            # Create surface plot
            fig = plot_ly(
              z = ~gene_expression,
              x = colnames(gene_expression),
              y = rownames(gene_expression),
              type = "surface"
            ) %>%
              layout(
                scene = list(
                  xaxis = list(title = ""),
                  yaxis = list(title = ""),
                  zaxis = list(title = "Log10 Level")
                ),
                title = "3D Surface of the 100 most variable genes"
              )
            }else{
                log_counts <- log10(counts + 1)
                gene_variances <- apply(log_counts, 1, var)
                top_genes <- names(sort(gene_variances, decreasing = TRUE))
                gene_expression <- as.matrix(log_counts[top_genes, ])
                # Create surface plot
                fig = plot_ly(
                  z = ~gene_expression,
                  x = colnames(gene_expression),
                  y = rownames(gene_expression),
                  type = "surface"
                ) %>%
                  layout(
                    scene = list(
                      xaxis = list(title = ""),
                      yaxis = list(title = ""),
                      zaxis = list(title = "Log10 Level")
                    ),
                    title = "3D Surface"
                  )
              }
            }
          else if (func_name == "violin")   {
            # Step 1: Log-transform counts
            log_counts <- log10(counts + 1)
            # Step 2: Add gene names as a column before melting
            log_counts$Gene <- rownames(log_counts)
            # Step 3: Melt into long format
            log_counts_long <- melt(log_counts, id.vars = "Gene", variable.name = "Sample", value.name = "Expression")
            violin_plot <- plot_ly(      
              data = log_counts_long,
              x = ~Sample,
              y = ~Expression,
              type = 'violin',
              split = ~Sample,
              box = list(
                visible = T
              ),
              meanline = list(
                visible = T
              )
            )%>%
              layout(
                title = "Violin Plot of log10(counts + 1)",
                xaxis = list(title = "Sample", tickangle = 90),
                yaxis = list(title = "log10(Expression + 1)")
              )
             fig = violin_plot
          } 
          else if (func_name == "corrheat"){
            top_counts <- cor(counts)
            # Step 3: Scale by row (z-score normalization per gene)
            scaled_counts <- top_counts
            # Step 4: Create Plotly heatmap
            heatmap_plot <- plot_ly(
              x = colnames(scaled_counts),
              y = rownames(scaled_counts),
              z = scaled_counts,
              type = "heatmap",
              colors = colorRampPalette(c("white","blue"))(500),
              showscale = TRUE,
            ) %>%
              layout(
                title = list(text = "Correlation Heatmap", font = list(size = 12)),
                xaxis = list(title = "", tickangle = 45),
                yaxis = list(title = "", tickfont = list(size = 8)),
                margin = list(l = 100, r = 50, b = 100, t = 50)
              )
            fig = heatmap_plot
          }
          else if (func_name == "volcano")  {
            res_df = counts
            print(res_df)
            res_df$significant <- res_df$FDR < 0.05
            res_df$gene <- rownames(res_df)
            volcano_plot <- plot_ly(
              data = res_df,
              x = ~log2FoldChange,
              y = ~-log10(FDR),
              type = 'scatter',
              mode = 'markers',
              color = ~significant,
              colors = c("gray", "red"),  # or customize your color scheme
              text = ~paste(
                "Gene: ", gene,
                "<br>log2FC: ", round(log2FoldChange, 2),
                "<br>-log10(p): ", round(-log10(FDR), 2)
              ),
              hoverinfo = 'text',
              marker = list(size = 6, opacity = 0.6)
            ) %>%
              layout(
                title = list(text = "Volcano Plot", font = list(size = 16)),
                xaxis = list(title = "log2 Fold Change"),
                yaxis = list(title = "-log10 FDR"),
                legend = list(title = list(text = "Significance"))
              )
            fig = volcano_plot
          }
          else if ((func_name == "cluster") || (func_name == "dendro")){
            x = counts
            x = log10(x+1)
            print("Clustering!")
            #x <- normalize.quantiles(as.matrix(x), copy = TRUE)
            colnames(x) <- colnames(counts)
            rownames(x) <- rownames(counts)
            
            clustering_used_method <- "ward.D"
            distance_used_method <- "euclidean"
            
            hclustfunc <- function(x) hclust(x, method = clustering_used_method)
            distfunc <- function(x) dist(x, method = distance_used_method)
            
            fit <- hclustfunc(distfunc(t(x)))
            
            # Convert hclust to dendrogram and then to dendrogram segments
            dend <- as.dendrogram(fit)
            p <- ggdendrogram(dend , rotate = FALSE, size = 2)
            fig=ggplotly(p)
          }
          else if (func_name == "density") {
            print("density")
            # Prepare the dataframe
            df_list <- lapply(seq_len(ncol(counts)), function(i) {
              data.frame(
                counts = counts[, i],
                idsample = colnames(counts)[i]
              )
            })
            
            df <- bind_rows(df_list)
            
            # Log-transform counts
            df$log_counts <- log(df$counts + 1)
            
            # Initialize plotly object
            p <- plot_ly()
            
            # Add a density curve per sample
            for (sample_name in unique(df$idsample)) {
              sample_data <- df %>% filter(idsample == sample_name)
              d <- density(sample_data$log_counts)
              
              p <- p %>%
                add_trace(
                  x = d$x,
                  y = d$y,
                  type = "scatter",
                  mode = "lines",
                  name = sample_name
                )
            }
            
            # Add layout
            p <- p %>%
              layout(
                title = "Log Gene Count Densities by Sample",
                xaxis = list(title = "log(count + 1)"),
                yaxis = list(title = "Density"),
                legend = list(title = list(text = "Sample"))
              )
            fig = p
          }
          else if (func_name =="pcacomp"){
            # Step 1: Remove genes (rows) with zero variance
            counts <-  counts[apply( counts, 1, var) != 0, ]
            pca_result <- prcomp(t(counts), scale. = TRUE)  # Transpose if samples are columns
            # Calculate variance explained
            var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
            # Create Plotly screeplot
            fig <- plot_ly(
              x = paste0("PC", 1:length(var_explained)),
              y = var_explained,
              type = "bar",
              marker = list(color = 'blue'),
              text = ~paste0(round(var_explained, 2), "%"),
              hoverinfo = 'text'
            ) %>%
              layout(
                title = "Principal Component Histogram",
                xaxis = list(title = "Principal Components"),
                yaxis = list(title = "Variance Explained (%)")
              )
            fig
          }
          else if (func_name == "MAplot")  {
            res_df = counts
            print(res_df)
            res_df$significant <- res_df$FDR < 0.05
            #plot(log2(res$baseMean + 1) , res$log2FoldChange ,col = "black", main="DESeq2 Fold Change Plot", xlab='Mean of Normalized Counts', ylab='log2FoldChange',pch=19,cex=0.3) 
            res_df$gene <- rownames(res_df)
            print(head(res_df))
            
            # Step 3: Create interactive MA plot
            ma_plot <- plot_ly(
              data = res_df,
              x = ~logCPM,
              y = ~log2FoldChange,
              type = 'scatter',
              mode = 'markers',
              color = ~significant,
              colors = c("gray", "red"),  # Customize colors if needed
              text = ~paste(
                "Gene: ", gene,
                "<br>logCPM: ", round(logCPM, 2),
                "<br>log2FC: ", round(log2FoldChange, 2)
              ),
              hoverinfo = 'text',
              marker = list(size = 6, opacity = 0.6)
            ) %>%
              layout(
                title = list(text =  "MA Plot", font = list(size = 16)),
                xaxis = list(title = "logCPM"),
                yaxis = list(title = "log2 Fold Change"),
                legend = list(title = list(text = "Significance"))
              )
            fig = ma_plot
          }
          else if (func_name == "dot")     {
            gene_symbol_list <- rownames(subset(counts, as.numeric(counts$FDR) < 0.05))
            entrez_ids <- mapIds(
              org.Hs.eg.db,
              keys = as.character(gene_symbol_list),
              column = "ENTREZID",
              keytype = "SYMBOL",
              multiVals = "first"
            )
            entrez_ids <- na.omit(entrez_ids)
            de <- as.character(entrez_ids)
            
            # KEGG enrichment
            x <- enrichKEGG(
              gene = de,
              organism = 'hsa',
              pAdjustMethod = "BH",
              qvalueCutoff = 0.001
            )
            
            # Convert to data frame
            x_df <- as.data.frame(x)
            
            # Keep only top N results (optional)
            x_df <- x_df %>% arrange(qvalue) #%>% head(20)
            
            # Create Plotly dot plot
            fig = plot_ly(
              data = x_df,
              x = ~GeneRatio,
              y = ~reorder(Description, -qvalue),
              type = 'scatter',
              mode = 'markers',
              marker = list(
                size = ~Count * 0.7,
                color = ~qvalue,
                colorscale = 'Magma',
                colorbar = list(title = "qvalue"),
                reversescale = TRUE,
                showscale = TRUE
              ),
              text = ~paste(
                "Pathway:", Description,
                "<br>Gene Ratio:", GeneRatio,
                "<br>Count:", Count,
                "<br>qvalue:", signif(qvalue, 4)
              ),
              hoverinfo = "text"
            ) %>%
              layout(
                title = "KEGG Pathway Enrichment (Interactive)",
                xaxis = list(title = "Gene Ratio"),
                yaxis = list(title = ""),
                margin = list(l = 150)  # for long pathway names
              )
          }
          else if (func_name == "KEGG") {
            # Step 1: Filter genes with padj < 0.05
            gene_symbol_list <- rownames(subset(counts, as.numeric(counts$FDR) < 0.05))
            entrez_ids <- mapIds(
              org.Hs.eg.db,
              keys = as.character(gene_symbol_list),
              column = "ENTREZID",
              keytype = "SYMBOL",
              multiVals = "first"
            )
            entrez_ids <- na.omit(entrez_ids)
            de <- as.character(entrez_ids)
            
            # Step 2: KEGG enrichment
            x <- enrichKEGG(
              gene = de,
              organism = 'hsa',
              pAdjustMethod = "BH",
              qvalueCutoff = 0.001
            )
            edox <- setReadable(x, "org.Hs.eg.db", "ENTREZID")
            
            # Step 3: Build edges and nodes (top 5 terms)
            top_n <- 10
            top_terms <- head(edox@result$Description, top_n)
            
            edges <- data.frame(from = character(), to = character(), stringsAsFactors = FALSE)
            nodes <- data.frame(id = character(), type = character(), stringsAsFactors = FALSE)
            
            for (term in top_terms) {
              genes <- strsplit(edox@result$geneID[edox@result$Description == term], "/")[[1]]
              edges <- rbind(edges, data.frame(from = term, to = genes))
              nodes <- rbind(nodes, data.frame(id = term, type = "pathway"))
              nodes <- rbind(nodes, data.frame(id = genes, type = "gene"))
            }
            
            nodes <- distinct(nodes)
            
            # Step 4: Graph layout
            g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
            lay <- layout_with_fr(g)
            
            # Assign layout coordinates to nodes
            node_coords <- data.frame(
              id = V(g)$name,
              x = lay[, 1],
              y = lay[, 2],
              type = V(g)$type
            )
            
            # Extract edges with coordinates
            edge_coords <- data.frame(
              x = numeric(),
              y = numeric(),
              xend = numeric(),
              yend = numeric()
            )
            
            for (e in 1:ecount(g)) {
              ends_ids <- ends(g, E(g)[e])
              from_coord <- node_coords[node_coords$id == ends_ids[1], c("x", "y")]
              to_coord <- node_coords[node_coords$id == ends_ids[2], c("x", "y")]
              edge_coords <- rbind(edge_coords, data.frame(
                x = from_coord$x,
                y = from_coord$y,
                xend = to_coord$x,
                yend = to_coord$y
              ))
            }
            
            # Step 5: Plot with Plotly
            p <- plot_ly(type = 'scatter', mode = 'lines')
            
            # Add edges
            for (i in 1:nrow(edge_coords)) {
              p <- add_segments(p,
                                x = edge_coords$x[i],
                                y = edge_coords$y[i],
                                xend = edge_coords$xend[i],
                                yend = edge_coords$yend[i],
                                line = list(color = 'blue', width = 0.5),
                                showlegend = FALSE
              )
            }
            
            # Add nodes
            p <- add_trace(p,
                           x = node_coords$x,
                           y = node_coords$y,
                           text = node_coords$id,
                           mode = 'markers+text',
                           textposition = 'top center',
                           marker = list(
                             size = 10,
                             color = ifelse(node_coords$type == "pathway", "black", "red")
                           ),
                           hoverinfo = "text",
                           showlegend = FALSE
            )
            
            p <- layout(p,
                        title = "Interactive KEGG cnetplot (Plotly)",
                        xaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                        yaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE)
            )
            
            fig = p
          }
          else if ((func_name == "net")|| (func_name == "integration")){
    
    # counts <- read.table("C:/Users/f.russo.ENGIBBC/Desktop/REDAC_submission/expression_file.txt", header = TRUE, sep = "\t", row.names = 1)
    # ppi_data  <- read.table("C:/Users/f.russo.ENGIBBC/Desktop/REDAC_submission/protein_protein_interaction_data.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    ppi_data  <- read.table(file ="9606.protein.physical.links.detailed.v12.0_with_gene_symbols_.txt", header = TRUE,sep = "\t")
    protein_info <- read.table(file ="9606.protein.info.v12.0.txt", header = TRUE,sep = "\t")
    colnames(ppi_data) <- c("protein_code_2", "protein_code_1", "experimental", "Protein1", "Protein2")
    # colnames(ppi_data) <- c("Protein1", "Protein2")
    #Adesso elimino i conteggi nulli o quelli inferiori a 10 su tutti i campioni
    #keep <- rowSums(counts > dim(counts)[2])
    #counts <- counts[keep,]
    gene_expr = counts
    interaction_threshold = 980
    ppi_data = subset(ppi_data, ppi_data$experimental > interaction_threshold)
    print(dim(ppi_data))
    # Ensure genes in PPI match those in the gene expression dataset
    valid_genes <- rownames(gene_expr)
    ppi_filtered <- ppi_data %>% filter(Protein1 %in% valid_genes & Protein2 %in% valid_genes)
    # Function to compute Pearson correlation between two interacting proteins
    compute_correlation2 <- function(protein1, protein2, gene_expr) {
      expr1 <- gene_expr[protein1, ]
      expr2 <- gene_expr[protein2, ]
      
      if (sd(expr1) == 0 || sd(expr2) == 0) {
        return(NA)  # Avoid zero variance issue
      }
      
      cor_value <- cor(as.numeric(expr1), as.numeric(expr2), method = "pearson", use = "complete.obs")
      cor_value 
    }
    
    ppi_filtered <- ppi_filtered[,3:5]
    ppi_filtered <- ppi_filtered[,c(2,3,1)]
    
    # Compute correlation for each interaction
    ppi_filtered <- ppi_filtered %>%
      rowwise() %>%
      mutate(Correlation = compute_correlation2(Protein1,Protein2, gene_expr) ) %>%
      drop_na()  # Remove NA values due to zero variance
    
    # Create a weighted graph
    ppi_graph <- graph_from_data_frame(ppi_filtered, directed = FALSE)
    
    # Generate layout
    layout <- layout_with_fr(ppi_graph)
    V(ppi_graph)$x <- layout[,1]
    V(ppi_graph)$y <- layout[,2]
    
    # Get edge and node data
    edges <- get.data.frame (ppi_graph)
    nodes <- data.frame(id = V(ppi_graph)$name,
                        x = V(ppi_graph)$x,
                        y = V(ppi_graph)$y)
    
    # Start plot
    fig <- plot_ly()
    
    # Add edge traces
    for (i in 1:nrow(edges)) {
      edge <- edges[i, ]
      v0 <- which(V(ppi_graph)$name == edge$from)
      v1 <- which(V(ppi_graph)$name == edge$to)
      
      fig <- fig %>% add_trace(
        type = "scatter",
        mode = "lines",  
        x = c(V(ppi_graph)$x[v0], V(ppi_graph)$x[v1]),
        y = c(V(ppi_graph)$y[v0], V(ppi_graph)$y[v1]),
        line = list(
          width = abs(edge$Correlation) * 2,
          color = ifelse(edge$Correlation > 0, "red", "blue"),
          opacity = 0.5
        ),
        hoverinfo = "text",
        text = paste0(edge$from, " - ", edge$to, "<br>Correlation: ", round(edge$Correlation, 2)),
        showlegend = FALSE
      )
    }
    
    # Add node trace
    fig <- fig %>% add_trace(
      type = "scatter",
      mode = "markers+text",
      x = nodes$x,
      y = nodes$y,
      text = nodes$id,
      textposition = "top center",
      marker = list(size = 2, color = 'skyblue'),
      hoverinfo = "text",
      showlegend = FALSE
    )
    
    # Final layout
    fig <- fig %>% layout(
      title = 'Gene Expression Weighted PPI Network',
      xaxis = list(title = "", showgrid = FALSE, zeroline = FALSE),
      yaxis = list(title = "", showgrid = FALSE, zeroline = FALSE)
    )
          }
        
      }
        n <- 100
        for (i in 1:n) {
          # Increment the progress bar, and update the detail text.
          incProgress(1/n, detail = paste("\n  ",i, "% completed!"))
          # Pause for 0.1 seconds to simulate a long computation.
          Sys.sleep(0.0005)
        }
      }) 
      fig
 }
})
  
  
  output$download_html_report <- downloadHandler(
    filename = function() {
      paste("DE_Report", Sys.Date(), ".html", sep = "")
    },
    content = function(file) {
      tempReport <- file.path("report.Rmd")
      writeLines(c(
        "---",
        "title: \"Report of the Analysis\"",
        "output: html_document",
        "params:",
        "  res: NA", 
        "  counts: NA", 
        "  treated_samples: NA", 
        "  wt_samples: NA", 
        "  which_results: NA",
        "---",
        "",
        "```{r setup, include=FALSE}",
        "knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)",
        "library(edgeR)",
        "library(ggplot2)",
        "```",
        "",
        "## Design and Dispersion",
        "```{r}",
        "Sys.time()",
        "",
        "counts <- params$counts",
        "res <- params$res",
        "#treated_samples <- colnames(counts)[as.numeric(params$treated_samples)]",
        "#wt_samples <-      colnames(counts)[as.numeric(params$wt_samples)]",
        "#which_results <- params$which_results",
        "",
        "# Build metadata",
        "#treated_condition <- data.frame(patients = treated_samples, condition = 'treated')",
        "#wt_condition <- data.frame(patients = wt_samples, condition = 'wt')",
        "#samples <- rbind(wt_condition, treated_condition)",
        "#rownames(samples) <- samples$patients",
        "",
        "# Subset and order count matrix",
        "#dataset <- counts[, samples$patients]",
        "#samples$condition <- factor(samples$condition, levels = c('wt', 'treated'))",
        "",
        "# Create DGEList",
        "#dge <- DGEList(counts = dataset, group = samples$condition)",
        "",
        "# Filter low counts",
        "#keep <- rowSums(counts) >= dim(counts)[2]",
        "#dge <- dge[keep, , keep.lib.sizes = FALSE]",
        "",
        "# Normalize and estimate dispersion",
        "#dge <- calcNormFactors(dge)",
        "#dge <- estimateDisp(dge)",
        "",
        "# Differential expression",
        "#et <- exactTest(dge, pair = c('wt', 'treated'))",
        "#res <- topTags(et, n = nrow(dge))$table",
        "#res <- res[order(res$FDR), ]",
        "#res$log2FoldChange <- round(res$logFC, 4)",
        "#res$logFC <- NULL",
        "",
        "# Filter DE results",
        "#if (which_results == 'down') {",
        "#  DEres <- res[res$log2FoldChange < 0, ]",
        "#} else if (which_results == 'up') {",
        "#  DEres <- res[res$log2FoldChange > 0, ]",
        "#} else {",
        "#  DEres <- res",
        "#}",
        "#DEres <- na.omit(DEres)",
        "",
        "# Plot",
        "#plotBCV(dge)",
        "```",
        "",
        "## Top 100 Differentially Expressed Genes",
        "```{r}",
        "knitr::kable(head(res, 100), caption = 'Top 100 Differentially Expressed Genes')",
        "",
        "# Plot Dendrogram",
        "x = counts",
        "x = log10(x+1)",
        "colnames(x) <- colnames(counts)",
        "rownames(x) <- rownames(counts)",
        "clustering_used_method <- 'ward.D'",
        "distance_used_method <- 'euclidean'",
        "hclustfunc <- function(x) hclust(x, method = clustering_used_method)",
        "distfunc <- function(x) dist(x, method = distance_used_method)",
        "fit <- hclustfunc(distfunc(t(x)))",
        "# Convert hclust to dendrogram and then to dendrogram segments",
        "dend <- as.dendrogram(fit)",
        "p <- ggdendrogram(dend , rotate = FALSE, size = 2)",
        "ggplotly(p)",
        "",
        "# Principal Component Histogram",
        "# Step 1: Remove genes (rows) with zero variance",
        "counts <-  counts[apply( counts, 1, var) != 0, ]",
        "pca_result <- prcomp(t(counts), scale. = TRUE)  # Transpose if samples are columns",
        "# Calculate variance explained",
        "var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100",
        "# Create Plotly screeplot",
        "scree_plot <- plot_ly(",
        "  x = paste0('PC', 1:length(var_explained)),",
        "  y = var_explained,",
        "  type = 'bar',",
        "  marker = list(color = 'blue'),",
        "  text = ~paste0(round(var_explained, 2), '%'),",
        "  hoverinfo = 'text'",
        ") %>%",
        "  layout(",
        "    title = 'Principal Component Histogram',",
        "    xaxis = list(title = 'Principal Components'),",
        "    yaxis = list(title = 'Variance Explained (%)')",
        "  )",
        "scree_plot",
        "",
        "# Principal Component Analysis",
        "x <- counts",
        "# Step 1: Remove genes (rows) with zero variance",
        "x <- x[apply(x, 1, var) != 0, ]",
        "# Step 2: Transpose and run PCA",
        "x_t <- t(x)",
        "pca_result <- prcomp(x_t, scale. = TRUE)",
        "# Step 3: Extract PCA scores",
        "pca_data <- as.data.frame(pca_result$x)",
        "pca_data$Sample <- rownames(pca_data)",
        "# Step 4: Get variance explained for labels",
        "variance_explained <- summary(pca_result)$importance[2, 1:2] * 100",
        "# Step 5: Plotly PCA plot",
        "pca_plot <- plot_ly(",
        "  data = pca_data,",
        "  x = ~PC1,",
        "  y = ~PC2,",
        "  type = 'scatter',",
        "  mode = 'markers+text',",
        "  text = ~Sample,",
        "  hoverinfo = 'text',",
        "  textposition = 'top right',",
        "  marker = list(size = 6, color = 'blue')",
        ") %>%",
        "  layout(",
        "    title = 'PCA: PC1 vs PC2',",
        "    xaxis = list(title = paste0('PC1 (', round(variance_explained[1], 1), '%)')),",
        "    yaxis = list(title = paste0('PC2 (', round(variance_explained[2], 1), '%)')),",
        "    showlegend = FALSE",
        "  )",
        "pca_plot",
        "",
        "# 3D Principal Component Analysis",
        "pca_result <- PCA(t(counts), graph = FALSE)  # Note the transpose 't()' if samples are columns",
        "# Extract the PCA scores (rows are now samples after transposing)",
        "pca_data <- data.frame(",
        "  Sample = colnames(counts),",
        "  PC1 = pca_result$ind$coord[, 1],",
        "  PC2 = pca_result$ind$coord[, 2],",
        "  PC3 = pca_result$ind$coord[, 3]",
        ")",
        "# Create interactive 3D plot with sample labels as hover text",
        "p <- plot_ly(",
        "  data = pca_data,",
        "  x = ~PC1, y = ~PC2, z = ~PC3,",
        "  type = 'scatter3d',",
        "  mode = 'markers+text',",
        "  text = ~Sample,",
        "  hoverinfo = 'text',",
        "  marker = list(size = 5, color = 'blue')",
        ")",
        "p <- p %>% layout(",
        "  title = '3D PCA Plot',",
        "  scene = list(",
        "    xaxis = list(title = 'PC1'),",
        "    yaxis = list(title = 'PC2'),",
        "    zaxis = list(title = 'PC3')",
        "  )",
        ")",
        "p",
        "",
        "# Violin Plot",
        "log_counts <- log10(counts + 1)",
        "# Step 2: Add gene names as a column before melting",
        "log_counts$Gene <- rownames(log_counts)",
        "# Step 3: Melt into long format",
        "log_counts_long <- melt(log_counts, id.vars = 'Gene', variable.name = 'Sample', value.name = 'Expression')",
        "violin_plot <- plot_ly(   ",   
        "  data = log_counts_long,",
        "  x = ~Sample,",
        "  y = ~Expression,",
        "  type = 'violin',",
        "  split = ~Sample,",
        "  box = list(",
        "    visible = T",
        "  ),",
        "  meanline = list(",
        "    visible = T",
        "  )",
        ")%>%",
        "  layout(",
        "    title = 'Violin Plot of log10(counts + 1)',",
        "    xaxis = list(title = 'Sample', tickangle = 90),",
        "    yaxis = list(title = 'log10(Expression + 1)')",
        "  )",
        "violin_plot",
        "",
        "# Densities",
        "df_list <- lapply(seq_len(ncol(counts)), function(i) {",
        "  data.frame(",
        "    counts = counts[, i],",
        "    idsample = colnames(counts)[i]",
        "  )",
        "})",
        "df <- bind_rows(df_list)",
        "# Log-transform counts",
        "df$log_counts <- log(df$counts + 1)",
        "# Initialize plotly object",
        "p <- plot_ly()",
        "# Add a density curve per sample",
        "for (sample_name in unique(df$idsample)) {",
        "  sample_data <- df %>% filter(idsample == sample_name)",
        "  d <- density(sample_data$log_counts)",
        "  p <- p %>%",
        "    add_trace(",
        "      x = d$x,",
        "      y = d$y,",
        "      type = 'scatter',",
        "      mode = 'lines',",
        "      name = sample_name",
        "    )",
        "}",
        "# Add layout",
        "p <- p %>%",
        "  layout(",
        "    title = 'Log Gene Count Densities by Sample',",
        "    xaxis = list(title = 'log(count + 1)'),",
        "    yaxis = list(title = 'Density'),",
        "    legend = list(title = list(text = 'Sample'))",
        "  )",
        "p",
        "",
        "# Heatmap",
        "if(dim(counts)[1]<1000){",
        "  # Compute log-transformed expression matrix",
        "  top_counts <- log(as.matrix(counts + 1))",
        "  # Scale by row (z-score normalization per gene)",
        "  scale_rows <- function(x) {",
        "    t(scale(t(x)))",
        "  }",
        "  scaled_counts <- scale_rows(top_counts)",
        "  # Create Plotly heatmap",
        "  heatmap_plot <- plot_ly(",
        "    y = colnames(scaled_counts),",
        "    x = rownames(scaled_counts),",
        "    z = t(scaled_counts),",
        "    type = 'heatmap',",
        "    colors = colorRampPalette(c('red', 'white', 'blue'))(500),",
        "    showscale = TRUE",
        "  ) %>%",
        "    layout(",
        "      title = list(text = 'Heatmap of the Most Variable 1000 Genes', font = list(size = 12)),",
        "      xaxis = list(title = '', tickangle = 90 ,tickfont = list(size = 10)),",
        "      yaxis = list(title = '', tickfont = list(size = 10)),",
        "      margin = list(l = 100, r = 50, b = 100, t = 50)",
        "    )",
        "  fig = heatmap_plot",
        "}else{",
        "log_counts <- log(as.matrix(counts + 1))",
        "# Select top 1000 most variable genes",
        "gene_variances <- apply(log_counts, 1, var)",
        "top_genes <- names(sort(gene_variances, decreasing = TRUE))[1:1000]",
        "top_counts <- log_counts[top_genes, ]",
        "# Scale by row (z-score normalization per gene)",
        "  scale_rows <- function(x) {",
        "    t(scale(t(x)))",
        "  }",
        "  scaled_counts <- scale_rows(top_counts)",
        "  # Step 4: Create Plotly heatmap",
        "  heatmap_plot <- plot_ly(",
        "    y = colnames(scaled_counts),",
        "    x = rownames(scaled_counts),",
        "    z = t(scaled_counts),",
        "    type = 'heatmap',",
        "    colors = colorRampPalette(c('red', 'white', 'blue'))(500),",
        "    showscale = TRUE",
        "  ) %>%",
        "    layout(",
        "      title = list(text = 'Heatmap of the Most Variable 1000 Genes', font = list(size = 12)),",
        "      xaxis = list(title = '', tickangle = 90 ,tickfont = list(size = 10)),",
        "      yaxis = list(title = '', tickfont = list(size = 10)),",
        "      margin = list(l = 100, r = 50, b = 100, t = 50)",
        "    )",
        "  fig = heatmap_plot",
        "} ",
        "fig",
        "",
        "# Volcano",
        "res_df <- res",
        "res_df$significant <- res_df$FDR < 0.05",
        "res_df$gene <- rownames(res_df)",
        "volcano_plot <- plot_ly(data = res_df, x = ~log2FoldChange, y = ~-log10(FDR),type = 'scatter',mode = 'markers', color = ~significant,colors = c('gray', 'red'),text = ~paste('Gene: ', gene,'<br>log2FC: '",
        ", round(log2FoldChange, 2),'<br>-log10(p): ', round(-log10(FDR), 2)),hoverinfo = 'text',marker = list(size = 6, opacity = 0.6)) %>% ", 
        "  layout(title = list(text = 'Volcano Plot', font = list(size = 16)),xaxis = list(title = 'log2 Fold Change'),yaxis = list(title = '-log10 FDR'),legend = list(title = list(text = 'Significance')))",
        "volcano_plot",
        "",
        "# Fold Change vs logCPM Plot",
        "ma_plot <- plot_ly(data = res_df, x = ~logCPM,y = ~log2FoldChange,type = 'scatter',mode = 'markers',color = ~significant,colors = c('gray', 'red'),text = ~paste(",
        "    'Gene: ', gene,'<br>logCPM: ', round(logCPM, 2),'<br>log2FC: ', round(log2FoldChange, 2)),hoverinfo = 'text',marker = list(size = 6, opacity = 0.6)) %>%",
        "  layout(title = list(text =  'Fold Change vs logCPM Plot', font = list(size = 16)),xaxis = list(title = 'logCPM'),yaxis = list(title = 'log2 Fold Change'),legend = list(title = list(text = 'Significance')))",
        "ma_plot",
        "",
        "Sys.time()",
        "",
        "sessionInfo()",
        "",
        "```"
      ), tempReport)
      
      
      params <- params_to_pass() #qui glieli passo
      
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv()))
    }
  )
  
  output$download_word_report <- downloadHandler(
    filename = function() {
      paste("DE_Report", Sys.Date(), ".docx", sep = "")
    },
    content = function(file) {
      tempReport <- file.path("report.Rmd")
      writeLines(c(
        "---",
        "title: \"Report of the Analysis\"",
        "output: word_document",
        "params:",
        "  res: NA", 
        "  counts: NA", 
        "  treated_samples: NA", 
        "  wt_samples: NA", 
        "  which_results: NA",
        "---",
        "",
        # "```{r setup, include=FALSE}",
        "```{r, fig.width=18, fig.height=12, out.width='200%', fig.dim = c(18, 16), fig.align='center'}",
        #"```{r, results='asis', echo=FALSE}",
        "knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)",
        "library(edgeR)",
        "library(ggplot2)",
        "```",
        "",
        "## Design and Dispersion",
        "```{r}",
        "Sys.time()",
        "",
        "counts <- params$counts",
        "res <- params$res",
        "#treated_samples <- colnames(counts)[as.numeric(params$treated_samples)]",
        "#wt_samples <-      colnames(counts)[as.numeric(params$wt_samples)]",
        "#which_results <- params$which_results",
        "",
        "# Build metadata",
        "#treated_condition <- data.frame(patients = treated_samples, condition = 'treated')",
        "#wt_condition <- data.frame(patients = wt_samples, condition = 'wt')",
        "#samples <- rbind(wt_condition, treated_condition)",
        "#rownames(samples) <- samples$patients",
        "",
        "# Subset and order count matrix",
        "#dataset <- counts[, samples$patients]",
        "#samples$condition <- factor(samples$condition, levels = c('wt', 'treated'))",
        "",
        "# Create DGEList",
        "#dge <- DGEList(counts = dataset, group = samples$condition)",
        "",
        "# Filter low counts",
        "#keep <- rowSums(counts) >= dim(counts)[2]",
        "#dge <- dge[keep, , keep.lib.sizes = FALSE]",
        "",
        "# Normalize and estimate dispersion",
        "#dge <- calcNormFactors(dge)",
        "#dge <- estimateDisp(dge)",
        "",
        "# Differential expression",
        "#et <- exactTest(dge, pair = c('wt', 'treated'))",
        "#res <- topTags(et, n = nrow(dge))$table",
        "#res <- res[order(res$FDR), ]",
        "#res$log2FoldChange <- round(res$logFC, 4)",
        "#res$logFC <- NULL",
        "",
        "# Filter DE results",
        "#if (which_results == 'down') {",
        "#  DEres <- res[res$log2FoldChange < 0, ]",
        "#} else if (which_results == 'up') {",
        "#  DEres <- res[res$log2FoldChange > 0, ]",
        "#} else {",
        "#  DEres <- res",
        "#}",
        "#DEres <- na.omit(DEres)",
        "",
        "# Plot",
        "#plotBCV(dge)",
        "```",
        "",
        "## Top 100 Differentially Expressed Genes",
        "```{r}",
        "knitr::kable(head(res, 100), caption = 'Top 100 Differentially Expressed Genes')",
        "",
        "# Plot Dendrogram",
        "x = counts",
        "x = log10(x+1)",
        "colnames(x) <- colnames(counts)",
        "rownames(x) <- rownames(counts)",
        "clustering_used_method <- 'ward.D'",
        "distance_used_method <- 'euclidean'",
        "hclustfunc <- function(x) hclust(x, method = clustering_used_method)",
        "distfunc <- function(x) dist(x, method = distance_used_method)",
        "fit <- hclustfunc(distfunc(t(x)))",
        "# Convert hclust to dendrogram and then to dendrogram segments",
        "dend <- as.dendrogram(fit)",
        "p <- ggdendrogram(dend , rotate = FALSE, size = 2)",
        "p",
        "",
        "# Principal Component Histogram",
        "# Step 1: Remove genes (rows) with zero variance",
        "counts <-  counts[apply( counts, 1, var) != 0, ]",
        "pca_result <- prcomp(t(counts), scale. = TRUE)  # Transpose if samples are columns",
        "# Calculate variance explained",
        "var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100",
        "# Create Plotly screeplot",
        "counts <- counts[apply(counts, 1, var) != 0, ]",
        "# Perform PCA",
        "pca_result <- prcomp(t(counts), scale. = TRUE)",
        "# Calculate variance explained",
        "var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100",
        "# Create a data frame for plotting",
        "scree_df <- data.frame(PC = paste0('PC', 1:length(var_explained)),Variance = var_explained)",
        "# Create ggplot scree plot",
        "ggplot(scree_df, aes(x = PC, y = Variance)) +",
        "  geom_bar(stat = 'identity', fill = 'blue') +",
        "  geom_text(aes(label = paste0(round(Variance, 2), '%')), vjust = -0.5, size = 3) +",
        "  labs(title = 'Principal Component Histogram',x = 'Principal Components', y = 'Variance Explained (%)') +",
        "  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))",
        "",
        "# Principal Component Analysis",
        "x <- counts",
        "# Step 1: Remove genes (rows) with zero variance",
        "x <- x[apply(x, 1, var) != 0, ]",
        "# Step 2: Transpose and run PCA",
        "x_t <- t(x)",
        "pca_result <- prcomp(x_t, scale. = TRUE)",
        "# Step 3: Extract PCA scores",
        "pca_data <- as.data.frame(pca_result$x)",
        "pca_data$Sample <- rownames(pca_data)",
        "# Step 4: Get variance explained for labels",
        "variance_explained <- summary(pca_result)$importance[2, 1:2] * 100",
        "# Step 5: Plotly PCA plot",
        "x_lab <- paste0('PC1 (', round(variance_explained[1], 1), '%)')",
        "y_lab <- paste0('PC2 (', round(variance_explained[2], 1), '%)')",
        "# Step 5: Create PCA scatter plot with ggplot2",
        "ggplot(pca_data, aes(x = PC1, y = PC2, label = Sample)) +",
          "geom_point(color = 'blue', size = 2) +",
          "geom_text(vjust = -0.5, hjust = 0.5, size = 2) +",
          "labs(title = 'PCA: PC1 vs PC2', x = x_lab, y = y_lab) +",
          "theme_minimal()",
        "",
        "# 3D Principal Component Analysis",
        "library(scatterplot3d)",
        "x <- counts",
        "# Step 1: Remove genes (rows) with zero variance",
        "x <- x[apply(x, 1, var) != 0, ]",
        "# Step 2: Transpose and run PCA",
        "x_t <- t(x)",
        "pca_result <- prcomp(x_t, scale. = TRUE)",
        "# Step 3: Extract PCA scores",
        "pca_data <- as.data.frame(pca_result$x)",
        "pca_data$Sample <- rownames(pca_data)",
        "# Step 4: Get variance explained for labels",
        "variance_explained <- summary(pca_result)$importance[2, 1:2] * 100",
        "s3d <- scatterplot3d(x = pca_data$PC1,y = pca_data$PC2,z = pca_result$ind$coord[, 3],pch = 19,color = 'blue',main = '3D PCA Plot',xlab = 'PC1', ylab = 'PC2', zlab = 'PC3')",
        "coords_2d <- s3d$xyz.convert(x = pca_data$PC1,y = pca_data$PC2,z = pca_result$ind$coord[, 3])",
        "text(x = coords_2d$x, y = coords_2d$y, labels = pca_data$Sample, pos = 3, cex = 0.3)",
        "",
        "# Violin Plot", 
        "data_list <- as.list(as.data.frame(log10(counts + 1)))",
        "vioplot(data_list,col = rainbow(length(data_list)),main = 'Violin Plot',las=2,cex.names=0.5)",
        "",
        "# Density Plot",
        "library(dplyr)",
        "df_list <- lapply(seq_len(ncol(counts)), function(i) {data.frame(counts = counts[, i],idsample = colnames(counts)[i])})",
        "df <- bind_rows(df_list)",
        "df$log_counts <- log10(df$counts + 1)",
        "ggplot(df, aes(x = log_counts, color = idsample)) +",
        "geom_density() +",
        "theme_minimal() +",
        "labs(title = 'Densities', x = 'log(count + 1)', y = 'Density')",
        "",
        "# Heatmap of the 1000 most variable genes",      
        "log_counts <- log10(counts + 1)",
        "gene_vars <- apply(counts, 1, var)",
        "top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:1000]",
        "top_counts <- log_counts[top_genes, ]",
        "scaled_counts <- t(scale(t(top_counts)))",
        "heatmap.2(t(scaled_counts),main='Heatmap',scale='col',trace='none', col = colorRampPalette(c('red', 'white', 'blue'))(500),Rowv=NA,cexRow=0.5,cexCol=0.2,dendrogram = 'none',keysize=2.2,margins=c(6,6),cex.main = 0.1)",
        "",
        "# Volcano",
        "res_df <- res",
        "res_df$significant <- res_df$FDR < 0.05",
        "res_df$gene <- rownames(res_df)",
        "res_df$neg_log10_FDR <- -log10(res_df$FDR)",
        "ggplot(res_df, aes(x = log2FoldChange, y = neg_log10_FDR)) +",
          "geom_point(aes(color = significant), size = 2.5, alpha = 0.6) +",
          "scale_color_manual(values = c('gray', 'red')) +",
          "theme_minimal() +",
          "labs(title = 'Volcano Plot',x = 'log2 Fold Change',y = '-log10 FDR',color = 'Significant') +",
          "theme(plot.title = element_text(size = 16, face = 'bold'),legend.title = element_text(size = 10))",
        "",
        "#Log2FoldChange vs logCPM Plot",
        "res_df <- res",
        "res_df$significant <- res_df$FDR < 0.05",
        "res_df$gene <- rownames(res_df)",
        "res_df$neg_log10_FDR <- -log10(res_df$FDR)",
        "ggplot(res_df, aes(x = logCPM, y = log2FoldChange)) +",
        "geom_point(aes(color = significant), size = 2.5, alpha = 0.6) +",
        "scale_color_manual(values = c('gray', 'red')) +",
        "theme_minimal() +",
        "labs(title = 'Log2FoldChange vs logCPM Plot',x = 'logCPM',y = 'log2FoldChange',color = 'Significant') +",
        "theme(plot.title = element_text(size = 16, face = 'bold'),legend.title = element_text(size = 10))",
        "",
        "Sys.time()",
        "",
        "sessionInfo()",
        "",
        "```"
      ), tempReport)
      params <- params_to_pass() #qui glieli passo
      rmarkdown::render(tempReport, output_file = file,params = params, envir = new.env(parent = globalenv()))
    }
  )
  
   params_to_pass <- eventReactive(input$run2, {
     req(input$file2, input$text2)
     #prendi i risultati già calcolati e non ricalcolarli
     result = dea_results2()$results
     # counts <- read.table("C:/Users/f.russo.ENGIBBC/Desktop/REDAC_submission/count_pippo.txt", header = TRUE, sep = "\t", row.names = 1)
     counts <- read_and_clean_colnames(input$file2$datapath, header = TRUE, sep = "\t", rownames = 1)
     #res = get_response_param(input$text2, input$file2)
     #treated_samples = res$treated_samples ##scrivere qui una funzione più snella
     #wt_samples      = res$wt_samples
     #which_results   = res$which_results
     list(res = result, counts = counts)
  })
   
   
   get_response_param <- function(prompt,inDataFile) {
     if (is.null(inDataFile) || (prompt=="")){return(NULL)}else{
       my_data <- read_and_clean_colnames(inDataFile$datapath, sep = "\t", header = TRUE, rownames = 1 )
       api_key <- "..."  # Replace with your API key
       url <- "https://api.together.xyz/v1/chat/completions"
       body <- list(
         model = "google/gemma-2-27b-it",
         messages = list(
           list(role = "system", content = "Respond only in JSON format. The JSON must have this structure: 
               {\"functiontoberun\":[\"analysis\"],
                \"treated\":[[\"col1\"],[\"col2\"]],\"wt\":[[\"col3\"],[\"col4\"]],
                \"regulated\":[\"down\",\"up\",\"de\"]}. Do not add explanations or other text."),
           list(role = "user", content = prompt)
         ),
         max_tokens = 500
       )
       response <- POST(
         url,
         add_headers(
           Authorization = paste("Bearer", api_key),
           `Content-Type` = "application/json"
         ),
         body = toJSON(body, auto_unbox = TRUE),
         encode = "json"
       )
       content <- httr::content(response, as = "text")
       print(content)
       parsed <- fromJSON(content)
       content = strsplit(as.character(parsed$choices$message$content),"\n" )[[1]][2] ## questa linea solo se usi Gemma LLM
       parsed <- fromJSON(content)
       treated_samples <- unlist(str_replace_all(parsed$treated, c("col" = "")))
       wt_samples      <- unlist(str_replace_all(parsed$wt, c("col" = "")))
       regulated       <- parsed$regulated
     }
     res = list(treated_samples = treated_samples, wt_samples = wt_samples, which_results = regulated)
     res
   }

})


  
  





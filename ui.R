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
source("definitions.R")
library(markdown)
library(reactome.db)
library(GO.db)
library(dplyr)
library(rentrez)
library(text2vec)
library(ReactomePA)
library(msigdbr)
library(rWikiPathways)
library(stopwords)
library(stringi)

options(shiny.maxRequestSize=100*2048^2)
options(repos = BiocManager::repositories())

# Define UI 
shinyUI(fluidPage(theme = shinytheme("united"),
                  tagList(
                    # CSS for navbar elements
                    tags$style(), 
                    # CSS for fonts
                    tags$style('h3 {font-weight: normal;}',
                               'h4 {font-weight: normal;}',
                               '* {font-family: Ubuntu;}'),
                    # CSS for all buttons (class .btn)
                    tags$style('.btn {color: #f1f518; 
                      background-color: #0d822c; 
                      border-color: #36447a;}
                .btn:hover {background-color: #734357;}'), 
                    # CSS for errors and validation messages
                    tags$style('.shiny-output-error-validation {
                 color: #e35300;
                 font-weight: bold;}'),
                    # CSS for individual elements: note the '#id' syntax
                    tags$style('#map-readme {color: #63071d;
                             background-color: #dce8e8;
                             border-color: #c84407;}
                #map-readme:hover {background-color: #e36a0075;}
                #trend-readme {color: #d1b6db;
                               background-color: #ff770050;
                               border-color: #c84407;}
                #trend-readme:hover {background-color: #e36a0075;}')
                  ), # end of style block
  titlePanel("REDAC: RNA-seq Expression Data Analysis Chatbot"),
  navbarPage("A Web App for analysing bulk RNA-seq data by asking questions written in English language (version: 2.3.9 released the 20/November/2025)", 
       tabPanel("Perform a Complete Analysis",
          sidebarLayout(
            sidebarPanel(
              helpText("You can use this chatbot by uploading your bulk RNA-seq raw count data (positive integers) in a tab separated format (as the one shown below)."),
              helpText(" "),
              helpText(HTML("Gene&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ID1&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ID2&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ID3&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ID4")) ,
              helpText(HTML("SEC24B-AS1&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;47	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;3&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;0&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;26")),
              helpText(HTML("A1BG0&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;410&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;3&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;14&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;4")),
              helpText(HTML("A1CF&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;192&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;202&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;156&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;63")),
              helpText(HTML("GGACT&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;28&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;23&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;17&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;15")),
              helpText(HTML("...  ...  ...")),
              helpText(" "),
              helpText("Then, write a request and click on 'Run Analysis!' button below."),
              helpText(" "),
              helpText("For more details, you can find the user manual for local installation and a complete user guide here: https://github.com/franruss/REDAC/blob/main/docs/REDAC_user_manual.pdf"),
              helpText(" "),
              helpText("If you need help using REDAC or if there is an error with the data you entered, please send an email to: francesco.russo AT cnr.it"),
              helpText(" "),
              helpText("Please note that errors are sometimes caused by the Shiny server (connection issues, temporarily unavailable server resources, etc.) on which REDAC is running, or by problems requesting LMMs via an API, which takes several seconds to execute and 
                       return to REDAC. In these cases, please retry the request and wait for a response."),
              helpText("---------------------------------------------------------------------------------- "),
              helpText(" "),
              helpText("YOU CAN TRY THIS SAMPLE FILE !!!"),
              helpText(" "),
              checkboxInput("use_example", "Check this box to try REDAC on a sample dataset!", value = FALSE),
              helpText(" "),
              helpText("This sample file is named 'Raw data RNAseq for submission with gene name 3 1 24.xlsx' and is available at https://zenodo.org/records/11057181"),
              helpText(" "),
              helpText("---------------------------------------------------------------------------------- "),
              fileInput('file2', "Please, upload your bulk RNA-seq raw count data (positive integers) in a tab separated format",accept=c('text/csv','text/comma-separated-values,text/plain','.csv')),
              textInput('text2', "Please, write your request (note: do not use special characters)","Example: perform an rnaseq analysis between treated 3,4 and wt 1,2 samples, to find up regulated genes",width="800px"),
              width = 30,
            ),
            # Show a tabset
            mainPanel(
              helpText(" "),
              tabsetPanel(tabPanel("The input file",DT::DTOutput('show_input_fun2'))),
              tags$hr(),
              helpText("__________________________________________________________________________________________________________________________"),
              helpText("<< DATA INSPECTION PLOTS >>"),
              tabsetPanel(type = "tabs",
                         tabPanel("Dendrogram", plotlyOutput("dendroPlot2", height = '400', width = '1200')),
                         tabPanel("PCA Components", plotlyOutput("pcaCompPlot2", height = '300', width = '1200')),
                         tabPanel("PCA", plotlyOutput("pcaPlot2", height = '400', width = '1200')),
                         tabPanel("PCA 3D", plotlyOutput("pca3DPlot2", height = '700', width = '1200')),
                         tabPanel("Violin Plot", plotlyOutput("violinPlot2", height = '400', width = '1600')),
                         tabPanel("Densities", plotlyOutput("densityPlot2", height = '400', width = '1200')),
                         tabPanel("Heatmap", plotlyOutput("heatmapPlot2", height = '700', width = '1600'))
              ),
              helpText(" "),
              tags$hr(),
              helpText("__________________________________________________________________________________________________________________________"),
              helpText(" "),
              actionButton("run2", "Run Analysis!"),
              helpText(" "),
              helpText(" "),
              helpText("__________________________________________________________________________________________________________________________"),
              helpText("<< ANSWER >>"),
              tabsetPanel(
                tabPanel("Llama's answer", verbatimTextOutput("resultsTable6"))
              ),
              helpText(" "),
              helpText("__________________________________________________________________________________________________________________________"),
              helpText("<< RESULT TABLE >>"),
              tabsetPanel(
                  tabPanel("Result Table", DT::DTOutput("resultsTable2"))
              ),
              downloadButton("download_results2", "Download the Result Table"),
              helpText(" "),
              helpText("__________________________________________________________________________________________________________________________"),
              helpText("<< RESULT INSPECTION PLOTS, DISCUSSIONS AND ALTERNATIVE CODES >>"),
              tabsetPanel(tabPanel("Volcano Plot", plotlyOutput("volcanoPlot2",height = '600', width = '1200')),
                          tabPanel("MA Plot", plotlyOutput("foldchangePlot2",height = '600', width = '1200')),
                          tabPanel("Code for Python developers (by Llama)", uiOutput("chat_output3")),
                          tabPanel("Alternative code for R developers (by Llama)", uiOutput("chat_output2")),
                          tabPanel("Analysis discussion (by Gemma)", uiOutput('chat_output_short_advice2'))
              ),
              helpText(" "),
              # helpText("__________________________________________________________________________________________________________________________"),
              # helpText("<<  >>"),
              # tabsetPanel( ),
              br(),
              br(),
              downloadButton("download_html_report", "Download the HTML Report"),
              br(),
              br(),
              downloadButton("download_word_report", "Download the Word Report"),
              br(),
              br(),
              br(),
              br(),
              tags$head(tags$style( HTML(".shiny-notification {background-color:yellow;position:fixed;top: 50%;left: 5%;right: 5%;}"))),
              tags$hr(),
              width = 30,
            )
          )
       ),
          
      tabPanel("Enrichment Analysis and Result Interpretation",
                   sidebarLayout(
                     sidebarPanel(
                       helpText("You can perform an enrichment analysis by uploading your edgeR result file in a tab separated format (as the one shown below)."),
                       helpText(" "),
                       helpText(HTML("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;logCPM&nbsp;&nbsp;&nbsp;&nbsp;PValue&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FDR&nbsp;&nbsp;&nbsp;log2FoldChange")) ,
                       helpText(HTML("SEC24B-AS1&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;3.86	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2.1e-12&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;...&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;...")),
                       helpText(HTML("A1BG0&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;3.04&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.3e-11&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;...&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;...")),
                       helpText(HTML("A1CF&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;6.07&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.1e-10&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;...&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;...")),
                       helpText(HTML("...  ...  ...")),
                       helpText(" "),
                       helpText(HTML("Then, choose an enrichment method and click on 'Enrich!' button below.")),
                       helpText(" "),
                       helpText("Finally, this chatbot can suggest an interpretation of your results via two LLMs: Gemma and Llama."),
                       helpText(" "),
                       helpText("You can find the user manual for local installation and a complete user guide here: https://github.com/franruss/REDAC/blob/main/docs/REDAC_user_manual.pdf"),
                       helpText(" "),
                       helpText("If you need help using REDAC or if there is an error with the data you entered, please send an email to: francesco.russo AT cnr.it"),
                       helpText(" "),
                       helpText(" "),
                       helpText("Please note that errors are sometimes caused by the Shiny server (connection issues, temporarily unavailable server resources, etc.) on which REDAC is running, or by problems requesting LMMs via an API, which takes several seconds to execute and 
                       return to REDAC. In these cases, please retry the request and wait for a response."),
                       helpText("---------------------------------------------------------------------------------- "),
                       helpText(" "),
                       helpText("YOU CAN TRY THIS EDGER RESULT SAMPLE FILE !!!"),
                       helpText(" "),
                       checkboxInput("use_example2", "Check this box to try REDAC on a edgeR result sample file!", value = FALSE),
                       helpText(" "),
                       helpText(" "),
                       helpText("---------------------------------------------------------------------------------- "),
                       fileInput('file3', "Please, upload your edgeR result file in a tab separated format",accept=c('text/csv','text/comma-separated-values,text/plain','.csv')), 
                       textInput('text7', "Please, write some keywords that can help both the LLMs understand the context and the scope of your analysis 
                                 (REDAC uses them to retrieve relevant articles from PubMed to build an RAG for the LLM interpretations):",
                                 "Example: Explain these results obtained by comparing... to test this hypothesis... to study these mechanisms... in this context... in order to understand... ",width="2000px"),
                       radioButtons('enrich_method', 'Select an enrichment method: ',
                                                              c(' KEGG    '='kegg',
                                                                ' Reactome               '='enrichPathway',
                                                                ' Disease Ontology       '='enrichDO',
                                                                ' Wiki Pathways          '='WikiPathways',
                                                                ' GO Biological Process  '='enrichGOBP',
                                                                ' GO Molecular Function  '='enrichGOMF',
                                                                ' GO Cellular Component  '='enrichGOCC',
                                                                ' Hallmark (MSD)         '='msigdbrH',
                                                                ' Positional (MSD)       '='msigdbrC1',
                                                                ' Curated (MSD)          '='msigdbrC2',
                                                                ' Regulatory target (MSD)'='msigdbrC3',
                                                                ' Computational (MSD)    '='msigdbrC4',
                                                                ' GO (MSD)               '='msigdbrC5',
                                                                ' Oncogenic (MSD)        '='msigdbrC6',
                                                                ' Immunologic (MSD)      '='msigdbrC7',
                                                                ' Cell type (MSD)        '='msigdbrC8'
                                                               )
                       ),
                       width = 26,
                     ),
                     # Show a tabset
                     mainPanel(
                       helpText(" "),
                       helpText("Data Tables"),
                       tabsetPanel(
                         tabPanel("The Input",DT::DTOutput('inputTableEnrich'))
                       ), 
                       tabsetPanel(
                         tabPanel("Answer", verbatimTextOutput("resultsTable10"))
                       ),
                       helpText(" "),
                       helpText(" "),
                       actionButton("run3", "Enrich!"),
                       helpText(" "),
                       helpText(" "),
                       tabsetPanel(  
                         tabPanel("Enrichment Results", DT::DTOutput("resultsTableEnrich")),
                         tabPanel("Dot Plot", plotlyOutput("generate_dotplot", height = '1800', width = '1300')),
                         #tabPanel("A possible interpretation (Gemma):", uiOutput('chat_output_interpretationGemma')),
                         tabPanel("A possible interpretation (by LlaMA powered by a PubMed-based RAG):", uiOutput('chat_output_interpretationRAGLlaMA')),
                         tabPanel("A possible interpretation (by Gemma powered by a PubMed-based RAG):", uiOutput('chat_output_interpretationRAGGemma'))
                       ),
                       downloadButton("download_results3", "Download Enrichment Results"),
                       helpText(" "),
                       #helpText("__________________________________________________________________________________________________________________________"),
                       #helpText("Interpretation of Results"),
                       #tabsetPanel(),
                       tags$head(tags$style( HTML(".shiny-notification {background-color:yellow;position:fixed;top: 50%;left: 5%;right: 5%;}"))),
                       tags$hr(),
                       width = 26,
                     )
                   )
      ),
      
      tabPanel("Plot Generation",
               sidebarLayout(
                 sidebarPanel(
                   helpText("You can use this Chatbot by uploading your data in a tab separated format."),
                   helpText("Write a request and click on 'Create Plot!' button below."),
                   helpText("You can create a plot on the count data within this list:"),
                   helpText("boxplot, violin, heatmap, correlation heatmap, pca, 3Dpca, dendrogram, density, pca components, network, surface."), 
                   helpText("Moreover, you can create a plot on a result file within this list:"),
                   helpText("volcano, maplot, dotplot, KEGGnet."),
                   helpText(" "),
                   helpText("You can find the user manual for local installation and a complete user guide here: https://github.com/franruss/REDAC/blob/main/docs/REDAC_user_manual.pdf"),
                   helpText(" "),
                   helpText("If you need help using REDAC or if there is an error with the data you entered, please send an email to: francesco.russo AT cnr.it"),
                   helpText(" "),
                   helpText(" "),
                   helpText("Please note that errors are sometimes caused by the Shiny server (connection issues, temporarily unavailable server resources, etc.) on which REDAC is running, or by problems requesting LMMs via an API, which takes several seconds to execute and 
                       return to REDAC. In these cases, please retry the request and wait for a response."),
                   helpText("---------------------------------------------------------------------------------- "),
                   helpText(" "),
                   helpText("YOU CAN EITHER TRY THIS COUNT SAMPLE FILE... "),
                   helpText(" "),
                   checkboxInput("use_example3", "Check this box to try these functions: boxplot, violin, heatmap, correlation heatmap, 
                                 pca, 3Dpca, dendrogram, density, pca components, network, surface.", value = FALSE),
                   helpText(" "),
                   helpText(" "),
                   helpText("---------------------------------------------------------------------------------- "),
                   helpText(" "),
                   helpText(" ...OR YOU CAN TRY THIS EDGER RESULT SAMPLE FILE !!!"),
                   helpText(" "),
                   checkboxInput("use_example4", "Check this box to try these functions: volcano, maplot, dotplot, KEGGnet. ", value = FALSE),
                   helpText(" "),
                   helpText(" "),
                   helpText("---------------------------------------------------------------------------------- "),
                   fileInput('file6', "Please, upload a file",accept=c('text/csv','text/comma-separated-values,text/plain','.csv')),
                   textInput('text6', "Please, write your request:","Example: create a heatmap",width="800px"),
                   helpText(" "),
                   actionButton("run6", "Create Plot!"),
                   #downloadButton("download_results2", "Download Results"),
                   width = 26,
                 ),
                 # Show a tabset
                 mainPanel( 
                     helpText(" "),
                     #tabsetPanel(tabPanel("Your Input", DT::DTOutput('show_input_fun3'))),
                     tabsetPanel(
                       tabPanel("Answer", verbatimTextOutput("resultsTable11"))
                     ),
                     helpText(" "),
                     helpText(" "),
                     helpText("__________________________________________________________________________________________________________________________"),
                     tabsetPanel(    
                           tabPanel("Default Plot", plotlyOutput("generate_plot3", height = '1800', width = '1300'))
                     ),
                     helpText(" "),
                     helpText("__________________________________________________________________________________________________________________________"),
                     tabsetPanel(tabPanel("Alternative code for Python developers (Llama)", uiOutput("chat_outputPythonDev")),
                                 tabPanel("Alternative code for R developers (Llama)", uiOutput("chat_outputRDev"))
                                ),
                   tags$head(tags$style( HTML(".shiny-notification {background-color:yellow;position:fixed;top: 50%;left: 5%;right: 5%;}"))),
                   tags$hr(),
                   width = 26,
                 )
               )
      ),
)))


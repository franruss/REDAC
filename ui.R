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
source("definitions.R")

options(shiny.maxRequestSize=100*2048^2)
options(repos = BiocManager::repositories())

# Define UI 
shinyUI(fluidPage(theme = shinytheme("united"),
  titlePanel("REDAC: RNA-seq Expression Data Analysis Chatbot"),
  navbarPage("A Web App for analysing bulk RNA-seq data by asking questions written in English language", 
       tabPanel("Perform a Complete Analysis",
          sidebarLayout(
            sidebarPanel(
              helpText("You can use this chatbot by uploading your bulk RNA-seq raw count data (positive integers) in a tab separetad format (as the one shown below)."),
              helpText(" "),
              helpText(HTML("Gene&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ID27&nbsp;&nbsp;&nbsp;&nbsp;ID33&nbsp;&nbsp;&nbsp;ID68&nbsp;&nbsp;&nbsp;ID70")) ,
              helpText(HTML("SEC24B-AS1&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;47	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;3&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;0&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;26")),
              helpText(HTML("A1BG0&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;410&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;3&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;14&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;4")),
              helpText(HTML("A1CF&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;192&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;202&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;156&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;63")),
              helpText(HTML("GGACT&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;28&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;23&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;17&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;15")),
              helpText(HTML("...  ...  ...")),
              helpText(" "),
              helpText("Then, write a request and click on 'Run DE Analysis!' button below."),
              helpText(" "),
              fileInput('file2', "Please, upload your bulk RNA-seq raw count data (positive integers) in a tab separetad format",accept=c('text/csv','text/comma-separated-values,text/plain','.csv')),
              textInput('text2', "Please, write your request (by Gemma)","Example: perform an rnaseq analysis between treated 3,4 and wt 1,2 samples, up regulated",width="800px"),
              width = 30,
            ),
            # Show a tabset
            mainPanel(
              helpText(" "),
              tabsetPanel(tabPanel("Your Input",DT::DTOutput('show_input_fun2'))),
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
              tags$hr(),
              helpText("<< RESULT TABLE >>"),
              tabsetPanel(
                  tabPanel("Result Table", DT::DTOutput("resultsTable2"))
              ),
              downloadButton("download_results2", "Download the Result Table"),
              helpText(" "),
              helpText("__________________________________________________________________________________________________________________________"),
              helpText("<< RESULT INSPECTION PLOTS, DISCUSSION AND ALTERNATIVE CODE>>"),
              tabsetPanel(tabPanel("Volcano Plot", plotlyOutput("volcanoPlot2",height = '600', width = '1200')),
                          tabPanel("MA Plot", plotlyOutput("foldchangePlot2",height = '600', width = '1200')),
                          tabPanel("Analysis discussion (by Gemma)", uiOutput('chat_output_short_advice2')),
                          tabPanel("Alternative code for R developers (by Llama)", uiOutput("chat_output2"))
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
                       helpText("You can perform an enrichment analysis by uploading your edgeR result file in a tab separetad format (as the one shown below)."),
                       helpText(" "),
                       helpText(HTML("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;logCPM&nbsp;&nbsp;&nbsp;&nbsp;PValue&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FDR&nbsp;&nbsp;&nbsp;log2FoldChange")) ,
                       helpText(HTML("SEC24B-AS1&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;3.86	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2.1e-12&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;...&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;...")),
                       helpText(HTML("A1BG0&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;3.04&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.3e-11&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;...&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;...")),
                       helpText(HTML("A1CF&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;6.07&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.1e-10&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;...&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;...")),
                       helpText(HTML("...  ...  ...")),
                       helpText(" "),
                       helpText("Finally, this chatbot can suggest an interpretation of your results via two LLMs, such as: Gemma and Llama."),
                       fileInput('file3', "Please, upload your edgeR result file in a tab separetad format",accept=c('text/csv','text/comma-separated-values,text/plain','.csv')),   
                       actionButton("run3", "Enrich!"),
                       
                       width = 26,
                     ),
                     # Show a tabset
                     mainPanel(
                       helpText(" "),
                       helpText("Data Tables"),
                       tabsetPanel(
                         tabPanel("Your Input",DT::DTOutput('inputTableEnrich')),
                         tabPanel("Enrichment Results", DT::DTOutput("resultsTableEnrich")),
                         tabPanel("Dot Plot", plotlyOutput("generate_dotplot", height = '1800', width = '1300')),
                         tabPanel("A possible interpretation (Gemma):", uiOutput('chat_output_interpretationGemma')),
                         tabPanel("Another possible interpretation (Llama):", uiOutput('chat_output_interpretationLlama'))
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
                   helpText("You can use this Chatbot by uploading your data in a tab separetad format."),
                   helpText("Write a request and click on 'Create Plot!' button below."),
                   helpText("You can create a plot on the count data within this list:"),
                   helpText("boxplot, violin, heatmap, correlation heatmap, pca, 3Dpca, dendrogram, density, pca components, network, surface."), 
                   helpText("Moreover, you can create a plot on a result file within this list:"),
                   helpText("volcano, maplot, dotplot, KEGGnet."),
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
                     #helpText("__________________________________________________________________________________________________________________________"),
                     tabsetPanel(    
                           tabPanel("Default Plot", plotlyOutput("generate_plot3", height = '1800', width = '1300'))
                     ),
                     helpText(" "),
                     helpText("__________________________________________________________________________________________________________________________"),
                     tabsetPanel(tabPanel("Alternative code for R developers (Llama)", uiOutput("chat_output"))),
                   tags$head(tags$style( HTML(".shiny-notification {background-color:yellow;position:fixed;top: 50%;left: 5%;right: 5%;}"))),
                   tags$hr(),
                   width = 26,
                 )
               )
      ),
)))


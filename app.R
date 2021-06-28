#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# Install any required packages.
packages <- c("shiny",
              "shinyjs",
              "shinydashboard",
              "shinyFiles",
              "DT",
              "tidyverse",
              "ggplot2",
              "reshape2",
              "biomaRt",
              "lattice",
              "ggrepel",
              "RColorBrewer")

install.packages(setdiff(packages, rownames(installed.packages())))


# Load packages.
library(shiny)
library(shinyjs)
library(shinydashboard)
library(shinyFiles)
library(DT)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(biomaRt)
library(lattice)
library(ggrepel)
library(RColorBrewer)


# Create parsing function.
parsing_fun <- function(text_input, symbol = c(",", " ", "\t", "\n", ";", "\\|")) {
    symbol_vec <- c()
    loc_vec <- c()
    for(sym in symbol){
        loc <- gregexpr(sym, text_input)[[1]][1]
        if(loc != -1){
            symbol_vec <- c(symbol_vec, sym)
            loc_vec <- c(loc_vec, loc)
        }
    }
    
    names(loc_vec) <- symbol_vec
    if(length(loc_vec) > 1){
        loc_vec <- loc_vec[order(loc_vec)]
    }
    
    return(names(loc_vec))
}

# Function to calculate top genes.
get_top_genes <- function(genes_df, cellex_df, annotLookup) {
    
    genes_output_df <- lapply(X = colnames(cellex_df)[-1], FUN = function(anno){
        es_nonzero <- cellex_df[c("gene", anno)] %>% filter_(paste0(anno, " > 0"))
        es_nonzero <- es_nonzero %>% arrange_(paste0("-", anno))
        es_nonzero <- es_nonzero %>% add_column(cellex_percentile = sort(1:nrow(es_nonzero), decreasing = TRUE)/nrow(es_nonzero))
        
        joined_df <- merge(es_nonzero, genes_df, by.x = "gene", by.y = "GENE")
        joined_df <- merge(joined_df, annotLookup, by.x = "gene", by.y = "ensembl_gene_id")
        joined_df <- joined_df %>% dplyr::select(-CHR, -START, -STOP, -NSNPS, -NPARAM, -N, -ZSTAT, -ZFITTED_BASE, -ZRESID_BASE, -uniprotswissprot)
        joined_df <- joined_df %>% 
            melt(id=c(1,3,4,5,6)) %>% 
            dplyr::rename(GENE = gene, 
                          GENE_HGNC = hgnc_symbol,
                          CELLTYPE_ANNOTATION = variable, 
                          CELLEX_VALUE = value,
                          CELLEX_PERCENTILE = cellex_percentile,
                          GENE_PERCENTILE = gene_percentile)
        joined_df <- joined_df %>% dplyr::relocate(GENE, GENE_HGNC, CELLTYPE_ANNOTATION, CELLEX_VALUE, P, CELLEX_PERCENTILE, GENE_PERCENTILE)
        joined_df <- joined_df %>% dplyr::arrange(GENE) %>% dplyr::distinct(GENE, .keep_all = TRUE)
        joined_df <- joined_df %>% arrange(-CELLEX_PERCENTILE, -GENE_PERCENTILE)
        return(joined_df)
    })
    
    genes_output_df <- bind_rows(genes_output_df)
    
    return(genes_output_df)
}


# Manhatten plot function.
manhatten.plot <- function(genes_df, annotLookup, gene_set = NULL, og_genes_df = FALSE, correction_method = c("bonferroni", "fdr", "none")[1], sig_thres = 0.05, return_just_genes_df = FALSE) {
    
    if(correction_method == "bonferroni") {
        
        sig <- sig_thres / nrow(genes_df)
        
        if(og_genes_df) {
            genes_df <- merge(genes_df, annotLookup, by.x = "GENE", by.y = "ensembl_gene_id")
        }
        
        genes_df <- genes_df %>% add_column(sign = genes_df$P < sig)
        
        data_cum <- genes_df %>% 
            group_by(CHR) %>% 
            summarise(max_bp = max(STOP)) %>% 
            mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
            dplyr::select(CHR, bp_add)
        
        genes_df <- genes_df %>% 
            inner_join(data_cum, by = "CHR") %>% 
            dplyr::mutate(bp_cum = STOP + bp_add)
        
        axis_set <- genes_df %>% 
            group_by(CHR) %>% 
            summarize(center = mean(bp_cum))
        
        ylim <- abs(floor(log10(min(genes_df$P)))) + 2
        
        
        manhplot <- ggplot(genes_df %>% filter(!sign), aes(x = bp_cum, y = -log10(P), 
                                                           color = as_factor(CHR), size = log10(NSNPS))) + # "P" or "p_adjust" or "NSNPS"
            geom_point(alpha = 0.7, size = genes_df %>% filter(!sign) %>% dplyr::select(NSNPS) %>% unlist(use.names = FALSE) %>% log10()) +
            scale_x_continuous(label = axis_set$CHR, breaks = axis_set$center) +
            scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
            scale_color_manual(values = rep(c("skyblue1", "lightskyblue1"), unique(length(axis_set$CHR)))) +
            scale_size_continuous(range = c(0.5,3)) +
            labs(x = "Chromosome", 
                 y = "-log10(p-value)") + 
            ggtitle(paste0(genes_df %>% filter(sign) %>% nrow(), " significant genes (p-value < ", sig, ")")) +
            theme_minimal() +
            theme( 
                legend.position = "none",
                panel.border = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                axis.title.x = element_text(size = 20),
                axis.title.y = element_text(size = 20),
                axis.text.x = element_text(angle = 90, size = 16, vjust = 0.5),
                axis.text.y = element_text(size = 16),
                plot.title = element_text(hjust = 0.5, size = 24)
            )
        
        #good_gene_names <- genes_df %>% filter(sign) %>% dplyr::select(hgnc_symbol) %>% unlist(use.names = FALSE)
        
        manhplot <- manhplot + 
            geom_point(data = genes_df %>% filter(sign), 
                       aes(x = bp_cum, y = -log10(P), size = log10(NSNPS)), # size = -log10(P)  /  -log10(NSNPS)
                       color = "darkgreen", 
                       alpha = 0.7) +
            geom_hline(yintercept = -log10(sig), color = "darkred", linetype = "dashed", size = 1)
        
        if(!is.null(gene_set)) {
            
            manhplot <- manhplot + # "P" or "p_adjust"
                geom_point(data = genes_df %>% filter(hgnc_symbol %in% gene_set), 
                           aes(x = bp_cum, y = -log10(P), size = log10(NSNPS)), # size = 8
                           color = "red2") +
                ggrepel::geom_text_repel(data = genes_df %>% filter(hgnc_symbol %in% gene_set),
                                         aes(x = bp_cum, y = -log10(P), label = hgnc_symbol), 
                                         size = 5, color = "red2", force = 0.75)
            
        }
        
    } else {
        
        if(og_genes_df) {
            genes_df <- merge(genes_df, annotLookup, by.x = "GENE", by.y = "ensembl_gene_id")
        }
        
        genes_df <- genes_df %>% add_column(p_adjust = p.adjust(genes_df$P, method = correction_method))
        
        data_cum <- genes_df %>% 
            group_by(CHR) %>% 
            summarise(max_bp = max(STOP)) %>% 
            mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
            dplyr::select(CHR, bp_add)
        
        genes_df <- genes_df %>% 
            inner_join(data_cum, by = "CHR") %>% 
            dplyr::mutate(bp_cum = STOP + bp_add)
        
        axis_set <- genes_df %>% 
            group_by(CHR) %>% 
            summarize(center = mean(bp_cum))
        
        ylim <- abs(floor(log10(min(genes_df$p_adjust)))) + 2 # "P" or "p_adjust"
        
        
        manhplot <- ggplot(genes_df %>% filter(p_adjust > sig_thres), aes(x = bp_cum, y = -log10(p_adjust), 
                                                                     color = as_factor(CHR), size = log10(NSNPS))) + # "P" or "p_adjust" or "NSNPS"
            geom_point(alpha = 0.7) +
            scale_x_continuous(label = axis_set$CHR, breaks = axis_set$center) +
            scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
            scale_color_manual(values = rep(c("skyblue1", "lightskyblue1"), unique(length(axis_set$CHR)))) +
            scale_size_continuous(range = c(0.5,3)) +
            labs(x = "Chromosome", 
                 y = "-log10(p-value)") + 
            ggtitle(paste0(genes_df %>% filter(p_adjust < sig_thres) %>% nrow(), " significant genes (p-value < ", sig_thres, ")")) +
            theme_minimal() +
            theme( 
                legend.position = "none",
                panel.border = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                axis.title.x = element_text(size = 20),
                axis.title.y = element_text(size = 20),
                axis.text.x = element_text(angle = 90, size = 16, vjust = 0.5),
                axis.text.y = element_text(size = 16),
                plot.title = element_text(hjust = 0.5, size = 24)
            )
        
        #good_gene_names <- genes_df %>% filter(sign) %>% dplyr::select(hgnc_symbol) %>% unlist(use.names = FALSE)
        
        manhplot <- manhplot + 
            geom_point(data = genes_df %>% filter(p_adjust < sig_thres), aes(x = bp_cum, y = -log10(p_adjust), size = log10(NSNPS)),
                       color = "darkgreen", alpha = 0.7) + # "P" or "p_adjust" or "NSNPS"
            geom_hline(yintercept = -log10(sig_thres), color = "darkred", linetype = "dashed", size = 1)
        
        
        if(!is.null(gene_set)) {
            
            manhplot <- manhplot + # "P" or "p_adjust"
                geom_point(data = genes_df %>% filter(hgnc_symbol %in% gene_set), 
                           aes(x = bp_cum, y = -log10(P), size = log10(NSNPS)),  # size = 8 or log10(NSNPS)
                           color = "red2") +
                ggrepel::geom_text_repel(data = genes_df %>% filter(hgnc_symbol %in% gene_set),
                                         aes(x = bp_cum, y = -log10(P), label = hgnc_symbol), 
                                         size = 5, color = "red2", force = 0.75)
            
        }
        
    }
    
    if(return_just_genes_df) {
        return(genes_df)
    } else {
        return(manhplot)
    }
    
}



# Define UI.
ui <- navbarPage("CELLECT-MAGMA ANALYSIS",
                 tabPanel("Run Analysis",
                          shinyjs::useShinyjs(),
                          sidebarLayout(fluid = FALSE,
                                        div(#style = "width: 500px;",
                                            sidebarPanel(width = 3,
                                                         div(id = "analysisParametersTitle",
                                                             h2("Select Analysis Parameters"),
                                                             hr()
                                                         ),
                                                         div(id = "analysisParametersPaths",
                                                             h3("File paths"),
                                                             textInput(inputId = "cellexPath",
                                                                       label = "Write path to CELLEX file:",
                                                                       placeholder = "Example: ~/CELLEX/outputFile.mu.csv.gz"),
                                                             textInput(inputId = "cellectPath",
                                                                       label = "Write path to CELLECT output folder:",
                                                                       placeholder = "Example: ~/CELLECT/OUTPUT/CELLECT-MAGMA"),
                                                             hr()
                                                         ),
                                                         div(id = "analysisParametersTopGenes",
                                                             h3("Top genes"),
                                                             markdown("The \"top genes\" from the analysis, by default, are selected based on being among the top 1,000 MAGMA genes (sorted by heritability value) and above the 90th percentile of CELLEX values.  
                                                                      _The user is, free to change these values to their liking, however I recommend leaving them as is._"),
                                                             numericInput(inputId = "topGenesNumber",
                                                                          label = "Number of MAGMA genes",
                                                                          value = 1000,
                                                                          min = 1,
                                                                          step = 1),
                                                             numericInput(inputId = "topGenesPercentile",
                                                                          label = "Percentile of CELLEX values",
                                                                          value = 0.9,
                                                                          min = 0.01,
                                                                          max = 1.0,
                                                                          step = 0.01),
                                                             hr()
                                                         ),
                                                         h3("(Optional) Investigate Gene Set"),
                                                         markdown("- If you want to study a set of genes more closely, you can upload it here.  
                                                                  - You can either upload a file with gene names (as a column), or paste the gene names into a text box.  
                                                                  - The tool accepts gene names in the formats of either Ensembl, HGNC, or Uniprot names."),
                                                         checkboxInput(inputId = "enableGeneSet",
                                                                       label = markdown("**Upload your own gene set to investigate**"),
                                                                       value = FALSE),
                                                         div(id = "analysisParametersGenes",
                                                             radioButtons(inputId = "inputType",
                                                                          label = "Input type:",
                                                                          choices = c("Text input" = "text",
                                                                                      "File input" = "file")),
                                                             radioButtons(inputId = "inputTypeHeader",
                                                                          label = "If the file has a header (column names):",
                                                                          choices = c("Yes" = TRUE,
                                                                                      "No" = FALSE),
                                                                          selected = TRUE),
                                                             numericInput(inputId = "inputTypeColumn",
                                                                          label = "What column the gene names are in, if stored in a column-based format (leave at default value if not applicable):",
                                                                          value = 1,
                                                                          min = 1,
                                                                          step = 1),
                                                             textAreaInput(inputId = "textInput",
                                                                           label = "Type gene names (each written on a new line, or separated by either commas, spaces or tabs):",
                                                                           placeholder = "Gene names",
                                                                           #width = "100%",
                                                                           resize = "vertical"),
                                                             fileInput(inputId = "fileInput",
                                                                       label = "Upload file with gene names:"),
                                                             selectInput(inputId = "geneFormat",
                                                                         label = "Gene/protein name format:",
                                                                         choices = c("ENSEMBL" = "ensembl",
                                                                                     "HGNC" = "gene",
                                                                                     "UNIPROT" = "uniprot")),
                                                             hr()
                                                         ),
                                                         div(id = "runAnalysisButton",
                                                             actionButton(inputId = "runAnalysis",
                                                                          label = "Run Analysis"),
                                                             actionButton(inputId = "resetAnalysis",
                                                                          label = "Reset Options")
                                                         )
                                            )
                                        ),
                                        div(#style = "width: 1000px;",
                                            mainPanel(width = 9,
                                                      hidden(
                                                          div(id = "tablesAndPlots",
                                                              tabsetPanel(
                                                                  tabPanel(title = "Plot 1 - Manhatten Plot of Genes",
                                                                           br(),
                                                                           br(),
                                                                           fluidRow(column(3,
                                                                                           selectInput(inputId = "plottab3.correct",
                                                                                                       label = "Correction method",
                                                                                                       choices = c("Bonferroni" = "bonferroni",
                                                                                                                   "FDR (False Discovery Rate)" = "fdr",
                                                                                                                   "None" = "none"),
                                                                                                       selected = "bonferroni",
                                                                                                       width = "100%")
                                                                                           ),
                                                                                    column(3,
                                                                                           numericInput(inputId = "plottab3.num",
                                                                                                        label = "Significance threshold",
                                                                                                        value = 0.05,
                                                                                                        min = 0.00,
                                                                                                        max = 1.00,
                                                                                                        step = 0.01)
                                                                                           )
                                                                                    ),
                                                                           br(),
                                                                           div(id = "plottab3",
                                                                               fluidRow(
                                                                                   column(width = 12,
                                                                                          plotOutput("plotManhatten",
                                                                                                     height = 600,
                                                                                                     click = "plot1_click",
                                                                                                     brush = brushOpts(
                                                                                                         id = "plot1_brush"
                                                                                                     )
                                                                                          )
                                                                                   )
                                                                               ),
                                                                               fluidRow(
                                                                                   column(width = 12,
                                                                                          h4("Selected points"),
                                                                                          verbatimTextOutput("brush_info_1")
                                                                                   )
                                                                               )
                                                                           )
                                                                  ),
                                                                  tabPanel(title = "Plot 2 - Heritability-Specificity",
                                                                           br(),
                                                                           br(),
                                                                           div(id = "plottab4",
                                                                               fluidRow(
                                                                                   column(width = 6,
                                                                                          plotOutput("plotHeritabilitySpecificity",
                                                                                                     height = 600,
                                                                                                     width = 600,
                                                                                                     click = "plot2_click",
                                                                                                     brush = brushOpts(
                                                                                                         id = "plot2_brush"
                                                                                                     )
                                                                                          )
                                                                                   )
                                                                               ),
                                                                               fluidRow(
                                                                                   column(width = 12,
                                                                                          h4("Selected points"),
                                                                                          verbatimTextOutput("brush_info_2")
                                                                                   )
                                                                               )
                                                                           )
                                                                  ),
                                                                  tabPanel(title = "Plot 3 - Heritability Distribution",
                                                                           br(),
                                                                           br(),
                                                                           div(id = "plottab5.3",
                                                                               markdown("### This page is only shown when a gene set is given as input.  ")
                                                                           ),
                                                                           div(id = "plottab5",
                                                                               markdown("**Mann-Whitney statistical test:**  "),
                                                                               plotOutput("plotHeritabilityDist"),
                                                                           ),
                                                                           br(),
                                                                           br(),
                                                                           div(id = "plottab5.2",
                                                                               markdown("**Mann-Whitney non-empirical statistical test:**  "),
                                                                               textOutput("w_test_pval"),
                                                                               # br(),
                                                                               # markdown("**Welch's T-test statistical test:**  "),
                                                                               # textOutput("t_test_pval"),
                                                                               # br(),
                                                                               # markdown("**Kolmogorov–Smirnov statistical test:**  "),
                                                                               # textOutput("ks_test_pval"),
                                                                               br(),
                                                                               br(),
                                                                               markdown("### Empirical results:  "),
                                                                               br(),
                                                                               markdown("**Mann-Whitney empirical statistical test:**  "),
                                                                               textOutput("w_test_pval_emp"),
                                                                               plotOutput("plotStatW")
                                                                               # br(),
                                                                               # markdown("**Welch's T-test empirical statistical test:**  "),
                                                                               # textOutput("t_test_pval_emp"),
                                                                               # plotOutput("plotStatT"),
                                                                               # br(),
                                                                               # markdown("**Kolmogorov–Smirnov empirical statistical test:**  "),
                                                                               # textOutput("ks_test_pval_emp"),
                                                                               # plotOutput("plotStatKS")
                                                                           )
                                                                  ),
                                                                  tabPanel(title = "Plot 4 - CELLEX Distribution",
                                                                           br(),
                                                                           br(),
                                                                           div(id = "plottab6.3",
                                                                               markdown("### This page is only shown when a gene set is given as input.  ")
                                                                               ),
                                                                           div(id = "plottab6",
                                                                               markdown("**Mann-Whitney statistical test:**  "),
                                                                               plotOutput("plotCellexDist")
                                                                           ),
                                                                           br(),
                                                                           br(),
                                                                           div(id = "plottab6.2",
                                                                               markdown("**Mann-Whitney non-empirical statistical test:**  "),
                                                                               textOutput("w_test_pval_cellex"),
                                                                               br(),
                                                                               br(),
                                                                               markdown("### Empirical results:  "),
                                                                               br(),
                                                                               markdown("**Mann-Whitney empirical statistical test:**  "),
                                                                               textOutput("w_test_pval_emp_cellex"),
                                                                               plotOutput("plotStatWcellex")
                                                                           )
                                                                  ),
                                                                  tabPanel(title = "Plot 5 - Cell-Type/Tissue Prioritization",
                                                                           br(),
                                                                           br(),
                                                                           div(id = "plottab1",
                                                                               plotOutput("plotPrioritization")
                                                                           )
                                                                  ),
                                                                  tabPanel(title = "Plot 6 - Heritability",
                                                                           br(),
                                                                           br(),
                                                                           div(id = "plottab2",
                                                                               plotOutput("plotHeritability")
                                                                           )
                                                                  ),
                                                                  tabPanel(title = "Table 1 - Cell-Type/Tissue Prioritization",
                                                                           br(),
                                                                           br(),
                                                                           div(id = "tabletab1",
                                                                               dataTableOutput("tablePrioritization")
                                                                           ),
                                                                           div(id = "tableTabSave1",
                                                                               shinySaveButton(id = "saveTable1",
                                                                                               label = "Save Table",
                                                                                               title = "Save file as ...", 
                                                                                               filetype = list(CSV="csv")),
                                                                           )
                                                                  ),
                                                                  tabPanel(title = "Table 2 - Top Genes",
                                                                           br(),
                                                                           br(),
                                                                           div(id = "tabletab2",
                                                                               dataTableOutput("tableTopGenes")
                                                                           ),
                                                                           div(id = "tableTabSave2",
                                                                               shinySaveButton(id = "saveTable2",
                                                                                               label = "Save Table",
                                                                                               title = "Save file as ...", 
                                                                                               filetype = list(CSV="csv")),
                                                                           )
                                                                  ),
                                                                  tabPanel(title = "Table 3 - CELLEX",
                                                                           br(),
                                                                           br(),
                                                                           div(id = "tabletab3",
                                                                               uiOutput("tabletab3Variables"),
                                                                               dataTableOutput("tableCellex")
                                                                           ),
                                                                           div(id = "tableTabSave3",
                                                                               shinySaveButton(id = "saveTable3",
                                                                                               label = "Save Table",
                                                                                               title = "Save file as ...", 
                                                                                               filetype = list(CSV="csv")),
                                                                           )
                                                                  ),
                                                                  tabPanel(title = "Table 4 - Heritability",
                                                                           br(),
                                                                           br(),
                                                                           div(id = "tabletab4",
                                                                               dataTableOutput("tableHeritability")
                                                                           ),
                                                                           div(id = "tableTabSave4",
                                                                               shinySaveButton(id = "saveTable4",
                                                                                               label = "Save Table",
                                                                                               title = "Save file as ...", 
                                                                                               filetype = list(CSV="csv")),
                                                                           )
                                                                  )
                                                              )
                                                            )
                                                      ),
                                                      div(id = "notTableAndPlots",
                                                          markdown("### Run analysis first to view output.")
                                                      ),
                                                      hidden(
                                                          div(id = "runningAnalysis",
                                                              markdown("### Running analysis - please wait."))
                                                      )
                                        )
                          )
                 )
             ),
             navbarMenu("Help",
                        tabPanel("Getting the required input"),
                        tabPanel("Running the web application"),
                        tabPanel("Explanation of figures and tables")
                    )
)



# Define server logic.
server <- function(input, output, session) {
    
    # Initialize values that need to be changed later on in the program.
    analysisVals <- reactiveValues(prioritization_df = NULL,
                                   cellex_df = NULL,
                                   cellex_df_max = NULL,
                                   heritability_df = NULL,
                                   heritability_df_new = NULL,
                                   genes_df = NULL,
                                   top_df = NULL,
                                   top_n_genes = NULL,
                                   cellex_subset = NULL,
                                   sumstats_names = NULL,
                                   annotLookup = NULL,
                                   gene_set = NULL,
                                   read_input = NULL,
                                   w_test = NULL,
                                   t_test = NULL,
                                   ks_test = NULL,
                                   w_test_emp = NULL,
                                   t_test_emp = NULL,
                                   ks_test_emp = NULL,
                                   w_test_cellex = NULL,
                                   w_test_emp_cellex = NULL)
    
    ## observe input being provided.
    shiny::observe({
        shinyjs::toggle(id = "textInput", condition = input$inputType == "text")
        shinyjs::toggle(id = "fileInput", condition = input$inputType == "file")
        shinyjs::toggle(id = "inputTypeHeader", condition = input$inputType == "file")
        shinyjs::toggle(id = "inputTypeColumn", condition = input$inputType == "file")
        shinyjs::toggle(id = "analysisParametersGenes", condition = input$enableGeneSet)
        shinyjs::toggle(id = "plottab5", condition = !is.null(analysisVals$gene_set))
        shinyjs::toggle(id = "plottab5.2", condition = !is.null(analysisVals$gene_set))
        shinyjs::toggle(id = "plottab5.3", condition = is.null(analysisVals$gene_set))
        shinyjs::toggle(id = "plottab6", condition = !is.null(analysisVals$gene_set))
        shinyjs::toggle(id = "plottab6.2", condition = !is.null(analysisVals$gene_set))
        shinyjs::toggle(id = "plottab6.3", condition = is.null(analysisVals$gene_set))
    })
    
    shiny::observeEvent(input$resetAnalysis, {
        shinyjs::reset("inputType")
        shinyjs::reset("inputTypeHeader")
        shinyjs::reset("inputTypeColumn")
        shinyjs::reset("textInput")
        shinyjs::reset("fileInput")
        shinyjs::reset("geneFormat")
        shinyjs::reset("cellexPath")
        shinyjs::reset("cellectPath")
        shinyjs::reset("enableGeneSet")
    })
    
    observeEvent(input$runAnalysis, {
        shinyjs::hide(id = "notTableAndPlots")
        shinyjs::hide(id = "tableAndPlots")
        shinyjs::show(id = "runningAnalysis")
        
        n_genes <- 1000
        magma_percentile_cutoff <- 0.9
        n_cores <- 8
        n_reps <- 1000 # 10000 or 1000
        adjust_method <- c("bonferroni", "fdr", "none")[3]
        
        redcols <- rev(brewer.pal(10, "YlOrRd"))
        pal <- colorRampPalette(redcols)
        palcols <- pal(10)
        
        
        
        if(is.null(analysisVals$annotLookup)) {
            mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl") # host = "uswest.ensembl.org", ensemblRedirect = FALSE
            analysisVals$annotLookup <- getBM(mart = mart,attributes = c('hgnc_symbol','ensembl_gene_id','uniprotswissprot'), uniqueRows = TRUE)
        }
        
        ### Parse input genes.
        if(input$enableGeneSet & str_trim(input$textInput) != "") {
            read_input <- str_trim(input$textInput)
            parsing <- parsing_fun(read_input)
            
            for(sym in parsing){
                read_input <- gsub(sym, "REPLACE", read_input)
            }
            read_input <- strsplit(read_input, "REPLACE")[[1]]
            
            for(sym in c("", " ", ",", "\n", "\t", ";", "\\|")){
                read_input <- read_input[read_input != sym]
            }
            
            analysisVals$read_input <- read_input
            
            try(
                if(input$geneFormat == "ensembl") {
                    alias_df <- analysisVals$annotLookup %>%
                        dplyr::filter(ensembl_gene_id %in% analysisVals$read_input) %>%
                        dplyr::arrange(desc(hgnc_symbol), desc(uniprotswissprot)) %>%
                        dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)
                    missing <- analysisVals$read_input[!(analysisVals$read_input %in% alias_df$ensembl_gene_id)]
                    alias_df <- alias_df %>% dplyr::add_row(data.frame(ensembl_gene_id = missing,
                                                                       hgnc_symbol = rep("", length(missing)),
                                                                       uniprotswissprot = rep("", length(missing))))
                    analysisVals$gene_set <- alias_df
                } else if(input$geneFormat == "gene") {
                    alias_df <- analysisVals$annotLookup %>%
                        dplyr::filter(hgnc_symbol %in% analysisVals$read_input) %>%
                        dplyr::arrange(desc(hgnc_symbol), desc(uniprotswissprot)) %>%
                        dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)
                    missing <- analysisVals$read_input[!(analysisVals$read_input %in% alias_df$hgnc_symbol)]
                    alias_df <- alias_df %>% dplyr::add_row(data.frame(ensembl_gene_id = rep("", length(missing)),
                                                                       hgnc_symbol = missing,
                                                                       uniprotswissprot = rep("", length(missing))))
                    analysisVals$gene_set <- alias_df
                } else if(input$geneFormat == "uniprot") {
                    alias_df <- analysisVals$annotLookup %>%
                        dplyr::filter(uniprotswissprot %in% analysisVals$read_input) %>%
                        dplyr::arrange(desc(hgnc_symbol), desc(uniprotswissprot)) %>%
                        dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)
                    missing <- analysisVals$read_input[!(analysisVals$read_input %in% alias_df$uniprotswissprot)]
                    alias_df <- alias_df %>% dplyr::add_row(data.frame(ensembl_gene_id = rep("", length(missing)),
                                                                       hgnc_symbol = rep("", length(missing)),
                                                                       uniprotswissprot = missing))
                    analysisVals$gene_set <- alias_df
                }
            )
            
            
            if(!is.null(analysisVals$gene_set)) {
                print(analysisVals$gene_set)
            }
            
            
        }
        
        ### Prioritization data and enrichment / heritability data.
        if(input$cellectPath != "" | ("results" %in% dir(input$cellectPath))) {
            
            # Load prioritization data.
            prioritization_path <- paste0(input$cellectPath, "/results/prioritization.csv")
            analysisVals$prioritization_df <- read_csv(prioritization_path)
            analysisVals$prioritization_df <- analysisVals$prioritization_df %>%
                dplyr::rename(GWAS_DATASET = gwas, 
                              EXPRESSION_DATASET = specificity_id,
                              CELLTYPE_ANNOTATION = annotation,
                              BETA = beta,
                              BETA_SE = beta_se,
                              P = pvalue) %>%
                dplyr::arrange(P)
            
            output$tablePrioritization <- renderDataTable({
                analysisVals$prioritization_df
            })
            
            # Get GWAS data names (if multiple, otherwise just a single name is returned).
            analysisVals$sumstats_names <- unique(analysisVals$prioritization_df$GWAS_DATASET) # sumstats_names
            
            # Load heritability data.
            analysisVals$heritability_df <- lapply(X = analysisVals$sumstats_names, FUN = function(sumstats_name) {
                heritability <- read_table2(paste0(input$cellectPath, "/precomputation/", sumstats_name, "/", sumstats_name, ".resid_correct_all.gsa.genes.out"), skip = 1)
                heritability <- heritability %>% add_column(P = 1-pnorm(heritability$ZSTAT))
                heritability <- heritability %>% arrange(P)
                heritability <- heritability %>% add_column(gene_percentile = sort(1:nrow(heritability), decreasing = TRUE)/nrow(heritability))
                return(heritability)
            })
            
            analysisVals$heritability_df <- bind_rows(analysisVals$heritability_df)
            
            # Modify the heritability results to only include one gene per location (multiple variants will point to the same gene).
            analysisVals$heritability_df_new <- merge(analysisVals$heritability_df, analysisVals$annotLookup, by.x = "GENE", by.y = "ensembl_gene_id")
            analysisVals$heritability_df_new <- analysisVals$heritability_df_new %>% dplyr::select(-uniprotswissprot)
            analysisVals$heritability_df_new <- analysisVals$heritability_df_new %>% dplyr::relocate(GENE, hgnc_symbol)
            analysisVals$heritability_df_new <- analysisVals$heritability_df_new %>% dplyr::arrange(GENE) %>% dplyr::distinct(hgnc_symbol, .keep_all = TRUE)
            analysisVals$heritability_df_new <- analysisVals$heritability_df_new %>% add_column(P_adjust = p.adjust(analysisVals$heritability_df_new$P, method = adjust_method))
            
            print(head(analysisVals$heritability_df_new))
            
            output$tableHeritability <- renderDataTable({
                # Use either the "old" or "new" heritability table.
                #analysisVals$heritability_df %>% dplyr::select(-gene_percentile) %>% dplyr::arrange(P)
                analysisVals$heritability_df_new %>% dplyr::select(-gene_percentile, -P_adjust) %>% dplyr::arrange(P)
            })
            
            output$plotPrioritization <- renderPlot({
                ggplot2::ggplot(data = analysisVals$prioritization_df %>% head(50), mapping = ggplot2::aes(x = reorder(CELLTYPE_ANNOTATION, P), y = -log10(P), fill = GWAS_DATASET)) +
                    ggplot2::geom_bar(stat = "identity") +
                    ggplot2::xlab("Annotation") +
                    ggplot2::ylab("-log10(p-value)") +
                    ggplot2::ggtitle(paste0("Barplot of prioritization results")) +
                    ggplot2::theme_bw() +
                    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
                                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 0, size = 12),
                                   axis.text.y = ggplot2::element_text(size = 12),
                                   axis.title.x = ggplot2::element_text(size = 16),
                                   axis.title.y = ggplot2::element_text(size = 16),
                                   legend.title=element_text(size=14),
                                   legend.text=element_text(size=12)
                    )
            },
            # width = 1000,
            # height = 600
            )
            
            # Create Heritability plot.
            output$plotHeritability <- renderPlot({
                ggplot2::ggplot(data = analysisVals$heritability_df_new %>% dplyr::arrange(P) %>% head(50), mapping = ggplot2::aes(x = reorder(hgnc_symbol, P), y = -log10(P), fill = as.factor(CHR))) +
                    ggplot2::geom_bar(stat = "identity") +
                    ggplot2::xlab("Gene") +
                    ggplot2::ylab("-log10(p-value)") +
                    ggplot2::labs(fill = "Chromosome") +
                    ggplot2::ggtitle(paste0("Barplot of gene heritability results")) +
                    ggplot2::theme_bw() +
                    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
                                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
                                   axis.text.y = ggplot2::element_text(size = 12),
                                   axis.title.x = ggplot2::element_text(size = 16),
                                   axis.title.y = ggplot2::element_text(size = 16),
                                   legend.title=element_text(size=14),
                                   legend.text=element_text(size=12)
                    )
            },
            # width = 1000,
            # height = 600
            )
            
            # Create Manhatten plot.
            output$plotManhatten <- renderPlot({
                if(is.null(analysisVals$gene_set)) {
                    manhatten.plot(analysisVals$heritability_df_new %>% dplyr::select(-P_adjust), analysisVals$annotLookup, correction_method = input$plottab3.correct, sig_thres = input$plottab3.num, og_genes_df = FALSE)
                } else if(!is.null(analysisVals$gene_set)) {
                    manhatten.plot(analysisVals$heritability_df_new %>% dplyr::select(-P_adjust), analysisVals$annotLookup, unique(analysisVals$gene_set$hgnc_symbol), correction_method = input$plottab3.correct, sig_thres = input$plottab3.num, og_genes_df = FALSE)
                }
                
            },
            # width = 1000,
            # height = 600
            )
            
            output$brush_info_1 <- renderPrint({
                if(is.null(analysisVals$gene_set)) {
                    brushedPoints(manhatten.plot(analysisVals$heritability_df_new, 
                                                 analysisVals$annotLookup,
                                                 correction_method = input$plottab3.correct,
                                                 sig_thres = input$plottab3.num,
                                                 og_genes_df = FALSE,
                                                 return_just_genes_df = TRUE) %>%
                                      dplyr::select(GENE, hgnc_symbol, CHR, START, STOP, NSNPS, ZSTAT, P, bp_cum) %>% 
                                      dplyr::mutate("-log10(P)"=-log10(P)) %>%
                                      dplyr::relocate(GENE, hgnc_symbol, P, "-log10(P)", ZSTAT, CHR, START, STOP, NSNPS, bp_cum) %>%
                                      dplyr::rename(GENE_HGNC = hgnc_symbol) %>%
                                      dplyr::arrange(P), 
                                  input$plot1_brush)
                } else if(!is.null(analysisVals$gene_set)) {
                    brushedPoints(manhatten.plot(analysisVals$heritability_df_new, 
                                                 analysisVals$annotLookup, 
                                                 unique(analysisVals$gene_set$hgnc_symbol),
                                                 correction_method = input$plottab3.correct,
                                                 sig_thres = input$plottab3.num,
                                                 og_genes_df = FALSE, 
                                                 return_just_genes_df = TRUE) %>%
                                      dplyr::select(GENE, hgnc_symbol, CHR, START, STOP, NSNPS, ZSTAT, P, bp_cum) %>% 
                                      dplyr::mutate("-log10(P)"=-log10(P)) %>%
                                      dplyr::relocate(GENE, hgnc_symbol, P, "-log10(P)", ZSTAT, CHR, START, STOP, NSNPS, bp_cum) %>%
                                      dplyr::rename(GENE_HGNC = hgnc_symbol) %>%
                                      dplyr::arrange(P), 
                                  input$plot1_brush)
                }
            })
            
            print("Manhattan plot done!")
            
            # Create distribution of heritability p-values.
            output$plotHeritabilityDist <- renderPlot({
                if(is.null(analysisVals$gene_set)) {
                    ggplot2::ggplot(data = analysisVals$heritability_df_new, mapping = ggplot2::aes(x = ZSTAT)) + # P or P_adjust or ZSTAT
                        ggplot2::geom_histogram(binwidth = 0.1) + # 0.01 for P values, 0.1 for Z-scores
                        # ggplot2::xlab("p-value") +
                        ggplot2::xlab("z-score") +
                        ggplot2::ylab("Count") +
                        ggplot2::ggtitle(paste0("Distribution of gene heritability results")) +
                        ggplot2::theme_bw() +
                        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
                                       axis.text.x = ggplot2::element_text(size = 12),
                                       axis.text.y = ggplot2::element_text(size = 12),
                                       axis.title.x = ggplot2::element_text(size = 16),
                                       axis.title.y = ggplot2::element_text(size = 16),
                                       legend.title=element_text(size=14),
                                       legend.text=element_text(size=12)
                        )
                    
                } else if(!is.null(analysisVals$gene_set)) {
                    ggplot2::ggplot(data = analysisVals$heritability_df_new, mapping = ggplot2::aes(x = ZSTAT)) + # P or P_adjust or ZSTAT
                        ggplot2::geom_histogram(binwidth = 0.1, fill = "skyblue1", color = "gray") + # 0.01 for P values, 0.1 for Z-scores
                        # ggplot2::geom_vline(xintercept = analysisVals$heritability_df_new %>%
                        #                                           dplyr::filter(hgnc_symbol %in% unique(analysisVals$gene_set$hgnc_symbol)) %>%
                        #                                           dplyr::select(ZSTAT) %>%
                        #                                           unlist(use.names = FALSE) %>%
                        #                                           abs(), color = "red") +
                        ggplot2::geom_segment(data = analysisVals$heritability_df_new %>%
                                                  dplyr::filter(hgnc_symbol %in% unique(analysisVals$gene_set$hgnc_symbol)) %>%
                                                  add_column(y_coord = rep(0, times = analysisVals$heritability_df_new %>%
                                                                               dplyr::filter(hgnc_symbol %in% unique(analysisVals$gene_set$hgnc_symbol)) %>%
                                                                               nrow())),
                                              aes(x = ZSTAT, y = y_coord, xend = ZSTAT, yend = -Inf),
                                              color = rep(palcols[1:5], times = analysisVals$heritability_df_new %>%
                                                              dplyr::filter(hgnc_symbol %in% unique(analysisVals$gene_set$hgnc_symbol)) %>%
                                                              nrow())[1:(analysisVals$heritability_df_new %>%
                                                                             dplyr::filter(hgnc_symbol %in% unique(analysisVals$gene_set$hgnc_symbol)) %>%
                                                                             nrow())]) +
                        ### Don't know whether to include this???
                        # ggplot2::geom_vline(xintercept = analysisVals$heritability_df_new %>%
                        #                       dplyr::filter(hgnc_symbol %in% unique(analysisVals$gene_set$hgnc_symbol)) %>%
                        #                       dplyr::select(ZSTAT) %>%
                        #                       unlist(use.names = FALSE) %>%
                        #                       abs() %>%
                        #                       quantile(0.95), 
                        #                     color = "darkred", linetype = "dashed") +
                        ggrepel::geom_text_repel(data = analysisVals$heritability_df_new %>%
                                                     dplyr::filter(hgnc_symbol %in% unique(analysisVals$gene_set$hgnc_symbol)),
                                                 aes(x = ZSTAT, y = 0, label = hgnc_symbol), 
                                                 color = rep(palcols[1:5], times = analysisVals$heritability_df_new %>%
                                                                 dplyr::filter(hgnc_symbol %in% unique(analysisVals$gene_set$hgnc_symbol)) %>%
                                                                 nrow())[1:(analysisVals$heritability_df_new %>%
                                                                                dplyr::filter(hgnc_symbol %in% unique(analysisVals$gene_set$hgnc_symbol)) %>%
                                                                                nrow())], size = 4) +
                        #ggplot2::xlab("p-value") +
                        ggplot2::xlab("z-score") +
                        ggplot2::ylab("Count") +
                        ggplot2::ggtitle(paste0("Distribution of gene heritability results")) +
                        ggplot2::theme_bw() +
                        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 24),
                                       axis.text.x = ggplot2::element_text(size = 16),
                                       axis.text.y = ggplot2::element_text(size = 16),
                                       axis.title.x = ggplot2::element_text(size = 20),
                                       axis.title.y = ggplot2::element_text(size = 20),
                                       legend.title=element_text(size=16),
                                       legend.text=element_text(size=14)
                        )
                }

                
                # if() {
                #     df <- analysisVals$heritability_df %>% filter()
                #     output_plot <- output_plot + 
                # }
                
                #output_plot
            },
            # width = 1000,
            # height = 600
            )
            
        }
        
        
        ### CELLEX data.
        if(input$cellexPath != "" | grepl("csv.gz$", input$cellexPath) | grepl("csv$", input$cellexPath)) {
            
            analysisVals$cellex_df <- read_csv(input$cellexPath)
            
            output$tabletab3Variables <- renderUI({
                selectInput("tab3Variables", "Variables:", choices = c("All", colnames(analysisVals$cellex_df)[-1]))
            })
            
            output$tableCellex <- renderDataTable({
                if(input$tab3Variables == "All") {
                    analysisVals$cellex_df
                } else {
                    analysisVals$cellex_df %>% dplyr::select("gene", input$tab3Variables)
                }
            })
            
            
        }
        
        if((input$cellexPath != "" | grepl("csv.gz$", input$cellexPath) | grepl("csv$", input$cellexPath)) & (input$cellectPath != "" | ("results" %in% dir(input$cellectPath))) & (input$enableGeneSet & str_trim(input$textInput) != "") ) {
            # Perform statistical significance tests (heritability).
            analysisVals$w_test <- wilcox.test(x = analysisVals$heritability_df_new %>% 
                                                   dplyr::filter(hgnc_symbol %in% unique(analysisVals$gene_set$hgnc_symbol)) %>%
                                                   dplyr::select(ZSTAT) %>%
                                                   unlist(use.names = FALSE),
                                               y = analysisVals$heritability_df_new %>% 
                                                   dplyr::filter(!(hgnc_symbol %in% unique(analysisVals$gene_set$hgnc_symbol))) %>%
                                                   dplyr::select(ZSTAT) %>%
                                                   unlist(use.names = FALSE), 
                                               alternative = "greater")
            # analysisVals$t_test <- t.test(x = analysisVals$heritability_df_new %>% 
            #                                   dplyr::filter(hgnc_symbol %in% unique(analysisVals$gene_set$hgnc_symbol)) %>%
            #                                   dplyr::select(ZSTAT) %>%
            #                                   unlist(use.names = FALSE) %>%
            #                                   abs(),
            #                               y = analysisVals$heritability_df_new %>% 
            #                                   dplyr::filter(!(hgnc_symbol %in% unique(analysisVals$gene_set$hgnc_symbol))) %>%
            #                                   dplyr::select(ZSTAT) %>%
            #                                   unlist(use.names = FALSE) %>%
            #                                   abs(), alternative = "greater")
            # analysisVals$ks_test <- ks.test(x = analysisVals$heritability_df_new %>% 
            #                                     dplyr::filter(hgnc_symbol %in% unique(analysisVals$gene_set$hgnc_symbol)) %>%
            #                                     dplyr::select(ZSTAT) %>%
            #                                     unlist(use.names = FALSE) %>%
            #                                     abs(),
            #                                 y = analysisVals$heritability_df_new %>% 
            #                                     dplyr::filter(!(hgnc_symbol %in% unique(analysisVals$gene_set$hgnc_symbol))) %>%
            #                                     dplyr::select(ZSTAT) %>%
            #                                     unlist(use.names = FALSE) %>%
            #                                     abs(), alternative = "less")
            
            output$w_test_pval <- renderPrint(
                writeLines(paste0("Statistic: ", analysisVals$w_test$statistic, "\n\nP-value: ", analysisVals$w_test$p.value))
            )
            
            # output$t_test_pval <- renderPrint(
            #     writeLines(paste0("Statistic: ", analysisVals$t_test$statistic, "\n\nP-value: ", analysisVals$t_test$p.value))
            # )
            # 
            # output$ks_test_pval <- renderPrint(
            #     writeLines(paste0("Statistic: ", analysisVals$ks_test$statistic, "\n\nP-value: ", analysisVals$ks_test$p.value))
            # )
            
            
            # Obtain dataframe with max values.
            analysisVals$cellex_df_max <- merge(analysisVals$cellex_df, analysisVals$annotLookup, by.x = "gene", by.y = "ensembl_gene_id")
            analysisVals$cellex_df_max <- analysisVals$cellex_df_max %>% dplyr::select(-uniprotswissprot)
            analysisVals$cellex_df_max <- analysisVals$cellex_df_max %>% dplyr::relocate(gene, hgnc_symbol)
            
            max_esmu <- apply(analysisVals$cellex_df[,-1], MARGIN = 1, FUN = function(x) {
                max(x)
            })
            max_esmu <- tibble(gene_names = analysisVals$cellex_df$gene, max_esmu = max_esmu)
            
            analysisVals$cellex_df_max <- analysisVals$cellex_df_max %>% 
                dplyr::select(gene, hgnc_symbol)
            analysisVals$cellex_df_max <- merge(analysisVals$cellex_df_max, max_esmu, by.x = "gene", by.y = "gene_names")
            analysisVals$cellex_df_max <- analysisVals$cellex_df_max %>% dplyr::arrange(gene) %>% dplyr::distinct(hgnc_symbol, .keep_all = TRUE)
            
            
            print("Max val df done!")
            
            
            # Perform statistical significance tests (cellex).
            analysisVals$w_test_cellex <- wilcox.test(x = analysisVals$cellex_df_max %>% 
                                                          dplyr::filter(hgnc_symbol %in% unique(analysisVals$gene_set$hgnc_symbol)) %>%
                                                          dplyr::select(max_esmu) %>%
                                                          unlist(use.names = FALSE) %>%
                                                          as.numeric(),
                                                      y = analysisVals$cellex_df_max %>% 
                                                          dplyr::filter(!(hgnc_symbol %in% unique(analysisVals$gene_set$hgnc_symbol))) %>%
                                                          dplyr::select(max_esmu) %>%
                                                          unlist(use.names = FALSE) %>%
                                                          as.numeric(), 
                                                      alternative = "greater")
            
            
            
            output$w_test_pval_cellex <- renderPrint(
                writeLines(paste0("Statistic: ", analysisVals$w_test_cellex$statistic, "\n\nP-value: ", analysisVals$w_test_cellex$p.value))
            )
            
            
            ### Perform statistical significance tests (empirical).
            set.seed(1234321)
            
            sample_list <- lapply(X = 1:n_reps, FUN = function(i) {
                sampled_names <- sample(x = analysisVals$heritability_df_new[!(analysisVals$heritability_df_new$hgnc_symbol %in% unique(analysisVals$gene_set$hgnc_symbol)),"hgnc_symbol"] %>% 
                                            unlist(use.names = FALSE) %>%
                                            unique(), 
                                        size = analysisVals$heritability_df_new[analysisVals$heritability_df_new$hgnc_symbol %in% unique(analysisVals$gene_set$hgnc_symbol),"hgnc_symbol"] %>% 
                                            unlist(use.names = FALSE) %>%
                                            unique() %>%
                                            length()
                )
                
                return(sampled_names)
            })
            print(head(sample_list))
            
            
            # W test (heritability).
            analysisVals$w_test_emp <- lapply(X = 1:n_reps, FUN = function(i){
                stat <- wilcox.test(x = analysisVals$heritability_df_new[analysisVals$heritability_df_new$hgnc_symbol %in% sample_list[[i]],"ZSTAT"] %>% unlist(use.names = FALSE),
                                    y = analysisVals$heritability_df_new[!(analysisVals$heritability_df_new$hgnc_symbol %in% sample_list[[i]]),"ZSTAT"] %>% unlist(use.names = FALSE), alternative = "greater")$statistic
                return(stat)
            })
            analysisVals$w_test_emp <- unlist(analysisVals$w_test_emp, use.names = FALSE)
            print(head(analysisVals$w_test_emp))
            
            
            
            # # T test.
            # analysisVals$t_test_emp <- lapply(X = 1:n_reps, FUN = function(i){
            #     stat <- t.test(x = analysisVals$heritability_df_new[analysisVals$heritability_df_new$hgnc_symbol %in% sample_list[[i]],"ZSTAT"] %>%
            #                             abs() %>% unlist(use.names = FALSE),
            #                         y = analysisVals$heritability_df_new[!(analysisVals$heritability_df_new$hgnc_symbol %in% sample_list[[i]]),"ZSTAT"] %>%
            #                             abs() %>% unlist(use.names = FALSE), alternative = "greater")$statistic
            #     return(stat)
            # })
            # analysisVals$t_test_emp <- unlist(analysisVals$t_test_emp, use.names = FALSE)
            # 
            # # KS test.
            # analysisVals$ks_test_emp <- lapply(X = 1:n_reps, FUN = function(i){
            #     stat <- ks.test(x = analysisVals$heritability_df_new[analysisVals$heritability_df_new$hgnc_symbol %in% sample_list[[i]],"ZSTAT"] %>%
            #                        abs() %>% unlist(use.names = FALSE),
            #                    y = analysisVals$heritability_df_new[!(analysisVals$heritability_df_new$hgnc_symbol %in% sample_list[[i]]),"ZSTAT"] %>%
            #                        abs() %>% unlist(use.names = FALSE), alternative = "less")$statistic
            #     return(stat)
            # })
            # analysisVals$ks_test_emp <- unlist(analysisVals$ks_test_emp, use.names = FALSE)
            
            # Calculate empirical p-values.
            # p_val_emp_df <- tibble(w_test = ( sum(analysisVals$w_test_emp > analysisVals$w_test$statistic) +1 ) / (n_reps+1),
            #                        t_test = ( sum(analysisVals$t_test_emp > analysisVals$t_test$statistic) +1 ) / (n_reps+1),
            #                        ks_test = ( sum(analysisVals$ks_test_emp > analysisVals$ks_test$statistic) +1 ) / (n_reps+1) )
            p_val_emp_df <- tibble(w_test = ( sum(analysisVals$w_test_emp > analysisVals$w_test$statistic) +1 ) / (n_reps+1) )
            
            output$w_test_pval_emp <- renderPrint(
                writeLines(paste0("Statistic: ", analysisVals$w_test$statistic, "\n\nP-value: ", p_val_emp_df$w_test))
            )
            
            # output$t_test_pval_emp <- renderPrint(
            #     writeLines(paste0("Statistic: ", analysisVals$t_test$statistic, "\n\nP-value: ", p_val_emp_df$t_test))
            # )
            # 
            # output$ks_test_pval_emp <- renderPrint(
            #     writeLines(paste0("Statistic: ", analysisVals$ks_test$statistic, "\n\nP-value: ", p_val_emp_df$ks_test))
            # )
            
            
            # Make plot results.
            output$plotStatW <- renderPlot({
                ggplot2::ggplot(data = tibble(stat = analysisVals$w_test_emp), mapping = ggplot2::aes(x = stat)) +
                    ggplot2::geom_histogram(fill = "skyblue1", color = "gray") + # 0.01 for P values, 0.1 for Z-scores
                    ggplot2::geom_vline(xintercept = analysisVals$w_test_emp %>% quantile(0.95), 
                                        color = "darkred", linetype = "dashed") +
                    ggplot2::geom_vline(xintercept = analysisVals$w_test$statistic, 
                                        color = "green") +
                    ggplot2::xlab("Statistic") +
                    ggplot2::ylab("Count") +
                    ggplot2::ggtitle(paste0("Distribution of gene heritability results")) +
                    ggplot2::theme_bw() +
                    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 24),
                                   axis.text.x = ggplot2::element_text(size = 16),
                                   axis.text.y = ggplot2::element_text(size = 16),
                                   axis.title.x = ggplot2::element_text(size = 20),
                                   axis.title.y = ggplot2::element_text(size = 20),
                                   legend.title=element_text(size=16),
                                   legend.text=element_text(size=14)
                    )
            },
            #width = 1000,
            #height = 600
            )
            
            # output$plotStatT <- renderPlot({
            #     ggplot2::ggplot(data = tibble(stat = analysisVals$t_test_emp), mapping = ggplot2::aes(x = stat)) +
            #         ggplot2::geom_histogram(fill = "skyblue1", color = "gray") + # 0.01 for P values, 0.1 for Z-scores
            #         ggplot2::geom_vline(xintercept = analysisVals$t_test_emp %>% quantile(0.95), 
            #                             color = "darkred", linetype = "dashed") +
            #         ggplot2::geom_vline(xintercept = analysisVals$t_test$statistic, 
            #                             color = "green") +
            #         ggplot2::xlab("Statistic") +
            #         ggplot2::ylab("Count") +
            #         ggplot2::ggtitle(paste0("Distribution of gene heritability results")) +
            #         ggplot2::theme_bw() +
            #         ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
            #                        axis.text.x = ggplot2::element_text(size = 12),
            #                        axis.text.y = ggplot2::element_text(size = 12),
            #                        axis.title.x = ggplot2::element_text(size = 16),
            #                        axis.title.y = ggplot2::element_text(size = 16),
            #                        legend.title=element_text(size=14),
            #                        legend.text=element_text(size=12)
            #         )
            # },
            # #width = 1000,
            # #height = 600
            # )
            
            # output$plotStatKS <- renderPlot({
            #     ggplot2::ggplot(data = tibble(stat = analysisVals$ks_test_emp), mapping = ggplot2::aes(x = stat)) +
            #         ggplot2::geom_histogram(fill = "skyblue1", color = "gray") + # 0.01 for P values, 0.1 for Z-scores
            #         ggplot2::geom_vline(xintercept = analysisVals$ks_test_emp %>% quantile(0.95), 
            #                             color = "darkred", linetype = "dashed") +
            #         ggplot2::geom_vline(xintercept = analysisVals$ks_test$statistic, 
            #                             color = "green") +
            #         ggplot2::xlab("Statistic") +
            #         ggplot2::ylab("Count") +
            #         ggplot2::ggtitle(paste0("Distribution of gene heritability results")) +
            #         ggplot2::theme_bw() +
            #         ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
            #                        axis.text.x = ggplot2::element_text(size = 12),
            #                        axis.text.y = ggplot2::element_text(size = 12),
            #                        axis.title.x = ggplot2::element_text(size = 16),
            #                        axis.title.y = ggplot2::element_text(size = 16),
            #                        legend.title=element_text(size=14),
            #                        legend.text=element_text(size=12)
            #         )
            # },
            # #width = 1000,
            # #height = 600
            # )
            
            
            # W test (cellex).
            analysisVals$w_test_emp_cellex <- lapply(X = 1:n_reps, FUN = function(i){
                stat <- wilcox.test(x = analysisVals$cellex_df_max[analysisVals$cellex_df_max$hgnc_symbol %in% sample_list[[i]],"max_esmu"] %>% unlist(use.names = FALSE),
                                    y = analysisVals$cellex_df_max[!(analysisVals$cellex_df_max$hgnc_symbol %in% sample_list[[i]]),"max_esmu"] %>% unlist(use.names = FALSE), alternative = "greater")$statistic
                return(stat)
            })
            analysisVals$w_test_emp_cellex <- unlist(analysisVals$w_test_emp_cellex, use.names = FALSE)
            print(head(analysisVals$w_test_emp_cellex))
            
            
            
            p_val_emp_df_cellex <- tibble(w_test = ( sum(analysisVals$w_test_emp_cellex > analysisVals$w_test_cellex$statistic) +1 ) / (n_reps+1) )
            
            output$w_test_pval_emp_cellex <- renderPrint(
                writeLines(paste0("Statistic: ", analysisVals$w_test_cellex$statistic, "\n\nP-value: ", p_val_emp_df_cellex$w_test))
            )
            
            # Make plot results.
            output$plotStatWcellex <- renderPlot({
                ggplot2::ggplot(data = tibble(stat = analysisVals$w_test_emp_cellex), mapping = ggplot2::aes(x = stat)) +
                    ggplot2::geom_histogram(fill = "skyblue1", color = "gray") + # 0.01 for P values, 0.1 for Z-scores
                    ggplot2::geom_vline(xintercept = analysisVals$w_test_emp_cellex %>% quantile(0.95), 
                                        color = "darkred", linetype = "dashed") +
                    ggplot2::geom_vline(xintercept = analysisVals$w_test_cellex$statistic, 
                                        color = "green") +
                    ggplot2::xlab("Statistic") +
                    ggplot2::ylab("Count") +
                    ggplot2::ggtitle(paste0("Distribution of gene heritability results")) +
                    ggplot2::theme_bw() +
                    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 24),
                                   axis.text.x = ggplot2::element_text(size = 16),
                                   axis.text.y = ggplot2::element_text(size = 16),
                                   axis.title.x = ggplot2::element_text(size = 20),
                                   axis.title.y = ggplot2::element_text(size = 20),
                                   legend.title=element_text(size=16),
                                   legend.text=element_text(size=14)
                    )
            },
            #width = 1000,
            #height = 600
            )
            
            # Create distribution of CELLEX ESµ values.
            output$plotCellexDist <- renderPlot({
                if(is.null(analysisVals$gene_set)) {
                    ggplot2::ggplot(data = analysisVals$cellex_df_max, mapping = ggplot2::aes(x = abs(max_esmu))) + # P or P_adjust or ZSTAT
                        ggplot2::geom_histogram(fill = "skyblue1", color = "gray") + # 0.01 for P values, 0.1 for Z-scores
                        # ggplot2::xlab("p-value") +
                        ggplot2::xlab("ESµ values") +
                        ggplot2::ylab("Count") +
                        ggplot2::ggtitle(paste0("Distribution of max per-gene ESµ values")) +
                        ggplot2::theme_bw() +
                        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 24),
                                       axis.text.x = ggplot2::element_text(size = 16),
                                       axis.text.y = ggplot2::element_text(size = 16),
                                       axis.title.x = ggplot2::element_text(size = 20),
                                       axis.title.y = ggplot2::element_text(size = 20),
                                       legend.title=element_text(size=16),
                                       legend.text=element_text(size=14)
                        )
                    
                } else if(!is.null(analysisVals$gene_set)) {
                    ggplot2::ggplot(data = analysisVals$cellex_df_max, mapping = ggplot2::aes(x = max_esmu)) + # P or P_adjust or ZSTAT
                        ggplot2::geom_histogram(fill = "skyblue1", color = "gray") + # 0.01 for P values, 0.1 for Z-scores
                        # ggplot2::geom_vline(xintercept = analysisVals$heritability_df_new %>%
                        #                                           dplyr::filter(hgnc_symbol %in% unique(analysisVals$gene_set$hgnc_symbol)) %>%
                        #                                           dplyr::select(ZSTAT) %>%
                        #                                           unlist(use.names = FALSE) %>%
                        #                                           abs(), color = "red") +
                        ggplot2::geom_segment(data = analysisVals$cellex_df_max %>%
                                                  dplyr::filter(hgnc_symbol %in% unique(analysisVals$gene_set$hgnc_symbol)) %>%
                                                  add_column(y_coord = rep(0, times = analysisVals$cellex_df_max %>%
                                                                               dplyr::filter(hgnc_symbol %in% unique(analysisVals$gene_set$hgnc_symbol)) %>%
                                                                               nrow())),
                                              aes(x = max_esmu, y = y_coord, xend = max_esmu, yend = -Inf),
                                              color = rep(palcols[1:5], times = analysisVals$cellex_df_max %>%
                                                              dplyr::filter(hgnc_symbol %in% unique(analysisVals$gene_set$hgnc_symbol)) %>%
                                                              nrow())[1:(analysisVals$cellex_df_max %>%
                                                                             dplyr::filter(hgnc_symbol %in% unique(analysisVals$gene_set$hgnc_symbol)) %>%
                                                                             nrow())]) +
                        # ggplot2::geom_vline(xintercept = analysisVals$cellex_df_max %>%
                        #                         dplyr::filter(hgnc_symbol %in% unique(analysisVals$gene_set$hgnc_symbol)) %>%
                        #                         dplyr::select(max_esmu) %>%
                        #                         unlist(use.names = FALSE) %>%
                        #                         quantile(0.95),
                        #                     color = "darkred", linetype = "dashed") +
                        ggrepel::geom_text_repel(data = analysisVals$cellex_df_max %>%
                                                     dplyr::filter(hgnc_symbol %in% unique(analysisVals$gene_set$hgnc_symbol)),
                                                 aes(x = max_esmu, y = 0, label = hgnc_symbol),
                                                 color = rep(palcols[1:5], times = analysisVals$cellex_df_max %>%
                                                                 dplyr::filter(hgnc_symbol %in% unique(analysisVals$gene_set$hgnc_symbol)) %>%
                                                                 nrow())[1:(analysisVals$cellex_df_max %>%
                                                                                dplyr::filter(hgnc_symbol %in% unique(analysisVals$gene_set$hgnc_symbol)) %>%
                                                                                nrow())], size = 4) +
                        #ggplot2::xlab("p-value") +
                        ggplot2::xlab("ESµ values") +
                        ggplot2::ylab("Count") +
                        ggplot2::ggtitle(paste0("Distribution of max per-gene ESµ values")) +
                        ggplot2::theme_bw() +
                        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 24),
                                       axis.text.x = ggplot2::element_text(size = 16),
                                       axis.text.y = ggplot2::element_text(size = 16),
                                       axis.title.x = ggplot2::element_text(size = 20),
                                       axis.title.y = ggplot2::element_text(size = 20),
                                       legend.title=element_text(size=16),
                                       legend.text=element_text(size=14)
                        )
                }
                
                
                # if() {
                #     df <- analysisVals$heritability_df %>% filter()
                #     output_plot <- output_plot +
                # }
                
                #output_plot
            },
            # width = 1000,
            # height = 600
            )
        }
        
        
        ### Calculate the top genes. # Perhaps change the if statement to look for no NULL values in the output.
        if(input$cellectPath != "" | ("precomputation" %in% dir(input$cellectPath)) & (input$cellexPath != "" | grepl("csv.gz$", input$cellexPath) | grepl("csv$", input$cellexPath))) {
            sumstats_name <- str_subset(dir(paste0(input$cellectPath, "/precomputation/")), "sumstats")
            sumstats_name <- substr(sumstats_name, start = 1, stop = nchar(sumstats_name) - 9)
            
            analysisVals$genes_output_df <- get_top_genes(analysisVals$heritability_df, analysisVals$cellex_df, analysisVals$annotLookup)
            analysisVals$top_n_genes <- head(analysisVals$heritability_df$GENE, 1000) # n_genes instead of 1000
            analysisVals$top_df <- analysisVals$genes_output_df %>% dplyr::filter(GENE %in% analysisVals$top_n_genes, CELLEX_PERCENTILE > 0.9)

            # Insert an if statement here regarding if custom gene set is uploaded to the tool.
            output$tableTopGenes <- renderDataTable({
                analysisVals$top_df %>% arrange(-CELLEX_PERCENTILE, -GENE_PERCENTILE)
            })
        }
        
        
        ### Create the big ESµ vs -log10(p-val) plot, using the top genes from earlier (also add option to display custom gene set).
        if(input$cellectPath != "" | ("results" %in% dir(input$cellectPath)) & (input$cellexPath != "" | grepl("csv.gz$", input$cellexPath) | grepl("csv$", input$cellexPath))) {
            output$plotHeritabilitySpecificity <- renderPlot({
                if(is.null(analysisVals$gene_set)) {
                    ggplot2::ggplot(data = analysisVals$genes_output_df %>%
                                        dplyr::arrange(GENE_HGNC, -CELLEX_VALUE) %>%
                                        dplyr::distinct(GENE_HGNC, .keep_all = TRUE),
                                    mapping = ggplot2::aes(x = -log10(P), y = CELLEX_VALUE)) + # analysisVals$genes_output_df or analysisVals$top_df
                        ggplot2::geom_point(alpha = 0.75, color = "skyblue1") +
                        ggplot2::geom_hline(yintercept = min(analysisVals$top_df$CELLEX_VALUE),
                                            color = "darkred",
                                            linetype = "dashed") +
                        ggplot2::geom_vline(xintercept = min(-log10(analysisVals$top_df$P)),
                                            color = "darkred",
                                            linetype = "dashed") +
                        ggplot2::xlab("-log10(p-value)") +
                        ggplot2::ylab("ESµ value") +
                        ggplot2::ggtitle(paste0("Plot of heritability and specificity")) +
                        ggplot2::theme_bw() +
                        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 24),
                                       axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 0, size = 16),
                                       axis.text.y = ggplot2::element_text(size = 16),
                                       axis.title.x = ggplot2::element_text(size = 20),
                                       axis.title.y = ggplot2::element_text(size = 20),
                                       legend.title=element_text(size=16),
                                       legend.text=element_text(size=14)
                        )
                } else if(!is.null(analysisVals$gene_set)) {
                    ggplot2::ggplot(data = analysisVals$genes_output_df %>%
                                        dplyr::arrange(GENE_HGNC, -CELLEX_VALUE) %>%
                                        dplyr::distinct(GENE_HGNC, .keep_all = TRUE), 
                                    mapping = ggplot2::aes(x = -log10(P), y = CELLEX_VALUE)) + # analysisVals$genes_output_df or analysisVals$top_df
                        ggplot2::geom_point(alpha = 0.75, color = "skyblue1") +
                        ggplot2::geom_hline(yintercept = min(analysisVals$top_df$CELLEX_VALUE),
                                            color = "darkred",
                                            linetype = "dashed") +
                        ggplot2::geom_vline(xintercept = -log10(max(analysisVals$top_df$P)),
                                            color = "darkred",
                                            linetype = "dashed") +
                        ggplot2::geom_point(data = analysisVals$genes_output_df %>% 
                                                filter(GENE_HGNC %in% unique(analysisVals$gene_set$hgnc_symbol)) %>%
                                                dplyr::arrange(GENE_HGNC, -CELLEX_VALUE) %>%
                                                dplyr::distinct(GENE_HGNC, .keep_all = TRUE), # analysisVals$genes_output_df or analysisVals$top_df
                                           aes(x = -log10(P), y = CELLEX_VALUE),
                                           color = "red2", size = 2) +
                        ggrepel::geom_text_repel(data = analysisVals$genes_output_df %>% 
                                                     filter(GENE_HGNC %in% unique(analysisVals$gene_set$hgnc_symbol)) %>%
                                                     dplyr::arrange(GENE_HGNC, -CELLEX_VALUE) %>%
                                                     dplyr::distinct(GENE_HGNC, .keep_all = TRUE), # analysisVals$genes_output_df or analysisVals$top_df
                                                 aes(x = -log10(P), y = CELLEX_VALUE, label = GENE_HGNC), 
                                                 color = "red2", size = 4) +
                         ggplot2::xlab("-log10(p-value)") +
                         ggplot2::ylab("ESµ value") +
                         ggplot2::ggtitle(paste0("Plot of heritability and specificity")) +
                         ggplot2::theme_bw() +
                         ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 24),
                                        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 0, size = 16),
                                        axis.text.y = ggplot2::element_text(size = 16),
                                        axis.title.x = ggplot2::element_text(size = 20),
                                        axis.title.y = ggplot2::element_text(size = 20),
                                        legend.title=element_text(size=16),
                                        legend.text=element_text(size=14)
                        )
                }

            })
            output$brush_info_2 <- renderPrint({
                brushedPoints(analysisVals$genes_output_df %>%
                                  dplyr::arrange(P, -CELLEX_VALUE) %>%
                                  dplyr::distinct(GENE_HGNC, .keep_all = TRUE) %>%
                                  dplyr::mutate("-log10(P)" = -log10(P)) %>%
                                  dplyr::relocate("-log10(P)", .after = P),
                              input$plot2_brush)
            })
        }
        
        
        
        shinyjs::hide(id = "runningAnalysis")
        shinyjs::show(id = "tablesAndPlots")
        
    })
    
    
    # "Observe events" for saving the tables as .csv files.
    observe({
        volumes <- c("UserFolder"="~")
        shinyFileSave(input, "saveTable1", roots=volumes, session=session)
        fileinfo <- parseSavePath(volumes, input$saveTable1)
        if (nrow(fileinfo) > 0) {
            write.csv(analysisVals$prioritization_df, as.character(fileinfo$datapath))
        }
    })
    
    observe({
        volumes <- c("UserFolder"="~")
        shinyFileSave(input, "saveTable2", roots=volumes, session=session)
        fileinfo <- parseSavePath(volumes, input$saveTable2)
        if (nrow(fileinfo) > 0) {
            write.csv(analysisVals$top_df %>% arrange(-CELLEX_PERCENTILE, -GENE_PERCENTILE), as.character(fileinfo$datapath))
        }
    })
    
    observe({
        volumes <- c("UserFolder"="~")
        shinyFileSave(input, "saveTable3", roots=volumes, session=session)
        fileinfo <- parseSavePath(volumes, input$saveTable3)
        if (nrow(fileinfo) > 0) {
            write.csv(analysisVals$cellex_df, as.character(fileinfo$datapath))
        }
    })
    
    observe({
        volumes <- c("UserFolder"="~")
        shinyFileSave(input, "saveTable4", roots=volumes, session=session)
        fileinfo <- parseSavePath(volumes, input$saveTable4)
        if (nrow(fileinfo) > 0) {
            write.csv(analysisVals$heritability_df_new %>% dplyr::select(-gene_percentile, -P_adjust) %>% dplyr::arrange(P), as.character(fileinfo$datapath))
        }
    })
    
    
    
}

# Run the application.
shinyApp(ui = ui, server = server)

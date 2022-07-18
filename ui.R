library(shinyWidgets)
library(shinythemes)
library(bslib)
library(DT)
library(shiny)
library(stats)
library(gplots)
library(graphics)
library(viridis)
library(utils)
library(shinycssloaders)
library(shinybusy)
library(shinyjs)
library(shinyalert)
library(stringi)
library(stringr)

ui <- fluidPage(
  useShinyjs(),
  add_busy_spinner(spin = "dots", position = "bottom-right", color = "#3E3F3A"),
  
  titlePanel(title=div(img(src="ODClogo.png", height = 80), "OutDeCo")),
  theme = bs_theme(version = 5, bootswatch = "sandstone", 
                  heading_font = font_google("Poppins"), 
                  base_font = font_collection(font_google("Roboto")),
                  success = "#325D88"),

  tags$style(HTML(".js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {background: #3E3F3A}")),
  #navbarPage is top menu bar
  navbarPage(title=NULL, collapsible = FALSE,

            ##################### HOME TAB #####################
            tabPanel(
              title="Home",
              icon = icon("home"),

              navlistPanel(
                id = "Header", selected = NULL, well = FALSE, fluid = TRUE, widths = c(3, 9),

                # about
                tabPanel(title="About",
                  h3("What is OutDeCo?"),
                  p("Outlier detection through co-expression. Using the tools on this website you can:"),
                  p("- Run a differential expression (DE) analysis "),
                  p("- Assess your DE results using gene co-expression properties"),
                  p("- Report a functional outlier assessment"),
                  p("- Run a network connectivity analysis of DE results within a gene co-expression network"),
                ),
                
                # differential analysis
                tabPanel(title="Differential Analysis",
                  h3("Differential Expression Analysis"),
                  p("Statistical analysis to discover quantitative changes in expression levels between experimental groups."),
                  h5("Methods:"),
                  h6(strong("wilcox")),
                  p("Compares means of two groups to analyse count data and test for differential expression."),
                  h6(strong("DESeq")),
                  p("Uses geometric normalisation to analyse count data and test for differential expression."),
                  h6(strong("edgeR")),
                  p("Uses weighted mean of log ratio to analyse count data and test for differential expression."),
                  br(),
                  em("Note: This tool is under construction")
                   
                ),

                # cluster genes
                tabPanel(title="Cluster Genes",
                  h3("Cluster Genes"),
                  p('Creates modules which are clusters of genes that are hightly co-expressed'),
                  h5("Plot Types"),
                  h6(strong("Heatmap")),  
                  p("up-regulated or down-regulated heatmap of genes"),
                  img(src="plot_coexpression_heatmap_down.png", height = 200),
                  h6(strong("Network")),  
                  p("up regulated or down-regulated network plot where nodes are genes and the weight of the edges corresponds to the co-expression levels between its endpoints"),
                  img(src="plot_network_down.png", height = 200),
                  h6(strong("Binarized heatmap")),  
                  p("up-regulated or down-regulated binary co-expression sub-network"),
                  br(),
                  em("Note: This tool is under construction")
                 ),

                # gene connectivity 
                tabPanel(title="Gene Connectivity",
                  h3("Gene Connectivity"),
                  p('Calculates node degrees to get a sense of the global and local connectivities of the gene'),
                  h5("Plot Types"),
                  h6(strong("Density")),  
                  p("up-regulated or down-regulated density plot of genes"),
                  img(src="plot_scatter_density_up.png", height = 200), 
                  h6(strong("Histogram")),
                  p("up-regulated or down-regulated histogram of genes"),
                  img(src="plot_scatter_hist_up.png", height = 200),
                  h6(strong("Subset by clusters")),  
                  img(src="plot_scatter_hist_down_colored.png", height = 200), 
                  br(),
                  em("Note: This tool is under construction"),
                   
                ),
                 
                # functional outliers
                tabPanel(title="Functional Outliers",
                  h3("Functional Outliers"),
                  p("Functional outliers are genes that have been identified to be potentially dysregulated. 
                  They are the genes that are Differentially Expressed but do not show local co-expression"),
                  p("Module Default Threshold: More than 6 genes"),
                  h5("Analysis Options"),
                  h6(strong("Coexpression Heatmap")),
                  p("up-regulated or down-regulated heatmap of genes detailing the outliers connectivity and expression"),
                  img(src="plot_coexpression_heatmap_down_filt.png", height = 200),
                  h6(strong("Network")),
                  p("up-regulated or down-regulated subnetwork plot detailing the outliers connectivity and expression"),
                  img(src="plot_network_down.png", height = 200),
                  h6(strong("Table")),
                  p("Genes Filtered and Genes Remaining"),
                  br(),
                  em("Note: This tool is under construction"),
                ),
                
                # GSEA
                tabPanel(
                  title="Gene Set Enrichment Analysis",
                  h3("Gene Set Enrichment Analysis (GSEA)"),
                  p("GSEA is a process of ranking genes by how statistically significant their differential gene expression is. 
                  This can remove false positives from the data. "),
                  h5("Analysis Options"),
                  h6(strong("Overlap")),
                  p("A map detailing the overlap"),
                  img(src="go_enrich.png", height = 200),
                  h6(strong("Ranking")),
                  img(src="go_enrich_ranked.png", height = 200),
                  br(),
                  em("Note: This tool is under construction")
                 ),
                 
              )
            ),

            ##################### RUN DE TAB #####################
            tabPanel(
              title="Run DE", 

              # side panel for upload options
              dropdown(
                tags$h3("Options"), 
                # side panel characteristics
                style = "jelly", icon = "OPTIONS",
                status = "primary", width = "300px", size = "sm",
               
                # title of sidepanel
                fluidPage(
                # inputs in the sidepanel
                  selectInput(
                    inputId = "DE_method",
                    label= "Choose DE Method",
                    choices = c("wilcox", "DESeq2", "edgeR"),
                    selected = NULL,
                  ),

                  fileInput("counts_file", label = "Upload counts"),

                  fileInput("labels_file", "Choose Labels File",
                    accept = c(
                    "text/csv",
                    "text/comma-separated-values,text/plain",
                    ".csv")
                  ),

                
                ), 
               ),


               br(),
               navlistPanel(
                widths = c(3, 9), well = FALSE,
                tabPanel(
                  title="View File",
                  uiOutput("UILabelsContent"),
                 ),
                tabPanel(
                  title="Plot DE",
                  dropdown(
                    # side panel characteristics
                    style = "minimal", icon = "OPTIONS",
                    status = "primary", width = "600px", size = "sm",
                  
                    h5(id = "case_selection", "Case/Control Selection"),
                    
                    selectInput(
                      inputId="select_column",
                      label= "Select label to group ",
                      choices = NULL #no choice before uploading
                    ),
                
                    selectInput(
                      inputId="select_case",
                      label= "Select case to analyse",
                      choices = NULL #no choice before column selected
                    ),
                    actionButton(inputId="run_DE", label = "Run DE"),
                  
                  ),

                    splitLayout(cellWidths = c("50%", "50%"), 
                    fluidPage(
                      textOutput("DE_V_text"),
                      plotOutput(outputId = "DEplot", height = "450px"), 
                    ),
                    fluidPage(
                      textOutput("DE_MA_text"),
                      plotOutput(outputId = "DEplot_average", height = "450px")
                    )
                          
                    )
                  
                 
                 ),
                 
              ),
            ),
            
            ##################### ASSESS DE TAB #####################
            tabPanel(
              title="Assess DE", 
              
              # options dropdown
              dropdown(

                # network selection
                tags$h4("Network Selection"),
                selectInput(
                  inputId = "network_type",
                  label = NULL,
                  choices = c("Blood", "Brain", "Generic"),
                  selected = "Generic"
                ),

                # gene list selection
                radioButtons(
                  inputId = "gene_list_selection",
                  label = tags$h4("Gene List Selection"),
                  choices = c("Upload Gene List", "Generate Gene List"),
                  selected = ""
                ),

                # generate gene list
                conditionalPanel(
                  condition = "input.gene_list_selection == 'Generate Gene List'", 
                  textInput(
                    inputId = 'chooseChrome', 
                    label = 'Choose Chromosome' , 
                    placeholder = "chrX"
                  ),
                  textInput(
                    inputId = 'chooseGeneNo', 
                    label = 'Choose Number of Genes',
                    placeholder = "100"
                  ),
                ),

                # gene list upload
                conditionalPanel(
                  condition = "input.gene_list_selection == 'Upload Gene List'", 
                  # upload file
                  fileInput(
                    inputId = "DEFile", 
                    label = "Choose Gene List File",
                    accept = c(".csv", ".tsv", ".txt")
                  ),
                  # div(style = "margin-top: -25px"),
                  # select delimiter (default is nothing until file is selected and handled in server side)
                  radioButtons(
                    inputId = 'sepButton', 
                    label = 'Delimiter Selector', 
                    choices = c(Default=''), 
                    selected = ''
                  ),
                ),

                
                
                # generate subnet button
                actionButton("generate_subnet", "Generate Subnetwork", 
                style="color: #fff; background-color: #3E3F3A; border-color: #20201F"),
      
                # side panel characteristics
                style = "jelly", icon = "OPTIONS",
                status = "primary", width = "300px", size = "sm",

               ),
              
              
              br(),
        
              navlistPanel(
                widths = c(3, 9), well = FALSE,
                
                # VIEW FILES
                tabPanel(
                  title="View Files",
                  tabsetPanel(
                    
                    # view subnetwork tab
                    tabPanel(
                      title="Subnetwork", 
                      br(),
                      tableOutput("subnetwork")
                    ),

                    # view file tab
                    tabPanel(
                      title="File",
                      br(),
                      uiOutput("UIDEContent")
                    ),
                  ),
                ),

                # CLUSTER GENES
                tabPanel(
                  title="Cluster Genes",
                  mainPanel(
                    h3("Cluster Genes"),

                    # options dropdown
                    dropdown(
                      inputId = "CG_dropdown",
                      style = "minimal", icon = "OPTIONS",
                      status = "primary", width = "300px", size = "sm",

                      # select plots
                      awesomeCheckboxGroup(
                        inputId = "clusterPlotOptions",
                        label = tags$h4("Select Plots"), 
                        choices = c("Network", "Heatmap", "Binarized Heatmap"),
                        status = ""
                      ),

                      # run button
                      actionButton(inputId = "CG_run", label = "Run",

                      # dropdown characteristics
                      style="color: #fff; background-color: #3E3F3A; border-color: #20201F"),
                    ),  

                    br(),

                    # error message
                    textOutput("CG_error"),
                    
                    tabsetPanel(

                      # plots tab
                      tabPanel(
                        br(),
                        title="Plots",
                        br(),

                        # network
                        conditionalPanel(
                          condition = "$.inArray('Network', input.clusterPlotOptions) > -1", 
                          h4("Network of Clustered Genes"), 
                          br(),
                          plotOutput(outputId = "network", height = "500px"),
                        ),
                        
                        # heatmap
                        conditionalPanel(
                          condition = "$.inArray('Heatmap', input.clusterPlotOptions) > -1", 
                          h4("Heatmap of Clustered Genes"),
                          br(),
                          plotOutput(outputId = "heatmap", height = "500px"),
                        ),

                        # binarized heatmap
                        conditionalPanel(
                          condition = "$.inArray('Binarized Heatmap', input.clusterPlotOptions) > -1", 
                          h4("Binarized Heatmap of Clustered Genes"), 
                          br(),
                          plotOutput(outputId = "Bheatmap", height = "500px"), 
                        ),
                      ),

                      # tables tab
                      tabPanel(
                        br(),
                        title="Tables",
                        br(),
                        
                        # clustering genes
                        h4("Clustering Genes"), 
                        br(),
                        fluidRow(
                          column( 11,
                                  dataTableOutput("CG_table"),
                          )
                        ),
                      )
                      
                    )

                  )
                ),
                
                # GENE CONNECTIVITY
                tabPanel(
                  title="Gene Connectivity",

                  mainPanel(
                    h3("Gene Connectivity"),
                    
                    # options dropdown
                    dropdown(
                      inputId = "GC_dropdown",
                      style = "minimal", icon = "OPTIONS",
                      status = "primary", width = "300px", size = "sm",

                      # select plots
                      awesomeCheckboxGroup(
                        inputId = "GCPlotOptions",
                        label = tags$h4("Select Plots"), 
                        choices = c("Density", "Histogram", "Clustered Density", "Clustered Histogram"),
                        status = ""
                      ),
                      
                      # filt_min slider
                      conditionalPanel(
                        condition = "$.inArray('Histogram', input.GCPlotOptions) > -1 || $.inArray('Clustered Histogram', input.GCPlotOptions) > -1" ,
                        sliderInput(
                          "xybreaks", 
                          label = "Number of breaks for histogram:",
                          min = 10, max = 150, value = 100, step = 10,
                        ),
                      ),
                      
                      br(),

                      # run button
                      actionButton(inputId = "GC_run", label = "Run", 

                      # dropdown characteristics
                      style="color: #fff; background-color: #3E3F3A; border-color: #20201F"),

                    ),
                    
                    br(),

                    # error message
                    textOutput("GC_error"),

                    # density
                    conditionalPanel(
                      br(),
                      condition = "$.inArray('Density', input.GCPlotOptions) > -1", 
                      h4("Density Plot of Gene Connectivity"), 
                      br(),
                      plotOutput(outputId = "GCdensityG", height = "500px",),
                      br(),
                    ),

                    # histogram
                    conditionalPanel(
                      br(),
                      condition = "$.inArray('Histogram', input.GCPlotOptions) > -1", 
                      h4("Histogram of Gene Connectivity"),
                      br(),
                      plotOutput(outputId = "GChistogramG", height = "500px",),
                      br(),
                    ),

                    # density (subset by clusters)
                    conditionalPanel(
                      br(),
                      condition = "$.inArray('Clustered Density', input.GCPlotOptions) > -1", 
                      h4("Density plot of Gene Connectivity subset by their clusters"), 
                      br(),
                      plotOutput(outputId = "GCdensitySubsetG", height = "500px",),
                      br(),
                    ),

                    # histogram (subset by clusters)
                    conditionalPanel(
                      br(),
                      condition = "$.inArray('Clustered Histogram', input.GCPlotOptions) > -1", 
                      h4("Histogram of Gene Connectivity subset by their clusters"), 
                      br(),
                      plotOutput(outputId = "GChistogramSubsetG", height = "500px",),
                      br(),
                    ),
                  )
                ),
                
                # FUNCTIONAL OUTLIERS
                tabPanel(
                  title="Functional Outliers",
                  
                  mainPanel(
                    h3("Functional Outliers"),

                    # options dropdown
                    dropdown(
                      inputId = "FO_dropdown",
                      style = "minimal", icon = "OPTIONS",
                      status = "primary", width = "300px", size = "sm",

                      # select plots
                      awesomeCheckboxGroup(
                        inputId = "FOPlotOptions",
                        label = tags$h4("Select Plots"),
                        choices = c("Network", "Heatmap"),
                        status = ""
                      ),

                      # select tables
                      awesomeCheckboxGroup(
                          inputId = "FO_table_options",
                          label = tags$h4("Select Tables"),
                          choices = c("Functional Outliers", "Genes in Module"),
                          status = ""
                      ),

                      # other options
                      tags$h4("Other"),
                      
                      # filt_min slider
                      sliderInput("filtmin", label = "Number of Genes to form Module",
                          min = 0, max = 20, value = 6, step = 1
                      ),

                      br(),
                      
                      # run button
                      actionButton(inputId = "FO_run", label = "Run", 

                      # dropdown characteristics
                      style="color: #fff; background-color: #3E3F3A; border-color: #20201F"),

                    ),

                    br(),
                    
                    # error message
                    textOutput("FO_error"),

                  ),
                  br(),
                  
                  tabsetPanel(

                    # plots tab
                    tabPanel(
                      br(),
                      title="Plots",
                      br(),
                      
                      # heatmap
                      conditionalPanel(
                        condition = "$.inArray('Network', input.FOPlotOptions) > -1", 
                        h4("Network"), 
                        plotOutput(outputId = "FO_network", height = "500px"),
                      ),

                      # network
                      conditionalPanel(
                        condition = "$.inArray('Heatmap', input.FOPlotOptions) > -1", 
                        h4("Heatmap"), 
                        plotOutput(outputId = "FO_heatmap", height = "500px"),
                        br(),
                      ),
                    ),

                    # tables tab
                    tabPanel(
                      br(),
                      title="Tables", 
                      br(),

                      # selected genes table output
                      conditionalPanel(
                        condition = "$.inArray('Genes in Module', input.FO_table_options) > -1", 
                        h4("Genes in Module"), 
                        br(),
                        fluidRow(
                          column( 11,
                                  dataTableOutput("genes_not_keep_table"),
                          )
                        ),
                      ),

                      br(),

                      # unselected genes table output
                      conditionalPanel(
                        condition = "$.inArray('Functional Outliers', input.FO_table_options) > -1", 
                        h4("Outliers"), 
                        br(),
                        fluidRow(
                          column( 11,
                                  dataTableOutput("genes_keep_table"),
                          )
                        ),
                      ),
                    )
                  ),
                ),

                # GENE SET ENRICHMENT ANALYSIS
                tabPanel(
                  title="Gene Set Enrichment Analysis",
                  mainPanel(
                    h3("Gene Set Enrichment Analysis"),
                    
                    # options dropdown
                    dropdown(
                      inputId = "GSEA_dropdown",
                      style = "minimal", icon = "OPTIONS",
                      status = "primary", width = "300px", size = "sm",
                      
                      # select GSEA type
                      radioButtons(
                        inputId = "GSEA_type",
                        label = tags$h4("GSEA Type"),
                        choices = c("Standard GSEA", "AUCs GSEA"),
                        selected = ""
                      ),
                      
                      # standard GSEA options
                      conditionalPanel(
                        condition = "input.GSEA_type == 'Standard GSEA'", 
                        awesomeCheckboxGroup(
                          inputId = "GSEA_std_PlotOptions",
                          label = tags$h4("Select Plots"), 
                          choices = c("P-value Heatmap"),
                          status = ""
                        ),
                        br(),
                        # run button
                        actionButton(
                          inputId = "GSEA_std_run",
                          label = "Run Standard GSEA",
                          style="color: #fff; background-color: #3E3F3A; border-color: #20201F"
                        ),
                      ),

                      # AUCs GSEA options
                      conditionalPanel(
                        condition = "input.GSEA_type == 'AUCs GSEA'", 
                        awesomeCheckboxGroup(
                          inputId = "GSEA_auc_PlotOptions",
                          label = tags$h4("Select Plots"), 
                          choices = c("P-value Heatmap", "AUROC Graph", "AUC Histogram", "Gene Set Size Histogram", "Scatter Plot (AUCs vs P-values)"),
                          status = ""
                        ),
                        br(),
                        # run button
                        actionButton(
                          inputId = "GSEA_auc_run",
                          label = "Run AUCs GSEA",
                          style="color: #fff; background-color: #3E3F3A; border-color: #20201F"
                        ),

                      ),
                      
                      # # filt_min slider
                      # conditionalPanel(
                      #   condition = "$.inArray('Histogram', input.GCPlotOptions) > -1 || $.inArray('Clustered Histogram', input.GCPlotOptions) > -1" ,
                      #   sliderInput(
                      #     "xybreaks", 
                      #     label = "Number of breaks for histogram:",
                      #     min = 10, max = 150, value = 100, step = 10,
                      #   ),
                      # ),
                      
                      br(),

                    ),
                    
                    br(),

                    # error message
                    textOutput("GSEA_error"),
                  ),
                  br(),
                  
                  # heatmap
                  conditionalPanel(
                    condition = "$.inArray('P-value Heatmap', input.GSEA_std_PlotOptions) > -1", 
                    h4("P-value Heatmap"), 
                    plotOutput(outputId = "GSEA_std_heatmap", height = "500px"),
                  ),

                  # AUROC Graph
                  conditionalPanel(
                    condition = "$.inArray('AUROC Graph', input.GSEA_auc_PlotOptions) > -1", 
                    h4("AUROC Graph"), 
                    plotOutput(outputId = "GSEA_AUROC", height = "500px"),
                    br(),
                  ),
                  
                ),

              
              ),
            )
  )
)
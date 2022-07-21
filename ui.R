library(OutDeCo)
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
library(EGAD)

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
  navbarPage(title=NULL, id="navpage", collapsible = FALSE,

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
                  fileInput("counts_file", label = "Choose Counts Data File"),
                  
                  radioButtons(
                    inputId = 'sepCountsButton', 
                    label = 'Delimiter Selector', 
                    choices = c(Comma=",", Semicolon=";", Tab="\t", Space=" "), 
                    selected = ''
                  ),

                  fileInput("labels_file", "Choose Labels File",
                    accept = c(
                    "text/csv",
                    "text/comma-separated-values,text/plain",
                    ".csv")
                  ),

                  radioButtons(
                      inputId = 'sepLabelsButton', 
                      label = 'Delimiter Selector', 
                      choices = c(Comma=",", Semicolon=";", Tab="\t", Space=" "),
                      selected = ''
                  ),

                ), 
              ),

              br(),

              navlistPanel(
                widths = c(3, 9), well = FALSE,

                # VIEW FILES
                tabPanel(
                  title="View File",

                  tabsetPanel(
                    id="counts_labels_tabset",
                    # Counts File tab
                    tabPanel(
                      title="Counts File",
                      dataTableOutput("UICountsContent")
                    ),
                    # Labels File tab
                    tabPanel(
                      title="Labels File",
                      dataTableOutput("UILabelContent")
                    ),
                  )
                ),

                # PLOT DE
                tabPanel(
                  title="Plot DE",
                  h4("Plot Differential Expression"),
                  p(id = "runDE_error", "Please upload counts and labels data in OPTIONS"),
                  
                  
                  dropdown(
                    inputId = "DE_options",
                    # side panel characteristics
                    style = "minimal", icon = "ANALYSIS OPTIONS",
                    status = "primary", width = "600px", size = "sm",

                    # select DE method
                    selectInput(
                      inputId = "DE_method",
                      label= tags$h5("Choose DE Method"),
                      choices = c("wilcox", "DESeq2", "edgeR"),
                      selected = NULL,
                      width = "600px",
                    ),
                    
                    # select case control method
                    radioButtons(
                      inputId = "case_control_method",
                      label = tags$h5("Case/Control Selection"),
                      choices = c("Choose Case by Label", "Choose Case/Controls individually"),
                      selected = ""
                    ),

                    # Choose Case by Label
                    conditionalPanel(condition = "input.case_control_method == 'Choose Case by Label'", 
                      selectInput(
                        inputId="select_column",
                        label= "Select label to group ",
                        choices = NULL, # no choice before uploading
                        width = "600px",
                      ),
                  
                      selectInput(
                        inputId="select_case",
                        label= "Select case to analyse",
                        choices = NULL, # no choice before column selected
                        width = "600px",
                      ),
                    ),

                    # Choose Case/Controls individually
                    conditionalPanel(condition = "input.case_control_method == 'Choose Case/Controls individually'", 
                      h6(strong("Select Cases")),
                      dataTableOutput("UILabelContentSelection"),   
                      h6(strong("Select Rows to ignore")),
                      dataTableOutput("UILabelContentRemoveSelection"),                   
                    ),
                    
                    # run DE button
                    actionButton(inputId="run_DE", label = "Run DE"),
                  
                  ),

                  splitLayout(
                    cellWidths = c("50%", "50%"), 

                    # Volcano Plot
                    fluidPage(
                      #textOutput("DE_V_text"),
                      h4(id="vol"," Volcano Plot", style="text-align: center;"),
                      column(12, plotOutput(outputId = "DEplot", height = "450px"), align = "center"), 
                    ),

                    # MA Plot
                    fluidPage(
                      #textOutput("DE_MA_text"),
                      h4(id="MA"," MA Plot", style="text-align: center;"),
                      column(12, plotOutput(outputId = "DEplot_average", height = "450px"), align = "center"),
                    )

                  ), 
                  
                  # assess DE button
                  actionButton(inputId="assess_run_de", label = "Assess DE Data"), 
                 
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
                actionButton("generate_subnet", "Generate Subnetwork",),
      
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
                    id="subnetwork_file_tabset",

                    # view file tab
                    tabPanel(
                      title="File",
                      uiOutput("UIDEContent"),
                    ),

                    # view subnetwork tab
                    tabPanel(
                      title="Subnetwork", 
                      br(),
                      tableOutput("subnetwork")
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
                        inputId = "CG_PlotOptions",
                        label = tags$h4("Select Plots"), 
                        choices = c("Network", "Heatmap", "Binarized Heatmap"),
                        status = ""
                      ),

                      # run button
                      actionButton(inputId = "CG_run", label = "Run",
                      style="color: #fff; background-color: #3E3F3A; border-color: #20201F"),
                    ),  

                    br(),
                    # error message
                    textOutput("CG_error"),
                    br(),

                    tabsetPanel(

                      # plots tab
                      tabPanel(
                        title="Plots",

                        # network
                        conditionalPanel(
                          condition = "$.inArray('Network', input.CG_PlotOptions) > -1", 
                          h5(id="CG_network_text", "Network of Clustered Genes"), 
                          br(),
                          plotOutput(outputId = "network", height = "500px"),
                        ),
                        
                        # heatmap
                        conditionalPanel(
                          condition = "$.inArray('Heatmap', input.CG_PlotOptions) > -1", 
                          h5(id="CG_heatmap_text", "Heatmap of Clustered Genes"),
                          br(),
                          plotOutput(outputId = "heatmap", height = "500px"),
                        ),

                        # binarized heatmap
                        conditionalPanel(
                          condition = "$.inArray('Binarized Heatmap', input.CG_PlotOptions) > -1", 
                          h5(id="CG_bheatmap_text", "Binarized Heatmap of Clustered Genes"), 
                          br(),
                          plotOutput(outputId = "Bheatmap", height = "500px"), 
                        ),

                        conditionalPanel(
                          condition = "$.inArray('Upregulated Network', input.CG_PlotOptions) > -1", 
                          h4("Network of Clustered, Upregulated Genes"), 
                          br(), 
                          plotOutput(outputId = "upregNetwork", height = "500px"), 
                        ),

                        conditionalPanel(
                          condition = "$.inArray('Upregulated Heatmap', input.CG_PlotOptions) > -1", 
                          h4("Heatmap of Clustered, Upregulated Genes"), 
                          br(), 
                          plotOutput(outputId = "upregHeatmap", height = "500px"),
                        ),

                        conditionalPanel(
                          condition = "$.inArray('Upregulated Binarized Heatmap', input.CG_PlotOptions) > -1", 
                          h4("Binarized Hatmap of Clustered, Upregulated Genes"), 
                          br(),
                          plotOutput(outputId = "upregbinHeatmap", height = "500px"), 
                        ), 

                        conditionalPanel(
                          condition = "$.inArray('Downregulated Network', input.CG_PlotOptions) > -1",
                          h4("Network of Clustered, Downregulated Genes"), 
                          br(), 
                          plotOutput(outputId = "downregNetwork", height = "500px"),
                        ), 

                        conditionalPanel(
                          condition = "$.inArray('Downregulated Heatmap', input.CG_PlotOptions) > -1", 
                          h4("Heatmap of Clustered, Downregulated Genes"), 
                          br(), 
                          plotOutput(outputId = "downregHeatmap", height = "500px"),
                        ), 

                        conditionalPanel(
                          condition = "$.inArray('Downregulated Binarized Heatmap', input.CG_PlotOptions) > -1", 
                          h4("Binarized Heatmap of Clustered, Downregulated Genes"), 
                          br(), 
                          plotOutput(outputId = "downregbinHeatmap", height = "500px"),    
                        ),

                      ),

                      # tables tab
                      tabPanel(
                        title="Tables",
                        
                        # clustering genes
                        h5(id="CG_table_text", "Clustering Genes"), 
                        br(),
                        fluidRow(
                          column(11,
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
                        inputId = "GC_PlotOptions",
                        label = tags$h4("Select Plots"), 
                        choices = c("Density", "Histogram", "Clustered Density", "Clustered Histogram"),
                        status = ""
                      ),
                      
                      # filt_min slider
                      conditionalPanel(
                        condition = "$.inArray('Histogram', input.GC_PlotOptions) > -1 || $.inArray('Clustered Histogram', input.GC_PlotOptions) > -1" ,
                        sliderInput(
                          inputId="xybreaks", 
                          label = "Number of breaks for histogram:",
                          min = 10, max = 150, value = 100, step = 10,
                        ),
                      ),
                      
                      br(),

                      # run button
                      actionButton(
                        inputId = "GC_run",
                        label = "Run", 
                        style="color: #fff; background-color: #3E3F3A; border-color: #20201F"
                      ),

                    ),
                    
                    br(),

                    # error message
                    textOutput("GC_error"),

                    # density
                    conditionalPanel(
                      br(),
                      condition = "$.inArray('Density', input.GC_PlotOptions) > -1", 
                      h5(id="GCdensityG_text", "Density Plot of Gene Connectivity"), 
                      br(),
                      plotOutput(outputId = "GCdensityG", height = "500px",),
                      br(),
                    ),

                    # histogram
                    conditionalPanel(
                      br(),
                      condition = "$.inArray('Histogram', input.GC_PlotOptions) > -1", 
                      h5(id="GChistogramG_text", "Histogram of Gene Connectivity"),
                      br(),
                      plotOutput(outputId = "GChistogramG", height = "500px",),
                      br(),
                    ),

                    # density (subset by clusters)
                    conditionalPanel(
                      br(),
                      condition = "$.inArray('Clustered Density', input.GC_PlotOptions) > -1", 
                      h5(id="GCdensitySubsetG_text", "Density plot of Gene Connectivity subset by their clusters"), 
                      br(),
                      plotOutput(outputId = "GCdensitySubsetG", height = "500px",),
                      br(),
                    ),

                    # histogram (subset by clusters)
                    conditionalPanel(
                      br(),
                      condition = "$.inArray('Clustered Histogram', input.GC_PlotOptions) > -1", 
                      h5(id="GChistogramSubsetG_text", "Histogram of Gene Connectivity subset by their clusters"), 
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
                      actionButton(
                        inputId = "FO_run",
                        label = "Run",
                        style="color: #fff; background-color: #3E3F3A; border-color: #20201F"
                      ),

                    ),

                    br(),
                    
                    # error message
                    textOutput("FO_error"),

                  ),
                  br(),
                  
                  tabsetPanel(

                    # plots tab
                    tabPanel(
                      title="Plots",
                      
                      # heatmap
                      conditionalPanel(
                        condition = "$.inArray('Network', input.FOPlotOptions) > -1", 
                        h5(id="FO_network_text", "Network"), 
                        plotOutput(outputId = "FO_network", height = "500px"),
                      ),

                      # network
                      conditionalPanel(
                        condition = "$.inArray('Heatmap', input.FOPlotOptions) > -1", 
                        h5(id="FO_heatmap_text", "Heatmap"), 
                        plotOutput(outputId = "FO_heatmap", height = "500px"),
                        br(),
                      ),
                    ),

                    # tables tab
                    tabPanel(
                      title="Tables", 

                      # selected genes table output
                      conditionalPanel(
                        condition = "$.inArray('Genes in Module', input.FO_table_options) > -1", 
                        h5(id="genes_not_keep_table_text", "Genes in Module"), 
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
                        h5(id="genes_keep_table_text", "Outliers"), 
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
                      awesomeCheckboxGroup(
                        inputId = "GSEA_type",
                        label = tags$h4("GSEA Type"),
                        choices = c("Standard GSEA"),
                        selected = ""
                      ),
                      
                      # standard GSEA options
                      conditionalPanel(
                        condition = "input.GSEA_type.includes('Standard GSEA')", 
                        awesomeCheckboxGroup(
                          inputId = "GSEA_std_PlotOptions",
                          label = tags$h4("Standard GSEA"), 
                          choices = c("P-value Heatmap"),
                          status = ""
                        ),
                      ),

                      # AUCs GSEA options
                      conditionalPanel(
                        condition = "input.GSEA_type.includes('AUCs GSEA')", 
                        awesomeCheckboxGroup(
                          inputId = "GSEA_auc_PlotOptions",
                          label = tags$h4("AUCs GSEA"), 
                          choices = c("AUROC Graph", "AUC Histogram", "Gene Set Size Histogram", "Scatter Plot (AUCs vs P-values)"),
                          status = ""
                        ),
                      ),
                      
                      br(),

                      # run button
                      actionButton(
                        inputId = "GSEA_run",
                        label = "Run",
                        style="color: #fff; background-color: #3E3F3A; border-color: #20201F"
                      ),

                    ),
                    
                    br(),

                    # error message
                    textOutput("GSEA_error"),
                  ),
                  br(),
                  
                  tabsetPanel(

                    # Standard GSEA tab
                    tabPanel(
                      title="Standard",
                      mainPanel(

                        # heatmap
                        conditionalPanel(
                          condition = "$.inArray('P-value Heatmap', input.GSEA_std_PlotOptions) > -1", 
                          h5(id="GSEA_heatmap_text", "P-value Heatmap"), 
                          plotOutput(outputId = "GSEA_heatmap", height = "500px"),
                        ),
                        br(),

                        # upregulated heatmap
                        conditionalPanel(
                          condition = "$.inArray('Upregulated P-value Heatmap', input.GSEA_std_PlotOptions) > -1", 
                          h5(id="GSEA_up_heatmap_text", "Upregulated P-value Heatmap"), 
                          plotOutput(outputId = "GSEA_up_heatmap", height = "500px"),
                        ),
                        br(),

                        # downregulated heatmap
                        conditionalPanel(
                          condition = "$.inArray('Downregulated P-value Heatmap', input.GSEA_std_PlotOptions) > -1", 
                          h5(id="GSEA_down_heatmap_text", "Downregulated P-value Heatmap"), 
                          plotOutput(outputId = "GSEA_down_heatmap", height = "500px"),
                        ),
                        br(),

                      )
                    ),

                    # AUCs GSEA tab
                    tabPanel(
                      title="AUC",
                      mainPanel(
                        # AUROC graph
                        conditionalPanel(
                          condition = "$.inArray('AUROC Graph', input.GSEA_auc_PlotOptions) > -1", 
                          h5(id="GSEA_auroc_text", "AUROC Graph"), 
                          plotOutput(outputId = "GSEA_auroc", height = "500px"),
                        ),
                        br(),

                      )
                    ),
                  )
                  

                  
                  
                ),

              
              ),
            )
  )
)
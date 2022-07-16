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

            #tabPanel is each tab in the navbarPage
            # home tab
            tabPanel(
              title="Home",
              icon = icon("home"),

              # navlistPanel is each tab on the side menu for each tabPanel
              navlistPanel(
                id = "Header", selected = NULL, well = FALSE, fluid = TRUE, widths = c(3, 9),
             
                tabPanel(title="About",
                  h3("What is OutDeCo?"),
                  p("Outlier detection through co-expression. Using the tools on this website you can:"),
                  p("- Run a differential expression (DE) analysis "),
                  p("- Assess your DE results using gene co-expression properties"),
                  p("- Report a functional outlier assessment"),
                  p("- Run a network connectivity analysis of DE results within a gene co-expression network"),
                ),
                 
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

            # run DE tab
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

                  # upload counts file
                  fileInput("counts_file", 
                            label = "Upload counts",
                            accept = c(".csv", ".tsv", ".txt")
                  ),

                  radioButtons(
                    inputId = 'sepCountsButton', 
                    label = 'Delimiter Selector', 
                    choices = c(Default=''), 
                    selected = ''
                  ),

                  fileInput(
                    inputId = "labels_file",
                    label = "Choose Labels File",
                    accept = c(".csv", ".tsv", ".txt")
                  ),

                  radioButtons(
                    inputId = 'sepLabelsButton', 
                    label = 'Delimiter Selector', 
                    choices = c(Default=''), 
                    selected = ''
                  )
                ), 
               ),


               br(),
               navlistPanel(
                widths = c(3, 9), well = FALSE,
                tabPanel(
                  title="View File",

                  tabsetPanel(
                    tabPanel(
                      title="Counts File",
                      uiOutput("UICountsContent")
                    ),
                    tabPanel(
                      title="Labels File",
                      dataTableOutput("UILabelContent")
                    )
                  )


                 ),
                tabPanel(
                  title="Plot DE",
                  dropdown(
                    # side panel characteristics
                    style = "minimal", icon = "OPTIONS",
                    status = "primary", width = "600px", size = "sm",
                  splitLayout(cellWidths = c("50%", "50%"), 
                  fluidPage(
                    h5(id = "case_selection", "Case Selection"),
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
                    br(),
                    br(), 
                    br(),
                  ),
                  fluidPage(
                    h5(id = "control_selection", "Control Selection"),
                    br(),
                    actionButton(inputId="run_DE", label = "Run DE"),
                  ),
                  ),
                  ),

                    splitLayout(cellWidths = c("50%", "50%"), 
                    plotOutput(outputId = "DEplot", height = "450px"), 
                    plotOutput(outputId = "DEplot_average", height = "450px"))
                          
                     
                  
                  #plotOutput(outputId = "DEplot", height = "500px"),
                  #plotOutput(outputId = "DEplot_average", height = "500px"),
                 ),
                 
              ),
            ),
            
            # Assess DE tab
            tabPanel(
              title="Assess DE", 
              
              dropdown(
                tags$h4("Network Selection"),
                selectInput(
                  inputId = "network_type",
                  label=NULL,
                  choices = c("Blood", "Brain", "Generic"),
                  selected = "Generic"
                ),

                # gene list
                radioButtons(
                  inputId = "gene_list_selection",
                  label = tags$h4("Gene List Selection"),
                  choices = c("Upload Gene List", "Generate Gene List"),
                  selected = ""
                ),
            
                # upload gene list
                conditionalPanel(
                  condition = "input.gene_list_selection == 'Upload Gene List'", 

                  # upload file
                  fileInput(
                    inputId = "DEFile", 
                    label = "Choose Gene List File",
                    accept = c(".csv", ".tsv", ".txt")
                  ),

                  # div(style = "margin-top: -25px"),
                  # button for selecting delimiter, default is nothing until file is selected and handled in server side
                  radioButtons(
                    inputId = 'sepDEButton', 
                    label = 'Delimiter Selector', 
                    choices = c(Default=''), 
                    selected = ''
                  ),
                ),

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
                tabPanel(
                  title="View Files",
                  
                  tabsetPanel(
                    tabPanel(
                      title="File",
                      uiOutput("UIDEContent")
                    ),
                    tabPanel(
                      title="Subnetwork", 
                      tableOutput("subnetwork"),
                      

                    )
                  ),
                ),

                tabPanel(
                  title="Cluster Genes",
                  mainPanel(
                    h3("Cluster Genes"),
                    dropdown(
                      inputId = "CG_dropdown",

                      # title of sidepanel
                      tags$h4("Cluster Genes Options"),

                      # inputs in the sidepanel
                      # side panel characteristics
                      style = "minimal", icon = "PLOT OPTIONS",
                      status = "primary", width = "300px", size = "sm",
                      
                      awesomeCheckboxGroup(
                        inputId = "clusterPlotOptions",
                        label = "Select Plots", 
                        choices = c("Network", "Heatmap", "Binarized Heatmap"),
                        status = ""
                          
                      ),
                      textOutput('error_msg'),
                      actionButton(inputId = "run", label = "Run", 
                      style="color: #fff; background-color: #3E3F3A; border-color: #20201F"),
              
                    ),  
                    br(),
                    
                    textOutput("CG_error"),
                    
                    conditionalPanel(
                      condition = "$.inArray('Network', input.clusterPlotOptions) > -1", 
                      textOutput("CNtext"), 
                      plotOutput(outputId = "network", height = "500px"),
                    ),

                    conditionalPanel(
                      condition = "$.inArray('Heatmap', input.clusterPlotOptions) > -1", 
                      textOutput("CHtext"),
                      plotOutput(outputId = "heatmap", height = "500px"),
                      br(),
                    ),

                    conditionalPanel(
                      condition = "$.inArray('Binarized Heatmap', input.clusterPlotOptions) > -1", 
                      textOutput("CHBtext"), 
                      plotOutput(outputId = "Bheatmap", height = "500px"), 
    
                    ),

                  )
                ),
                
                tabPanel(
                  title="Gene Connectivity",
                  h3("Gene Connectivity"),
                  mainPanel(
                    
                    dropdown(
                      inputId = "GC_dropdown",
                      # title of sidepanel
                      tags$h4("Gene Connectivity Options"),

                      # inputs in the sidepanel
                      # side panel characteristics
                      style = "minimal", icon = "PLOT OPTIONS",
                      status = "primary", width = "300px", size = "sm",
                     
                      awesomeCheckboxGroup(
                        inputId = "GCPlotOptions",
                        label = "Select Plots", 
                        choices = c("Density", "Histogram", "Clustered Density", "Clustered Histogram"),
                        status = ""
                        
                      ),

                      # ADD XYBREAKS SLIDER FOR HISTOGRAM 
                      conditionalPanel(
                        condition = "$.inArray('Histogram', input.GCPlotOptions) > -1 || $.inArray('Clustered Histogram', input.GCPlotOptions) > -1" ,
                        sliderInput("xybreaks", label = "Number of breaks for histogram:",
                          min = 10, max = 150, value = 100, step = 10
                        ),
                      ),
                      
                      br(),
                      actionButton(inputId = "runGC", label = "Run", 
                      style="color: #fff; background-color: #3E3F3A; border-color: #20201F"),
              
                    ),
                    
                    br(),
                    textOutput("GC_error"),

                    # Run Gene Connectivity

                    # Density Plot Selected
                    conditionalPanel(
                      br(),
                      condition = "$.inArray('Density', input.GCPlotOptions) > -1", 
                      textOutput("GCdensityGtext"), 
                      plotOutput(outputId = "GCdensityG", height = "500px",),
                      br(),
                    ),

                    # Histogram Selected
                    conditionalPanel(
                      br(),
                      condition = "$.inArray('Histogram', input.GCPlotOptions) > -1", 
                      textOutput("GChistogramGtext"), 
                      plotOutput(outputId = "GChistogramG", height = "500px",),
                      br(),
                    ),

                    # Density Selected - subset by clusters
                    conditionalPanel(
                      br(),
                      condition = "$.inArray('Clustered Density', input.GCPlotOptions) > -1", 
                      textOutput("GCdensitySubsetGtext"), 
                      plotOutput(outputId = "GCdensitySubsetG", height = "500px",),
                      br(),

                    ),

                    # Histogram Selected - subset by clusters
                    conditionalPanel(
                      br(),
                      condition = "$.inArray('Clustered Histogram', input.GCPlotOptions) > -1", 
                      textOutput("GChistogramSubsetGtext"), 
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

                    dropdown(
                      inputId = "FO_dropdown", 

                      
                      awesomeCheckboxGroup(
                        inputId = "FOPlotOptions",
                        label = tags$h4("Plot Options"),
                        choices = c("Network", "Heatmap"),
                        status = ""
                      ),

                      awesomeCheckboxGroup(
                          inputId = "FO_table_options",
                          label = tags$h4("Table Options"),
                          choices = c("Selected Genes", "Unselected Genes"),
                          status = ""
                      ),
                      
                      sliderInput("filtmin", label = "Number of Genes to form Module",
                          min = 0, max = 20, value = 6, step = 2
                      ),

                      br(),
                      actionButton(inputId = "runFO", label = "Run", 
                      style="color: #fff; background-color: #3E3F3A; border-color: #20201F"),

                      style = "minimal", icon = "FUNCTIONAL OUTLIERS OPTIONS",
                      status = "primary", width = "300px", size = "sm",

                    ), 
                    br(),
                    
                    textOutput("FO_error"),

                  ),
                  br(),
                  tabsetPanel(

                    # plots
                    tabPanel(
                      br(),
                      title="Plots",
                      br(),
                      
                      # heatmap output
                      conditionalPanel(
                        condition = "$.inArray('Network', input.FOPlotOptions) > -1", 
                        textOutput("FO_network_text"), 
                        plotOutput(outputId = "FO_network", height = "500px"),
                      ),

                      # network output
                      conditionalPanel(
                        condition = "$.inArray('Heatmap', input.FOPlotOptions) > -1", 
                        textOutput("FO_heatmap_text"), 
                        plotOutput(outputId = "FO_heatmap", height = "500px"),
                        br(),
                      ),
                    ),

                    # tables
                    tabPanel(
                      br(),
                      title="Tables", 
                      br(),

                      # selected genes table output
                      conditionalPanel(
                        condition = "$.inArray('Selected Genes', input.FO_table_options) > -1", 
                        textOutput("FO_selected_text"), 
                        fluidRow(
                          column( 11,
                                  dataTableOutput("genes_selected_table"),
                          )
                        ),
                      ),

                      br(),

                      # unselected genes table output
                      conditionalPanel(
                        condition = "$.inArray('Unselected Genes', input.FO_table_options) > -1", 
                        textOutput("FO_unselected_text"), 
                        fluidRow(
                          column( 11,
                                  dataTableOutput("genes_unselected_table"),
                          )
                        ),
                      ),
                    )
                  ),
                ),

                
                tabPanel(
                  title="Gene Set Enrichment Analysis",
                  "GSE Page",
                ),
               ),
             )
  ),
  
)

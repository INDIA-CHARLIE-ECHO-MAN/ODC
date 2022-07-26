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
library(OutDeCo)
library(EGAD)


ui <- fluidPage(
  
  ######################################################################
  #                                                                    #
  #                               THEME                                #
  #                                                                    #
  ######################################################################

  useShinyjs(),
  chooseSliderSkin("Flat",  color = "#3E3F3A"),
  add_busy_spinner(spin = "dots", position = "bottom-right", color = "#3E3F3A"),
  
  titlePanel(title=div(img(src="ODClogo.png", height = 80), "OutDeCo")),

  theme = bs_theme(version = 5, bootswatch = "sandstone", 
                  heading_font = font_google("Poppins"), 
                  base_font = font_collection(font_google("Roboto")),
                  success = "#325D88"),

  tags$style(HTML(".js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {background: #3E3F3A}")),
  




  navbarPage(
    title = NULL, id = "navpage", collapsible = FALSE,

            ######################################################################
            #                                                                    #
            #                             HOME TAB                               #
            #                                                                    #
            ######################################################################

            tabPanel(title = "Home",
              icon = icon("home"),

              navlistPanel(
                id = "Header", selected = NULL, well = FALSE, fluid = TRUE, widths = c(3, 9),

                ########################### ABOUT ###########################

                tabPanel(title = "About",
                  h3("What is OutDeCo?"),
                  p("Outlier detection through co-expression. Using the tools on this website you can:"),
                  p("- Run a differential expression (DE) analysis "),
                  p("- Assess your DE results using gene co-expression properties"),
                  p("- Report a functional outlier assessment"),
                  p("- Run a network connectivity analysis of DE results within a gene co-expression network"),
                ),
                
                ########################### DIFFERENTIAL ANALYSIS ###########################

                tabPanel(title = "Differential Expression Analysis",
                  h3("Differential Expression Analysis"),
                  p("Statistical analysis to discover quantitative changes in expression levels between experimental groups."),
                  h5("Methods:"),
                  h6(strong("wilcox")),
                  p("Compares means of two groups to analyse count data and test for differential expression."),
                  h6(strong("DESeq")),
                  p("Uses geometric normalisation to analyse count data and test for differential expression."),
                  h6(strong("edgeR")),
                  p("Uses weighted mean of log ratio to analyse count data and test for differential expression."),
                ),

                ########################### CLUSTER GENES ###########################

                tabPanel(title = "Cluster Genes",
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
                ),

                ########################### GENE CONNECTIVITY ###########################

                tabPanel(title = "Gene Connectivity",
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
                ),
                 
                ########################### FUNCTIONAL OUTLIERS ###########################

                tabPanel(title = "Functional Outliers",
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
                ),
                
                ########################### GSEA ###########################

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
                ),
                 
              )
            ),

            ######################################################################
            #                                                                    #
            #                            RUN DE TAB                              #
            #                                                                    #
            ######################################################################
            
            tabPanel(title = "Run DE", 

              # options dropdown
              dropdown(

                # dropdown characteristics
                style = "jelly", icon = "FILE OPTIONS",
                status = "primary", width = "300px", size = "sm",

                # data selection
                tags$h3("Data Selection"),

                # counts data file
                fileInput(inputId = "DE_run_counts_file", label = "Choose Counts Data File"),
                
                radioButtons( # delimiter selector
                  inputId = 'DE_run_counts_sep', 
                  label = 'Delimiter Selector',
                  choices = c(Comma = ",", Semicolon = ";", Tab = "\t", Space = " "), 
                  selected = ''
                ),

                # labels file
                fileInput(inputId = "DE_run_labels_file", label = "Choose Labels File",
                  accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")
                ),

                radioButtons( # delimiter selector
                  inputId = 'DE_run_labels_sep', 
                  label = 'Delimiter Selector', 
                  choices = c(Comma = ",", Semicolon = ";", Tab = "\t", Space = " "),
                  selected = ''
                ),

              ),

              br(),

              navlistPanel(
                widths = c(3, 9), well = FALSE,
                
                ########################### VIEW FILE ###########################

                tabPanel(title = "View File",

                  tabsetPanel(
                    id = "DE_run_VF_tabset",

                    # counts file
                    tabPanel(
                      title = "Counts File",
                      dataTableOutput("DE_run_counts_table")
                    ),

                    # labels file
                    tabPanel(
                      title = "Labels File",
                      dataTableOutput("UILabelContent")
                    ),

                  )
                ),

                ########################### PLOT DE ###########################

                tabPanel(title = "Plot DE",

                  h4("Plot Differential Expression"),
                  p(id = "DE_plot_error", "Please upload counts and labels data in FILE OPTIONS", style = "color:red"),

                  # options dropdown
                  dropdown(

                    inputId = "DE_run_PD_options",

                    # dropdown characteristics
                    style = "minimal", icon = "OPTIONS",
                    status = "primary", width = "600px", size = "sm",

                    # DE method
                    selectInput(
                      inputId = "DE_run_PD_options_method",
                      label = tags$h5("Choose DE Method"),
                      choices = c("wilcox", "DESeq2", "edgeR"),
                      selected = NULL,
                      width = "600px",
                    ),
                    
                    # case control method
                    radioButtons(
                      inputId = "DE_run_PD_options_control_case",
                      label = tags$h5("Case/Control Selection"),
                      choices = c("Choose Case by Label", "Choose Case/Controls individually"),
                      selected = "",
                    ),


                    conditionalPanel(
                      condition = "input.DE_run_PD_options_control_case == 'Choose Case by Label'",

                      selectInput(
                        inputId = "DE_run_PD_options_column",
                        label = "Select label to group",
                        choices = NULL, # no choice before uploading
                        width = "600px",
                      ),
                  
                      selectInput(
                        inputId = "DE_run_PD_options_case",
                        label= "Select case to analyse",
                        choices = NULL, # no choice before column selected
                        width = "600px",
                      ),

                    ),

                    conditionalPanel(
                      condition = "input.DE_run_PD_options_case == 'Choose Case/Controls individually'",

                      h6(strong("Select Cases")),
                      dataTableOutput("UILabelContentSelection"),   
                      h6(strong("Select Conditions")),
                      dataTableOutput("UILabelContentRemoveSelection"),                   
                    ),
                    
                    # run DE button
                    actionButton(inputId = "DE_run", label = "Run DE"),
                  
                  ),

                  splitLayout(cellWidths = c("50%", "50%"), 

                    # volcano plot
                    fluidPage(
                      h4(id = "DE_vol_text", label = "Volcano Plot", style = "text-align: center;"),
                      column(12, plotOutput(outputId = "DE_vol_plot", height = "450px"), align = "center"), 
                    ),
                    
                    # MA plot
                    fluidPage(
                      h4(id = "DE_MA_text", label = "MA Plot", style = "text-align: center;"),
                      column(12, plotOutput(outputId = "DE_MA_plot", height = "450px"), align = "center"),
                    )

                  ),

                  # run button
                  actionButton(inputId = "DE_run_assess", label = "Assess DE Data"), 
                  
                  
                 
                ),
                 
              ),
            ),
            
            ######################################################################
            #                                                                    #
            #                           ASSESS DE TAB                            #
            #                                                                    #
            ######################################################################

            tabPanel(title = "Assess DE", 

              # options dropdown
              dropdown(

                # dropdown characteristics
                style = "jelly", icon = "NETWORK OPTIONS",
                status = "primary", width = "300px", size = "sm",

                # network selection
                selectInput(
                  inputId = "DE_options_network",
                  label = tags$h4("Network Selection"),
                  choices = c("Blood", "Brain", "Generic"),
                  selected = "Generic"
                ),

                # gene list selection
                radioButtons(
                  inputId = "DE_options_data",
                  label = tags$h4("DE Data Selection"),
                  choices = c("Use DE Results"),
                  selected = "Use DE Results"
                ),
                
                # generate subnet button
                actionButton("DE_generate_subnet", "Generate Subnetwork",),

              ),
              
              br(),
        
              navlistPanel(
                widths = c(3, 9), well = FALSE,
                
                ########################### VIEW FILES ###########################

                tabPanel(title = "View Files",
                  
                  tabsetPanel(
                    id = "DE_VF_tabset",

                    # DE data
                    tabPanel(title = "DE Data",
                      br(),
                      dataTableOutput("DE_VF_table"),
                    ),

                    # DE subnetwork
                    tabPanel(title = "Subnetwork", 
                      br(),
                      tableOutput("DE_VF_subnetwork")
                    ),
                   
                  ),
                ),

                ########################### CLUSTER GENES ###########################

                tabPanel(title = "Cluster Genes",

                  mainPanel(
                    h3("Cluster Genes"),

                    # options dropdown
                    dropdown(
                      inputId = "DE_CG_options",

                      # dropdown characteristics
                      style = "minimal", icon = "OPTIONS",
                      status = "primary", width = "300px", size = "sm",

                      # upregulated
                      awesomeCheckboxGroup(
                        inputId = "DE_CG_options_up", 
                        label = tags$h4("Upregulated"),
                        choices = c("Network", "Heatmap", "Binarized Heatmap"),
                        status = ""
                      ), 
                      
                      # downregulated
                      awesomeCheckboxGroup(
                        inputId = "DE_CG_options_down", 
                        label = tags$h4("Downregulated"), 
                        choices = c("Network", "Heatmap", "Binarized Heatmap"), 
                        status = ""
                      ),

                      # run button
                      actionButton(inputId = "DE_CG_run", label = "Run",),
                    ),  
                    
                    # error message
                    p(id = "DE_CG_error", "Please generate a subnetwork in OPTIONS.", style = "color:red"),
                    br(),
                    
                    tabsetPanel(

                      # upregulated tab
                      tabPanel(title="Upregulated",

                        # upregulated - network
                        conditionalPanel(
                          condition = "$.inArray('Network', input.DE_CG_options_up) > -1", 
                          h5(id = "DE_CG_up_network_text", "Network of Clustered, Upregulated Genes"), 
                          plotOutput(outputId = "DE_CG_up_network_plot", height = "500px"), 
                          br(), 
                        ),

                        # upregulated - heatmap
                        conditionalPanel(
                          condition = "$.inArray('Heatmap', input.DE_CG_options_up) > -1", 
                          h5(id = "DE_CG_up_heatmap_text","Heatmap of Clustered, Upregulated Genes"), 
                          plotOutput(outputId = "DE_CG_up_heatmap_plot", height = "500px"),
                          br(), 
                        ),

                        # upregulated - binarized heatmap
                        conditionalPanel(
                          condition = "$.inArray('Binarized Heatmap', input.DE_CG_options_up) > -1", 
                          h5(id = "DE_CG_up_bheatmap_text","Binarized Heatmap of Clustered, Upregulated Genes"), 
                          plotOutput(outputId = "DE_CG_up_bheatmap_plot", height = "500px"),
                          br(),
                        ), 

                      ),

                      # downregulated tab
                      tabPanel(title="Downregulated",

                        # downregulated - network
                        conditionalPanel(
                          condition = "$.inArray('Network', input.DE_CG_options_down) > -1",
                          h5(id = "DE_CG_down_network_text", "Network of Clustered, Downregulated Genes"), 
                          plotOutput(outputId = "DE_CG_down_network_plot", height = "500px"),
                          br(), 
                        ), 

                        # downregulated - heatmap
                        conditionalPanel(
                          condition = "$.inArray('Heatmap', input.DE_CG_options_down) > -1", 
                          h5(id = "DE_CG_down_heatmap_text", "Heatmap of Clustered, Downregulated Genes"), 
                          plotOutput(outputId = "DE_CG_down_heatmap_plot", height = "500px"),
                          br(), 
                        ), 

                        # downregulated - binarized heatmap
                        conditionalPanel(
                          condition = "$.inArray('Binarized Heatmap', input.DE_CG_options_down) > -1", 
                          h5(id = "DE_CG_down_bheatmap_text", "Binarized Heatmap of Clustered, Downregulated Genes"), 
                          plotOutput(outputId = "DE_CG_down_bheatmap_plot", height = "500px"),
                          br(), 
                        ),

                      ),

                    )
                  )
                ),
                
                ########################### GENE CONNECTIVITY ###########################

                tabPanel(title = "Gene Connectivity",

                  mainPanel(
                    h3("Gene Connectivity"),
                    
                    # options dropdown
                    dropdown(
                      inputId = "DE_GC_options",
                      style = "minimal", icon = "OPTIONS",
                      status = "primary", width = "300px", size = "sm",

                      # select plots
                    
                      awesomeCheckboxGroup(
                        inputId = "DE_GC_options_up", 
                        label = tags$h4("Upregulated"),
                        choices = c("Density", "Histogram", "Clustered Density", "Clustered Histogram"), 
                        status = ""
                      ), 

                      awesomeCheckboxGroup(
                        inputId = "DE_GC_options_down", 
                        label = tags$h4("Downregulated"),
                        choices =  c("Density", "Histogram", "Clustered Density", "Clustered Histogram"),
                        status = ""

                      ),
                      
                      # histogram breaks slider
                      conditionalPanel(
                        condition = "$.inArray('Histogram', input.DE_GC_options_up) > -1 || $.inArray('Clustered Histogram', input.DE_GC_options_up) > -1 || $.inArray('Histogram', input.DE_GC_options_down) > -1 || $.inArray('Clustered Histogram', input.DE_GC_options_down) > -1" ,
                        sliderInput(
                          inputId = "DE_GC_options_xybreaks", 
                          label = "Number of breaks for histogram:",
                          min = 10, max = 100, value = 100, step = 10,
                        ),
                      ),

                      # run button
                      actionButton(inputId = "DE_GC_run", label = "Run"),

                    ),

                    # error message
                    p(id = "DE_GC_error", "Please generate a subnetwork in OPTIONS.", style = "color:red"),
                    br(),

                    tabsetPanel(

                      # upregulated tab
                      tabPanel(title="Upregulated",

                        # upregulated - density
                        conditionalPanel(
                          condition = "$.inArray('Density', input.DE_GC_options_up) > -1", 
                          h5(id = "DE_GC_up_density_text", "Density Plot of Upregulated Gene Connectivity"), 
                          plotOutput(outputId = "DE_GC_up_density_plot", height = "500px",),
                          br(),
                        ),

                        # upregulated - histogram
                        conditionalPanel(
                          condition = "$.inArray('Histogram', input.DE_GC_options_up) > -1", 
                          h5(id = "DE_GC_up_hist_text", "Histogram of Upregulated Gene Connectivity"),
                          plotOutput(outputId = "DE_GC_up_hist_plot", height = "500px",),
                          br(),
                        ),

                        # upregulated - density (subset by clusters)
                        conditionalPanel(
                          condition = "$.inArray('Clustered Density', input.DE_GC_options_up) > -1", 
                          h5(id = "DE_GC_up_density_subset_text", "Density plot of Upregulated Gene Connectivity subset by their clusters"), 
                          plotOutput(outputId = "DE_GC_up_density_subset_plot", height = "500px",),
                          br(),
                        ),

                        # upregulated - histogram (subset by clusters)
                        conditionalPanel(
                          condition = "$.inArray('Clustered Histogram', input.DE_GC_options_up) > -1", 
                          h5(id = "DE_GC_up_hist_subset_text", "Histogram of Upregulated Gene Connectivity subset by their clusters"), 
                          plotOutput(outputId = "DE_GC_up_hist_subset_plot", height = "500px",),
                          br(),
                        ),
                        
                      ),

                      # downregulated tab
                      tabPanel(title="Downregulated",
                        
                        # downregulated - density 
                        conditionalPanel(
                          condition = "$.inArray('Density', input.DE_GC_options_down) > -1", 
                          h5(id = "DE_GC_down_density_text", "Density Plot of Downregulated Gene Connectivity"), 
                          plotOutput(outputId = "DE_GC_down_density_plot", height = "500px",),
                          br(),
                        ),

                        # downregulated - histogram
                        conditionalPanel(
                          condition = "$.inArray('Histogram', input.DE_GC_options_down) > -1", 
                          h5(id = "DE_GC_down_hist_text", "Histogram of Downregulated Gene Connectivity"),
                          plotOutput(outputId = "DE_GC_down_hist_plot", height = "500px",),
                          br(),
                        ),

                        # downregulated - density (subset by clusters)
                        conditionalPanel(
                          condition = "$.inArray('Clustered Density', input.DE_GC_options_down) > -1", 
                          h5(id = "DE_GC_down_density_subset_text", "Density plot of Downregulated Gene Connectivity subset by their clusters"), 
                          plotOutput(outputId = "DE_GC_down_density_subset_plot", height = "500px",),
                          br(),
                        ),

                        # downregulated - histogram (subset by clusters)
                        conditionalPanel(
                          condition = "$.inArray('Clustered Histogram', input.DE_GC_options_down) > -1", 
                          h5(id = "DE_GC_down_hist_subset_text", "Histogram of Dowregulated Gene Connectivity subset by their clusters"), 
                          plotOutput(outputId = "DE_GC_down_hist_subset_plot", height = "500px",),
                          br(),
                        ),

                      ),
                    ),
                  ),
                ),
                
                ########################### FUNCTIONAL OUTLIERS ###########################

                tabPanel(title = "Functional Outliers",
                  
                  mainPanel(
                    h3("Functional Outliers"),

                    # options dropdown
                    dropdown(
                      inputId = "DE_FO_options",
                      style = "minimal", icon = "OPTIONS",
                      status = "primary", width = "300px", size = "sm",

                      # select plots
                      awesomeCheckboxGroup(
                        inputId = "DE_FO_options_plots", 
                        label = tags$h4("Select Plots"),
                        choices = c("Upregulated Network", "Upregulated Heatmap", "Downregulated Network", "Downregulated Heatmap"), 
                        status = ""
                      ),

                      # filt_min slider
                      tags$h4("Other"),
                      sliderInput(
                        inputId = "DE_FO_filtmin",
                        label = "Number of Genes to form Module",
                        min = 0, max = 20, value = 6, step = 1
                      ),
                      
                      # run button
                      actionButton(inputId = "DE_FO_run", label = "Run"),

                    ),
                    
                    # error message
                    p(id = "DE_FO_error", "Please generate a subnetwork in OPTIONS.", style = "color:red"),
                    br(),
                    

                    tabsetPanel(

                      # plots tab
                      tabPanel(title = "Plots",                      

                        # upregulated - network
                        conditionalPanel(
                          condition = "$.inArray('Upregulated Network', input.DE_FO_options_plots) > -1", 
                          h5(id = "DE_FO_up_network_text", "Upregulated Network"), 
                          plotOutput(outputId = "DE_FO_up_network_plot", height = "500px"), 
                        ), 

                        # upregulated - heatmap
                        conditionalPanel(
                          condition = "$.inArray('Upregulated Heatmap', input.DE_FO_options_plots) > -1",
                          h5(id = "DE_FO_up_heatmap_text", "Upregulated Heatmap"), 
                          plotOutput(outputId = "DE_FO_up_heatmap_plot", height = "500px"),
                        ), 

                        # downregulated - network
                        conditionalPanel(
                          condition = "$.inArray('Downregulated Network', input.DE_FO_options_plots) > -1", 
                          h5(id = "DE_FO_down_network_text", "Downregulated Network"), 
                          plotOutput(outputId = "DE_FO_down_network_plot", height = "500px"),
                        ), 

                        # upregulated - heatmap
                        conditionalPanel(
                          condition = "$.inArray('Downregulated Heatmap', input.DE_FO_options_plots) > -1", 
                          h5(id = "DE_FO_down_heatmap_text", "Downregulated Heatmap"), 
                          plotOutput(outputId = "DE_FO_down_heatmap_plot", height = "500px"),
                        ),

                      ),

                      # tables tab
                      tabPanel(title = "Tables", 

                        # selected genes table output
                        # conditionalPanel(
                        #   condition = "$.inArray('Genes in Module', input.GL_FO_options_tables) > -1", 
                        #   h5(id = "GL_FO_in_text", "Genes in Module"), 
                        #   br(),
                        #   fluidRow(
                        #     column( 11,
                        #             dataTableOutput("GL_FO_in_table"),
                        #     )
                        #   ),
                        # ),

                        # br(),

                        # # unselected genes table output
                        # conditionalPanel(
                        #   condition = "$.inArray('Functional Outliers', input.GL_FO_options_tables) > -1", 
                        #   h5(id = "GL_FO_out_text", "Outliers"), 
                        #   br(),
                        #   fluidRow(
                        #     column( 11,
                        #             dataTableOutput("GL_FO_out_table"),
                        #     )
                        #   ),
                        # ),
                      ),
                    ),
                  ),
                  
                  
                ),

                ########################### GSEA ###########################

                tabPanel(title = "Gene Set Enrichment Analysis",
                  
                  mainPanel(
                    h3("Gene Set Enrichment Analysis"),
                    
                    # options dropdown
                    dropdown(
                      inputId = "DE_GSEA_options",
                      style = "minimal", icon = "OPTIONS",
                      status = "primary", width = "300px", size = "sm",
                      
                      # select GSEA type
                      awesomeCheckboxGroup(
                        inputId = "GSEA_type",
                        label = tags$h4("GSEA Type"),
                        choices = c("Standard GSEA", "AUCs GSEA"),
                        selected = ""
                      ),
                      
                      # standard GSEA options
                      conditionalPanel(
                        condition = "input.GSEA_type.includes('Standard GSEA')", 
                        awesomeCheckboxGroup(
                          inputId = "GSEA_std_PlotOptions",
                          label = tags$h4("Standard GSEA"), 
                          choices = c("Upregulated P-value Heatmap", "Downregulated P-value Heatmap"),
                          status = ""
                        ),
                      ),

                      # run button
                      actionButton(inputId = "DE_GSEA_run", label = "Run"),

                    ),

                    # error message
                    p(id = "DE_GSEA_error", "Please generate a subnetwork in OPTIONS.", style = "color:red"),
                    br(),

                    tabsetPanel(

                      # Standard GSEA tab
                      tabPanel(title="Standard",
                      
                        mainPanel(

                          # upregulated heatmap
                          conditionalPanel(
                            condition = "$.inArray('Upregulated P-value Heatmap', input.GSEA_std_PlotOptions) > -1", 
                            h5(id = "DE_GSEA_up_heatmap_text", "Upregulated P-value Heatmap"), 
                            br(),
                            plotOutput(outputId = "GSEA_up_heatmap", height = "500px"),
                          ),

                          # downregulated heatmap
                          conditionalPanel(
                            condition = "$.inArray('Downregulated P-value Heatmap', input.GSEA_std_PlotOptions) > -1", 
                            h5(id = "DE_GSEA_down_heatmap_text", "Downregulated P-value Heatmap"), 
                            br(),
                            plotOutput(outputId = "GSEA_down_heatmap", height = "500px"),
                          ),

                        )
                      ),

                      # AUCs GSEA tab
                      tabPanel(title = "AUC",

                        # AUCs graphs
                        conditionalPanel(
                          condition = "$.inArray('AUCs GSEA', input.GSEA_type) > -1", 
                          plotOutput(outputId = "GSEA_auc", height = "500px"),
                        ),
                      ),
                    )
                  ),

                ),


               ),
            ),

            ######################################################################
            #                                                                    #
            #                       ASSESS GENE LIST TAB                         #
            #                                                                    #
            ######################################################################
             
             tabPanel(title = "Assess Gene List",

              # options dropdown
              dropdown(

                # dropdown characteristics
                style = "jelly", icon = "NETWORK OPTIONS",
                status = "primary", width = "300px", size = "sm",

                # network selection
                tags$h4("Network Selection"),
                selectInput(
                  inputId = "GL_options_network",
                  label = NULL,
                  choices = c("Blood", "Brain", "Generic"),
                  selected = "Generic"
                ),

                # gene list selection
                radioButtons(
                  inputId = "GL_options_gene_list",
                  label = tags$h4("Gene List Selection"),
                  choices = c("Upload Gene List", "Generate Gene List"),
                  selected = ""
                ),

                # generate gene list
                conditionalPanel(
                  condition = "input.GL_options_gene_list == 'Generate Gene List'", 

                  # choose chromosome
                  textInput(
                    inputId = 'GL_chrome', 
                    label = 'Choose Chromosome' , 
                    placeholder = "chrX"
                  ),

                  # choose number of genes
                  textInput(
                    inputId = 'GL_genes_no', 
                    label = 'Choose Number of Genes',
                    placeholder = "100"
                  ),
                ),

                # upload gene list
                conditionalPanel(
                  condition = "input.GL_options_gene_list == 'Upload Gene List'",

                  # upload gene list
                  fileInput(
                    inputId = "GL_gene_list", 
                    label = "Choose Gene List File",
                    accept = c(".csv", ".tsv", ".txt")
                  ),

                  radioButtons( # select delimiter
                    inputId = 'GL_gene_list_sep', 
                    label = 'Delimiter Selector', 
                    choices = c(Default=''), 
                    selected = ''
                  ),

                  radioButtons( # select gene list type
                    inputId = 'GL_gene_list_type', 
                    label = 'Gene List Type', 
                    choices = c("Gene Names", "Entrez Id"), 
                    selected = ''
                  ),
                ),

                # generate subnet button
                actionButton("generate_subnet", "Generate Subnetwork"),
      

              ),

              br(),
        
              navlistPanel(
                widths = c(3, 9), well = FALSE,

                ########################### VIEW FILES ###########################

                tabPanel(title = "View Files",

                  tabsetPanel(
                    id = "GL_VF_tabset",

                     # view file tab
                    tabPanel(title = "File",
                      uiOutput("UIDEContent"),
                    ),
                    
                    # view subnetwork tab
                    tabPanel(title = "Subnetwork", 
                      tableOutput("subnetwork")
                    ),

                   
                  ),
                ),

                ########################### CLUSTER GENES ###########################

                tabPanel(title = "Cluster Genes",

                  mainPanel(
                    h3("Cluster Genes"),

                    # options dropdown
                    dropdown(

                      inputId = "GL_CG_options",

                      # dropdown characteristics
                      style = "minimal", icon = "OPTIONS",
                      status = "primary", width = "300px", size = "sm",

                      # select plots
                      awesomeCheckboxGroup(
                        inputId = "GL_CG_options_plots",
                        label = tags$h4("Select Plots"), 
                        choices = c("Network", "Heatmap", "Binarized Heatmap"),
                        status = ""
                      ),

                      # run button
                      actionButton(inputId = "GL_CG_run", label = "Run",),
                    ),  

                    # error message
                    p(id = "GL_CG_error", "Please generate a subnetwork in OPTIONS.", style = "color:red"),
                    br(),

                    tabsetPanel(

                      # plots tab
                      tabPanel(title = "Plots",
                        br(),

                        # network
                        conditionalPanel(
                          condition = "$.inArray('Network', input.GL_CG_options_plots) > -1", 
                          h5(id = "GL_CG_network_text", "Network of Clustered Genes"), 
                          plotOutput(outputId = "GL_CG_network_plot", height = "500px"),
                          br(),
                        ),
                        
                        # heatmap
                        conditionalPanel(
                          condition = "$.inArray('Heatmap', input.GL_CG_options_plots) > -1", 
                          h5(id = "GL_CG_heatmap_text", "Heatmap of Clustered Genes"),
                          plotOutput(outputId = "GL_CG_heatmap_plot", height = "500px"),
                          br(),
                        ),

                        # binarized heatmap
                        conditionalPanel(
                          condition = "$.inArray('Binarized Heatmap', input.GL_CG_options_plots) > -1", 
                          h5(id = "GL_CG_bheatmap_text", "Binarized Heatmap of Clustered Genes"), 
                          plotOutput(outputId = "GL_CG_bheatmap_plot", height = "500px"), 
                          br(),
                        ),

                      
                      ),

                      # tables tab
                      tabPanel(title = "Tables",
                        br(),

                        # cluster genes table
                        fluidRow(column(11, dataTableOutput("GL_CG_table_plot"),)),
                      ),
                      
                    ), 

                  ), 
                ), 

                ########################### GENE CONNECTIVITY ###########################

                tabPanel(title = "Gene Connectivity", 

                  mainPanel(
                    h3("Gene Connectivity"),
                    
                    # options dropdown
                    dropdown(
                      inputId = "GL_GC_options",

                      # dropdown characteristics
                      style = "minimal", icon = "OPTIONS",
                      status = "primary", width = "300px", size = "sm",

                      # select plots
                      awesomeCheckboxGroup(
                        inputId = "GL_GC_options_plots",
                        label = tags$h4("Select Plots"), 
                        choices = c("Density", "Histogram", "Clustered Density", "Clustered Histogram"),
                        status = ""
                      ),
                      
                      # breaks slider
                      conditionalPanel(
                        condition = "$.inArray('Histogram', input.GL_GC_options_plots) > -1 || $.inArray('Clustered Histogram', input.GL_GC_options_plots) > -1" ,
                        tags$h4("Other"),
                        sliderInput(
                          inputId = "GL_GC_options_xybreaks", 
                          label = "Number of breaks for histogram:",
                          min = 10, max = 100, value = 100, step = 10,
                        ),
                      ),
                      
                      # run button
                      actionButton(inputId = "GL_GC_run", label = "Run"),

                    ),


                    # error message
                    p(id = "GL_GC_error", "Please generate a subnetwork in OPTIONS.", style = "color:red"),
                    br(),

                    # density
                    conditionalPanel(
                      condition = "$.inArray('Density', input.GL_GC_options_plots) > -1", 
                      h5(id = "GL_GC_density_text", "Density Plot of Gene Connectivity"),
                      plotOutput(outputId = "GL_GC_density_plot", height = "500px",),
                      br(),
                    ),

                    # histogram
                    conditionalPanel(
                      condition = "$.inArray('Histogram', input.GL_GC_options_plots) > -1", 
                      h5(id = "GL_GC_hist_text", "Histogram of Gene Connectivity"),
                      plotOutput(outputId = "GL_GC_hist_plot", height = "500px",),
                      br(),
                    ),

                    # density (subset by clusters)
                    conditionalPanel(
                      condition = "$.inArray('Clustered Density', input.GL_GC_options_plots) > -1", 
                      h5(id = "GL_GC_density_subset_text", "Density plot of Gene Connectivity subset by their clusters"), 
                      plotOutput(outputId = "GL_GC_density_subset_plot", height = "500px",),
                      br(),
                    ),

                    # histogram (subset by clusters)
                    conditionalPanel(
                      condition = "$.inArray('Clustered Histogram', input.GL_GC_options_plots) > -1", 
                      h5(id = "GL_GC_hist_subset_text", "Histogram of Gene Connectivity subset by their clusters"), 
                      plotOutput(outputId = "GL_GC_hist_subset_plot", height = "500px",),
                      br(),
                    ),

                  ),
                ), 

                ########################### FUNCTIONAL OUTLIERS ###########################

                tabPanel(title = "Functional Outliers",

                  mainPanel(
                    h3("Functional Outliers"),

                    # options dropdown
                    dropdown(
                      inputId = "GL_FO_options",

                      # dropdown characteristics
                      style = "minimal", icon = "OPTIONS",
                      status = "primary", width = "300px", size = "sm",

                      # select plots
                      awesomeCheckboxGroup(
                        inputId = "GL_FO_options_plots",
                        label = tags$h4("Select Plots"),
                        choices = c("Network", "Heatmap"),
                        status = ""
                      ),

                      # select tables
                      awesomeCheckboxGroup(
                          inputId = "GL_FO_options_tables",
                          label = tags$h4("Select Tables"),
                          choices = c("Functional Outliers", "Genes in Module"),
                          status = ""
                      ),

                      # other options
                      tags$h4("Other"),
                      
                      # filt_min slider
                      sliderInput("GL_FO_filtmin", label = "Number of Genes to form Module",
                          min = 0, max = 20, value = 6, step = 1
                      ),

                      # run button
                      actionButton(inputId = "GL_FO_run", label = "Run"),

                    ),
                    
                    # error message
                    p(id = "GL_FO_error", "Please generate a subnetwork in OPTIONS.", style = "color:red"),

                  ), 

                  br(), 

                  
                  tabsetPanel(

                    # plots tab
                    tabPanel(title="Plots",
                      br(),

                      # heatmap
                      conditionalPanel(
                        condition = "$.inArray('Network', input.GL_FO_options_plots) > -1", 
                        h5(id = "GL_FO_network_text", "Network"), 
                        plotOutput(outputId = "GL_FO_network_plot", height = "500px"),
                        br(),
                      ),

                      # network
                      conditionalPanel(
                        condition = "$.inArray('Heatmap', input.GL_FO_options_plots) > -1", 
                        h5(id = "GL_FO_heatmap_text", "Heatmap"), 
                        plotOutput(outputId = "GL_FO_heatmap_plot", height = "500px"),
                        br(),
                        
                      ),

                    ),

                    # tables tab
                    tabPanel(title="Tables", 
                      br(),

                      # selected genes table output
                      conditionalPanel(
                        condition = "$.inArray('Genes in Module', input.GL_FO_options_tables) > -1", 
                        h5(id = "GL_FO_in_text", "Genes in Module"), 
                        fluidRow(column(11, dataTableOutput("GL_FO_in_table"))),
                        br(),
                      ),

                      # unselected genes table output
                      conditionalPanel(
                        condition = "$.inArray('Functional Outliers', input.GL_FO_options_tables) > -1", 
                        h5(id = "GL_FO_out_text", "Outliers"), 
                        fluidRow(column(11, dataTableOutput("GL_FO_out_table"))),
                        br(),
                      ),
                    )
                  ),


                ),

                ########################### GSEA ###########################

                tabPanel(title = "Gene Set Enrichment Analysis",

                  mainPanel(
                    h3("Gene Set Enrichment Analysis"),
                    
                    # options dropdown
                    dropdown(
                      inputId = "GL_GSEA_options",

                      # dropdown characteristics
                      style = "minimal", icon = "OPTIONS",
                      status = "primary", width = "300px", size = "sm",
                      
                      # standard GSEA options
                      awesomeCheckboxGroup(
                        inputId = "GL_GSEA_options_std",
                        label = tags$h4("Standard GSEA"), 
                        choices = c("P-value Heatmap"),
                        status = ""
                      ),
                      
                      # run button
                      actionButton(inputId = "GL_GSEA_run", label = "Run"),

                    ),

                    # error message
                    p(id = "GL_GSEA_error", "Please generate a subnetwork in OPTIONS.", style = "color:red"),
                    br(),
                  
                    # heatmap
                    conditionalPanel(
                      condition = "$.inArray('P-value Heatmap', input.GL_GSEA_options_std) > -1", 
                      h5(id = "GL_GSEA_heatmap_text", "P-value Heatmap"), 
                      plotOutput(outputId = "GL_GSEA_heatmap_plot", height = "500px"),
                      br(),
                    ),
                  ),
                ),

              ),

            ),

  )
)
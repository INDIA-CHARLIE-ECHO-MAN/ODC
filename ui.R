library(shinyWidgets)
library(shinythemes)
library(bslib)
library(DT)
library(shiny)
library(rhdf5)
ui <- fluidPage(
 
  titlePanel(title=div(img(src="ODClogo.png", height = 80), "OutDeCo")),
  theme = bs_theme(version = 3, bootswatch = "sandstone"),
  
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
                  h6("Outlier detection through co-expression. Using the tools on this website you can:"),
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

                # title of sidepanel
                 tags$h3("Options"),

                 # inputs in the sidepanel
                fileInput("file1", "Choose DE File",
                  accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                   ".csv")
                ),

                # side panel characteristics
                style = "jelly", icon = "FILE UPLOAD",
                status = "success", width = "300px", size = "sm",
               ),

               navlistPanel(
                 tabPanel(
                   title="wilcox",
                   "wilcox placeholder",
                 ),
                 
                 tabPanel(
                   title="DESeq",
                   "DESeq placeholder",
                 ),
                 
                 tabPanel(
                   title="edgeR",
                   "edgeR Placeholder",
                 ),
               ),
             ),
            
            # Assess DE tab
             tabPanel(
              title="Assess DE",
              dropdown(

                # title of sidepanel
                 tags$h3("Options"),

                 # inputs in the sidepanel
                fileInput("DEFile", "Choose DE File",
                  accept = c(
                   ".csv",
                   ".tsv",
                   ".txt"
                   )
                ),

                # button for selecting delimiter, default is nothing until file is selected and handled in server side
                radioButtons(inputId = 'sepButton', label = 'Delimiter Selector', choices = c(Default=''), selected = ''),
    
                pickerInput(
                            inputId = "networkSelect",
                            label = "Network Selection",
                            choices = c("Blood", "Brain", "Generic"),
                ),
                textInput(inputId = 'chooseChrome', label = 'Choose Chromosome' , value = "", width = NULL, placeholder = "chrX"),
                textInput(inputId = 'chooseGeneNo', label = 'Choose Number of Genes' , value = "", width = NULL, placeholder = "100"),
                actionButton("subnets", "Generate Subnetwork"),

                # side panel characteristics
                style = "jelly", icon = "FILE UPLOAD",
                status = "success", width = "300px", size = "sm",

               ),

        
              navlistPanel(

                tabPanel(
                  title="View File", "View File Page", 
 
                  tabsetPanel(
                    tabPanel(
                      title="Sub Network", 
                      tableOutput("subnetwork")
                    ),
                    tabPanel(
                      title="Plot 3"
                    )
                  ),
                ),

                tabPanel(
                  title="Cluster Genes",
                  "Cluster genes Page",

                  # Navigation Bar for types of plots inside cluster
                  tabsetPanel(
                    tabPanel(
                      title="Cluster Options"
                    ),
                    tabPanel(
                      title="Plot 3"
                    )
                  ),
                ),
                
                tabPanel(
                  title="Gene Connectivity",
                  "Gene Connectivity Page",
                ),
                 
                 tabPanel(
                   title="Functional Outliers",
                   "Functional Outliers Page",
                 ),
                 
                 tabPanel(
                   title="Gene Set Enrichment Analysis",
                   "GSE Page",
                 ),
               ),
             )
  ),
)

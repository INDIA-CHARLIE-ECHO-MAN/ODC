library(shinyWidgets)
library(DT)
library(shiny)
ui <- fluidPage(
  titlePanel(title=div(img(src="ODClogo.png", height = 50), "OutDeCo")),
  
  #navbarPage is top menu bar
  navbarPage("",

            #tabPanel is each tab in the navbarPage
            # home tab
             tabPanel(
               title="Home",

              # navlistPanel is each tab on the side menu for each tabPanel
               navlistPanel(
                 tabPanel(
                   title="About",
                   "What is OutDeCo?"       
                 ),
                 
                 tabPanel(
                   title="Differential Analysis",
                   "information about DE",
                 ),
                 
                 tabPanel(
                   title="Cluster Genes",
                   "information about Cluster genes",
                 ),
                 
                 tabPanel(
                   title="Gene Connectivity",
                   "information about Gene Connectivity",
                 ),
                 
                 tabPanel(
                   title="Functional Outliers",
                   "information about Functional Outliers",
                 ),
                 
                 tabPanel(
                   title="Gene Set Enrichment Analysis",
                   "information about Enrichment Analysis",
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
                style = "gradient", icon = icon("cog"),
                status = "primary", width = "300px",
                animate = animateOptions(
                enter = animations$fading_entrances$fadeInLeftBig,
                exit = animations$fading_exits$fadeOutLeftBig
                )
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
            navbarMenu("Assess DE Menu",
                        
              tabPanel(
                title="Cluster Genes",
                "Cluster genes Page",
                sidebarPanel(
                position = c("left"),
                width = c(4),
                # inputs in the sidepanel
                fileInput("DEFile", "Choose DE File",
                  accept = c(
                  ".csv",
                  ".tsv",
                  ".txt"
              )
              ),
                             
              radioButtons(inputId = 'sep', label = 'Separator', choices = c(Comma=',',Semicolon=';',Tab='\t', Space=''), selected = ','),
                          
              ),
                           
                           mainPanel(
                           # Navigation Bar for types of plots inside cluster
                           tabsetPanel(
                             tabPanel(
                               title="View file",
                               mainPanel(
                                 uiOutput("UIDEContent") 
                               )
                               
                             ),
                             tabPanel(
                               title="Plot 2"
                             ),
                             tabPanel(
                               title="Plot 3"
                             )
                           ),
              )
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
             
  ),
)

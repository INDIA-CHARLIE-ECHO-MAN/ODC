library(shinyWidgets)
library(DT)
library(shiny)
ui <- fluidPage(
  titlePanel(title=div(img(src="ODClogo.png", height = 50), "OutDeCo")),
  
  #navbarPage is top menu bar
  navbarPage("",

            #tabPanel is each tab in the navbarPage          
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
                  title="Cluster Genes",
                  "Cluster genes Page",

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
                ),
               ),
             )
  ),
)

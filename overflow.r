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

server <- function(input, output, session) {
  output$distPlot <- renderPlot({
    hist(rnorm(input$obs), col = 'darkgray', border = 'white')
  })
  
  observe({
    # DEFile from fileInput() function
    ServerDEFile <- req(input$DEFile)
    
    # extensions tool for format validation
    extDEFile <- tools::file_ext(ServerDEFile$datapath)
    if(is.null(input$DEFile)){return()
    }else{
      if (extDEFile == "txt") {
        label = paste("Delimiters for", extDEFile, "file")
        choice <-c(Comma=",", Semicolon=";", Tab="\t", Space=" ")
      }else if (extDEFile == "tsv") {
        label = paste("Delimiter: Tab")
        choice <- (Tab="\t")
      }else {
        label = paste("Delimiter: Comma")
        choice <- (Comma=",")
      }
      updateRadioButtons(session, "sepButton", label = label, choices = choice)
    }
  })
  
  # reactive converts the upload file into a reactive expression known as data
  data <- reactive({

    # DEFile from fileInput() function
    ServerDEFile <- input$DEFile

    # extensions tool for format validation
    extDEFile <- tools::file_ext(ServerDEFile$datapath)

    # file format checking
    req(ServerDEFile)
    validate(need(extDEFile == c("csv", "tsv", "txt"), "Please upload a csv, tsv or txt file."))

    # convert data into file format
    if(is.null(extDEFile)){return()}

    read.table(file=ServerDEFile$datapath, sep=input$sepButton)
  })

  # creates reactive table called DEFileContent
  output$DEFileContent <- renderTable({
    if(is.null(data())){return ()}
    data()
  })

  # handles rendering of reactive object on tb on ui
  output$UIDEContent <- renderUI({
    tableOutput("DEFileContent")
  })
}

shinyApp(ui = ui, server = server)
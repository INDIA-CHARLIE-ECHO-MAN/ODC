server <- function(input, output, session) {
  output$distPlot <- renderPlot({
    hist(rnorm(input$obs), col = 'darkgray', border = 'white')
  })
  # reactive converts the upload file into a reactive expression known as data
  data <- eventReactive(input$DEFile,{

  # DEFile from fileInput() function
  ServerDEFile <- input$DEFile

  # extensions tool for format validation
  extDEFile <- tools::file_ext(ServerDEFile$datapath)

  # file format checking
  req(ServerDEFile)
  validate(need(extDEFile == c("csv", "tsv", "txt"), "Please upload a csv, tsv or txt file."))

  # convert data into file format
  if(is.null(extDEFile)){return()} 

  if (extDEFile == "txt") {
    choice <-c(Comma=",", Semicolon=";", Tab="\t", Space=" ")
    updateRadioButtons(session, "sepButton",
                      label = paste("Delimiters for", extDEFile, "file"),
                      choices = choice,
                      )
  }
  else if (extDEFile == "tsv") {
    choice <- (Tab="\t")
    updateRadioButtons(session, "sepButton",
              label = paste("Delimiter: Tab"),
              choices = choice
              )
  }
  else {
    choice <- (Comma=",")
    updateRadioButtons(session, "sepButton",
              label = paste("Delimiter: Comma"),
              choices = choice
              )
  }

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

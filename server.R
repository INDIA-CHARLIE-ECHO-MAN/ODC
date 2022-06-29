server <- function(input, output, session) {
  output$distPlot <- renderPlot({
    hist(rnorm(input$obs), col = 'darkgray', border = 'white')
  })
  
  output$upregheatmap <- renderPlot({
    gene_list <- sample(  EGAD::attr.human$name[EGAD::attr.human$chr=="chrX"], 100 )
    network_type <- 'generic'
    sub_nets <- subset_network_hdf5_gene_list(gene_list, network_type, dir="./ODC_backend1/")
    sub_net <- sub_nets$sub_net
    node_degrees <-  sub_nets$node_degrees
    medK <-  as.numeric(sub_nets$median)
    clust_net <- list() 
    clust_net[["genes"]]  <- cluster_coexp( sub_net$genes, medK = medK, flag_plot = FALSE )
    plot_coexpression_heatmap( sub_net$gene, clust_net$gene)

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
    # validate(need(extDEFile == c("csv", "tsv", "txt"), "Please upload a csv, tsv or txt file."))

    # convert data into file format
    if(is.null(extDEFile)){return()}

    read.table(file=ServerDEFile$datapath, sep=input$sepButton, header=TRUE)
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
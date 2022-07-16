source("./src/plot_functions.R", local = TRUE)
source("./src/cluster_coexp.R", local = TRUE)
source("./src/subset_network_hdf5.R", local = TRUE)
source("./src/calc_DE.R", local = TRUE)

server <- function(input, output, session) {
  # Removing elements that are not functional without subnetwork
  # hide(id = "runGC")
  # hide(id = "run")
  hide(id = "GC_dropdown")
  hide(id = "cluster_dropdown")
  hide(id = "CG_dropdown")
  hide(id = "FO_dropdown")

  # Run DE hide
  #hide(id = "select_column")
  #hide(id = "select_case")  

  #hide(id ="run_DE")
  
  #################### Assess DE upload DE data ###################
  # reactive converts the upload file into a reactive expression known as data
  DEData <- reactive({

    # DEFile from fileInput() function
    ServerDEFile <- input$DEFile

    # extensions tool for format validation
    extDEFile <- tools::file_ext(ServerDEFile$datapath)

    # file format checking
    req(ServerDEFile)
    validate(need(extDEFile == c("csv", "tsv", "txt"), "Please upload a csv, tsv or txt file."))

    # convert data into file format
    if (is.null(extDEFile)) {
      return ()
    }

    read.table(file=ServerDEFile$datapath, sep=input$septxtButton, header=TRUE, nrows=5)
  })

  # creates reactive table called DEFileContent
  output$DEFileContent <- renderTable({
    if (is.null(DEData())) {
      return ()
    }
    DEData()
  })

  # handles rendering of reactive object on tb on ui
  output$UIDEContent <- renderUI({
    tableOutput("DEFileContent")
  })

  labelsData <- reactive({
    # Labels File from fileInput() function
    server_labels_file <- input$labels_file

    # extensions tool for format validation
    ext_labels_file <- tools::file_ext(server_labels_file$datapath)

    # file format checking
    req(server_labels_file)
    validate(need(ext_labels_file == c("csv", "tsv", "txt"), "Please upload a csv, tsv or txt file."))
  
    # convert data into file format
    if (is.null(ext_labels_file)) {
      return ()
    }
    read.table(file=server_labels_file$datapath, sep=input$sepDEButton, header=TRUE)
  })

  observe({
      # DEFile from fileInput() function
      ServerDEFile <- req(input$DEFile)
      
      # extensions tool for format validation
      extDEFile <- tools::file_ext(ServerDEFile$datapath)
      if (is.null(input$DEFile)) {
        return ()
      } else{
        if (extDEFile == "txt") {
          label = paste("Delimiters for", extDEFile, "file")
          choice <-c(Comma=",", Semicolon=";", Tab="\t", Space=" ")
        } else if (extDEFile == "tsv") {
          label = paste("Delimiter: Tab")
          choice <- (Tab="\t")
        } else {
          label = paste("Delimiter: Comma")
          choice <- (Comma=",")
        }
        updateRadioButtons(session, "septxtButton", label = label, choices = choice)
      }
    })
  ########################### RUN DE UPLOAD COUNTS DATA #######################
  # Make countsData
  countsData <- reactive({
    ServerCountsFile <- input$counts_file
    extDEFile <- tools::file_ext(ServerDEFile$datapath)
    req(ServerDEFile)
    validate(need(extDEFile == c("csv", "tsv", "txt"), "Please upload a csv, tsv or txt file."))
    if (is.null(extDEFile)) {
      return ()
    }
    read.table(file=ServerDEFile$datapath, sep=input$septxtButton, header=TRUE, nrows=5)
  })
  
  # creates reactive table called CountsFileContent
  output$CountsFileContent <- renderTable({
    if (is.null(countsData())) {
      return ()
    }
    countsData()
  })

  # handles rendering of reactive object on tb on ui
  output$UIRawContent <- renderUI({
    tableOutput("CountsFileContent")
  })

  observe({
    # counts_file from fileInput() function
    ServerCountsFile <- req(input$counts_file)
    
    # extensions tool for format validation
    extCountsFile <- tools::file_ext(ServerCountsFile$datapath)
    if (is.null(input$counts_file)) {
      return ()
    } else{
      if (extCountsFile == "txt") {
        label = paste("Delimiters for", extCountsFile, "file")
        choice <-c(Comma=",", Semicolon=";", Tab="\t", Space=" ")
      } else if (extCountsFile == "tsv") {
        label = paste("Delimiter: Tab")
        choice <- (Tab="\t")
      } else {
        label = paste("Delimiter: Comma")
        choice <- (Comma=",")
      }
      updateRadioButtons(session, "septxtButton", label = label, choices = choice)
      }
    })

  # Jaz's countsData
#  countsData <- reactive({
#    if ( is.null(input$counts_file)) return(NULL)
#    inFile <- input$counts_file
#    file <- inFile$datapath
    # load the file into new environment and get it from there
#    e = new.env()
#    name <- load(file, envir = e)
#    data <- e[[name]]
#  })
  ##################### RUN DE UPLOAD LABELS DATA ###########################

  # Plot the data
  observeEvent(input$labels_file, {
      show(id = "select_column")
      show(id = "select_case")

      show(id = "run_DE")
      options <- names(labelsData())
      updateSelectInput(session, inputId="select_column","Select column to group", choices = options[2:length(options)], selected = NULL)
      })

  observeEvent(input$select_column, {
      var <- labelsData()[[input$select_column]]
      lvl <- levels(as.factor(var))
      updateSelectInput(session, inputId="select_case", "Select case to analyse", choices = lvl, selected = NULL)

  })

  # creates reactive table called labelsFileContent
  output$labelsFileContent <- renderTable({
    if (is.null(labelsData())) {
      return ()
    }
  })

  # Output labels file
  output$UILabelsContent <- renderUI({
    tableOutput("labelsFileContent")
  })

  observe({
    # labels_file from fileInput() function
    ServerLabelsFile <- req(input$labels_file)
    
    # extensions tool for format validation
    extLabelsFile <- tools::file_ext(ServerLabelsFile$datapath)
    if (is.null(input$labels_file)) {
      return ()
    } else{
      if (extLabelsFile == "txt") {
        label = paste("Delimiters for", extLabelsFile, "file")
        choice <-c(Comma=",", Semicolon=";", Tab="\t", Space=" ")
      } else if (extLabelsFile == "tsv") {
        label = paste("Delimiter: Tab")
        choice <- (Tab="\t")
      } else {
        label = paste("Delimiter: Comma")
        choice <- (Comma=",")
      }
      updateRadioButtons(session, "sepLabelButton", label = label, choices = choice)
      }
    })

  ###################################################################

  # reactive converts the upload file into a reactive expression known as data
  LabelData <- reactive({

    # RawFile from fileInput() function
    ServerLabelFile <- input$labels_file

    # extensions tool for format validation
    extLabelFile <- tools::file_ext(ServerLabelFile$datapath)

    # file format checking
    req(ServerLabelFile)
     validate(need(extLabelFile == c("csv", "tsv", "txt"), "Please upload a csv, tsv or txt file."))

    # convert data into file format
    if (is.null(extLabelFile)) {
      return ()
    }

    read.table(file=ServerLabelFile$datapath, sep=input$septxtButton, header=TRUE, nrows=5)
  })

  # creates reactive table called LabelFileContent
  output$LabelFileContent <- renderTable({
    if (is.null(LabelData())) {
      return ()
    }
    LabelData()
  })

  # handles rendering DT table of labels file
  output$UILabelContent <- renderDataTable(
    LabelData()    
  )

  # Run DE
  observeEvent(input$run_DE, {
    labels <- labelsData()
    counts_data <- countsData()
  
    # var <- labelsData()[[input$select_column]]

    var <- input$select_column
    case <- input$select_case

    # Format labels$var
    labels_var <- labels[[paste0(var)]]

    #Initialise the variables of the chosen column to all be 1
    groups <- rep(1, length(labels_var))
    # Pick the case, relabel as 2
    groups[labels_var == case] = 2   

    #Uncomment for Sex/Gender - doesn't work with status
    #groups[labels$Family == 1] <- 0
    #groups[labels$Relationship == "prb"] <- 0

    filt = groups != 0 
    deg = calc_DE(counts_data[,filt], groups[filt], input$DE_method) 
    output$DEplot <- renderPlot(
            {plot( deg$degs$log2_fc, -log10(deg$degs$pvals),  
            pch=19, bty="n", 
            xlab="log2 FC", ylab="-log10 p-vals" )},
            width = 450,
            height = 450
    )
    output$DEplot_average <- renderPlot(
            {plot( log2(deg$degs$mean_cpm),  deg$degs$log2_fc,  
            pch=19, bty="n", 
            ylab="log2 FC", xlab="Average expression (log2 CPM + 1)")},
            width = 450,
            height = 450
    )
    }
  )

  # Add the Run buttons 
  observeEvent(
    input$generate_subnet,
    {
      # cluster genes
      show(id = "CG_dropdown")
      show(id = "run")
      hide(id = "CG_error")

      # gene connectivity
      show(id = "GC_dropdown")
      show(id = "runGC")
      hide(id = "GC_error")

      # functional outliers
      show(id = "FO_dropdown")
      show(id = "runFO")
      hide(id = "FO_error")

    }
  )



  sn <- reactiveValues(sub_nets = NULL)
  # sub_nets
  sn <- reactiveValues(
    sub_nets = NULL,
  )

  # generate sub_nets
  observeEvent(
    input$generate_subnet,
    {
      # determine gene_list
      if (input$gene_list_selection == "Generate Gene List") {
        validate(
          need(input$chooseChrome, 'Please enter a Chromosme between 1-22, X, Y'),
          need(input$chooseGeneNo != "" && input$chooseGeneNo > 0, 'Please enter a valid Gene Number'),
        )
        gene_list <- sample(EGAD::attr.human$name[EGAD::attr.human$chr==input$chooseChrome], input$chooseGeneNo)
      } else {
        gene_list <- read.delim(file = input$DEFile$datapath, header = FALSE, sep = "\n", dec = ".")[,1]
      }
      # generate sub_nets & store under sn
      sn$sub_nets <- subset_network_hdf5_gene_list(gene_list, tolower(input$network_type), dir="../networks/")
    }
  )


  # clust_net
  clust_net <- reactive({
    sub_net <- sn$sub_nets$sub_net
    node_degrees <- sn$sub_nets$node_degrees
    medK <- as.numeric(sn$sub_nets$median)

    clust_net <- list() 
    clust_net[["genes"]] <- cluster_coexp( sub_net$genes, medK = medK, flag_plot = FALSE )
    return(clust_net)

  })
  
  
  # Output of subnetowrk table
  observeEvent(
    input$generate_subnet, 
    {output$subnetwork <- renderTable(sn$sub_nets)}
  )
  
  







  # CLUSTER GENES

  observeEvent(
    {input$run},
    {
      sub_net <- sn$sub_nets$sub_net
      node_degrees <- sn$sub_nets$node_degrees
      medK <- as.numeric(sn$sub_nets$median)
     

      # network output
      output$CNtext = renderText("Network of Clustered Genes")
      output$network <- renderPlot(
        {plot_network( sub_net$gene, clust_net()$gene, medK )},
        width = 500,
        height = 500
      )


      # heatmap output
      output$CHtext = renderText("Heatmap of Clustered Genes")
      output$heatmap <- renderPlot(
        {plot_coexpression_heatmap( sub_net$gene, clust_net()$gene, flag_plot_bin = FALSE)},
        width = 500,
        height = 500
      )


      # binarized heatmap output
      output$CHBtext = renderText("Binarized Heatmap of Clustered Genes")
      output$Bheatmap <- renderPlot(
        {plot_coexpression_heatmap( sub_net$gene, clust_net()$gene )},
        width = 500,
        height = 500
      )
      

      # clustering genes table output
      output$CG_table_text = renderText({
        # input$runGC
        # req(input$runGC) # to prevent print at first lauch
        # isolate(print("Filtered Genes"))
        "Clustering Genes"
      })
      output$CG_table <- renderDT(
        {EGAD::attr.human[match(clust_net()$genes$clusters$genes[!genes_keep],EGAD::attr.human$name[EGAD::attr.human$chr==input$chooseChrome], input$chooseGeneNo),]},
        # options=list(columnDefs = list(list(visible=FALSE, targets=c(0,1,2,3))))
      )


    }
  )



  # GENE CONNECTIVITY

  observeEvent(
    {input$runGC},
    {
      sub_net <- sn$sub_nets$sub_net
      node_degrees <- sn$sub_nets$node_degrees  
      medK <- as.numeric(sn$sub_nets$median)
      m <- match(clust_net()$genes$clusters$genes, rownames(sub_net$genes))


      # density output
      output$GCdensityGtext = renderText("Density plot of Gene Connectivity")
      output$GCdensityG <- renderPlot(
        {plot_scatter(node_degrees$genes[,1]/node_degrees$n_genes_total, 
                    node_degrees$genes[,2]/node_degrees$n_genes, 
                    xlab="Global node degree", 
                    ylab="Local node degree", flag= "density")},
        width = 500,
        height = 500
      )


      # histogram output
      output$GChistogramGtext = renderText("Histogram of Gene Connectivity")
      output$GChistogramG <- renderPlot(
        {plot_scatter(node_degrees$genes[,1]/node_degrees$n_genes_total, 
                    node_degrees$genes[,2]/node_degrees$n_genes, 
                    xybreaks = input$xybreaks,
                    xlab="Global node degree", 
                    ylab="Local node degree", flag= "hist")},
        width = 500,
        height = 500
      )


      # density output - subset by clusters
      output$GCdensitySubsetGtext = renderText("Density plot of Gene Connectivity subset by their clusters")
      output$GCdensitySubsetG <- renderPlot(
        { plot_scatter(node_degrees$genes[m,1]/node_degrees$n_genes_total, 
                      node_degrees$genes[m,2]/node_degrees$n_genes, 
                      xlab="Global node degree", 
                      ylab="Local node degree", 
                      clusters = clust_net()$genes$clusters, flag = "density")},
        width = 500,
        height = 500
      )


      # histogram output - subset by clusters
      output$GChistogramSubsetGtext = renderText("Histogram of Gene Connectivity subset by their clusters")
      output$GChistogramSubsetG <- renderPlot(
        { plot_scatter(node_degrees$genes[m,1]/node_degrees$n_genes_total, 
                      node_degrees$genes[m,2]/node_degrees$n_genes, 
                      xybreaks = input$xybreaks,
                      xlab="Global node degree", 
                      ylab="Local node degree", 
                      clusters = clust_net()$genes$clusters, flag = "hist")},
        width = 500,
        height = 500
      )


    }
  )


  # FUNCTIONAL OUTLIERS

  observeEvent(
    {input$runFO},
    {
      sub_net <- sn$sub_nets$sub_net
      filt_min <- input$filtmin
      clust_size <- plyr::count(clust_net()$genes$clusters$labels)
      clust_keep <- clust_size[clust_size[,2] < filt_min ,1]
      genes_keep <- !is.na(match(clust_net()$genes$clusters$labels, clust_keep))
      medK <- as.numeric(sn$sub_nets$median)


      # heatmap output
      output$FO_heatmap_text = renderText("Heatmap")
      output$FO_heatmap <- renderPlot(
        {plot_coexpression_heatmap(sub_net$genes, clust_net()$genes, filt=TRUE, flag_plot_bin=FALSE)},
        width = 500,
        height = 500
      )

      # network output
      output$FO_network_text = renderText("Network")
      output$FO_network <- renderPlot(
        {plot_network(1-sub_net$genes, clust_net()$genes, 1 - medK)},
        width = 500,
        height = 500
      )

      # unselected genes table output
      output$FO_unselected_text = renderText("Unselected Genes")
      output$genes_unselected_table <- renderDT(
        {EGAD::attr.human[match(clust_net()$genes$clusters$genes[!genes_keep],EGAD::attr.human$name[EGAD::attr.human$chr==input$chooseChrome], input$chooseGeneNo),]},
        # options=list(columnDefs = list(list(visible=FALSE, targets=c(0,1,2,3))))
      )


      # selected genes table output
      output$FO_selected_text = renderText("Selected Genes")
      output$genes_selected_table <- renderDataTable(
        {EGAD::attr.human[match(clust_net()$genes$clusters$genes[genes_keep],EGAD::attr.human$name[EGAD::attr.human$chr==input$chooseChrome], input$chooseGeneNo),]},
        # options=list(columnDefs = list(list(visible=FALSE, targets=c(0,1,2,3))))
      )
      
    }
  )

  # ERROR MESSAGES

  output$CG_error <- renderText({
    print("Please upload/generate a gene list in OPTIONS")
  })

  output$GC_error <- renderText({
    print("Please upload/generate a gene list in OPTIONS")
  }) 

  output$FO_error <- renderText({
    print("Please upload/generate a gene list in OPTIONS")
  }) 

}



source("./src/plot_functions.R", local = TRUE)
source("./src/cluster_coexp.R", local = TRUE)
source("./src/subset_network_hdf5.R", local = TRUE)
source("./src/calc_DE.R", local = TRUE)
source("./src/gene_set_enrichment.R", local = TRUE)

# Warnings silenced for wilcox
options(warn=-1)
defaultW <- getOption("warn")

server <- function(input, output, session) {

  # Removing elements that are not functional without subnetwork
  # hide(id = "GC_run")
  # hide(id = "CG_run")
  hide(id = "GC_dropdown")
  hide(id = "CG_dropdown")
  hide(id = "GC_dropdown")
  hide(id = "FO_dropdown")
  hide(id = "GSEA_dropdown")
  hide(id = "assess_run_de")


  hide(id = "sepLabelsButton")
  hide(id = "sepCountsButton")

  #Run Differential Expression
  hide(id="DE_options")
  hide(id="run_DE")



  # PLOTS/TABLES HEADERS
  # Run DE
  hide(id = "vol")
  hide(id = "MA")
  # Clustering 
  hide(id="CG_network_text")
  hide(id="CG_heatmap_text")
  hide(id="CG_bheatmap_text")
  hide(id="CG_table_text")
  # Gene Connectivity
  hide(id="GCdensityG_text")
  hide(id="GChistogramG_text")
  hide(id="GCdensitySubsetG_text")
  hide(id="GChistogramSubsetG_text")
  # Functional Outliers
  hide(id="FO_network_text")
  hide(id="FO_heatmap_text")
  hide(id="genes_not_keep_table_text")
  hide(id="genes_keep_table_text")
  # GSEA
  hide(id="GSEA_heatmap_text")
  hide(id="GSEA_up_heatmap_text")
  hide(id="GSEA_down_heatmap_text")
  hide(id="GSEA_auroc_text")



  
  ##########################################################################################
  #                                                                                        #
  #                                    RUN DE                                              #
  #                                                                                        #
  ##########################################################################################
  
  
  #####################  UPLOAD COUNTS DATA ###########################

  # Make countsData
  countsData <- reactive({
    ServerCountsFile <- input$counts_file
    extCountsFile <- tools::file_ext(ServerCountsFile$datapath)
    req(ServerCountsFile)
    validate(need(extCountsFile == c("csv", "tsv", "txt"), "Please upload a csv, tsv or txt file."))
    if (is.null(extCountsFile)) {
      return ()
    }
    if (extCountsFile == "csv") {
      read.table(file=ServerCountsFile$datapath, sep=input$sepCountsButton, header=TRUE, row.names = 1)
    } else {
      read.table(file=ServerCountsFile$datapath, sep=input$sepCountsButton, header=TRUE) 
    }
    
    
  })

  observe({
     # counts_file from fileInput() function
    ServerCountsFile <- req(input$counts_file)
    
    # extensions tool for format validation
    extCountsFile <- tools::file_ext(ServerCountsFile$datapath)
    if (is.null(input$counts_file)) {
      return ()
    } else {
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
      updateRadioButtons(session, "sepCountsButton", label = label, choices = choice)
      }

      #print(counts_data[:])
    })

    output$UICountsContent <- renderDataTable(
        countsData(), options = list(
          pageLength = 25
        )
      )
      
    observeEvent(input$counts_file, {
      show(id = "sepCountsButton")
    })

  ########################### UPLOAD LABELS DATA ###########################

  # Make labelsData
  labelsData <- reactive({
    ServerLabelsFile <- input$labels_file
    extLabelsFile <- tools::file_ext(ServerLabelsFile$datapath)
    req(ServerLabelsFile)
    validate(need(extLabelsFile == c("csv", "tsv", "txt"), "Please upload a csv, tsv or txt file."))
    if (is.null(extLabelsFile)) {
      return ()
    }
    if (extLabelsFile == "csv") {
      read.table(file=ServerLabelsFile$datapath, sep=input$sepLabelsButton, header=TRUE, row.names = 1)
    } else {
      read.table(file=ServerLabelsFile$datapath, sep=input$sepLabelsButton, header=TRUE)
    }
    
    
  })

  observeEvent(input$labels_file, {
    show(id = "sepLabelsButton")
  })



  observe({
     # labels_file from fileInput() function
    ServerLabelsFile <- req(input$labels_file)
    
  #   # extensions tool for format validation
    extLabelsFile <- tools::file_ext(ServerLabelsFile$datapath)
    if (is.null(input$labels_file)) {
      return ()
    } else {
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
      updateRadioButtons(session, "sepLabelsButton", label = label, choices = choice)
      }
    })


  # handles rendering DT table of labels file

  output$UILabelContent <- renderDataTable(
    labelsData(), options = list(
      pageLength = 100
    )
  )

  # rendering DT table for RUN DE options (to select cases)
  output$UILabelContentSelection <- renderDataTable(
    labelsData(), options = list(
      pageLength = 100
    )
  )

    # rendering DT table for RUN DE options (to remove cases)
  output$UILabelContentRemoveSelection <- renderDataTable(
    labelsData(), options = list(
      pageLength = 100
    ) 
  )

  # Plot the data
  observeEvent(input$case_control_method, {
      options <- names(labelsData())
      updateSelectInput(session, 
        inputId="select_column",
        "Select column to group", 
        choices = options[2:length(options)], 
        selected = NULL
      )
      show(id="run_DE")
  })

  # Group by label option 
  observeEvent(input$select_column, {
      # Update the case selection with levels of selected column 
      var <- labelsData()[[input$select_column]]
      lvl <- levels(as.factor(var))
      updateSelectInput(session, 
        inputId="select_case", 
        "Select case to analyse", 
        choices = lvl, 
        selected = NULL
      )
  })

  # countsData <- reactive({
  #   if ( is.null(input$counts_file)) return(NULL)
  #   inFile <- input$counts_file
  #   file <- inFile$datapath
  #   # load the file into new environment and get it from there
  #   e = new.env()
  #   name <- load(file, envir = e)
  #   data <- e[[name]]
  # })

  case_selected <- reactive({
    input$UILabelContentSelection_rows_selected
  })

  remove_selected <- reactive({
    input$UILabelContentRemoveSelection_rows_selected
  })

  # Switch to labels tab if labels file is uploaded
  observeEvent(input$labels_file, {
    updateTabsetPanel(session, "counts_labels_tabset", selected = "Labels File")
  })
  # Switch to counts tab if counts file is uploaded
  observeEvent(input$counts_file, {
    updateTabsetPanel(session, "counts_labels_tabset", selected = "Counts File")
  })




  ########################### RUN DE ###########################
  observe({
    if (!is.null(input$labels_file) && !is.null(input$labels_file)) {
      show(id="DE_options")
      hide(id="runDE_error")
    }
  })
  
  de <- reactiveValues(
    deg_output = NULL, 
  )


  observeEvent(
    input$run_DE, 
    {
      labels <- labelsData()
      counts_data <- countsData()
      deg <- NULL

      # var <- labelsData()[[input$select_column]]
      if (input$case_control_method == "Choose Case by Label") {
        var <- input$select_column
        case <- input$select_case

        # Format labels$var
        labels_var <- labels[[paste0(var)]]

        #Initialise the variables of the chosen column to all be 1
        groups <- rep(1, length(labels_var))
        
        # Pick the case, relabel as 2
        groups[labels_var == case] = 2   

        filt = groups != 0 
        deg <- calc_DE(counts_data[,filt], groups[filt], input$DE_method) 
        de$deg_output <- deg

      } else {
        cases <- case_selected()
        cases_removed <- remove_selected()
        

        
        # Initalise all values to 1
        groups <- rep(1, nrow(labels))
        check <- rep(1, nrow(labels))

        for (c in cases) {
          groups[c] = 2
        }
        
        for (d in cases_removed) {
          groups[d] = 0
        }
        # No cases have been selected
        print(groups)
        print(check)
        if (groups == check) {
          shinyalert(title = "Invalid Input", text = "Please select cases to assess", type = "error")
        } else {
          deg <- calc_DE(counts_data, groups, input$DE_method) 
          de$deg_output <- deg
        }
        

      }
      
      if (!is.null(deg)) {

        # Volcano Plot
        show(id="vol")
        output$DEplot <- renderPlot(
                {plot(deg$degs$log2_fc, -log10(deg$degs$pvals),  
                pch=19, bty="n", 
                xlab="log2 FC", ylab="-log10 p-vals" )},
                width = 450,
                height = 450
        )

        # MA Plot
        show(id="MA")
        #output$DE_MA_text = renderText("MA Plot")
        output$DEplot_average <- renderPlot(
                {plot( log2(deg$degs$mean_cpm),  deg$degs$log2_fc,  
                pch=19, bty="n", 
                ylab="log2 FC", xlab="Average expression (log2 CPM + 1)")},
                width = 450,
                height = 450
        )
      }
      show(id = "assess_run_de")
    }
  )

  observeEvent(input$assess_run_de, { 
    updateTabsetPanel(session, inputId="navpage", selected="Assess DE")
    updateTabsetPanel(session, "subnetwork_file_tabset", selected = "Subnetwork")
    updateRadioButtons(session, inputId="gene_list_selection", choices=c("Upload Gene List", "Generate Gene List", "Use DE results"), selected = "Use DE results")
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
        updateRadioButtons(session, "sepButton", label = label, choices = choice)
      }
    })
  
  
  ##########################################################################################
  #                                                                                        #
  #                                        ASSESS DE                                       #
  #                                                                                        #
  ##########################################################################################

  # sub_nets
  sn <- reactiveValues(
    sub_nets = NULL,
  )

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

    read.table(file=ServerDEFile$datapath, sep=input$sepButton, header=TRUE)
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


  # generate sub_nets

  observeEvent(
    input$generate_subnet, 
    {
      # subnetwork from DE results 
      if (input$gene_list_selection == "Use DE results") { 
          sn$sub_nets <- subset_network_hdf5(de$deg_output$degs, tolower(input$network_type), dir="../networks/")
          updateAwesomeCheckboxGroup(session, inputId="CG_PlotOptions", choices=c("Upregulated Network", "Upregulated Heatmap", "Upregulated Binarized Heatmap", "Downregulated Network", "Downregulated Heatmap", "Downregulated Binarized Heatmap"))
          updateAwesomeCheckboxGroup(session, inputId="GSEA_type", choices=c("Standard GSEA", "AUCs GSEA"))
          updateAwesomeCheckboxGroup(session, inputId="GSEA_std_PlotOptions", choices=c("Upregulated P-value Heatmap", "Downregulated P-value Heatmap"))
      } 
      
      # subnetwork from gene_list 
      else { 
        
        # generate gene_list
        if (input$gene_list_selection == "Generate Gene List") {

          chromNum <- strtoi(substring(input$chooseChrome,4))

          # error: invalid chromosome format
          if (!str_detect(input$chooseChrome, "chr[XY]") && !str_detect(input$chooseChrome, "chr[0-9]")) {
            shinyalert(title = "Invalid Input", text = "Please enter a Chromosome between 1 - 22, X, Y", type = "error")
            gene_list <- NULL
          } 
          
          # error: invalid chromosome number
          else if (str_detect(substring(input$chooseChrome,4), "[0-9]") && (chromNum < 1 || chromNum > 22)) {
            shinyalert(title = "Invalid Input", text = "Please enter a Chromosome between 1 - 22, X, Y", type = "error")
            gene_list <- NULL
          } 
          
          # error: invalid number of genes
          else if (input$chooseGeneNo == "" || input$chooseGeneNo <= 0) { 
            shinyalert(title = "Invalid Input", text = "Please enter a valid number of Genes", type = "error")
            gene_list <- NULL
          } 
          
          # generate gene_list from EGAD
          else { 
            gene_list <- sample( EGAD::attr.human$name[EGAD::attr.human$chr==input$chooseChrome], input$chooseGeneNo,)
          } 
        }
        
        # upload gene_list
        else {
          gene_list <- read.delim(file = input$DEFile$datapath, header = FALSE, sep = "\n", dec = ".")[,1]
        }

        # generate sub_nets
        if (!is.null(gene_list)) { 
          sn$sub_nets <- subset_network_hdf5_gene_list(gene_list, tolower(input$network_type), dir="../networks/")
        } 
      }
    }
  )
  
  # view files tab selection
  observeEvent(input$gene_list_selection, {
    if (input$gene_list_selection == "Upload Gene List") {
      updateTabsetPanel(session, "subnetwork_file_tabset", selected = "File")
    } else {
      updateTabsetPanel(session, "subnetwork_file_tabset", selected = "Subnetwork")
    }
  })





  



  # Add the Run buttons 
  observeEvent(
    input$generate_subnet,
    {
      # cluster genes
      show(id = "CG_dropdown")
      show(id = "CG_run")
      hide(id = "CG_error")

      # gene connectivity
      show(id = "GC_dropdown")
      show(id = "GC_run")
      hide(id = "GC_error")

      # functional outliers
      show(id = "FO_dropdown")
      show(id = "FO_run")
      hide(id = "FO_error")
      if (input$gene_list_selection == "Generate Gene List") {
        updateSliderInput(session, inputId = "xybreaks", min = 10, max = 150, value = input$chooseGeneNo, step = 10)
      }
      
      # GSEA
      show(id = "GSEA_dropdown")
      show(id = "GSEA_run")
      hide(id = "GSEA_error")

    }
  )



  # clust_net
  clust_net <- reactive({

    sub_net <- sn$sub_nets$sub_net
    node_degrees <- sn$sub_nets$node_degrees
    medK <- as.numeric(sn$sub_nets$median)

    clust_net <- list() 
    # For DE data 
    if (input$gene_list_selection == "Use DE results") {
      deg_sig <- sn$sub_nets$deg_sig
      fc_sig  <- sn$sub_nets$fc_sig
      clust_net[["down"]]  <- cluster_coexp(sub_net$down, medK = medK, flag_plot = FALSE)
      clust_net[["up"]]  <- cluster_coexp( sub_net$up, medK = medK, flag_plot = FALSE)

    } 
    # For gene list
    else {  
      clust_net[["genes"]] <- cluster_coexp(sub_net$genes, medK = medK, flag_plot = FALSE)
    }
    return(clust_net)
  })
  
  
  # Output of subnetwork table
  observeEvent(
    input$generate_subnet, 
    {output$subnetwork <- renderTable(sn$sub_nets)}
  )

  # output$subnetwork <- renderTable({
  #   sub_nets()
  # })
  
  






  ##################### CLUSTER GENES #####################

  observeEvent(
    {input$run},
    {
      sub_net <- sn$sub_nets$sub_net
      node_degrees <- sn$sub_nets$node_degrees
      medK <- as.numeric(sn$sub_nets$median)

      # network output
      show(id="CG_network_text")
      output$network <- renderPlot(
        {plot_network(sub_net$genes, clust_net()$genes, medK)},
        width = 500,
        height = 500
      )


      # heatmap output
      show(id="CG_heatmap_text")
      output$heatmap <- renderPlot(
        {plot_coexpression_heatmap(sub_net$genes, clust_net()$genes, flag_plot_bin = FALSE)},
        width = 500,
        height = 500
      )


      # binarized heatmap output
      show(id="CG_bheatmap_text")
      output$Bheatmap <- renderPlot(
        {plot_coexpression_heatmap(sub_net$genes, clust_net()$genes)},
        width = 500,
        height = 500
      )

      # upregulated network 
      output$upregNetwork <- renderPlot(
        {plot_network(sub_net$up, clust_net()$up, medK)}, 
        width = 500, 
        height = 500 
      )

      # upregulated heatmap 
      output$upregHeatmap <- renderPlot(
        {plot_coexpression_heatmap(sub_net$up, clust_net()$up, flag_plot_bin = FALSE)}, 
        width = 500,
        height = 500
      )

      # upregulated binarized heatmap 
      output$upregbinHeatmap <- renderPlot(
        {plot_coexpression_heatmap(sub_net$up, clust_net()$up)}, 
        width = 500, 
        height = 500
      )

      output$downregNetwork <- renderPlot(
        {plot_network(sub_net$down, clust_net()$down, medK)},
        width = 500, 
        height = 500
      )

      output$downregHeatmap <- renderPlot(
        {plot_coexpression_heatmap(sub_net$down, clust_net()$down, flag_plot_bin = FALSE)}, 
        width = 500, 
        height = 500 
      )

      output$downregbinHeatmap <- renderPlot(
        {plot_coexpression_heatmap(sub_net$down, clust_net()$down)}, 
        width = 500, 
        height = 500
      )

      # clustering genes table output
      show(id="CG_table_text")
      output$CG_table <- renderDataTable(
        {EGAD::attr.human[match(clust_net()$genes$clusters$genes,EGAD::attr.human$name[EGAD::attr.human$chr==input$chooseChrome],input$chooseGeneNo),]},
        # options=list(columnDefs = list(list(visible=FALSE, targets=c(0,1,2,3))))
      )


    }
  )



  ##################### GENE CONNECTIVITY #####################

  observeEvent(
    {input$GC_run},
    {
      sub_net <- sn$sub_nets$sub_net
      node_degrees <- sn$sub_nets$node_degrees  
      medK <- as.numeric(sn$sub_nets$median)
      m <- match(clust_net()$genes$clusters$genes, rownames(sub_net$genes))


      # density output
      show(id="GCdensityG_text")
      output$GCdensityG <- renderPlot(
        {plot_scatter(node_degrees$genes[,1]/node_degrees$n_genes_total, 
                    node_degrees$genes[,2]/node_degrees$n_genes, 
                    xlab="Global node degree", 
                    ylab="Local node degree", flag= "density")},
        width = 500,
        height = 500
      )


      # histogram output
      show(id="GChistogramG_text")
      output$GChistogramG <- renderPlot(
        {plot_scatter(node_degrees$genes[,1]/node_degrees$n_genes_total, 
                    node_degrees$genes[,2]/node_degrees$n_genes, 
                    xybreaks = input$xybreaks,
                    xlab="Global node degree", 
                    ylab="Local node degree", flag= "hist")},
        width = 500,
        height = 500
      )

      show(id="GCdensitySubsetG_text")
      # density output - subset by clusters
      output$GCdensitySubsetG <- renderPlot(
        {plot_scatter(node_degrees$genes[m,1]/node_degrees$n_genes_total, 
                      node_degrees$genes[m,2]/node_degrees$n_genes, 
                      xlab="Global node degree", 
                      ylab="Local node degree", 
                      clusters = clust_net()$genes$clusters, flag = "density")},
        width = 500,
        height = 500
      )


      # histogram output - subset by clusters
      show(id="GChistogramSubsetG_text")
      output$GChistogramSubsetG <- renderPlot(
        {plot_scatter(node_degrees$genes[m,1]/node_degrees$n_genes_total, 
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


  ##################### FUNCTIONAL OUTLIERS #####################

  observeEvent(
    {input$FO_run},
    {
      sub_net <- sn$sub_nets$sub_net
      node_degrees <- sn$sub_nets$node_degrees  
      medK <- as.numeric(sn$sub_nets$median)

      sub_net <- sn$sub_nets$sub_net
      filt_min <- input$filtmin
      clust_size <- plyr::count(clust_net()$genes$clusters$labels)
      clust_keep <- clust_size[clust_size[,2] < filt_min ,1]
      genes_keep <- !is.na(match(clust_net()$genes$clusters$labels, clust_keep))
      medK <- as.numeric(sn$sub_nets$median)


      # heatmap output
      show(id="FO_heatmap_text")
      output$FO_heatmap <- renderPlot(
        {plot_coexpression_heatmap(sub_net$genes, clust_net()$genes, filt = TRUE, flag_plot_bin = FALSE)},
        width = 500,
        height = 500
      )

      # network output
      show(id="FO_network_text")
      output$FO_network <- renderPlot(
        {plot_network(1-sub_net$genes, clust_net()$genes, 1 - medK)},
        width = 500,
        height = 500
      )

      # genes in module table output
      show(id="genes_not_keep_table_text")
      output$genes_not_keep_table <- renderDataTable(
        {EGAD::attr.human[match(clust_net()$genes$clusters$genes[!genes_keep],EGAD::attr.human$name[EGAD::attr.human$chr==input$chooseChrome], input$chooseGeneNo),]},
        # options=list(columnDefs = list(list(visible=FALSE, targets=c(0,1,2,3))))
      )


      # functional outliers table output
      show(id="genes_keep_table_text")
      output$genes_keep_table <- renderDataTable(
        {EGAD::attr.human[match(clust_net()$genes$clusters$genes[genes_keep],EGAD::attr.human$name[EGAD::attr.human$chr==input$chooseChrome], input$chooseGeneNo),]},
        # options=list(columnDefs = list(list(visible=FALSE, targets=c(0,1,2,3))))
      )
      
    }
  )


  ##################### GSEA #####################  
  observeEvent(
    {input$GSEA_run},
    {
      
      # Standard GSEA
      if ("Standard GSEA" %in% input$GSEA_type) {
        data(go_slim)
        data(go_voc)
        
        # heatmap
        show(id="GSEA_heatmap_text")
        output$GSEA_heatmap <- renderPlot(
          {
            filt <- colSums(go_slim) < 5000 & colSums(go_slim) >= 10
            gene_list <- clust_net()$genes$clusters$genes[clust_net()$genes$order]
            go_enrich <- gene_set_enrichment(gene_list, go_slim[filt,], go_voc) 
            plot_gene_set_enrichment(go_enrich, gene_list, go_slim[filt,])
          },
          width = 500,
          height = 500
        )

        # Use DE results
        if (input$gene_list_selection == "Use DE results") {
          # upregulated heatmap
          show(id="GSEA_up_heatmap_text")
          output$GSEA_up_heatmap <- renderPlot(
            {
              filt <- colSums(go_slim_entrez) < 5000 & colSums(go_slim_entrez) >= 10
              gene_list <- clust_net()$up$clusters$genes[clust_net()$up$order]
              go_enrich <- gene_set_enrichment(gene_list, go_slim_entrez[filt,], go_voc) 
              plot_gene_set_enrichment( go_enrich, gene_list, go_slim_entrez[filt,]) 
            },
            width = 500,
            height = 500
          )
          
          # downregulated heatmap
          show(id="GSEA_down_heatmap_text")
          output$GSEA_down_heatmap <- renderPlot(
            {
              filt <- colSums(go_slim_entrez) < 5000 & colSums(go_slim_entrez) >= 10
              gene_list <- clust_net()$down$clusters$genes[clust_net()$down$order]
              go_enrich <- gene_set_enrichment(gene_list, go_slim_entrez[filt,], go_voc) 
              plot_gene_set_enrichment(go_enrich, gene_list, go_slim_entrez[filt,]) 
            },
            width = 500,
            height = 500
          )
        }
        
      }

      # AUCs GSEA
      if ("AUCs GSEA" %in% input$GSEA_type) {
        
        data(go_slim_entrez)
      
        gene_rankings <- order(log10(deg$degs$pvals), abs(deg$degs$log2_fc)) 
        names(gene_rankings) <- rownames(deg$degs)
        gene_rankings_rev <- rank(max(gene_rankings) - gene_rankings) 
        
        m <- match(rownames(go_slim_entrez), names(gene_rankings_rev))
        f.g = !is.na(m)
        f.r = m[f.g]
        gene_sets = go_slim_entrez[f.g,]
        gene_rankings_rev = rank(gene_rankings_rev[f.r])

        gene_set_aucs <- gene_set_enrichment_aucs(gene_sets, gene_rankings_rev) 

        # AUROC graph
        show(id="GSEA_auroc_text")
        output$GSEA_auroc <- renderPlot(
          {plot_gene_set_enrichment_ranked(gene_set_aucs, gene_rankings_rev, gene_list, go_slim_entrez)},
          width = 500,
          height = 500
        )
      }
    }
  )


  ##################### ERROR MESSAGES #####################

  output$CG_error <- renderText({
    print("Please generate a subnetwork in OPTIONS")
  })

  output$GC_error <- renderText({
    print("Please generate a subnetwork in OPTIONS")
  }) 

  output$FO_error <- renderText({
    print("Please generate a subnetwork in OPTIONS")
  })

  output$GSEA_error <- renderText({
    print("Please generate a subnetwork in OPTIONS")
  }) 

}



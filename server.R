source("./src/plot_functions.R", local = TRUE)
source("./src/cluster_coexp.R", local = TRUE)
source("./src/subset_network_hdf5.R", local = TRUE)
source("./src/calc_DE.R", local = TRUE)

# Warnings silenced for wilcox
options(warn=-1)
defaultW <- getOption("warn")


server <- function(input, output, session) {

  # Removing elements that are not functional without subnetwork
  hide(id = "GL_GC_options")
  hide(id = "DE_GC_options")
  hide(id = "CG_dropdown")
  hide(id = "DE_CG_options")
  hide(id = "GL_FO_options")
  hide(id = "DE_FO_options")
  hide(id = "DE_GSEA_dropdown")
  hide(id = "GL_GSEA_options")
  hide(id = "DE_run_assess")


  hide(id = "DE_run_labels_sep")
  hide(id = "DE_run_counts_sep")

  #Run Differential Expression
  hide(id="DE_run_PD_options")
  hide(id="DE_run")



  # PLOTS/TABLES HEADERS
  # Run DE
  hide(id = "DE_vol_text")
  hide(id = "DE_MA_text")
  # Clustering 
  hide(id="GL_CG_network_text")
  hide(id="GL_CG_heatmap_text")
  hide(id="GL_CG_bheatmap_text")
  hide(id="GL_CG_table_text")
  hide(id="DE_CG_up_network_text")
  hide(id="DE_CG_up_heatmap_text")
  hide(id="DE_CG_up_bheatmap_text")
  hide(id="DE_CG_down_network_text")
  hide(id="DE_CG_down_heatmap_text")
  hide(id="DE_CG_down_bheatmap_text")


  # Gene Connectivity
  hide(id="GL_GC_density_text")
  hide(id="GL_GC_hist_text")
  hide(id="GL_GC_density_subset_text")
  hide(id="GL_GC_hist_subset_text")
  hide(id="DE_GC_up_density_text")
  hide(id="DE_GC_up_hist_text")
  hide(id="DE_GC_up_density_subset_text")
  hide(id="DE_GC_up_hist_subset_text")
  hide(id="GL_GC_density_plot_downreg_text")
  hide(id="GL_GC_hist_plot_downreg_text")
  hide(id="GL_GC_density_subset_plot_downreg_text")
  hide(id="GL_GC_hist_subset_plot_downreg_text")

  # Functional Outliers
  hide(id="GL_FO_network_text")
  hide(id="GL_FO_heatmap_text")
  hide(id="DE_FO_up_network_text")
  hide(id="DE_FO_up_heatmap_text")
  hide(id="DE_FO_down_network_text")
  hide(id="DE_FO_down_heatmap_plot")
  hide(id="GL_FO_in_text")
  hide(id="GL_FO_out_text")

  # GSEA
  hide(id="GL_GSEA_heatmap_text")
  hide(id="GSEA_up_heatmap_text")
  hide(id="GSEA_down_heatmap_text")


  ######################################################################
  #                                                                    #
  #                               RUN DE                               #
  #                                                                    #
  ######################################################################
  
  
  #####################  UPLOAD COUNTS DATA ###########################

  # Make countsData
  countsData <- reactive({
    ServerCountsFile <- input$DE_run_counts_file
    extCountsFile <- tools::file_ext(ServerCountsFile$datapath)
    req(ServerCountsFile)
    validate(need(extCountsFile == c("csv", "tsv", "txt"), "Please upload a csv, tsv or txt file."))
    if (is.null(extCountsFile)) {
      return ()
    }
    if (extCountsFile == "csv") {
      read.table(file=ServerCountsFile$datapath, sep=input$DE_run_counts_sep, header=TRUE, row.names = 1)
    } else {
      read.table(file=ServerCountsFile$datapath, sep=input$DE_run_counts_sep, header=TRUE) 
    }
    
    
  })

  observe({
     # DE_run_counts_file from fileInput() function
    ServerCountsFile <- req(input$DE_run_counts_file)
    
  #   # extensions tool for format validation
    extCountsFile <- tools::file_ext(ServerCountsFile$datapath)
    if (is.null(input$DE_run_counts_file)) {
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
      updateRadioButtons(session, "DE_run_counts_sep", label = label, choices = choice)
      }

      #print(counts_data[:])
    })

    output$UICountsContent <- renderDataTable(
        countsData(), options = list(
          pageLength = 25
        )
      )
      
    observeEvent(input$DE_run_counts_file, {
      show(id = "DE_run_counts_sep")
    })

  ########################### UPLOAD LABELS DATA ###########################
  # Make labelsData
  labelsData <- reactive({
    ServerLabelsFile <- input$DE_run_labels_file
    extLabelsFile <- tools::file_ext(ServerLabelsFile$datapath)
    req(ServerLabelsFile)
    validate(need(extLabelsFile == c("csv", "tsv", "txt"), "Please upload a csv, tsv or txt file."))
    if (is.null(extLabelsFile)) {
      return ()
    }
    if (extLabelsFile == "csv") {
      read.table(file=ServerLabelsFile$datapath, sep=input$DE_run_labels_sep, header=TRUE, row.names = 1)
    } else {
      read.table(file=ServerLabelsFile$datapath, sep=input$DE_run_labels_sep, header=TRUE)
    }
    
    
  })

  observeEvent(input$DE_run_labels_file, {
    show(id = "DE_run_labels_sep")
  })



  observe({
     # DE_run_labels_file from fileInput() function
    ServerLabelsFile <- req(input$DE_run_labels_file)
    
  #   # extensions tool for format validation
    extLabelsFile <- tools::file_ext(ServerLabelsFile$datapath)
    if (is.null(input$DE_run_labels_file)) {
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
      updateRadioButtons(session, "DE_run_labels_sep", label = label, choices = choice)
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
  observeEvent(input$DE_run_PD_options_control_case, {
      options <- names(labelsData())
      updateSelectInput(session, 
        inputId = "DE_run_PD_options_column",
        "Select column to group", 
        choices = options[2:length(options)], 
        selected = NULL
      )
      show(id="DE_run")
  })

  # Group by label option 
  observeEvent(input$DE_run_PD_options_column, {
      # Update the case selection with levels of selected column 
      var <- labelsData()[[input$DE_run_PD_options_column]]
      lvl <- levels(as.factor(var))
      updateSelectInput(session, 
        inputId="DE_run_PD_options_case", 
        "Select case to analyse", 
        choices = lvl, 
        selected = NULL
      )
  })

  # countsData <- reactive({
  #   if ( is.null(input$DE_run_counts_file)) return(NULL)
  #   inFile <- input$DE_run_counts_file
  #   file <- inFile$datapath
  #   # load the file into new environment and get it from there
  #   e = new.env()
  #   name <- load(file, envir = e)
  #   data <- e[[name]]
  # })

  case_selected <- reactive({
    input$UILabelContentSelection_rows_selected
  })

  conditions_selected <- reactive({
    input$UILabelContentRemoveSelection_rows_selected
  })

  
    
  # Switch to labels tab if labels file is uploaded
  observeEvent(input$DE_run_labels_file, {
    updateTabsetPanel(session, "DE_run_VF_tabset", selected = "Labels File")
  })
  # Switch to counts tab if counts file is uploaded
  observeEvent(input$DE_run_counts_file, {
    updateTabsetPanel(session, "DE_run_VF_tabset", selected = "Counts File")
  })




  ########################### RUN DE ###########################
  observe({
    if (!is.null(input$DE_run_labels_file) && !is.null(input$DE_run_labels_file)) {
      show(id="DE_run_PD_options")
      hide(id="DE_plot_error")
    }
  })
  
  de <- reactiveValues(
    deg_output = NULL, 
  )


  observeEvent(input$DE_run, {
    labels <- labelsData()
    counts_data <- countsData()
    deg <- NULL

    # var <- labelsData()[[input$DE_run_PD_options_column]]
    if (input$DE_run_PD_options_control_case == "Choose Case by Label") {
      var <- input$DE_run_PD_options_column
      case <- input$DE_run_PD_options_case

      # Format labels$var
      labels_var <- labels[[paste0(var)]]

      #Initialise the variables of the chosen column to all be 1
      groups <- rep(1, length(labels_var))
      
      # Pick the case, relabel as 2
      groups[labels_var == case] = 2   


      filt = groups != 0 
      deg <- calc_DE(counts_data[,filt], groups[filt], input$DE_run_PD_options_method) 
      de$deg_output <- deg

    } else {
      cases <- case_selected()
      conditions <- conditions_selected()
      
      #Initalise all values to 1
      groups <- rep(0, nrow(labels))

      for (c in cases) {
        groups[c] = 2
      }
      
      for (d in conditions) {
        groups[d] = 1
      }
      # No cases have been selected
      filt = groups != 0 
      if (is.null(cases)) {
        shinyalert(title = "Invalid Input", text = "Please select cases to assess", type = "error")
      } else if (is.null(conditions)) {
        shinyalert(title = "Invalid Input", text = "Please select conditions to assess", type = "error")
      # Valid input - cases and control selected
      } else {
        deg <- calc_DE(counts_data[,filt], groups[filt], input$DE_run_PD_options_method) 
        de$deg_output <- deg
      }
      

    }
    

    
    if (!is.null(deg)) {
      show(id="DE_vol_text")

      # Volcano Plot
      output$DE_vol_plot <- renderPlot(
        {plot(deg$degs$log2_fc, -log10(deg$degs$pvals), pch=19, bty="n", xlab="log2 FC", ylab="-log10 p-vals")},
        width = 450,
        height = 450,
      )

      #MA Plot
      show(id="DE_MA_text")
      #output$DE_MA_text = renderText("MA Plot")
      output$DE_MA_plot <- renderPlot(
        {plot(log2(deg$degs$mean_cpm), deg$degs$log2_fc, pch=19, bty="n", ylab="log2 FC", xlab="Average expression (log2 CPM + 1)")},
        width = 450,
        height = 450,
      )
    show(id = "DE_run_assess")
    }
    
    }
  )

  observeEvent(input$DE_run_assess, { 
    updateTabsetPanel(session, inputId="navpage", selected="Assess DE")
    updateTabsetPanel(session, "DE_view_file_tabset", selected = "Subnetwork")
    sn$sub_nets <- NULL
    output$DE_VF_table <- renderDataTable(
        {de$deg_output$degs},
    )
  })

  observe({
      # DEFile from fileInput() function
      ServerDEFile <- req(input$GL_gene_list)
      
      # extensions tool for format validation
      extDEFile <- tools::file_ext(ServerDEFile$datapath)
      if (is.null(input$GL_gene_list)) {
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
        updateRadioButtons(session, "GL_gene_list_sep", label = label, choices = choice)
      }
    })
  
  ##########################################################################################
  #                                                                                        #
  #                                    ASSESS DE DATA                                      #
  #                                                                                        #
  ##########################################################################################

  # sub_nets
  sn <- reactiveValues(
    sub_nets_DE = NULL,
    sub_nets = NULL, 
  )

  observeEvent(input$DE_generate_subnet, {
    if (is.null(de$deg_output)) {
      shinyalert(title = "Invalid Input", text = "Please first Run DE", type = "error")
      updateTabsetPanel(session, inputId="navpage", selected="Run DE")
    } else {
      sn$sub_nets_DE <- subset_network_hdf5(de$deg_output$degs, tolower(input$GL_options_network), dir="../networks/")
      show(id = "DE_CG_options")
      hide(id = "DE_CG_error")
      show(id = "DE_GC_options")
      hide(id = "DE_GC_error")
      show(id = "DE_FO_options")
      hide(id = "DE_FO_error")
      # GSEA
      show(id = "DE_GSEA_dropdown")
      hide(id = "DE_GSEA_error")
    }
  })

  observeEvent(
    input$DE_generate_subnet, 
    {output$DE_VF_subnetwork <- renderTable(sn$sub_nets_DE)}
  )

  clust_net_DE <- reactive({

    sub_net <- sn$sub_nets_DE$sub_net
    node_degrees <- sn$sub_nets_DE$node_degrees
    medK <- as.numeric(sn$sub_nets_DE$median)

    clust_net_DE <- list() 
    
    # For DE data 
    deg_sig <- sn$sub_nets_DE$deg_sig
    fc_sig  <- sn$sub_nets_DE$fc_sig
    clust_net_DE[["down"]]  <- cluster_coexp(sub_net$down, medK = medK, flag_plot = FALSE)
    clust_net_DE[["up"]]  <- cluster_coexp( sub_net$up, medK = medK, flag_plot = FALSE)

    return(clust_net_DE)
  })


  ################################ CLUSTER GENES ########################################

  observeEvent(
    {input$DE_CG_run},
    {
      sub_net <- sn$sub_nets_DE$sub_net
      node_degrees <- sn$sub_nets_DE$node_degrees
      medK <- as.numeric(sn$sub_nets_DE$median)

      
      # upregulated network 
      show(id="DE_CG_up_network_text")
      output$DE_CG_up_network_plot <- renderPlot(
        {plot_network(sub_net$up, clust_net_DE()$up, medK)}, 
        width = 500, 
        height = 500 
      )

      # upregulated heatmap 
      show(id="DE_CG_up_heatmap_text")
      output$DE_CG_up_heatmap_plot <- renderPlot(
        {plot_coexpression_heatmap(sub_net$up, clust_net_DE()$up, flag_plot_bin = FALSE)}, 
        width = 500,
        height = 500
      )

      # upregulated binarized heatmap 
      show(id="DE_CG_up_bheatmap_text")
      output$DE_CG_up_bheatmap_plot <- renderPlot(
        {plot_coexpression_heatmap(sub_net$up, clust_net_DE()$up)}, 
        width = 500, 
        height = 500
      )
      
      # downregulated network 
      show(id="DE_CG_down_network_text")
      output$DE_CG_down_network_plot <- renderPlot(
        {plot_network(sub_net$down, clust_net_DE()$down, medK)},
        width = 500, 
        height = 500
      )

      # downregulated heatmap
      show(id="DE_CG_down_heatmap_text")
      output$DE_CG_down_heatmap_plot <- renderPlot(
        {plot_coexpression_heatmap(sub_net$down, clust_net_DE()$down, flag_plot_bin = FALSE)}, 
        width = 500, 
        height = 500 
      )

      # downregulated binarized heatmap
      show(id="DE_CG_down_bheatmap_text")
      output$DE_CG_down_bheatmap_plot <- renderPlot(
        {plot_coexpression_heatmap(sub_net$down, clust_net_DE()$down)}, 
        width = 500, 
        height = 500
      )

      # clustering genes table output
      # show(id="GL_CG_table_text")
      # output$GL_CG_table_plot <- renderDataTable(
      #   {EGAD::attr.human[match(clust_net()$genes$clusters$genes,EGAD::attr.human$name[EGAD::attr.human$chr==input$GL_chrome],input$GL_genes_no),]},
      #   # options=list(columnDefs = list(list(visible=FALSE, targets=c(0,1,2,3))))
      # )


    }
  )

  ################################ GENE CONNECTIVITY ######################################

  observeEvent(
    {input$DE_GC_run},
    {
      
      sub_net <- sn$sub_nets_DE$sub_net
      node_degrees <- sn$sub_nets_DE$node_degrees  
      medK <- as.numeric(sn$sub_nets_DE$median)
      

      # density - upreg
      show(id="DE_GC_up_density_text")
      output$DE_GC_up_density_plot <- renderPlot(
        {plot_scatter(node_degrees$up[,1]/node_degrees$n_genes_total, 
                    node_degrees$up[,2]/node_degrees$n_genes_up, 
                    xlab="Global node degree", 
                    ylab="Local node degree", flag= "density")},
        width = 500,
        height = 500
      )

      # histogram - upreg 
      show(id="DE_GC_up_hist_text")
      output$DE_GC_up_hist_plot <- renderPlot(
        {plot_scatter(node_degrees$up[,1]/node_degrees$n_genes_total, 
                    node_degrees$up[,2]/node_degrees$n_genes_up, 
                    xybreaks = input$DE_GC_options_xybreaks,
                    xlab="Global node degree", 
                    ylab="Local node degree", flag= "hist")},
        width = 500,
        height = 500
      )

      # density by subsets - upreg 
      show(id="DE_GC_up_density_subset_text")
      output$DE_GC_up_density_subset_plot <- renderPlot(
        {m <- match(clust_net_DE()$up$clusters$genes, rownames(sub_net$up))
         plot_scatter(node_degrees$up[m,1]/node_degrees$n_genes_total, 
                      node_degrees$up[m,2]/node_degrees$n_genes_up, 
                      xlab="Global node degree", 
                      ylab="Local node degree", 
                      clusters = clust_net_DE()$up$clusters, flag = "density")},
        width = 500,
        height = 500
      )

      # histogram by subsets - upreg 
      show(id="DE_GC_up_hist_subset_text")
      output$DE_GC_up_hist_subset_plot <- renderPlot(
        {m <- match(clust_net_DE()$up$clusters$genes, rownames(sub_net$up))
         plot_scatter(node_degrees$up[m,1]/node_degrees$n_genes_total, 
                      node_degrees$up[m,2]/node_degrees$n_genes_up, 
                      xybreaks = input$DE_GC_options_xybreaks,
                      xlab="Global node degree", 
                      ylab="Local node degree", 
                      clusters = clust_net_DE()$up$clusters, flag = "hist")},
        width = 500,
        height = 500
      )

      # density - downreg
      show(id="GL_GC_density_plot_downreg_text")
      output$DE_GC_down_density_plot <- renderPlot(
        {plot_scatter(node_degrees$down[,1]/node_degrees$n_genes_total, 
                    node_degrees$down[,2]/node_degrees$n_genes_down, 
                    xlab="Global node degree", 
                    ylab="Local node degree", flag= "density")},
        width = 500,
        height = 500
      )

      # histogram - downreg 
      show(id="GL_GC_hist_plot_downreg_text")
      output$DE_GC_down_hist_plot <- renderPlot(
        {plot_scatter(node_degrees$down[,1]/node_degrees$n_genes_total, 
                    node_degrees$down[,2]/node_degrees$n_genes_down, 
                    xybreaks = input$DE_GC_options_xybreaks,
                    xlab="Global node degree", 
                    ylab="Local node degree", flag= "hist")},
        width = 500,
        height = 500
      )

      # density by subsets - downreg 
      show(id="GL_GC_density_subset_plot_downreg_text")
      output$DE_GC_down_density_subset_plot <- renderPlot(
        {m <- match(clust_net_DE()$down$clusters$genes, rownames(sub_net$down))
         plot_scatter(node_degrees$down[m,1]/node_degrees$n_genes_total, 
                      node_degrees$down[m,2]/node_degrees$n_genes_down, 
                      xlab="Global node degree", 
                      ylab="Local node degree", 
                      clusters = clust_net_DE()$down$clusters, flag = "density")},
        width = 500,
        height = 500
      )

      # histogram by subsets - downreg 
      show(id="GL_GC_hist_subset_plot_downreg_text")
      output$DE_GC_down_hist_subset_plot <- renderPlot(
        {m <- match(clust_net_DE()$down$clusters$genes, rownames(sub_net$down))
         plot_scatter(node_degrees$down[m,1]/node_degrees$n_genes_total, 
                      node_degrees$down[m,2]/node_degrees$n_genes_down, 
                      xybreaks = input$DE_GC_options_xybreaks,
                      xlab="Global node degree", 
                      ylab="Local node degree", 
                      clusters = clust_net_DE()$down$clusters, flag = "hist")},
        width = 500,
        height = 500
      )

    }
  )

  ################################ FUNCTIONAL OUTLIERS ######################################
  observeEvent(
    {input$DE_FO_run},
    {

      sub_net <- sn$sub_nets_DE$sub_net
      node_degrees <- sn$sub_nets_DE$node_degrees  
      medK <- as.numeric(sn$sub_nets_DE$median)

      filt_min <- input$DE_FO_filtmin

      show(id="DE_FO_up_heatmap_text")
      output$DE_FO_up_heatmap_plot <- renderPlot(
        {plot_coexpression_heatmap(sub_net$up, clust_net_DE()$up, filt = TRUE, flag_plot_bin = FALSE)}, 
        width = 500,
        height = 500 
      )

      show(id="DE_FO_up_network_text")
      output$DE_FO_up_network_plot <- renderPlot(
        {plot_network(1-sub_net$up, clust_net_DE()$up, 1 - medK)}, 
        width = 500, 
        height = 500
      )

      show(id="DE_FO_down_heatmap_plot")
      output$DE_FO_down_heatmap_plot <- renderPlot(
        {plot_coexpression_heatmap(sub_net$down, clust_net_DE()$down, filt = TRUE, flag_plot_bin = FALSE)}, 
        width = 500, 
        height = 500 
      )

      show(id="DE_FO_down_network_text")
      output$DE_FO_down_heatmap_plot <- renderPlot(
        {plot_network(1 - sub_net$down, clust_net_DE()$down, 1 - medK)}, 
        width = 500, 
        height = 500
      )

      # # genes in module table output
      # show(id="GLGL_FO_in_text")
      # output$GL_FO_in_table <- renderDataTable(
      #   { clust_size <- plyr::count(clust_net()$genes$clusters$labels)
      #     clust_keep <- clust_size[clust_size[,2] < filt_min ,1]
      #     genes_keep <- !is.na(match(clust_net()$genes$clusters$labels, clust_keep))
      #     EGAD::attr.human[match(clust_net()$genes$clusters$genes[!genes_keep],EGAD::attr.human$name[EGAD::attr.human$chr==input$GL_chrome], input$GL_genes_no),]},
      #   # options=list(columnDefs = list(list(visible=FALSE, targets=c(0,1,2,3))))
      # )


      # # functional outliers table output
      # show(id="GL_FO_out_text")
      # output$GL_FO_out_table <- renderDataTable(
      #   { clust_size <- plyr::count(clust_net()$genes$clusters$labels)
      #     clust_keep <- clust_size[clust_size[,2] < filt_min ,1]
      #     genes_keep <- !is.na(match(clust_net()$genes$clusters$labels, clust_keep))
      #     EGAD::attr.human[match(clust_net()$genes$clusters$genes[genes_keep],EGAD::attr.human$name[EGAD::attr.human$chr==input$GL_chrome], input$GL_genes_no),]},
      #   # options=list(columnDefs = list(list(visible=FALSE, targets=c(0,1,2,3))))
      # )
      
    }
  )
  

  ##########################################################################################
  #                                                                                        #
  #                                  ASSESS GENE LIST                                      #
  #                                                                                        #
  ##########################################################################################

  # reactive converts the upload file into a reactive expression known as data
  DEData <- reactive({

    # DEFile from fileInput() function
    ServerDEFile <- input$GL_gene_list

    # extensions tool for format validation
    extDEFile <- tools::file_ext(ServerDEFile$datapath)

    # file format checking
    req(ServerDEFile)
     validate(need(extDEFile == c("csv", "tsv", "txt"), "Please upload a csv, tsv or txt file."))

    # convert data into file format
    if (is.null(extDEFile)) {
      return ()
    }

    read.table(file=ServerDEFile$datapath, sep=input$GL_gene_list_sep, header=TRUE)
  })

  # creates reactive table called DEFileContent
  output$GL_gene_listContent <- renderTable({
    if (is.null(DEData())) {
      return ()
    }
    DEData()
  })

  # handles rendering of reactive object on tb on ui
  output$UIDEContent <- renderUI({
    tableOutput("DEFileContent")
  })

  observeEvent(input$GL_options_gene_list, {
    if (input$GL_options_gene_list == "Upload Gene List") {
      updateTabsetPanel(session, "GL_VF_tabset", selected = "File")
    } else {
      updateTabsetPanel(session, "GL_VF_tabset", selected = "Subnetwork")
    }
  })

  observeEvent(input$generate_subnet, {
    gene_list <- NULL
    if (is.null(input$GL_options_gene_list)) {
      shinyalert(title = "Invalid Input", text = "Please choose a gene list method", type = "error")
    } else {
      # GENERATE GENE LIST
      if (input$GL_options_gene_list == "Generate Gene List") {

        if (str_detect(input$GL_chrome, "chr[XY]") == FALSE && str_detect(input$GL_chrome, "chr[0-9]") == FALSE) {
          shinyalert(title = "Invalid Input", text = "Please enter a Chromosome between 1 - 22, X, Y", type = "error")
          gene_list <- NULL

        } else if (str_detect(substring(input$GL_chrome,4), "[0-9]")) {
          if (strtoi(substring(input$GL_chrome,4)) < 1 || strtoi(substring(input$GL_chrome,4)) > 22) {
            shinyalert(title = "Invalid Input", text = "Please enter a Chromosome between 1 - 22, X, Y", type = "error")
            gene_list <- NULL

          } else {
            gene_list <- sample( EGAD::attr.human$name[EGAD::attr.human$chr==input$GL_chrome], input$GL_genes_no,)
            print(gene_list)
          }

        } else if (input$GL_genes_no == "" || input$GL_genes_no < 0) { 
          shinyalert(title = "Invalid Input", text = "Please enter a valid number of Genes", type = "error")
          gene_list <- NULL

        } else { 
          gene_list <- sample( EGAD::attr.human$name[EGAD::attr.human$chr==input$GL_chrome], input$GL_genes_no,)

        }
          
      } else {
        # Invalid Input - user hasn't uploaded file
        if (is.null(input$GL_gene_list)) {
          shinyalert(title = "Invalid Input", text = "Please upload a gene list file", type = "error")
        } else {
          gene_list <- read.delim(file = input$GL_gene_list$datapath, header = FALSE, sep = "\n", dec = ".")[,1]
        }

      }
      
      # Valid Input
      if (!is.null(gene_list)) { 
        sn$sub_nets <- subset_network_hdf5_gene_list(gene_list, tolower(input$GL_options_network), dir="../networks/")
        show(id = "CG_dropdown")
        hide(id = "GL_CG_error")
        show(id = "GL_GC_options")
        hide(id = "GL_GC_error")
        show(id = "GL_FO_options")
        hide(id = "GL_FO_error")
        # GSEA
        show(id = "GL_GSEA_options")
        hide(id = "GL_GSEA_error")
        # Clear data
        output$GL_CG_network_plot <- NULL
        output$GL_CG_network_text <- NULL
        output$GL_genes_no <- NULL
        output$GL_CG_table_plot <- NULL
        output$GL_GC_density_plot <- NULL
        output$GL_GC_hist_plot <- NULL
        output$GL_GC_density_subset_plot<- NULL
        output$GL_GC_hist_subset_plot <- NULL
        output$GL_FO_heatmap_plot <- NULL
        output$GL_FO_network_plot <- NULL
        output$GL_FO_in_table <- NULL
        output$GL_FO_out_table <- NULL
        # Reset Checkboxes
        updateAwesomeCheckboxGroup(
          inputId = "clusterPlotOptions_genelist",
          choices = c("Network", "Heatmap", "Binarized Heatmap"),
          status = ""
        )
        updateAwesomeCheckboxGroup(
          inputId = "GL_GC_options_plots",
          choices = c("Density", "Histogram", "Clustered Density", "Clustered Histogram"),
          status = ""
        )
        updateAwesomeCheckboxGroup(
          inputId = "GL_FO_options_plots",
          choices = c("Network", "Heatmap"),
          status = ""
        )
        updateAwesomeCheckboxGroup(
          inputId = "GL_FO_options_tables",
          choices = c("Functional Outliers", "Genes in Module"),
          status = ""
        )
      } 
    }
  })

  observeEvent(
    input$generate_subnet, 
    {output$subnetwork <- renderTable(sn$sub_nets)}
  )

  clust_net <- reactive({

    sub_net <- sn$sub_nets$sub_net
    node_degrees <- sn$sub_nets$node_degrees
    medK <- as.numeric(sn$sub_nets$median)

    clust_net <- list() 
    clust_net[["genes"]] <- cluster_coexp(sub_net$genes, medK = medK, flag_plot = FALSE)
    
    return(clust_net)

  })

  ##################### CLUSTER GENES #####################

  observeEvent(
    {input$GL_CG_run},
    {
      sub_net <- sn$sub_nets$sub_net
      node_degrees <- sn$sub_nets$node_degrees
      medK <- as.numeric(sn$sub_nets$median)

      # network output
      show(id="GL_CG_network_text")
      output$GL_CG_network_plot <- renderPlot(
        {plot_network(sub_net$genes, clust_net()$genes, medK)},
        width = 500,
        height = 500
      )


      # heatmap output
      show(id="GL_CG_heatmap_text")
      output$GL_CG_heatmap_plot <- renderPlot(
        {plot_coexpression_heatmap(sub_net$genes, clust_net()$genes, flag_plot_bin = FALSE)},
        width = 500,
        height = 500
      )


      # binarized heatmap output
      show(id="GL_CG_bheatmap_text")
      output$GL_CG_bheatmap_plot <- renderPlot(
        {plot_coexpression_heatmap(sub_net$genes, clust_net()$genes)},
        width = 500,
        height = 500
      )

      # clustering genes table output
      show(id="GL_CG_table_text")
      output$GL_CG_table_plot <- renderDataTable(
        {EGAD::attr.human[match(clust_net()$genes$clusters$genes,EGAD::attr.human$name[EGAD::attr.human$chr==input$GL_chrome],input$GL_genes_no),]},
        # options=list(columnDefs = list(list(visible=FALSE, targets=c(0,1,2,3))))
      )

    }
  )


  ##################### GENE CONNECTIVITY #####################

  observeEvent(
    {input$GL_GC_run},
    {
      
      sub_net <- sn$sub_nets$sub_net
      node_degrees <- sn$sub_nets$node_degrees  
      medK <- as.numeric(sn$sub_nets$median)
      m <- match(clust_net()$genes$clusters$genes, rownames(sub_net$genes))

      # density output
      show(id="GL_GC_density_text")
      output$GL_GC_density_plot <- renderPlot(
        {plot_scatter(node_degrees$genes[,1]/node_degrees$n_genes_total, 
                    node_degrees$genes[,2]/node_degrees$n_genes, 
                    xlab="Global node degree", 
                    ylab="Local node degree", flag= "density")},
        width = 500,
        height = 500
      )


      # histogram output
      show(id="GL_GC_hist_text")
      output$GL_GC_hist_plot <- renderPlot(
        {plot_scatter(node_degrees$genes[,1]/node_degrees$n_genes_total, 
                    node_degrees$genes[,2]/node_degrees$n_genes, 
                    xybreaks = input$GL_GC_options_xybreaks,
                    xlab="Global node degree", 
                    ylab="Local node degree", flag= "hist")},
        width = 500,
        height = 500
      )

      show(id="GL_GC_density_subset_text")
      # density output - subset by clusters
      output$GL_GC_density_subset_plot <- renderPlot(
        {plot_scatter(node_degrees$genes[m,1]/node_degrees$n_genes_total, 
                      node_degrees$genes[m,2]/node_degrees$n_genes, 
                      xlab="Global node degree", 
                      ylab="Local node degree", 
                      clusters = clust_net()$genes$clusters, flag = "density")},
        width = 500,
        height = 500
      )


      # histogram output - subset by clusters
      show(id="GL_GC_hist_subset_text")
      output$GL_GC_hist_subset_plot <- renderPlot(
        {plot_scatter(node_degrees$genes[m,1]/node_degrees$n_genes_total, 
                      node_degrees$genes[m,2]/node_degrees$n_genes, 
                      xybreaks = input$GL_GC_options_xybreaks,
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
    {input$GL_FO_run},
    {

      sub_net <- sn$sub_nets$sub_net
      node_degrees <- sn$sub_nets$node_degrees  
      medK <- as.numeric(sn$sub_nets$median)

      filt_min <- input$GL_FO_filtmin

      clust_size <- plyr::count(clust_net()$genes$clusters$labels)
      clust_keep <- clust_size[clust_size[,2] < filt_min ,1]
      genes_keep <- !is.na(match(clust_net()$genes$clusters$labels, clust_keep))
     

      # heatmap output
      show(id="GL_FO_heatmap_text")
      output$GL_FO_heatmap_plot <- renderPlot(
        {plot_coexpression_heatmap(sub_net$genes, clust_net()$genes, filt = TRUE, flag_plot_bin = FALSE)},
        width = 500,
        height = 500
      )

      # network output
      show(id="GL_FO_network_text")
      output$GL_FO_network_plot <- renderPlot(
        {plot_network(1-sub_net$genes, clust_net()$genes, 1 - medK)},
        width = 500,
        height = 500
      )

      # genes in module table output
      show(id="GLGL_FO_in_text")
      output$GL_FO_in_table <- renderDataTable(
        {EGAD::attr.human[match(clust_net()$genes$clusters$genes[!genes_keep],EGAD::attr.human$name[EGAD::attr.human$chr==input$GL_chrome], input$GL_genes_no),]},
        # options=list(columnDefs = list(list(visible=FALSE, targets=c(0,1,2,3))))
      )


      # functional outliers table output
      show(id="GL_FO_out_text")
      output$GL_FO_out_table <- renderDataTable(
        {EGAD::attr.human[match(clust_net()$genes$clusters$genes[genes_keep],EGAD::attr.human$name[EGAD::attr.human$chr==input$GL_chrome], input$GL_genes_no),]},
        # options=list(columnDefs = list(list(visible=FALSE, targets=c(0,1,2,3))))
      )
      
    }
  )



  ##################### GSEA #####################  


  observeEvent(
    {input$DE_GSEA_run},
    {
      # Standard GSEA
      if ("Standard GSEA" %in% input$GSEA_type) {
        data(go_slim)
        data(go_voc)
        
        # upregulated heatmap
        show(id="GSEA_up_heatmap_text")
        output$GSEA_up_heatmap <- renderPlot(
          {
            filt <- colSums(go_slim_entrez) < 5000 & colSums(go_slim_entrez) >= 10
            gene_list <- clust_net_DE()$up$clusters$genes[clust_net_DE()$up$order]
            go_enrich <- gene_set_enrichment(gene_list, go_slim_entrez[filt,], go_voc) 
            plot_gene_set_enrichment(go_enrich, gene_list, go_slim_entrez[filt,]) 
          },
          width = 500,
          height = 500
        )
        
        # downregulated heatmap
        show(id="GSEA_down_heatmap_text")
        output$GSEA_down_heatmap <- renderPlot(
          {
            filt <- colSums(go_slim_entrez) < 5000 & colSums(go_slim_entrez) >= 10
            gene_list <- clust_net_DE()$down$clusters$genes[clust_net_DE()$down$order]
            go_enrich <- gene_set_enrichment(gene_list, go_slim_entrez[filt,], go_voc) 
            plot_gene_set_enrichment(go_enrich, gene_list, go_slim_entrez[filt,]) 
          },
          width = 500,
          height = 500
        )
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

        # AUC graphs
        output$GSEA_auc <- renderPlot(
          {plot_gene_set_enrichment_ranked(gene_set_aucs, gene_rankings_rev, gene_list, go_slim_entrez)},
          width = 1000,
          height = 1000
        )
      }
    }
  )

  # Gene List GSEA
  observeEvent(
    {input$GL_GSEA_run},
    {
      data(go_slim)
      data(go_voc)

      # heatmap
      show(id="GL_GSEA_heatmap_text")
      output$GL_GSEA_heatmap_plot <- renderPlot(
        {
          filt <- colSums(go_slim) < 5000 & colSums(go_slim) >= 10
          gene_list <- clust_net()$genes$clusters$genes[clust_net()$genes$order]
          go_enrich <- gene_set_enrichment(gene_list, go_slim[filt,], go_voc)
          plot_gene_set_enrichment(go_enrich, gene_list, go_slim[filt,])
        },
        width = 500,
        height = 500
      )
    }
  )

}



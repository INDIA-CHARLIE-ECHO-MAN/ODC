source("./src/plot_functions.R", local = TRUE)
source("./src/cluster_coexp.R", local = TRUE)
source("./src/subset_network_hdf5.R", local = TRUE)
source("./src/calc_DE.R", local = TRUE)

# Warnings silenced for wilcox
options(warn=-1)
defaultW <- getOption("warn")


server <- function(input, output, session) {

  # Removing elements that are not functional without subnetwork
  hide(id = "GC_dropdown")
  hide(id = "GC_dropdown_DE")
  hide(id = "CG_dropdown")
  hide(id = "CG_dropdown_DE")
  hide(id = "FO_dropdown")
  hide(id = "FO_dropdown_DE")
  hide(id = "DE_GSEA_dropdown")
  hide(id = "GL_GSEA_dropdown")
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
  hide(id="CGupreg_network_text")
  hide(id="CGupreg_heatmap_text")
  hide(id="CGupreg_bheatmap_text")
  hide(id="CGdownreg_network_text")
  hide(id="CGdownreg_heatmap_text")
  hide(id="CGdownreg_bheatmap_text")
  hide(id="CG_table_text_upreg")
  hide(id="CG_table_text_downreg")


  # Gene Connectivity
  hide(id="GCdensityG_text")
  hide(id="GChistogramG_text")
  hide(id="GCdensitySubsetG_text")
  hide(id="GChistogramSubsetG_text")
  hide(id="GCdensityG_upreg_text")
  hide(id="GChistogramG_upreg_text")
  hide(id="GCdensitySubsetG_upreg_text")
  hide(id="GChistogramSubsetG_upreg_text")
  hide(id="GCdensityG_downreg_text")
  hide(id="GChistogramG_downreg_text")
  hide(id="GCdensitySubsetG_downreg_text")
  hide(id="GChistogramSubsetG_downreg_text")

  # Functional Outliers
  hide(id="FO_network_text")
  hide(id="FO_heatmap_text")
  hide(id="FOnetwork_upreg_text")
  hide(id="FOheatmap_upreg_text")
  hide(id="FOnetwork_downreg_text")
  hide(id="FOheatmap_downreg_text")
  hide(id="genes_not_keep_table_text")
  hide(id="genes_keep_table_text")
  hide(id="genes_not_keep_table_text_upreg")
  hide(id="genes_keep_table_text_upreg")
  hide(id="genes_not_keep_table_text_downreg")
  hide(id="genes_keep_table_text_downreg")

  # GSEA
  hide(id="GSEA_heatmap_text")
  hide(id="GSEA_up_heatmap_text")
  hide(id="GSEA_down_heatmap_text")


  #Download Buttons
  hide(id="volcano_download")
  hide(id= "MA_download")
  hide(id="DE_table_download")
  hide(id="CG_up_network_download")
  hide(id= "CG_up_heatmap_download")
  hide(id="CG_up_bheatmap_download")
  hide(id="CG_down_network_download")
  hide(id="CG_down_heatmap_download")
  hide(id="CG_down_bheatmap_download")
  hide(id="GC_up_density_download")
  hide(id="GC_up_histogram_download")
  hide(id="GC_up_densitySubset_download")
  hide(id="GC_up_histSubset_download")
  hide(id="GC_down_density_download")
  hide(id="GC_down_hist_download")
  hide(id="GC_down_densitySubset_download")
  hide(id="GC_down_histSubset_download")
  hide(id="FO_up_network_download")
  hide(id="FO_up_heatmap_download")
  hide(id="FO_down_network_download")
  hide(id="FO_down_heatmap_download")
  hide(id="GSEA_up_heatmap_download")
  hide(id="GSEA_down_heatmap_download")
  hide(id="GSEAauc_download")
  hide(id="CG_network_download")
  hide(id="CG_heatmap_download")
  hide(id="CG_bheatmap_download")
  hide(id="table_genelist_download")
  hide(id="GC_density_download")
  hide(id="GC_histogram_download")
  hide(id="GC_densitySubset_download")
  hide(id="GC_histogramSubset_download")
  hide(id="FO_network_download")
  hide(id="FO_heatmap_download")
  hide(id="genes_not_keep_table_download")
  hide(id="genes_keep_table_table_download")
  hide(id="GSEA_heatmap_download")
  hide(id = "DE_subnet_download")
  hide(id = "GL_subnet_download")

  hide(id="CG_table_upreg_download")
  hide(id="CG_table_downreg_download")
  hide(id= "genes_not_keep_table_upreg_download")
  hide(id="genes_keep_table_upreg_download")
  hide(id="genes_not_keep_table_downreg_download")
  hide(id="genes_keep_table_downreg_download")


  #Hide return buttons
  hide(id="SN_return")
  hide(id="CG_return")
  hide(id="GC_return")
  hide(id="FO_return")
  hide(id="GSEA_return")
  hide(id="SN_return_GL")
  hide(id="CG_return_GL")
  hide(id="GC_return_GL")
  hide(id="FO_return_GL")
  hide(id="GSEA_return_GL")

  #Download Table Separator
  observe({
    if (input$download_table_format ==  ".csv" || input$ASSESS_DE_download_table_format  ==  ".csv" || input$ASSESS_GENE_LIST_download_table_format ==  ".csv") {
      separator <<- ","
    } else if (input$download_table_format ==  ".tsv" || input$ASSESS_DE_download_table_format ==  ".tsv"  || input$ASSESS_GENE_LIST_download_table_format ==  ".tsv" ) {
      separator <<- "\t"
    } else if (input$download_table_format ==  ".txt" || input$ASSESS_DE_download_table_format ==  ".txt" || input$ASSESS_GENE_LIST_download_table_format==  ".txt" ) {
      separator <<- " "
    }
  })
  ##########################################################################################
  #                                                                                        #
  #                                    HELP                                                #
  #                                                                                        #
  ##########################################################################################
  observeEvent(input$help, {
   
    if (input$navpage == "Run DE") {
      updateTabsetPanel(session, "Header", selected = "Differential Expression Analysis")
    } else if (input$navpage == "Assess DE") {
      if (input$assessDE_navList == "View Files") {
        updateTabsetPanel(session, "Header", selected = "Generate Subnetwork")
        show(id="SN_return")
        
      } else if (input$assessDE_navList == "Cluster Genes") {
        updateTabsetPanel(session, "Header", selected = "Cluster Genes")
        show(id="CG_return")
        
      } else if (input$assessDE_navList == "Gene Connectivity") {
        updateTabsetPanel(session, "Header", selected = "Gene Connectivity")
        show(id="GC_return")
        
      } else if (input$assessDE_navList == "Functional Outliers") {
        updateTabsetPanel(session, "Header", selected = "Functional Outliers")
        show(id="FO_return")
        
      } else if (input$assessDE_navList == "Gene Set Enrichment Analysis") {
        updateTabsetPanel(session, "Header", selected = "Gene Set Enrichment Analysis")
        show(id="GSEA_return")
        
      }
      

      
  } else if (input$navpage == "Assess Gene List") {
      if (input$assessGL_navList == "View Files") {
        updateTabsetPanel(session, "Header", selected = "Generate Subnetwork")
        show(id="SN_return_GL")
        
      } else if (input$assessGL_navList == "Cluster Genes") {
        updateTabsetPanel(session, "Header", selected = "Cluster Genes")
        show(id="CG_return_GL")
        
      } else if (input$assessGL_navList == "Gene Connectivity") {
        updateTabsetPanel(session, "Header", selected = "Gene Connectivity")
        show(id="GC_return_GL")
        
      } else if (input$assessGL_navList == "Functional Outliers") {
        updateTabsetPanel(session, "Header", selected = "Functional Outliers")
        show(id="FO_return_GL")
        
      } else if (input$assessGL_navList == "Gene Set Enrichment Analysis") {
        updateTabsetPanel(session, "Header", selected = "Gene Set Enrichment Analysis")
        show(id="GSEA_return_GL")
        
      }

    # Move to the home page


    }
  updateTabsetPanel(session, "navpage", selected = "Home")
  })

  observeEvent(input$SN_return, {
          updateTabsetPanel(session, "navpage", selected = "Assess DE")
          updateTabsetPanel(session, "assessDE_navList", selected = "View Files")
          hide(id="SN_return")
  })
  observeEvent(input$DE_return, {
        updateTabsetPanel(session, "navpage", selected = "Run DE")
  })
  
  observeEvent(input$CG_return, {
          updateTabsetPanel(session, "navpage", selected = "Assess DE")
          updateTabsetPanel(session, "assessDE_navList", selected = "Cluster Genes")
          hide(id="CG_return")
  })
  observeEvent(input$GC_return, {
          updateTabsetPanel(session, "navpage", selected = "Assess DE")
          updateTabsetPanel(session, "assessDE_navList", selected = "Gene Connectivity")
          hide(id="GC_return")
  })
  observeEvent(input$FO_return, {
          updateTabsetPanel(session, "navpage", selected = "Assess DE")
          updateTabsetPanel(session, "assessDE_navList", selected = "Functional Outliers")
          hide(id="FO_return")
  })
  observeEvent(input$GSEA_return, {
          updateTabsetPanel(session, "navpage", selected = "Assess DE")
          updateTabsetPanel(session, "assessDE_navList", selected = "Gene Set Enrichment Analysis")
          hide(id="GSEA_return")
  })

  observeEvent(input$SN_return_GL, {
          updateTabsetPanel(session, "navpage", selected = "Assess Gene List")
          updateTabsetPanel(session, "assessGL_navList", selected = "View Files")
          hide(id="SN_return_GL")
  })
  observeEvent(input$CG_return_GL, {
          updateTabsetPanel(session, "navpage", selected = "Assess Gene List")
          updateTabsetPanel(session, "assessGL_navList", selected = "Cluster Genes")
          hide(id="CG_return_GL")
  })

  observeEvent(input$GC_return_GL, {
          updateTabsetPanel(session, "navpage", selected = "Assess Gene List")
          updateTabsetPanel(session, "assessGL_navList", selected = "Gene Connectivity")
          hide(id="GC_return_GL")
  })
  observeEvent(input$FO_return_GL, {
          updateTabsetPanel(session, "navpage", selected = "Assess Gene List")
          updateTabsetPanel(session, "assessGL_navList", selected = "Functional Outliers")
          hide(id="FO_return_GL")
  })

  observeEvent(input$GSEA_return_GL, {
          updateTabsetPanel(session, "navpage", selected = "Assess Gene List")
          updateTabsetPanel(session, "assessGL_navList", selected = "Gene Set Enrichment Analysis")
          hide(id="GSEA_return_GL")
  })

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
    
  #   # extensions tool for format validation
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

  conditions_selected <- reactive({
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


  observeEvent(input$run_DE, {
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
        shinyalert(title = "Missing Selection", text = "Please select cases to assess", type = "error")
      } else if (is.null(conditions)) {
        shinyalert(title = "Missing Selection", text = "Please select conditions to assess", type = "error")
      # Valid input - cases and control selected
      } else {
        deg <- calc_DE(counts_data[,filt], groups[filt], input$DE_method) 
        de$deg_output <- deg
      }
      

    }
    

    # Volcano Plot
    if (!is.null(deg)) {
      show(id="vol")
      show(id="volcano_download")
      vol <- function(){plot( deg$degs$log2_fc, -log10(deg$degs$pvals),  
              pch=19, bty="n", 
              xlab="log2 FC", ylab="-log10 p-vals" )}
      output$DEplot <- renderPlot(
              {vol()},
              width = 450,
              height = 450
      )

      #MA Plot
      show(id="MA")
      show(id= "MA_download")
      MA <- function(){plot( log2(deg$degs$mean_cpm),  deg$degs$log2_fc,  
              pch=19, bty="n", 
              ylab="log2 FC", xlab="Average expression (log2 CPM + 1)")}
      output$DEplot_average <- renderPlot(
              {MA()},
              width = 450,
              height = 450
      )
      DE_table <<- function(){de$deg_output$degs}
      output$DE_table <- renderDataTable(
          {DE_table()},
      )
      show(id="DE_table_download")

      # Download
      output$DE_table_download <- downloadHandler(
        filename = function() {
          paste("DE_data", input$download_table_format, sep="")
        },
        content = function(file) {
          write.table(DE_table(), file, row.names = TRUE, sep = separator, col.names = TRUE)
        }
      )
    show(id = "assess_run_de")

    }
    
     #------------------ DOWNLOAD ----------------------#
    #Download plots    
    output$volcano_download <- downloadHandler(
      filename = function() {
        paste("volcano_plot", input$download_format)
      },
      content = function(file) {
        if (input$download_format == ".png") {
          png(file, width=1500, height=1500)
        } else if (input$download_format == ".pdf") {
          pdf(file)
        }
        vol()
        dev.off()
      }
    )

    output$MA_download <- downloadHandler(
      filename = function() {
        paste("MA_plot", input$download_format)
      },
      content = function(file) {
        if (input$download_format == ".png") {
          png(file, width=1500, height=1500)
        } else if (input$download_format == ".pdf") {
          pdf(file)
        }
        MA()
        dev.off()
      }
    )
    }
  )

  
  ##########################################################################################
  #                                                                                        #
  #                                    ASSESS DE DATA                                      #
  #                                                                                        #
  ##########################################################################################
  # USE DE data from Run DE
  observeEvent(input$assess_run_de, { 
    updateTabsetPanel(session, inputId="navpage", selected="Assess DE")

    updateRadioButtons(session,
                  inputId = "DE_data_selection",
                  choices = c("Use DE Results", "Upload DE Data"),
                  selected = "Use DE Results"
    )
    
    output$AssessDE_table <- renderDataTable(
        {DE_table()},
    )



  })

  
  # USE DE data from Run DE
    DE <- reactive({
    ServerDEFile <- input$DE_file
    extDEFile <- tools::file_ext(ServerDEFile$datapath)
    req(ServerDEFile)
    validate(need(extDEFile == c("csv", "tsv", "txt"), "Please upload a csv, tsv or txt file."))
    if (is.null(extDEFile)) {
      return ()
    }
    if (extDEFile == "csv") {
      read.table(file=ServerDEFile$datapath, sep=input$sepDEButton, header=TRUE, row.names = 1)
    } else {
      read.table(file=ServerDEFile$datapath, sep=input$sepDEButton, header=TRUE)
    }
    
    
  })



  # Upload own DE data
  observeEvent(input$DE_file, {updateTabsetPanel(session, "subnetwork_file_tabset_DE", selected = "DE Data")}) 
  observe({
     # DE_file from fileInput() function
    ServerDEFile <- req(input$DE_file)
    
  #   # extensions tool for format validation
    extDEFile <- tools::file_ext(ServerDEFile$datapath)
    if (is.null(input$DE_file)) {
      return ()
    } else {
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
      updateRadioButtons(session, "sepDEButton", label = label, choices = choice)
      }
    })


  # handles rendering DT table of DE file

  output$UIDE_loaded_Content <- renderDataTable(
    DE(), options = list(
      pageLength = 10
    )
  )

  

#################################################

  # sub_nets
  sn <- reactiveValues(
    sub_nets_DE = NULL,
    sub_nets = NULL, 
  )

  # network upload UI output
  output$select.folder <-
    renderUI(expr = selectInput(inputId = 'folder.name',
                                label = 'Network Name',
                                choices = list.dirs(path = "../networks",
                                                    full.names = FALSE,
                                                    recursive = FALSE)))

  # network upload select network
  network_type <- reactive({
    if (!is.null(input$folder.name)) {
      return(input$folder.name)
    }
  })

  network_path <- reactive({
    attachDir <- paste0("../networks/", network_type())
    path <- paste0(attachDir, "/")
    return(path)
  }) 

  observeEvent(input$generate_subnet_DE, {
    success <- FALSE
    updateTabsetPanel(session, "assessDE_navList", selected = "View Files")
    updateTabsetPanel(session, "subnetwork_file_tabset_DE", selected = "Subnetwork")
    if (input$DE_data_selection == "Use DE Results" && is.null(de$deg_output)) {
      shinyalert(title = "Error", text = "Please first run DE", type = "error")
      updateTabsetPanel(session, inputId="navpage", selected="Run DE")
    } else if (input$DE_data_selection == "Upload DE Data" && is.null(input$DE_file)) {
      shinyalert(title = "Missing DE Data File", text = "Please upload a DE Data File", type = "error")
    #Valid file for Upload DE Data 
    } else if (input$DE_data_selection == "Upload DE Data") {
      if (input$is_occr == "Yes") {
        occr <- paste0(network_type(), ".occr")
        err_genes <- paste0(occr, ".genes.h5")
        err_median <- paste0(occr, ".med.h5")
        err_net <- paste0(occr, ".net.h5")
        genes <- paste0(network_path(), occr, ".genes.h5")
        median <- paste0(network_path(), occr, ".med.h5")
        net <- paste0(network_path(), occr, ".net.h5")
        if (!file.exists(genes)) {
          errorMess <- paste("Please ensure", err_genes, "exists in", network_type(), "folder")
          shinyalert(title = "Missing Network File", text = errorMess, type = "error")
        }
        else if (!file.exists(median)) {
          errorMess <- paste("Please ensure", err_median, "exists in", network_type(), "folder")
          shinyalert(title = "Missing Network File", text = errorMess, type = "error")
        }
        else if (!file.exists(net)) {
          errorMess <- paste("Please ensure", err_net, "exists in", network_type(), "folder")
          shinyalert(title = "Missing Network File", text = errorMess, type = "error")
        }
        else {
          sn$sub_nets_DE <- subset_network_hdf5(DE(), tolower(network_type()), dir=network_path())
          success <- TRUE
        }
      } else {
        # standard networks
        err_genes <- paste0(network_type(), ".genes.h5")
        err_median <- paste0(network_type(), ".med.h5")
        err_net <- paste0(network_type(), ".net.h5")
        genes <- paste0(network_path(), network_type(), ".genes.h5")
        median <- paste0(network_path(), network_type(), ".med.h5")
        net <- paste0(network_path(), network_type(), ".net.h5")
        if (!file.exists(genes)) {
          errorMess <- paste("Please ensure", err_genes, "exists in", network_type(), "folder")
          shinyalert(title = "Missing Network File", text = errorMess, type = "error")
        }
        else if (!file.exists(median)) {
          errorMess <- paste("Please ensure", err_median, "exists in", network_type(), "folder")
          shinyalert(title = "Missing Network File", text = errorMess, type = "error")
        }
        else if (!file.exists(net)) {
          errorMess <- paste("Please ensure", err_net, "exists in", network_type(), "folder")
          shinyalert(title = "Missing Network File", text = errorMess, type = "error")
        }
        else {
          sn$sub_nets_DE <- subset_network_hdf5(DE(), tolower(network_type()), dir=network_path(), flag_occr = FALSE)
          success <- TRUE
        }
      }
      
    # Valid for Use DE results
    } else if (input$DE_data_selection == "Use DE Results") {
      if (input$is_occr == "Yes") {
        occr <- paste0(network_type(), ".occr")
        err_genes <- paste0(occr, ".genes.h5")
        err_median <- paste0(occr, ".med.h5")
        err_net <- paste0(occr, ".net.h5")
        genes <- paste0(network_path(), occr, ".genes.h5")
        median <- paste0(network_path(), occr, ".med.h5")
        net <- paste0(network_path(), occr, ".net.h5")
        if (!file.exists(genes)) {
          errorMess <- paste("Please ensure", err_genes, "exists in", network_type(), "folder")
          shinyalert(title = "Missing Network File", text = errorMess, type = "error")
        }
        else if (!file.exists(median)) {
          errorMess <- paste("Please ensure", err_median, "exists in", network_type(), "folder")
          shinyalert(title = "Missing Network File", text = errorMess, type = "error")
        }
        else if (!file.exists(net)) {
          errorMess <- paste("Please ensure", err_net, "exists in", network_type(), "folder")
          shinyalert(title = "Missing Network File", text = errorMess, type = "error")
        }
        else {
          sn$sub_nets_DE <- subset_network_hdf5(de$deg_output$degs, tolower(network_type()), dir=network_path())
          success <- TRUE
        }
      } else {
        # standard networks
        err_genes <- paste0(network_type(), ".genes.h5")
        err_median <- paste0(network_type(), ".med.h5")
        err_net <- paste0(network_type(), ".net.h5")
        genes <- paste0(network_path(), network_type(), ".genes.h5")
        median <- paste0(network_path(), network_type(), ".med.h5")
        net <- paste0(network_path(), network_type(), ".net.h5")
        if (!file.exists(genes)) {
          errorMess <- paste("Please ensure", err_genes, "exists in", network_type(), "folder")
          shinyalert(title = "Missing Network File", text = errorMess, type = "error")
        }
        else if (!file.exists(median)) {
          errorMess <- paste("Please ensure", err_median, "exists in", network_type(), "folder")
          shinyalert(title = "Missing Network File", text = errorMess, type = "error")
        }
        else if (!file.exists(net)) {
          errorMess <- paste("Please ensure", err_net, "exists in", network_type(), "folder")
          shinyalert(title = "Missing Network File", text = errorMess, type = "error")
        }
        else {
          sn$sub_nets_DE <- subset_network_hdf5(de$deg_output$degs, tolower(network_type()), dir=network_path(), flag_occr = FALSE)
          success <- TRUE
        }


      }
    }

    if (success == TRUE) {
      show(id = "CG_dropdown_DE")
      hide(id = "CG_error_DE")
      show(id = "GC_dropdown_DE")
      hide(id = "GC_error_DE")
      show(id = "FO_dropdown_DE")
      hide(id = "FO_error_DE")
      # GSEA
      show(id = "DE_GSEA_dropdown")
      hide(id = "DE_GSEA_error")

      #Update Checkbox inputs
      updateAwesomeCheckboxGroup(inputId = "clusterPlotOptions_upreg", choices = c("Network", "Heatmap", "Binarized Heatmap"),status = "")
      updateAwesomeCheckboxGroup(inputId = "clusterPlotOptions_downreg", choices = c("Network", "Heatmap", "Binarized Heatmap"), status = "")
      updateAwesomeCheckboxGroup(inputId = "GCPlotOptions_upreg", choices = c("Density", "Histogram", "Clustered Density", "Clustered Histogram"), status = "")
      updateAwesomeCheckboxGroup(inputId = "GCPlotOptions_downreg", choices =  c("Density", "Histogram", "Clustered Density", "Clustered Histogram"), status = "")
      updateAwesomeCheckboxGroup(inputId = "FO_upreg_options", choices = c("Network", "Heatmap", "Genes in Module", "Functional Outliers"), status = "")
      updateAwesomeCheckboxGroup(inputId = "GSEA_type",choices = c("Standard GSEA", "AUCs GSEA"),status = "",)
      updateAwesomeCheckboxGroup(inputId = "GSEA_std_PlotOptions",choices = c("Upregulated P-value Heatmap", "Compare upregulated P-value to Clustered Heatmap", "Downregulated P-value Heatmap", "Compare downregulated P-value to Clustered Heatmap"),status = "")
      # Clear plots
      output$upregNetwork <- NULL
      output$upregHeatmap <- NULL
      output$upregbinHeatmap <- NULL
      output$downregNetwork <- NULL
      output$downregHeatmap <- NULL
      output$downregbinHeatmap <- NULL
      output$GCdensityGupreg <- NULL
      output$GChistogramGupreg <- NULL
      output$GCdensitySubsetGupreg <- NULL
      output$GChistogramSubsetGupreg <- NULL
      output$GCdensityGdownreg <- NULL
      output$GChistogramGdownreg <- NULL
      output$GCdensitySubsetGdownreg <- NULL
      output$GChistogramSubsetGdownreg <- NULL
      output$FOnetwork_upreg <- NULL
      output$FOheatmap_upreg <- NULL
      output$FOnetwork_downreg <- NULL
      output$FOheatmap_downreg <- NULL
      output$GSEA_up_heatmap <- NULL
      output$GSEA_down_heatmap <- NULL
      output$GSEA_auc <- NULL
      DE_subnet_table <- function(){sn$sub_nets_DE}
      output$subnetwork_DE <- renderTable(sn$sub_nets_DE)
      show(id = "DE_subnet_download")
            
      output$DE_subnet_download <- downloadHandler(
      filename = function() {
        paste("DE_subnetwork", input$ASSESS_DE_download_table_format, sep="")
      },
      content = function(file) {
        write.table(DE_subnet_table(), file, row.names = FALSE, sep = separator, col.names = TRUE)
      }
    )

    }


      
    

  })
  
    



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
    {input$runCGDE},
    {
      sub_net <- sn$sub_nets_DE$sub_net
      node_degrees <- sn$sub_nets_DE$node_degrees
      medK <- as.numeric(sn$sub_nets_DE$median)

      
      # upregulated network 
      show(id="CGupreg_network_text")
      show(id="CG_up_network_download")
      CG_up_network <- function(){plot_network(sub_net$up, clust_net_DE()$up, medK)}
      output$upregNetwork <- renderPlot(
        {CG_up_network()}, 
        width = 500, 
        height = 500 
      )

      # upregulated heatmap 
      show(id="CGupreg_heatmap_text")
      show(id= "CG_up_heatmap_download")
      CG_up_heatmap <- function(){plot_coexpression_heatmap(sub_net$up, clust_net_DE()$up, flag_plot_bin = FALSE)}
      output$upregHeatmap <- renderPlot(
        {CG_up_heatmap()}, 
        width = 500,
        height = 500
      )

      # upregulated binarized heatmap 
      show(id="CGupreg_bheatmap_text")
      show(id="CG_up_bheatmap_download")
      CG_up_bheatmap <- function(){plot_coexpression_heatmap(sub_net$up, clust_net_DE()$up)}
      output$upregbinHeatmap <- renderPlot(
        {CG_up_bheatmap()}, 
        width = 500, 
        height = 500
      )
      
      # downregulated network 
      show(id="CGdownreg_network_text")
      show(id="CG_down_network_download")
      CG_down_network <- function(){plot_network(sub_net$down, clust_net_DE()$down, medK)}
      output$downregNetwork <- renderPlot(
        {CG_down_network()},
        width = 500, 
        height = 500
      )

      # downregulated heatmap
      show(id="CGdownreg_heatmap_text")
      show(id="CG_down_heatmap_download")
      CG_down_heatmap <- function(){plot_coexpression_heatmap(sub_net$down, clust_net_DE()$down, flag_plot_bin = FALSE)}
      output$downregHeatmap <- renderPlot(
        {CG_down_heatmap()}, 
        width = 500, 
        height = 500 
      )

      # downregulated binarized heatmap
      show(id="CGdownreg_bheatmap_text")
      show(id="CG_down_bheatmap_download")
      CG_down_bheatmap <- function(){plot_coexpression_heatmap(sub_net$down, clust_net_DE()$down)}
      output$downregbinHeatmap <- renderPlot(
        {CG_down_bheatmap()}, 
        width = 500, 
        height = 500
      )

      # clustering genes table output
      show(id="CG_table_text_upreg")
      show(id="CG_table_upreg_download")
      CG_table_upreg <- function(){EGAD::attr.human[match(clust_net_DE()$up$clusters$genes, EGAD::attr.human$entrezID),]}
      output$CG_table_upreg <- renderDataTable(
        {CG_table_upreg()},
        # options=list(columnDefs = list(list(visible=FALSE, targets=c(0,1,2,3))))
      )

      show(id="CG_table_text_downreg")
      show(id="CG_table_downreg_download")
      CG_table_downreg <- function(){EGAD::attr.human[match(clust_net_DE()$down$clusters$genes, EGAD::attr.human$entrezID),]}
      output$CG_table_downreg <- renderDataTable(
        {CG_table_downreg()},
        # options=list(columnDefs = list(list(visible=FALSE, targets=c(0,1,2,3))))
      )
      
      #------------------ DOWNLOAD ----------------------#
      #Download plots    
      output$CG_up_network_download <- downloadHandler(
        filename = function() {
          paste("plot_coexpression_network_up", input$ASSESS_DE_download_format)
        },
        content = function(file) {
          if (input$ASSESS_DE_download_format == ".png") {
            png(file, width=1500, height=1500)
          } else if (input$ASSESS_DE_download_format == ".pdf") {
            pdf(file)
          }
          CG_up_network()
          dev.off()
        }
      )

      output$CG_up_heatmap_download <- downloadHandler(
        filename = function() {
          paste("plot_coexpression_heatmap_up", input$ASSESS_DE_download_format)
        },
        content = function(file) {
          if (input$ASSESS_DE_download_format == ".png") {
            png(file, width=1500, height=1500)
          } else if (input$ASSESS_DE_download_format == ".pdf") {
            pdf(file)
          }
          CG_up_heatmap()
          dev.off()
        }
      )

      output$CG_up_bheatmap_download <- downloadHandler(
        filename = function() {
          paste("plot_coexpression_binarized_heatmap_up", input$ASSESS_DE_download_format)
        },
        content = function(file) {
          if (input$ASSESS_DE_download_format == ".png") {
            png(file, width=1500, height=1500)
          } else if (input$ASSESS_DE_download_format == ".pdf") {
            pdf(file)
          }
          CG_up_bheatmap()
          dev.off()
        }
      )

      output$CG_down_network_download <- downloadHandler(
        filename = function() {
          paste("plot_coexpression_network_down", input$ASSESS_DE_download_format)
        },
        content = function(file) {
          if (input$ASSESS_DE_download_format == ".png") {
            png(file, width=1500, height=1500)
          } else if (input$ASSESS_DE_download_format == ".pdf") {
            pdf(file)
          }
          CG_down_network()
          dev.off()
        }
      )

      output$CG_down_heatmap_download <- downloadHandler(
        filename = function() {
          paste("plot_coexpression_heatmap_down", input$ASSESS_DE_download_format)
        },
        content = function(file) {
          if (input$ASSESS_DE_download_format == ".png") {
            png(file, width=1500, height=1500)
          } else if (input$ASSESS_DE_download_format == ".pdf") {
            pdf(file)
          }
          CG_down_heatmap()
          dev.off()
        }
      )

      output$CG_down_bheatmap_download <- downloadHandler(
        filename = function() {
          paste("plot_coexpression_binarized_heatmap_down", input$ASSESS_DE_download_format)
        },
        content = function(file) {
          if (input$ASSESS_DE_download_format == ".png") {
            png(file, width=1500, height=1500)
          } else if (input$ASSESS_DE_download_format == ".pdf") {
            pdf(file)
          }
          CG_down_bheatmap()
          dev.off()
        }
      )

    output$CG_table_upreg_download <- downloadHandler(
      filename = function() {
        paste("clustered_up", input$ASSESS_DE_download_table_format, sep="")
      },
      content = function(file) {
        write.table(CG_table_upreg(), file, row.names = FALSE, sep = separator, col.names = TRUE)
      }
    )

    output$CG_table_downreg_download <- downloadHandler(
      filename = function() {
        paste("clustered_down", input$ASSESS_DE_download_table_format, sep="")
      },
      content = function(file) {
        write.table(CG_table_downreg(), file, row.names = FALSE, sep = separator, col.names = TRUE)
      }
    )

    }
  )

  ################################ GENE CONNECTIVITY ######################################

  observeEvent(
    {input$runGCDE},
    {
      
      sub_net <- sn$sub_nets_DE$sub_net
      node_degrees <- sn$sub_nets_DE$node_degrees  
      medK <- as.numeric(sn$sub_nets_DE$median)
      

      # density - upreg
      show(id="GCdensityG_upreg_text")
      show(id="GC_up_density_download")
      GC_up_density <- function(){plot_scatter(node_degrees$up[,1]/node_degrees$n_genes_total, 
                    node_degrees$up[,2]/node_degrees$n_genes_up, 
                    xlab="Global node degree", 
                    ylab="Local node degree", flag= "density")}
      output$GCdensityGupreg <- renderPlot(
        {GC_up_density()},
        width = 500,
        height = 500
      )

      # histogram - upreg 
      show(id="GChistogramG_upreg_text")
      show(id="GC_up_histogram_download")
      GC_up_histogram <- function(){plot_scatter(node_degrees$up[,1]/node_degrees$n_genes_total, 
                    node_degrees$up[,2]/node_degrees$n_genes_up, 
                    xybreaks = input$xybreaks_DE,
                    xlab="Global node degree", 
                    ylab="Local node degree", flag= "hist")}
      output$GChistogramGupreg <- renderPlot(
        {GC_up_histogram()},
        width = 500,
        height = 500
      )

      # density by subsets - upreg 
      show(id="GCdensitySubsetG_upreg_text")
      show(id="GC_up_densitySubset_download")
      GC_up_densitySubset <- function(){m <- match(clust_net_DE()$up$clusters$genes, rownames(sub_net$up))
         plot_scatter(node_degrees$up[m,1]/node_degrees$n_genes_total, 
                      node_degrees$up[m,2]/node_degrees$n_genes_up, 
                      xlab="Global node degree", 
                      ylab="Local node degree", 
                      clusters = clust_net_DE()$up$clusters, flag = "density")}
      output$GCdensitySubsetGupreg <- renderPlot(
        {GC_up_densitySubset()},
        width = 500,
        height = 500
      )

      # histogram by subsets - upreg 
      show(id="GChistogramSubsetG_upreg_text")
      show(id="GC_up_histSubset_download")
      GC_up_histSubset <- function(){m <- match(clust_net_DE()$up$clusters$genes, rownames(sub_net$up))
         plot_scatter(node_degrees$up[m,1]/node_degrees$n_genes_total, 
                      node_degrees$up[m,2]/node_degrees$n_genes_up, 
                      xybreaks = input$xybreaks_DE,
                      xlab="Global node degree", 
                      ylab="Local node degree", 
                      clusters = clust_net_DE()$up$clusters, flag = "hist")}
      output$GChistogramSubsetGupreg <- renderPlot(
        {GC_up_histSubset()},
        width = 500,
        height = 500
      )

      # density - downreg
      show(id="GCdensityG_downreg_text")
      show(id="GC_down_density_download")
      GC_down_density <- function(){plot_scatter(node_degrees$down[,1]/node_degrees$n_genes_total, 
                    node_degrees$down[,2]/node_degrees$n_genes_down, 
                    xlab="Global node degree", 
                    ylab="Local node degree", flag= "density")}
      output$GCdensityGdownreg <- renderPlot(
        {GC_down_density()},
        width = 500,
        height = 500
      )

      # histogram - downreg 
      show(id="GChistogramG_downreg_text")
      show(id="GC_down_hist_download")
      GC_down_hist <- function(){plot_scatter(node_degrees$down[,1]/node_degrees$n_genes_total, 
                    node_degrees$down[,2]/node_degrees$n_genes_down, 
                    xybreaks = input$xybreaks_DE,
                    xlab="Global node degree", 
                    ylab="Local node degree", flag= "hist")}
      output$GChistogramGdownreg <- renderPlot(
        {GC_down_hist()},
        width = 500,
        height = 500
      )

      # density by subsets - downreg 
      show(id="GCdensitySubsetG_downreg_text")
      show(id="GC_down_densitySubset_download")
      GC_down_densitySubset <- function(){m <- match(clust_net_DE()$down$clusters$genes, rownames(sub_net$down))
         plot_scatter(node_degrees$down[m,1]/node_degrees$n_genes_total, 
                      node_degrees$down[m,2]/node_degrees$n_genes_down, 
                      xlab="Global node degree", 
                      ylab="Local node degree", 
                      clusters = clust_net_DE()$down$clusters, flag = "density")}
      output$GCdensitySubsetGdownreg <- renderPlot(
      {GC_down_densitySubset()},
        width = 500,
        height = 500
      )

      # histogram by subsets - downreg 
      show(id="GChistogramSubsetG_downreg_text")
      show(id="GC_down_histSubset_download")
      GC_down_histSubset <- function(){m <- match(clust_net_DE()$down$clusters$genes, rownames(sub_net$down))
         plot_scatter(node_degrees$down[m,1]/node_degrees$n_genes_total, 
                      node_degrees$down[m,2]/node_degrees$n_genes_down, 
                      xybreaks = input$xybreaks_DE,
                      xlab="Global node degree", 
                      ylab="Local node degree", 
                      clusters = clust_net_DE()$down$clusters, flag = "hist")}
      output$GChistogramSubsetGdownreg <- renderPlot(
        {GC_down_histSubset()},
        width = 500,
        height = 500
      )

      #------------------ DOWNLOAD ----------------------#
      #Download plots    
      output$GC_up_density_download <- downloadHandler(
        filename = function() {
          paste("plot_scatter_density_up", input$ASSESS_DE_download_format)
        },
        content = function(file) {
          if (input$ASSESS_DE_download_format == ".png") {
            png(file, width=1500, height=1500)
          } else if (input$ASSESS_DE_download_format == ".pdf") {
            pdf(file)
          }
          GC_up_density()
          dev.off()
        }
      )

      output$GC_up_histogram_download <- downloadHandler(
        filename = function() {
          paste("plot_scatter_hist_up", input$ASSESS_DE_download_format)
        },
        content = function(file) {
          if (input$ASSESS_DE_download_format == ".png") {
            png(file, width=1500, height=1500)
          } else if (input$ASSESS_DE_download_format == ".pdf") {
            pdf(file)
          }
          GC_up_histogram()
          dev.off()
        }
      )

      output$GC_up_densitySubset_download <- downloadHandler(
        filename = function() {
          paste("plot_scatter_density_up_colored", input$ASSESS_DE_download_format)
        },
        content = function(file) {
          if (input$ASSESS_DE_download_format == ".png") {
            png(file, width=1500, height=1500)
          } else if (input$ASSESS_DE_download_format == ".pdf") {
            pdf(file)
          }
          GC_up_densitySubset()
          dev.off()
        }
      )

      output$GC_up_histSubset_download <- downloadHandler(
        filename = function() {
          paste("plot_scatter_hist_up_colored", input$ASSESS_DE_download_format)
        },
        content = function(file) {
          if (input$ASSESS_DE_download_format == ".png") {
            png(file, width=1500, height=1500)
          } else if (input$ASSESS_DE_download_format == ".pdf") {
            pdf(file)
          }
          GC_up_histSubset()
          dev.off()
        }
      )

      output$GC_down_density_download <- downloadHandler(
        filename = function() {
          paste("plot_scatter_density_down", input$ASSESS_DE_download_format)
        },
        content = function(file) {
          if (input$ASSESS_DE_download_format == ".png") {
            png(file, width=1500, height=1500)
          } else if (input$ASSESS_DE_download_format == ".pdf") {
            pdf(file)
          }
          GC_down_density()
          dev.off()
        }
      )

      output$GC_down_hist_download <- downloadHandler(
        filename = function() {
          paste("plot_scatter_hist_down", input$ASSESS_DE_download_format)
        },
        content = function(file) {
          if (input$ASSESS_DE_download_format == ".png") {
            png(file, width=1500, height=1500)
          } else if (input$ASSESS_DE_download_format == ".pdf") {
            pdf(file)
          }
          GC_down_hist()
          dev.off()
        }
      )

      output$GC_down_densitySubset_download <- downloadHandler(
        filename = function() {
          paste("plot_scatter_density_down_colored", input$ASSESS_DE_download_format)
        },
        content = function(file) {
          if (input$ASSESS_DE_download_format == ".png") {
            png(file, width=1500, height=1500)
          } else if (input$ASSESS_DE_download_format == ".pdf") {
            pdf(file)
          }
          GC_down_densitySubset()
          dev.off()
        }
      )

      output$GC_down_histSubset_download <- downloadHandler(
        filename = function() {
          paste("plot_scatter_hist_down_colored", input$ASSESS_DE_download_format)
        },
        content = function(file) {
          if (input$ASSESS_DE_download_format == ".png") {
            png(file, width=1500, height=1500)
          } else if (input$ASSESS_DE_download_format == ".pdf") {
            pdf(file)
          }
          GC_down_histSubset()
          dev.off()
        }
      )



    }
  )

  ################################ FUNCTIONAL OUTLIERS ######################################
  observeEvent(
    {input$runFODE},
    {

      sub_net <- sn$sub_nets_DE$sub_net
      node_degrees <- sn$sub_nets_DE$node_degrees  
      medK <- as.numeric(sn$sub_nets_DE$median)

      filt_min <- input$filtmin_DE

      show(id="FOnetwork_upreg_text")
      show(id="FO_up_network_download")
      FO_up_network <- function(){plot_network(1-sub_net$up, clust_net_DE()$up, 1 - medK)}
      output$FOnetwork_upreg <- renderPlot(
        {FO_up_network()}, 
        width = 500, 
        height = 500
      )

      show(id="FOheatmap_upreg_text")
      show(id="FO_up_heatmap_download")
      FO_up_heatmap <- function(){plot_coexpression_heatmap(sub_net$up, clust_net_DE()$up, filt = TRUE, flag_plot_bin = FALSE)}
      output$FOheatmap_upreg <- renderPlot(
        {FO_up_heatmap()}, 
        width = 500,
        height = 500 
      )

      show(id="FOnetwork_downreg_text")
      show(id="FO_down_network_download")
      FO_down_network <- function(){plot_network(1 - sub_net$down, clust_net_DE()$down, 1 - medK)}
      output$FOnetwork_downreg <- renderPlot(
        {FO_down_network()}, 
        width = 500, 
        height = 500
      )

      show(id="FOheatmap_downreg_text")
      show(id="FO_down_heatmap_download")
      FO_down_heatmap <- function(){plot_coexpression_heatmap(sub_net$down, clust_net_DE()$down, filt = TRUE, flag_plot_bin = FALSE)}
      output$FOheatmap_downreg <- renderPlot(
        {FO_down_heatmap()}, 
        width = 500, 
        height = 500 
      )


      # # genes in module table output
      show(id="genes_not_keep_table_text_upreg")
      show(id= "genes_not_keep_table_upreg_download")
      genes_not_keep_table_upreg <- function(){
          clust_size <- plyr::count(clust_net_DE()$up$clusters$labels)
          clust_keep <- clust_size[clust_size[,2] < filt_min ,1]
          genes_keep <- !is.na(match(clust_net_DE()$up$clusters$labels, clust_keep))
          EGAD::attr.human[match(clust_net_DE()$up$clusters$genes[!genes_keep], EGAD::attr.human$entrezID),]
          }

      output$genes_not_keep_table_upreg <- renderDataTable(
        {genes_not_keep_table_upreg()},
        # options=list(columnDefs = list(list(visible=FALSE, targets=c(0,1,2,3))))
      )

      # # functional outliers table output
      show(id="genes_keep_table_text_upreg")
      show(id="genes_keep_table_upreg_download")
      genes_keep_table_upreg <- function(){
        clust_size <- plyr::count(clust_net_DE()$up$clusters$labels)
        clust_keep <- clust_size[clust_size[,2] < filt_min ,1]
        genes_keep <- !is.na(match(clust_net_DE()$up$clusters$labels, clust_keep))
        EGAD::attr.human[match(clust_net_DE()$up$clusters$genes[genes_keep], EGAD::attr.human$entrezID),]
      }
      output$genes_keep_table_upreg <- renderDataTable(
        {genes_keep_table_upreg()},
        # options=list(columnDefs = list(list(visible=FALSE, targets=c(0,1,2,3))))
      )

      show(id="genes_not_keep_table_text_downreg")
      show(id="genes_not_keep_table_downreg_download")
      genes_not_keep_table_downreg <- function(){
        clust_size <- plyr::count(clust_net_DE()$down$clusters$labels)
        clust_keep <- clust_size[clust_size[,2] < filt_min ,1]
        genes_keep <- !is.na(match(clust_net_DE()$down$clusters$labels, clust_keep))
        EGAD::attr.human[match(clust_net_DE()$down$clusters$genes[!genes_keep], EGAD::attr.human$entrezID),]
      }
      output$genes_not_keep_table_downreg <- renderDataTable(
        {genes_not_keep_table_downreg()},
        # options=list(columnDefs = list(list(visible=FALSE, targets=c(0,1,2,3))))
      )

      show(id="genes_keep_table_text_downreg")
      show(id="genes_keep_table_downreg_download")
      genes_keep_table_downreg <- function(){
        clust_size <- plyr::count(clust_net_DE()$down$clusters$labels)
        clust_keep <- clust_size[clust_size[,2] < filt_min ,1]
        genes_keep <- !is.na(match(clust_net_DE()$down$clusters$labels, clust_keep))
        EGAD::attr.human[match(clust_net_DE()$down$clusters$genes[genes_keep], EGAD::attr.human$entrezID),]
      }
      output$genes_keep_table_downreg <- renderDataTable(
        {genes_keep_table_downreg()},
        # options=list(columnDefs = list(list(visible=FALSE, targets=c(0,1,2,3))))
      )

      #------------------ DOWNLOAD ----------------------#
      #Download plots    
      output$FO_up_network_download <- downloadHandler(
        filename = function() {
          paste("plot_coexpression_network_up_filtered", input$ASSESS_DE_download_format)
        },
        content = function(file) {
          if (input$ASSESS_DE_download_format == ".png") {
            png(file, width=1500, height=1500)
          } else if (input$ASSESS_DE_download_format == ".pdf") {
            pdf(file)
          }
          FO_up_network()
          dev.off()
        }
      ) 
      output$FO_up_heatmap_download <- downloadHandler(
        filename = function() {
          paste("plot_coexpression_heatmap_up_filtered", input$ASSESS_DE_download_format)
        },
        content = function(file) {
          if (input$ASSESS_DE_download_format == ".png") {
            png(file, width=1500, height=1500)
          } else if (input$ASSESS_DE_download_format == ".pdf") {
            pdf(file)
          }
          FO_up_heatmap()
          dev.off()
        }
      ) 
      output$FO_down_network_download <- downloadHandler(
        filename = function() {
          paste("plot_coexpression_network_down_filtered", input$ASSESS_DE_download_format)
        },
        content = function(file) {
          if (input$ASSESS_DE_download_format == ".png") {
            png(file, width=1500, height=1500)
          } else if (input$ASSESS_DE_download_format == ".pdf") {
            pdf(file)
          }
          FO_down_network()
          dev.off()
        }
      ) 
      output$FO_down_heatmap_download <- downloadHandler(
        filename = function() {
          paste("plot_coexpression_heatmap_down_filtered", input$ASSESS_DE_download_format)
        },
        content = function(file) {
          if (input$ASSESS_DE_download_format == ".png") {
            png(file, width=1500, height=1500)
          } else if (input$ASSESS_DE_download_format == ".pdf") {
            pdf(file)
          }
          FO_down_heatmap()
          dev.off()
        }
      )

      # Table Download
      output$genes_not_keep_table_upreg_download <- downloadHandler(
        filename = function() {
          paste("genes_not_keep_table_up", input$ASSESS_DE_download_table_format, sep="")
        },
        content = function(file) {
          write.table(genes_not_keep_table_upreg(), file, row.names = FALSE, sep = separator, col.names = TRUE)
        }
      ) 

      output$genes_keep_table_upreg_download <- downloadHandler(
        filename = function() {
          paste("genes_keep_table_up", input$ASSESS_DE_download_table_format, sep="")
        },
        content = function(file) {
          write.table(genes_keep_table_upreg(), file, row.names = FALSE, sep = separator, col.names = TRUE)
        }
      ) 

      output$genes_not_keep_table_downreg_download <- downloadHandler(
        filename = function() {
          paste("genes_not_keep_table_down", input$ASSESS_DE_download_table_format, sep="")
        },
        content = function(file) {
          write.table(genes_not_keep_table_downreg(), file, row.names = FALSE, sep = separator, col.names = TRUE)
        }
      ) 

      output$genes_keep_table_downreg_download <- downloadHandler(
        filename = function() {
          paste("genes_keep_table_down", input$ASSESS_DE_download_table_format, sep="")
        },
        content = function(file) {
          write.table(genes_keep_table_downreg(), file, row.names = FALSE, sep = separator, col.names = TRUE)
        }
      ) 

    }
  )

  ################################ GENE SET ENRICHMENT ANALYSIS ######################################  
  observeEvent(
    {input$DE_GSEA_run},
    {
      # Standard GSEA
      if ("Standard GSEA" %in% input$GSEA_type) {
        data(go_slim)
        data(go_voc)
        
        # upregulated heatmap
        show(id="GSEA_up_heatmap_text")
        show(id="GSEA_up_heatmap_download")
        GSEA_up_heatmap <- function(){
          filt <- colSums(go_slim_entrez) < 5000 & colSums(go_slim_entrez) >= 10
          gene_list <- clust_net_DE()$up$clusters$genes[clust_net_DE()$up$order]
          go_enrich <- gene_set_enrichment(gene_list, go_slim_entrez[filt,], go_voc) 
          plot_gene_set_enrichment(go_enrich, gene_list, go_slim_entrez[filt,]) 
        }
        output$GSEA_up_heatmap <- renderPlot(
          {GSEA_up_heatmap()},
          width = 500,
          height = 500
        )

        # upregulated compare
        # DE
        output$DE_UP_GSEA_heatmap_compare <- renderPlot{
          {GSEA_up_heatmap()},
          width = 500, 
          height = 500
        }
        
        # GSEA
        

        # downregulated heatmap
        show(id="GSEA_down_heatmap_text")
        show(id="GSEA_down_heatmap_download")
        GSEA_down_heatmap <- function(){filt <- colSums(go_slim_entrez) < 5000 & colSums(go_slim_entrez) >= 10
            gene_list <- clust_net_DE()$down$clusters$genes[clust_net_DE()$down$order]
            go_enrich <- gene_set_enrichment(gene_list, go_slim_entrez[filt,], go_voc) 
            plot_gene_set_enrichment(go_enrich, gene_list, go_slim_entrez[filt,])}
        output$GSEA_down_heatmap <- renderPlot(
          {GSEA_down_heatmap()},
          width = 500,
          height = 500
        )

        # downregulated compare
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
        show(id="GSEAauc_download")
        GSEAauc <- function(){plot_gene_set_enrichment_ranked(gene_set_aucs, gene_rankings_rev, gene_list, go_slim_entrez)}
        output$GSEA_auc <- renderPlot(
          {GSEAauc()},
          width = 500,
          height = 500
        )
      }


    #------------------ DOWNLOAD ----------------------#
    #Download plots    
    output$GSEA_up_heatmap_download <- downloadHandler(
      filename = function() {
        paste("GSEA_up", input$ASSESS_DE_download_format)
      },
      content = function(file) {
        if (input$ASSESS_DE_download_format == ".png") {
          png(file, width=1500, height=1500)
        } else if (input$ASSESS_DE_download_format == ".pdf") {
          pdf(file)
        }
        GSEA_up_heatmap()
        dev.off()
      }
    )
    output$GSEA_down_heatmap_download <- downloadHandler(
      filename = function() {
        paste("GSEA_down", input$ASSESS_DE_download_format)
      },
      content = function(file) {
        if (input$ASSESS_DE_download_format == ".png") {
          png(file, width=1500, height=1500)
        } else if (input$ASSESS_DE_download_format == ".pdf") {
          pdf(file)
        }
        GSEA_down_heatmap()
        dev.off()
      }
    )
    output$GSEAauc_download <- downloadHandler(
      filename = function() {
        paste("GSEA_auc", input$ASSESS_DE_download_format)
      },
      content = function(file) {
        if (input$ASSESS_DE_download_format == ".png") {
          png(file, width=1500, height=1500)
        } else if (input$ASSESS_DE_download_format == ".pdf") {
          pdf(file)
        }
        GSEAauc()
        dev.off()
      }
    )

    }
  )
  

  ##########################################################################################
  #                                                                                        #
  #                                  ASSESS GENE LIST                                      #
  #                                                                                        #
  ##########################################################################################

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

  observeEvent(input$gene_list_selection, {
    if (input$gene_list_selection == "Upload Gene List") {
      updateTabsetPanel(session, "subnetwork_file_tabset", selected = "File")
    } else {
      updateTabsetPanel(session, "subnetwork_file_tabset", selected = "Subnetwork")
    }
  })

  # network upload UI output
  output$select.folder_gene_list <-
    renderUI(expr = selectInput(inputId = 'folder.name_gene_list',
                                label = 'Network Name',
                                choices = list.dirs(path = "../networks",
                                                    full.names = FALSE,
                                                    recursive = FALSE)))

  # network upload select network
  network_type_gene_list <- reactive({
    if (!is.null(input$folder.name_gene_list)) {
      return(input$folder.name_gene_list)
    }
  })

  network_path_gene_list <- reactive({
    attachDir <- paste0("../networks/", network_type_gene_list())
    path <- paste0(attachDir, "/")
    return(path)
  }) 

  observeEvent(input$generate_subnet, {
    success <- FALSE
    gene_list <- NULL
    # updateTabsetPanel()
    updateTabsetPanel(session, "assessGL_navList", selected = "View Files")
    updateTabsetPanel(session, "subnetwork_file_tabset", selected = "Subnetwork")
    if (is.null(input$gene_list_selection)) {
      shinyalert(title = "Missing Selection", text = "Please choose a gene list method", type = "error")
    } else {
      # GENERATE GENE LIST
      if (input$gene_list_selection == "Generate Gene List") {

        if (str_detect(input$chooseChrome, "chr[XY]") == FALSE && str_detect(input$chooseChrome, "chr[0-9]") == FALSE) {
          shinyalert(title = "Invalid Input", text = "Please enter a Chromosome between 1 - 22, X, Y", type = "error")
          gene_list <- NULL

        } else if (str_detect(substring(input$chooseChrome,4), "[0-9]")) {
          if (strtoi(substring(input$chooseChrome,4)) < 1 || strtoi(substring(input$chooseChrome,4)) > 22) {
            shinyalert(title = "Invalid Input", text = "Please enter a Chromosome between 1 - 22, X, Y", type = "error")
            gene_list <- NULL

          } else {
            gene_list <- sample( EGAD::attr.human$name[EGAD::attr.human$chr==input$chooseChrome], input$chooseGeneNo,)
          }

        } else if (input$chooseGeneNo == "" || input$chooseGeneNo < 0) { 
          shinyalert(title = "Invalid Input", text = "Please enter a valid number of Genes", type = "error")
          gene_list <- NULL

        } else { 
          gene_list <- sample( EGAD::attr.human$name[EGAD::attr.human$chr==input$chooseChrome], input$chooseGeneNo,)

        }
          
      } else {
        # Invalid Input - user hasn't uploaded file
        if (is.null(input$DEFile)) {
          shinyalert(title = "Invalid Input", text = "Please upload a gene list file", type = "error")
        } else if (is.null(input$GL_gene_list_type)) {
          shinyalert(title = "Missing Selection", text = "Please select gene list type", type = "error")
        }
        else {
          gene_list <- read.delim(file = input$DEFile$datapath, header = FALSE, sep = "\n", dec = ".")[,1]
        }

      }
      #Valid Input
      if (!is.null(gene_list)) { 

        # occr network
        if (input$is_occr_gene_list == "Yes") {
          occr <- paste0(network_type_gene_list(), ".occr")
          err_genes <- paste0(occr, ".genes.h5")
          err_median <- paste0(occr, ".med.h5")
          err_net <- paste0(occr, ".net.h5")
          genes <- paste0(network_path_gene_list(), occr, ".genes.h5")
          median <- paste0(network_path_gene_list(), occr, ".med.h5")
          net <- paste0(network_path_gene_list(), occr, ".net.h5")
          if (!file.exists(genes)) {
            errorMess <- paste("Please ensure", err_genes, "exists in", network_type_gene_list(), "folder")
            shinyalert(title = "Missing Network File", text = errorMess, type = "error")
          }
          else if (!file.exists(median)) {
            errorMess <- paste("Please ensure", err_median, "exists in", network_type_gene_list(), "folder")
            shinyalert(title = "Missing Network File", text = errorMess, type = "error")
          }
          else if (!file.exists(net)) {
            errorMess <- paste("Please ensure", err_net, "exists in", network_type_gene_list(), "folder")
            shinyalert(title = "Missing Network File", text = errorMess, type = "error")
          }
          else {
            sn$sub_nets <- subset_network_hdf5_gene_list(gene_list, tolower(network_type_gene_list()), dir=network_path_gene_list())
            success <- TRUE
          }
        } else {
          # standard network
          err_genes <- paste0(network_type_gene_list(), ".genes.h5")
          err_median <- paste0(network_type_gene_list(), ".med.h5")
          err_net <- paste0(network_type_gene_list(), ".net.h5")
          genes <- paste0(network_path_gene_list(), network_type_gene_list(), ".genes.h5")
          median <- paste0(network_path_gene_list(), network_type_gene_list(), ".med.h5")
          net <- paste0(network_path_gene_list(), network_type_gene_list(), ".net.h5")
          if (!file.exists(genes)) {
            errorMess <- paste("Please ensure", err_genes, "exists in", network_type_gene_list(), "folder")
            shinyalert(title = "Missing Network File", text = errorMess, type = "error")
          }
          else if (!file.exists(median)) {
            errorMess <- paste("Please ensure", err_median, "exists in", network_type_gene_list(), "folder")
            shinyalert(title = "Missing Network File", text = errorMess, type = "error")
          }
          else if (!file.exists(net)) {
            errorMess <- paste("Please ensure", err_net, "exists in", network_type_gene_list(), "folder")
            shinyalert(title = "Missing Network File", text = errorMess, type = "error")
          }
          else {
            sn$sub_nets <- subset_network_hdf5_gene_list(gene_list, tolower(network_type_gene_list()), dir=network_path_gene_list(), flag_occr = FALSE)
            success <- TRUE
          }
        }
        if (success == TRUE) {
            show(id = "CG_dropdown")
            hide(id = "CG_error")
            show(id = "GC_dropdown")
            hide(id = "GC_error")
            show(id = "FO_dropdown")
            hide(id = "FO_error")
            show(id = "GL_GSEA_dropdown")
            hide(id = "GL_GSEA_error")
            # Clear data
            output$network <- NULL
            output$heatmap <- NULL
            output$Bheatmap <- NULL
            output$CG_table <- NULL
            output$GCdensityG <- NULL
            output$GChistogramG <- NULL
            output$GCdensitySubsetG<- NULL
            output$GChistogramSubsetG <- NULL
            output$FO_heatmap <- NULL
            output$FO_network <- NULL
            output$genes_not_keep_table <- NULL
            output$genes_keep_table <- NULL
            # Reset Checkboxes
            updateAwesomeCheckboxGroup(
              inputId = "clusterPlotOptions_genelist",
              choices = c("Network", "Heatmap", "Binarized Heatmap"),
              status = ""
            )
            updateAwesomeCheckboxGroup(
              inputId = "GCPlotOptions_genelist",
              choices = c("Density", "Histogram", "Clustered Density", "Clustered Histogram"),
              status = ""
            )
            updateAwesomeCheckboxGroup(
              inputId = "FOPlotOptions_genelist",
              choices = c("Network", "Heatmap"),
              status = ""
            )
            updateAwesomeCheckboxGroup(
              inputId = "FO_table_options",
              choices = c("Functional Outliers", "Genes in Module"),
              status = ""
            )
          
          output$subnetwork <- renderTable(sn$sub_nets)

          show(id="GL_subnet_download")
          GL_subnet_table <- function(){sn$sub_nets}
          output$GL_subnet_download <- downloadHandler(
            filename = function() {
              paste("gene_list_subnetwork", input$ASSESS_GENE_LIST_download_table_format, sep="")
            },
            content = function(file) {
              write.table(GL_subnet_table(), file, row.names = FALSE, sep = separator, col.names = TRUE)
            }
          )
        }
      } 
    }
  })


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
    {input$run},
    {
      sub_net <- sn$sub_nets$sub_net
      node_degrees <- sn$sub_nets$node_degrees
      medK <- as.numeric(sn$sub_nets$median)

      # network output
      show(id="CG_network_text")
      show(id="CG_network_download")
      CG_network <- function(){plot_network(sub_net$genes, clust_net()$genes, medK)}
      output$network <- renderPlot(
        {CG_network()},
        width = 500,
        height = 500
      )



      # heatmap output
      show(id="CG_heatmap_text")
      show(id="CG_heatmap_download")
      CG_heatmap <- function(){plot_coexpression_heatmap(sub_net$genes, clust_net()$genes, flag_plot_bin = FALSE)}
      output$heatmap <- renderPlot(
        {CG_heatmap()},
        width = 500,
        height = 500
      )


      # binarized heatmap output
      show(id="CG_bheatmap_text")
      show(id="CG_bheatmap_download")
      CG_bheatmap <- function(){plot_coexpression_heatmap(sub_net$genes, clust_net()$genes)}
      output$Bheatmap <- renderPlot(
        {CG_bheatmap()},
        width = 500,
        height = 500
      )


      # clustering genes table output
      show(id="CG_table_text")
      show(id="table_genelist_download")
      CG_table <- function(){{
        is_ENSG <- grep('^ENSG', head(clust_net()$genes$clusters$genes))
          is_entrez <- grep("^\\d+", head (clust_net()$genes$clusters$genes))

          if ( length(is_ENSG) > 0 ) { 
            EGAD::attr.human[match(clust_net()$genes$clusters$genes, EGAD::attr.human$ensemblID),]

          } else { 
            if ( length(is_entrez) > 0 ) {
              EGAD::attr.human[match(clust_net()$genes$clusters$genes, EGAD::attr.human$entrezID),]
            } else { 
              EGAD::attr.human[match(clust_net()$genes$clusters$genes, EGAD::attr.human$name),]
            }

          } 
      }}
      output$CG_table <- renderDataTable(
        {CG_table()},
        # options=list(columnDefs = list(list(visible=FALSE, targets=c(0,1,2,3))))
      )
    

    #------------------ DOWNLOAD CG ----------------------#
    #Download plots    
      output$CG_network_download <- downloadHandler(
        filename = function() {
          paste("plot_coexpression_network_up", input$ASSESS_GENE_LIST_download_format)
        },
        content = function(file) {
          if (input$ASSESS_GENE_LIST_download_format == ".png") {
            png(file, width=3000, height=3000)
          } else if (input$ASSESS_GENE_LIST_download_format == ".pdf") {
            pdf(file)
          }
          CG_network()
          dev.off()
        }
      )
      output$CG_heatmap_download <- downloadHandler(
  
        filename = function() {
          paste("plot_coexpression_heatmap_genes", input$ASSESS_GENE_LIST_download_format)
        },
        content = function(file) {
          if (input$ASSESS_GENE_LIST_download_format == ".png") {
            png(file, width=3000, height=3000)
          } else if (input$ASSESS_GENE_LIST_download_format == ".pdf") {
            pdf(file)
          }
          CG_heatmap()
          dev.off()
        }
      )

      output$CG_bheatmap_download <- downloadHandler(
        filename = function() {
          paste("plot_coexpression_binarized_heatmap_up", input$ASSESS_GENE_LIST_download_format)
        },
        content = function(file) {
          if (input$ASSESS_GENE_LIST_download_format == ".png") {
            png(file, width=3000, height=3000)
          } else if (input$ASSESS_GENE_LIST_download_format == ".pdf") {
            pdf(file)
          }
          CG_bheatmap()
          dev.off()
        }
      )

    output$table_genelist_download <- downloadHandler(
      filename = function() {
        paste("clustered", input$ASSESS_GENE_LIST_download_table_format, sep="")
      },
      content = function(file) {
        write.table(CG_table(), file, row.names = TRUE, sep = separator, col.names = TRUE)
      }
    )
    }
    
  )


  ##################### GENE CONNECTIVITY #####################

  observeEvent(
    {input$runGC},
    {
      
      sub_net <- sn$sub_nets$sub_net
      node_degrees <- sn$sub_nets$node_degrees  
      medK <- as.numeric(sn$sub_nets$median)
      m <- match(clust_net()$genes$clusters$genes, rownames(sub_net$genes))

      # density output
      show(id="GCdensityG_text")
      show(id="GC_density_download")
      GC_densityG <- function(){plot_scatter(node_degrees$genes[,1]/node_degrees$n_genes_total, 
                    node_degrees$genes[,2]/node_degrees$n_genes, 
                    xlab="Global node degree", 
                    ylab="Local node degree", flag= "density")}
      output$GCdensityG <- renderPlot(
        {GC_densityG()},
        width = 500,
        height = 500
      )


      # histogram output
      show(id="GChistogramG_text")
      show(id="GC_histogram_download")
      GC_histogramG <- function(){plot_scatter(node_degrees$genes[,1]/node_degrees$n_genes_total, 
                    node_degrees$genes[,2]/node_degrees$n_genes, 
                    xybreaks = input$xybreaks,
                    xlab="Global node degree", 
                    ylab="Local node degree", flag= "hist")}
      output$GChistogramG <- renderPlot(
        {GC_histogramG()},
        width = 500,
        height = 500
      )

      show(id="GCdensitySubsetG_text")
      show(id="GC_densitySubset_download")
      # density output - subset by clusters
      GC_densitySubsetG <- function(){plot_scatter(node_degrees$genes[m,1]/node_degrees$n_genes_total, 
                      node_degrees$genes[m,2]/node_degrees$n_genes, 
                      xlab="Global node degree", 
                      ylab="Local node degree", 
                      clusters = clust_net()$genes$clusters, flag = "density")}
      output$GCdensitySubsetG <- renderPlot(
        {GC_densitySubsetG()},
        width = 500,
        height = 500
      )


      # histogram output - subset by clusters
      show(id="GChistogramSubsetG_text")
      show(id="GC_histogramSubset_download")
      GC_histogramSubsetG <- function(){plot_scatter(node_degrees$genes[m,1]/node_degrees$n_genes_total, 
                      node_degrees$genes[m,2]/node_degrees$n_genes, 
                      xybreaks = input$xybreaks,
                      xlab="Global node degree", 
                      ylab="Local node degree", 
                      clusters = clust_net()$genes$clusters, flag = "hist")}
      output$GChistogramSubsetG <- renderPlot(
        {GC_histogramSubsetG()},
        width = 500,
        height = 500
      )

    
    #------------------ DOWNLOAD GC ----------------------#
    #Download plots    
      output$GC_density_download <- downloadHandler(
        filename = function() {
          paste("plot_scatter_density_genes", input$ASSESS_GENE_LIST_download_format)
        },
        content = function(file) {
          if (input$ASSESS_GENE_LIST_download_format == ".png") {
            png(file, width=3000, height=3000)
          } else if (input$ASSESS_GENE_LIST_download_format == ".pdf") {
            pdf(file)
          }
          GC_densityG()
          dev.off()
        }
      )

      output$GC_histogram_download <- downloadHandler(
        filename = function() {
          paste("plot_scatter_hist_genes", input$ASSESS_GENE_LIST_download_format)
        },
        content = function(file) {
          if (input$ASSESS_GENE_LIST_download_format == ".png") {
            png(file, width=3000, height=3000)
          } else if (input$ASSESS_GENE_LIST_download_format == ".pdf") {
            pdf(file)
          }
          GC_histogramG()
          dev.off()
        }
      )

      output$GC_densitySubset_download <- downloadHandler(
        filename = function() {
          paste("plot_scatter_density_genes_coloured", input$ASSESS_GENE_LIST_download_format)
        },
        content = function(file) {
          if (input$ASSESS_GENE_LIST_download_format == ".png") {
            png(file, width=3000, height=3000)
          } else if (input$ASSESS_GENE_LIST_download_format == ".pdf") {
            pdf(file)
          }
          GC_densitySubsetG()
          dev.off()
        }
      )

      output$GC_histogramSubset_download <- downloadHandler(
        filename = function() {
          paste("plot_scatter_hist_genes_coloured", input$ASSESS_GENE_LIST_download_format)
        },
        content = function(file) {
          if (input$ASSESS_GENE_LIST_download_format == ".png") {
            png(file, width=3000, height=3000)
          } else if (input$ASSESS_GENE_LIST_download_format == ".pdf") {
            pdf(file)
          } 
          GC_histogramSubsetG()
          dev.off()
        }
      )
    }
  )


  ##################### FUNCTIONAL OUTLIERS #####################

  observeEvent(
    {input$runFO},
    {

      sub_net <- sn$sub_nets$sub_net
      node_degrees <- sn$sub_nets$node_degrees  
      medK <- as.numeric(sn$sub_nets$median)

      filt_min <- input$filtmin

      clust_size <- plyr::count(clust_net()$genes$clusters$labels)
      clust_keep <- clust_size[clust_size[,2] < filt_min ,1]
      genes_keep <- !is.na(match(clust_net()$genes$clusters$labels, clust_keep))
     

      # heatmap output
      show(id="FO_heatmap_text")
      show(id="FO_network_download")
      FO_heatmap <- function(){plot_coexpression_heatmap(sub_net$genes, clust_net()$genes, filt = TRUE, flag_plot_bin = FALSE, filt_min = filt_min)}
      output$FO_heatmap <- renderPlot(
        {FO_heatmap()},
        width = 500,
        height = 500
      )

      # network output
      show(id="FO_network_text")
      show(id="FO_heatmap_download")
      FO_network <- function(){plot_network(1-sub_net$genes, clust_net()$genes, 1 - medK, filt_min = filt_min)}
      output$FO_network <- renderPlot(
        {FO_network()},
        width = 500,
        height = 500
      )

      # genes in module table output
      show(id="genes_not_keep_table_text")
      show(id="genes_not_keep_table_download")
      FO_genes_not_keep <- function(){
        is_ENSG <- grep('^ENSG', head(clust_net()$genes$clusters$genes))
          is_entrez <- grep("^\\d+", head (clust_net()$genes$clusters$genes))

          if ( length(is_ENSG) > 0 ) { 
            EGAD::attr.human[match(clust_net()$genes$clusters$genes[!genes_keep],EGAD::attr.human$ensemblID),]

          } else { 
            if ( length(is_entrez) > 0 ) {
              EGAD::attr.human[match(clust_net()$genes$clusters$genes[!genes_keep],EGAD::attr.human$entrezID),]
            } else { 
              EGAD::attr.human[match(clust_net()$genes$clusters$genes[!genes_keep],EGAD::attr.human$name),]
            }

          }
      }
      output$genes_not_keep_table <- renderDataTable(
        {FO_genes_not_keep()},
      )

      # functional outliers table output
      show(id="genes_keep_table_text")
      show(id="genes_keep_table_table_download")
      FO_genes_keep <- function(){
        is_ENSG <- grep('^ENSG', head(clust_net()$genes$clusters$genes))
          is_entrez <- grep("^\\d+", head (clust_net()$genes$clusters$genes))

          if ( length(is_ENSG) > 0 ) { 
            EGAD::attr.human[match(clust_net()$genes$clusters$genes[genes_keep],EGAD::attr.human$ensemblID),]

          } else { 
            if ( length(is_entrez) > 0 ) {
              EGAD::attr.human[match(clust_net()$genes$clusters$genes[genes_keep],EGAD::attr.human$entrezID),]
            } else { 
              EGAD::attr.human[match(clust_net()$genes$clusters$genes[genes_keep],EGAD::attr.human$name),]
            }

          }
      }
      output$genes_keep_table <- renderDataTable(
        {FO_genes_keep()},
      )

     #------------------ DOWNLOAD FO ----------------------#
    #Download plots    
      output$FO_network_download <- downloadHandler(
        filename = function() {
          paste("plot_scatter_density_genes", input$ASSESS_GENE_LIST_download_format)
        },
        content = function(file) {
          if (input$ASSESS_GENE_LIST_download_format == ".png") {
            png(file, width=3000, height=3000)
          } else if (input$ASSESS_GENE_LIST_download_format == ".pdf") {
            pdf(file)
          }
          FO_network()
          dev.off()
        }
      )

      output$FO_heatmap_download <- downloadHandler(
        filename = function() {
          paste("plot_scatter_hist_genes", input$ASSESS_GENE_LIST_download_format)
        },
        content = function(file) {
          if (input$ASSESS_GENE_LIST_download_format == ".png") {
            png(file, width=3000, height=3000)
          } else if (input$ASSESS_GENE_LIST_download_format == ".pdf") {
            pdf(file)
          }
          FO_heatmap()
          dev.off()
        }
      )

      output$genes_not_keep_table_download <- downloadHandler(
      filename = function() {
        paste("genes_not_keep", input$ASSESS_GENE_LIST_download_table_format, sep="")
      },
      content = function(file) {
        write.table(FO_genes_not_keep(), file, row.names = TRUE, sep = separator, col.names = TRUE)
      }
    )

    output$genes_keep_table_table_download <- downloadHandler(
      filename = function() {
        paste("genes_keep", input$ASSESS_GENE_LIST_download_table_format, sep="")
      },
      content = function(file) {
        write.table(FO_genes_keep(), file, row.names = TRUE, sep = separator, col.names = TRUE)
      }
    )
      
    }
  )

  # Gene List GSEA
  observeEvent(
    {input$GL_GSEA_run},
    {
      # heatmap
      show(id = "GL_GSEA_heatmap_text")

      data(go_voc)
      gene_list <- clust_net()$genes$clusters$genes[clust_net()$genes$order]
      
      if (input$gene_list_selection == "Generate Gene List" || input$GL_gene_list_type == "Gene Names") {
        data(go_slim)
        output$GL_GSEA_heatmap_plot <- renderPlot(
          {
            filt <- colSums(go_slim) < 5000 & colSums(go_slim) >= 10
            go_enrich <- gene_set_enrichment(gene_list, go_slim[filt,], go_voc)
            plot_gene_set_enrichment(go_enrich, gene_list, go_slim[filt,])
          },
          width = 500,
          height = 500
        )
        output$GL_GSEA_heatmap_compare <- renderPlot(
          {
            filt <- colSums(go_slim) < 5000 & colSums(go_slim) >= 10
            go_enrich <- gene_set_enrichment(gene_list, go_slim[filt,], go_voc)
            plot_gene_set_enrichment(go_enrich, gene_list, go_slim[filt,])
          },
          width = 500,
          height = 500
        )
        output$GL_GSEA_heatmap_cluster_compare <- renderPlot(
          {
            sub_net <- sn$sub_nets$sub_net
            node_degrees <- sn$sub_nets$node_degrees
            medK <- as.numeric(sn$sub_nets$median)
            plot_coexpression_heatmap(sub_net$genes, clust_net()$genes, flag_plot_bin = FALSE)
          },
            width = 500,
            height = 500
        )
      } else {
        data(go_slim_entrez)
        output$GL_GSEA_heatmap_plot <- renderPlot(
          {
            filt <- colSums(go_slim_entrez) < 5000 & colSums(go_slim_entrez) >= 10
            go_enrich <- gene_set_enrichment(gene_list, go_slim_entrez[filt,], go_voc)
            plot_gene_set_enrichment(go_enrich, gene_list, go_slim_entrez[filt,])
          },
          width = 500,
          height = 500
        )
        output$GL_GSEA_heatmap_compare <- renderPlot(
          {
            filt <- colSums(go_slim_entrez) < 5000 & colSums(go_slim_entrez) >= 10
            go_enrich <- gene_set_enrichment(gene_list, go_slim_entrez[filt,], go_voc)
            plot_gene_set_enrichment(go_enrich, gene_list, go_slim_entrez[filt,])
          },
          width = 500,
          height = 500
        )
        output$GL_GSEA_heatmap_cluster_compare <- renderPlot(
          {
            sub_net <- sn$sub_nets$sub_net
            node_degrees <- sn$sub_nets$node_degrees
            medK <- as.numeric(sn$sub_nets$median)
            plot_coexpression_heatmap(sub_net$genes, clust_net()$genes, flag_plot_bin = FALSE)
          },
            width = 500,
            height = 500
        )
      }
      
      

      # heatmap
      show(id="GSEA_heatmap_text")
      show(id="GSEA_heatmap_download")
      GSEA_heatmap <- function(){
          filt <- colSums(go_slim) < 5000 & colSums(go_slim) >= 10
          gene_list <- clust_net()$genes$clusters$genes[clust_net()$genes$order]
          go_enrich <- gene_set_enrichment(gene_list, go_slim[filt,], go_voc)
          plot_gene_set_enrichment(go_enrich, gene_list, go_slim[filt,])}
      output$GL_GSEA_heatmap <- renderPlot(
        {GSEA_heatmap()},
        width = 500,
        height = 500
      )
      
      #Download plots    
      output$GSEA_heatmap_download <- downloadHandler(
        filename = function() {
          paste("GSEA_heatmap", input$ASSESS_GENE_LIST_download_format)
        },
        content = function(file) {
          if (input$ASSESS_GENE_LIST_download_format == ".png") {
            png(file, width=1500, height=1500)
          } else if (input$ASSESS_GENE_LIST_download_format == ".pdf") {
            pdf(file)
          }
          GSEA_heatmap()
          dev.off()
        }
      )
    
    
    }
  )



  ##################### ERROR MESSAGES #####################

  output$CG_error <- renderText({
    print("Please generate a subnetwork in NETWORK OPTIONS")
  })

  output$CG_error_DE <- renderText({
    print("Please generate a subnetwork in NETWORK OPTIONS")
  })

  output$GC_error <- renderText({
    print("Please generate a subnetwork in NETWORK OPTIONS")
  }) 

  output$GC_error_DE <- renderText({
    print("Please generate a subnetwork in NETWORK OPTIONS")
  })

  output$FO_error <- renderText({
    print("Please generate a subnetwork in NETWORK OPTIONS")
  }) 

  output$FO_error_DE <- renderText({
    print("Please generate a subnetwork in NETWORK OPTIONS")
  })

  output$DE_GSEA_error <- renderText({
    print("Please generate a subnetwork in NETWORK OPTIONS")
  }) 

  output$GL_GSEA_error <- renderText({
    print("Please generate a subnetwork in NETWORK OPTIONS")
  }) 

}



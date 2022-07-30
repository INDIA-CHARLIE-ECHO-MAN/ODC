library(shinyWidgets)
library(bslib)
library(DT)
library(shiny)
library(stats)
library(gplots)
library(graphics)
library(viridis)
library(utils)
library(shinybusy)
library(shinyjs)
library(shinyalert)
library(stringi)
library(stringr)
library(ggplot2)
library(ggplotify)
library(OutDeCo)
library(EGAD)

ui <- fluidPage(

  ######################################################################
  #                                                                    #
  #                               STYLE                                #
  #                                                                    #
  ######################################################################

  div(style="display:inline-block; float:right", circleButton("help", icon=icon("question"), status="default", size="default")),

  useShinyjs(),
  chooseSliderSkin("Flat",  color="#325D88"),
  add_busy_spinner(spin="dots", position="bottom-right", color="#3E3F3A"),
  
  titlePanel(title=div(img(src="ODClogo.png", height=80), "OutDeCo")),
 
  theme=bs_theme(version=5, bootswatch="sandstone", 
                  heading_font=font_google("Poppins"), 
                  base_font=font_collection(font_google("Roboto")),
                  success="#325D88"),

  tags$style(HTML(".js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {background: #3E3F3A}")),
  # Style for Download Button
  tags$head(tags$style(" .download_style{color: #1D89FF;} .download_style{font-family: Poppins}")),
  # Style for all font
  tags$head(tags$style(HTML('
        body, label, input, button, select { 
          font-family: "Poppins"; 
        }'))),
  
  
  navbarPage(title=NULL, id="navpage", collapsible=FALSE, 
            
            ######################################################################
            #                                                                    #
            #                             HOME TAB                               #
            #                                                                    #
            ######################################################################

            tabPanel(title="Home",
              icon=icon("home"),
  
              navlistPanel(
                id="Header", selected=NULL, well=FALSE, fluid=TRUE, widths=c(3, 9),

                ########################### ABOUT ###########################

                tabPanel(title="About",
                  h3("What is OutDeCo?"),
                  br(),
                  p(strong("Outlier Detection through Co-expression"), "is an R package developed by Dr Sara Ballouz."),
                  p("The purpose of this package is to assess genes - more specifically differentially expressed genes, with respect to their co-expression properties."),
                  p("Using the tools on this website you can:"),
                  tags$li("Run a differential expression (DE) analysis."),
                  tags$li("Assess your DE results."),
                  tags$li("Assess an uploaded/generated gene list."),
                  br(),
                  tags$a(href="https://github.com/ballouzlab/OutDeCo_lite", "Github source code"),
                  br(),
                  br(),
                  h5(strong("Differential Expression"), "vs", strong("Gene List"), "?"),
                  tags$li("Assess DE will have analysis options for up and down regulated genes"),
                  tags$li("Assess Gene List will have analysis options for genes"),

                ),
                
                ########################### DIFFERENTIAL ANALYSIS ###########################

                tabPanel(title="Differential Expression Analysis",
                  div(style="display:inline-block; float:right", circleButton("DE_return", icon=icon("angle-right"), status="default", size="default")),
                  h3("Differential Expression Analysis"),
                  br(),
                  p(strong("Differential Expression Analysis"), "is a statistical analysis to discover quantitative changes in expression levels between experimental groups."),
                  br(),

                  mainPanel(
                    tabsetPanel(

                      tabPanel(title="User Guide",
                        br(),
                        h5(strong("Files/Data")),
                        splitLayout(cellWidths=c("30%", "70%"), 
                            fluidPage(
                              img(src="run_de_file_options.jpg", height=200),
                            ),
                            fluidPage(
                              h6(strong("Counts")),

                              tags$div("A matrix in either csv, tsv, txt format",tags$br(),
                              "- Columns contain individual samples  labelled by sample ID", tags$br(),
                              "- Rows contain genes                labelled by entrez ID ", ),
                            
                              h6(strong("Labels")),
                              tags$div("A matrix in either csv, tsv, txt format",tags$br(),
                              "- Columns  labelled by category/label", tags$br(),
                              "- Rows  labelled by run", ),
                            ),
                          ), 
                        br(),
                        br(),

                        h5(strong("DE Methods")),
                        img(src="DE_method.png", height=150),
                        h6(strong("wilcox")),
                        p("Compares means of two groups to analyse count data and test for differential expression."),
                        h6(strong("DESeq")),
                        p("Uses geometric normalisation to analyse count data and test for differential expression."),
                        h6(strong("edgeR")),
                        p("Uses weighted mean of log ratio to analyse count data and test for differential expression."),
                        br(),
                        br(),

                        h5(strong("Case/Conditions Selection:")),
                        img(src="choose_label.jpg", height=250),
                        h6(strong("By label")),
                        tags$div("- Cases are the rows where the chosen column (e.g. Relationship) has row equal to the chosen case (e.g. fthr)",tags$br(),
                              "- Condition are the remaining rows", tags$br(), ),
                        
                        h6(strong("Individually")),
                        p("Cases and Conditions are selected by clicking the associated row on the labels table"),
                      ),

                      tabPanel(title="Output",
                        br(),
                        h5(strong("Plot Types")),
                        h6(strong("Volcano Plot")),
                        img(src="volcano_plot.png", height=300),
                        br(),
                        br(),
                        h6(strong("MA Plot")),
                        img(src="MA_plot.png", height=300),
                      ),

                    )
                  ),
                  
                
                  

                ),

                ########################### GENERATE SUBNETWORK ###########################
                
                tabPanel(title="Generate Subnetwork",
                  div(style="display:inline-block; float:right", circleButton(inputId="SN_return", icon=">DE", status="default", size="default")),
                  div(style="display:inline-block; float:right", circleButton(inputId="SN_return_GL", icon=">GL", status="default", size="default")),
                  h3("Generating a Subnetwork"),
                  
                  br(),
                  tags$div(strong("Generating a subnetwork"), "with a chosen network produces a matrix of weights between genes.",tags$br(),
                          "This subnetwork can then be clustered, analysed for gene connectivity, functional outliers and gene set enrichment analysis.", p("")),

                  p(em("A subnetwork must be generated before assessing Differentially Expressed Data or a Gene List.")),
                  br(),

                  mainPanel(
                    tabsetPanel(

                      tabPanel(title="User Guide",
                        br(),
                        h5(strong("Step 1")),
                        p("Select ",  img(src="network_options.png", height=40), "Button."),
                        br(),

                        h5(strong("Step 2")),
                        p("Ensure the", code("networks"), "folder is in parent directory of the package."),
                        p("Inside", code("networks"), "are subfolders titled with the name of the network and contain network files in HDF5 (.h5) format."), 
                        br(),

                        
                        div(style="position:relative; left:calc(6%);", 
                          fluidPage(
                            p("Network Files are named after the folder and consist of..."),
                            tags$li("matrix of network itself (suffix", code(".net.h5"), ")"), 
                            tags$li("the genes (suffix", code(".genes.h5"),")"), 
                            tags$li("the median value of their network (suffix", code(".med.h5"), ")"), 
                          )
                        ),
                        
                        br(),
                        
                        h5(strong("Step 3")),
                        p("Select whether you wish to use", strong("occur"), "or", strong("standard"), "network."),

                        div(style="position:relative; left:calc(6%);", 
                          fluidPage(
                            tags$li(strong("Occr"), " - a tally of how often we see genes co-expressed"),
                            tags$li(strong("Standard"), "- an average of the weights/edges of coexpression values"),
                          )
                        ),
                        
                        br(),
                        
                        h5(strong("Step 4")),
                        p("Select ", strong("data.")), 

                        div(style="position:relative; left:calc(6%);", 
                          fluidPage(
                            h6(strong("Assess DE")),
                            tags$li(em("Use DE Results"), "- Takes the data from RUN DE"),  
                            tags$li(em("Upload DE Data"), "- A matrix in either csv, tsv, txt format. Columns required for", code("mean_cmp, lof2_fc, pvals, padj")), 
                            h6(strong("Assess Gene List:")),
                            tags$li(em("Upload Gene List"), "- A txt file containing a list of genes. Genes can be represented as: geneIDs, ensemblIDs or entrezIDs"), 
                            tags$li(em("Generate Gene List"), "- Requires a human chromosome and the number of genes from that chromosome."),
                          )
                        ),
                        br(),
                        
                        
                        h5(strong("Step 5")),
                        p("5) Click ", img(src="generate_subnetwork.png", height=40), "."),

                      ),

                      tabPanel(title="Output",
                        br(),
                        h5(strong("Subnetwork Table")),
                        br(),
                        img(src="subnetwork_table.png", height=100),
                      ),

                    ),
                  
                  br(), 
                  ),
                
                ),

                ########################### CLUSTER GENES ###########################

                tabPanel(title="Cluster Genes",
                  div(style="display:inline-block; float:right", circleButton(inputId="CG_return", icon=">DE", status="default", size="default")),
                  div(style="display:inline-block; float:right", circleButton(inputId="CG_return_GL", icon=">GL", status="default", size="default")),
                  h3("Cluster Genes"),
                  br(),
                  p("The", strong("Cluster Genes"), "function generates modules which are clusters of genes that are hightly co-expressed."),
                  br(),


                  mainPanel(
                    tabsetPanel(

                      tabPanel(title="User Guide",
                        br(),
                        # tags$div("If you have not generated a subnetwork, the following error will appear.", img(src="nosub_error.jpg"), "Make sure you have generated a subnetwork before continuing to these steps."),
                        # br(),

                        h5(strong("Step 1")),
                        p("Select the", img(src="options.jpg", height=40), "dropdown. The following options exist for clustering genes:"),
                        
                        div(style="position:relative; left:calc(6%);", 
                          fluidPage(
                            splitLayout(cellWidths=c("50%", "50%"), 
                              fluidPage(
                                p("For Assessing DE Data,", br(), "the following options exist:"),
                                h6(strong("Upregulated")), 
                                tags$li("Network"), 
                                tags$li("Heatmap"), 
                                tags$li("Binarized Heatmap"),
                                h6(strong("Downregulated")),
                                tags$li("Network"), 
                                tags$li("Heatmap"), 
                                tags$li("Binarized Heatmap"),
                                
                              ), 
                              fluidPage(
                                p("For Assessing gene lists,", br(), "the following options exist:"),
                                h6(strong("Genes")), 
                                tags$li("Network"), 
                                tags$li("Heatmap"), 
                                tags$li("Binarized Heatmap"),
                              ), 
                            ),
                          ),
                        ),
                        
                        br(), 
                        
                        h5(strong("Step 2")),
                        p("Select the plots you want to view by checking the checkbox associated with it."), 
                        br(),

                        h5(strong("Step 3")),
                        p("Click ", img(src="run_button.jpg"), "."),
                      ),

                      tabPanel(title="Output",
                        br(), 
                        h5(strong("Plots")), 
                        h6(strong("Heatmap")),  
                        p("Heatmap showing the co-expression of Genes. The trees on the top and left axis of the map represents the clusters and the colour key demonstrates the level of co-expression for each gene pair."),
                        img(src="heatmap_of_genes_CG.jpg", height=400),
                        br(),

                        h6(strong("Network")),  
                        p("Network plot where nodes are genes and the weight of the edges corresponds to the co-expression levels between its endpoints."),
                        img(src="network_of_genes_CG.jpg", height=400),
                        br(),

                        h6(strong("Binarized heatmap")),  
                        p("Binary Heatmap showing co-expression of genes. All values in the heatmap are either 0 or 1 instead of ranging between 0 and 1, Yellow=1, there is a high level of co-expression between the genes, Purple=0, co-expression level of genes is 0."),
                        img(src="binarized_heatmap_CG.jpg", height=400),
                        br(), 
                        br(),

                        h5(strong("Tables")), 
                        p("Gene details of the genes that were involved in clustering. Each row is a different gene and the columns are different fields pertaining to the gene."),
                        img(src="gene_deetz_clustering.jpg", height=400), 
                      ),
                    ),
                  ),

                ),

                ########################### GENE CONNECTIVITY ###########################

                tabPanel(title="Gene Connectivity",
                  div(style="display:inline-block; float:right", circleButton(inputId="GC_return", icon=">DE", status="default", size="default")),
                  div(style="display:inline-block; float:right", circleButton(inputId="GC_return_GL", icon=">GL", status="default", size="default")),
                  h3("Gene Connectivity"),
                  br(),
                  p("The", strong("Gene Connectivity"), "function calculates node degrees to get a sense of the global and local connectivities of the gene."),
                  br(),


                  mainPanel(
                    tabsetPanel(

                      tabPanel(title="User Guide",
                        br(),
                        # tags$div("If you have not generated a subnetwork, the following error will appear.", img(src="nosub_error.jpg"), "Make sure you have generated a subnetwork before continuing to these steps."),
                        # br(),

                        h5(strong("Step 1")),
                        p("Select the",  img(src="options.jpg", height=40), "dropdown."),

                        p("Options for visualisation of gene connectivity include:", strong("Density,"), strong("Clustered Density"), "and", strong("Clustered Histogram."), "When any histogram is chosen, a slider appears for the user to change the number of breaks in the histogram."), 
                        
                        div(style="position:relative; left:calc(6%);", 
                          fluidPage(
                            splitLayout(cellWidths=c("50%", "50%"), 
                              fluidPage(
                                p("For Assessing DE Data,", br(), "the user has the following options"),
                                h6(strong("Upregulated")), 
                                tags$li("Density"), 
                                tags$li("Histogram"), 
                                tags$li("Clustered Density"),
                                tags$li("Clustered Histogram"),
                                h6(strong("Downregulated")),
                                tags$li("Density"), 
                                tags$li("Histogram"), 
                                tags$li("Clustered Density"),
                                tags$li("Clustered Histogram"),
                              ), 
                              fluidPage(
                                p("For Assessing a Gene List,", br(), "the user has the following options"),
                                h6(strong("Genes")), 
                                tags$li("Density"), 
                                tags$li("Histogram"), 
                                tags$li("Clustered Density"),
                                tags$li("Clustered Histogram"),
                              ),
                            ),
                          ),
                        ),
                        
                        br(), 
                        
                        h5(strong("Step 2")),
                        p("Select the plots you want to view by checking the checkbox associated with it.", p("")),
                        br(),

                        h5(strong("Step 3")),
                        p("Click ", img(src="run_button.jpg"), "."),
                      ),

                      tabPanel(title="Output",
                        br(), 
                        h5(strong("Plots")),
                        h6(strong("Density")),  
                        p("Density plot of gene connectivity."),
                        img(src="plot_scatter_density_up.png", height=400), 
                        br(),

                        h6(strong("Histogram")),
                        p("Histogram of gene connectivity."),
                        img(src="plot_scatter_hist_up.png", height=400),
                        br(),

                        h6(strong("Subset by clusters")),  
                        p("Gene connectivity can be subset by their clusters as well."),
                        img(src="plot_scatter_hist_down_colored.png", height=400), 
                      ),

                    ),
                  ),

              
                   
                ),
                 
                ########################### FUNCTIONAL OUTLIERS ###########################

                tabPanel(title="Functional Outliers",
                  div(style="display:inline-block; float:right", circleButton(inputId="FO_return", icon=">DE", status="default", size="default")),
                  div(style="display:inline-block; float:right", circleButton(inputId="FO_return_GL", icon=">GL", status="default", size="default")),
                  h3("Functional Outliers"),
                  br(),
                  p("The", strong("Functional Outliers"), "function outputs genes that have been identified to be potentially dysregulated. 
                  They are the genes that are Differentially Expressed but do not show local co-expression."),
                  p(em("Module Default Threshold: More than 6 genes.")),
                  br(),


                  mainPanel(
                    tabsetPanel(

                      tabPanel(title="User Guide",
                        br(),
                        # p("If you have not generated a subnetwork, the following error will appear.", img(src="nosub_error.jpg"), "Make sure you have generated a subnetwork before continuing to these steps."),
                        
                        h5(strong("Step 1")),
                        p("Select the", img(src="options.jpg", height=40), "dropdown."),
                        p("Options for visualising functional outliers include: ", strong("Network, Heatmap, Genes in Modules,"), "and", strong("Functional Outliers."), "A slider will also appear, giving the user the option of choosing how many genes are required to form a module (default=6)."),
                        
                        div(style="position:relative; left:calc(6%);", 
                          fluidPage(
                            splitLayout(cellWidths=c("50%", "50%"),
                              fluidPage(
                                p("For Assessing DE Data,", br(), "the user has the following options"),
                                h6(strong("Upregulated")), 
                                tags$li("Network"), 
                                tags$li("Heatmap"), 
                                tags$li("Genes in Module"),
                                tags$li("Functional Outliers"),
                                h6(strong("Downregulated")),
                                tags$li("Network"), 
                                tags$li("Heatmap"), 
                                tags$li("Genes in Module"),
                                tags$li("Functional Outliers"),
                              ),
                              fluidPage(
                                p("For Assessing a Gene List,", br(), "the user has the following options"),
                                h6(strong("Genes")), 
                                tags$li("Network"), 
                                tags$li("Heatmap"), 
                                tags$li("Genes in Module"),
                                tags$li("Functional Outliers"),
                              ),
                            ),
                          ),
                        ),
                        
                        br(), 
                        
                        h5(strong("Step 2")),
                        p("Select the plots you want to view by checking the checkbox associated with it."),
                        br(),

                        h5(strong("Step 3")),
                        p("Click ", img(src="run_button.jpg"), "."),
                      ),

                      tabPanel(title="Output",
                        br(),
                        h5(strong("Plots")),
                        h6(strong("Coexpression Heatmap")),
                        p("Heatmap of genes detailing the outliers connectivity and expression."),
                        img(src="plot_coexpression_heatmap_down_filt.png", height=400),
                        br(),

                        h6(strong("Network")),
                        p("Network plot detailing the outliers connectivity and expression."),
                        img(src="plot_network_down.png", height=400),
                        br(),
                        br(),

                        h5(strong("Tables")),
                        h6(strong("Genes in Module")),
                        p("These are the gene details of the genes that were filtered away and formed modules with other genes."),
                        img(src="gene_deetz_notkeep_FO.jpg", height=400),
                        br(),

                        h6(strong("Functional Outliers")),
                        p("These are the gene details of the genes that were remaining once the other genes formed modules also known as Functional Outliers."),
                        img(src="genes_keep_deetz_FO.jpg", height=400),
                      ),

                    ),
                  ),

                ),
                
                ########################### GSEA ###########################

                tabPanel(
                  div(style="display:inline-block; float:right", circleButton(inputId="GSEA_return", icon=">DE", status="default", size="default")),
                  div(style="display:inline-block; float:right", circleButton(inputId="GSEA_return_GL", icon=">GL", status="default", size="default")),
                  title="Gene Set Enrichment Analysis",
                  h3("Gene Set Enrichment Analysis (GSEA)"),
                  br(),
                  p(strong("GSEA"), "is a process of ranking genes by how statistically significant their differential gene expression is. 
                  This can remove false positives from the data."),
                  br(),

                  mainPanel(
                    tabsetPanel(

                      tabPanel(title="User Guide",
                        br(),
                        # p("If you have not generated a subnetwork, the following error will appear.", img(src="nosub_error.jpg"), "Make sure you have generated a subnetwork before continuing to these steps."),
                        
                        h5(strong("Step 1")),
                        p("Select the",  img(src="options.jpg", height=40), "dropdown."),
                        p("There are two main types of Gene Set Enrichment Analysis that can occur. Standard GSEA and AUCS GSEA. AUCS GSEA can only be run on DE Data. The following options for GSEA include:"),
                        
                        div(style="position:relative; left:calc(6%);", 
                          fluidPage(
                            splitLayout(cellWidths=c("50%", "50%"),
                              fluidPage( 
                                p("For Assessing DE Data,", br(), "the user has the following options"),
                                h6(strong("GSEA Type")), 
                                tags$li("Standard GSEA"), 
                                tags$li("AUCs GSEA"), 
                                h6(strong("Standard GSEA")),
                                em("(appears when Standard GSEA is checked)"),
                                tags$li("Upregulated P-value Heatmap"), 
                                tags$li("Downregulated P-value Heatmap"), 
                              ),
                              fluidPage(
                                p("For Assessing a Gene List,", br(), "the user has the following options"),
                                h6(strong("Standard GSEA")),
                                tags$li("P-value Heatmap"), 
                              ), 
                            ),
                          ),
                        ),                        
                        br(), 

                        h5(strong("Step 2")),
                        p("Select the plots you want to view by checking the checkbox associated with it."), 
                        br(),

                        h5(strong("Step 3")),
                        p("Click ", img(src="run_button.jpg"), "."),
                      ),

                      tabPanel(title="Output",
                        br(),
                        h5(strong("Standard GSEA")),
                        p("P-value heatmap of genes."),
                        img(src="go_enrich.png", height=400),
                        br(),
                        br(),

                        h5(strong("Ranking")),
                        tags$li("AUROC graph"),
                        tags$li("Scatter plot (AUCs vs P-values)"),
                        tags$li("Gene set size histogram"),
                        tags$li("AUC histogram"),
                        img(src="go_enrich_ranked.png", height=400),
                      ),

                    ),
                  ),

                 ),
                 
              )
            ),

            ######################################################################
            #                                                                    #
            #                            RUN DE TAB                              #
            #                                                                    #
            ######################################################################

            tabPanel(title="Run DE", 

              ########################### OPTIONS ###########################

              dropdown(
                tags$h3("Options"), 
                # side panel characteristics
                style="jelly", icon="FILE OPTIONS",
                status="primary", width="335px", size="sm",
               
                # title of sidepanel
                fluidPage(
                # inputs in the sidepanel
                  fileInput("counts_file", label="Choose Counts Data File"),
                  
                  radioButtons(
                    inputId='sepCountsButton', 
                    label='Delimiter Selector', 
                    choices=c(Comma=",", Semicolon=";", Tab="\t", Space=" "), 
                    selected=''
                  ),

                  fileInput("labels_file", "Choose Labels File",
                    accept=c(
                    "text/csv",
                    "text/comma-separated-values,text/plain",
                    ".csv")
                  ),

                  radioButtons(
                      inputId='sepLabelsButton', 
                      label='Delimiter Selector', 
                      choices=c(Comma=",", Semicolon=";", Tab="\t", Space=" "),
                      selected=''
                  ),
                ), 
              ),
              
              div(
                style="display:inline-block; float:right", 
                dropdown(
                  right=TRUE,
                  tags$h4("Download Options"),
                  inputId="download_dropdown",

                  style="minimal", icon="DOWNLOAD OPTIONS",
                  status="primary", width="300px", size="sm",
                  
                  selectInput(
                    inputId="download_format",
                    label= tags$h6("Choose Plot Download Format"),
                    choices=c(".png", ".pdf"),
                    selected=".png",
                    width="600px",
                  ),

                  selectInput(
                    inputId="download_table_format",
                    label= tags$h6("Choose Table Download Format"),
                    choices=c(".csv", ".tsv", ".txt"),
                    selected=".png",
                    width="600px",
                  ),                   
                   
                ),
              ),

              br(),

              navlistPanel(
                widths=c(3, 9), well=FALSE,

                ########################### VIEW FILES ###########################

                tabPanel(title="View File",

                  tabsetPanel(
                    id="counts_labels_tabset",

                    tabPanel(
                      title="Counts File",
                      dataTableOutput("UICountsContent")
                    ),

                    tabPanel(
                      title="Labels File",
                      dataTableOutput("UILabelContent")
                    ),

                  )
                ),

                ########################### RUN DE ###########################

                tabPanel(title="Run Differential Expression",
                  h3("Plot Differential Expression"),

                  p(id="runDE_error", "Please upload counts and labels in FILE OPTIONS.", style="color:red"),

                  dropdown(
                    inputId="DE_options",
                    # side panel characteristics
                    style="minimal", icon="OPTIONS",
                    status="primary", width="600px", size="sm",

                    selectInput(
                      inputId="DE_method",
                      label= tags$h5("Choose DE Method"),
                      choices=c("wilcox", "DESeq2", "edgeR"),
                      selected=NULL,
                      width="600px",
                    ),
                  
                    
                    radioButtons(
                      inputId="case_control_method",
                      label=tags$h5("Case/Control Selection"),
                      choices=c("Choose Case by Label", "Choose Case/Controls individually"),
                      selected=""
                    ),
                   
                    conditionalPanel(condition="input.case_control_method == 'Choose Case by Label'", 
                      selectInput(
                        inputId="select_column",
                        label= "Select label to group ",
                        choices=NULL, #no choice before uploading
                        width="600px",
                      ),
                  
                      selectInput(
                        inputId="select_case",
                        label= "Select case to analyse",
                        choices=NULL, #no choice before column selected
                        width="600px",
                      ),

                    ),

                    conditionalPanel(condition="input.case_control_method == 'Choose Case/Controls individually'", 
                      h6(strong("Select Cases")),
                      dataTableOutput("UILabelContentSelection"),   
                      h6(strong("Select Conditions")),
                      dataTableOutput("UILabelContentRemoveSelection"),                   
                    ),
                    
                    actionButton(inputId="run_DE", label="Run DE"),
                  
                  ),

                  br(),

                  tabsetPanel(

                    tabPanel(title="Plots",
                      splitLayout(
                        cellWidths=c("50%", "50%"), 
                        fluidPage(
                          #textOutput("DE_V_text"),
                          h4(id="vol"," Volcano Plot", style="text-align: center;"),
                          column(12, plotOutput(outputId="DEplot", height="450px"), align="center"), 
                          div(style="margin-left: 375px;", downloadLink("volcano_download", label="Download", class="download_style")),
                          
                        ),
                        fluidPage(
                          #textOutput("DE_MA_text"),
                          h4(id="MA"," MA Plot", style="text-align: center;"),
                          column(12, plotOutput(outputId="DEplot_average", height="450px"), align="center"),
                          div(style="margin-left: 375px;", downloadLink("MA_download", label="Download", class="download_style")),
                          
                        ),
                      ), 
                    ),
                    
                    tabPanel(title="Table",
                      dataTableOutput("DE_table"),
                      downloadLink("DE_table_download", label="Download", class="download_style"),
                    ),

                  ),

                  br(),

                  actionButton(inputId="assess_run_de", label="Assess DE Data"), 

                ),

              ),
            ),
            
            ######################################################################
            #                                                                    #
            #                           ASSESS DE TAB                            #
            #                                                                    #
            ######################################################################

            tabPanel(title="Assess DE", 
              
              ########################### OPTIONS ###########################

              dropdown(

                # network selection
                tags$h4("Network Selection"),
                uiOutput("select.folder"),

                # occr network selection
                radioButtons(
                  inputId="is_occr",
                  label="Use occr network?",
                  choices=c("Yes", "No"),
                  selected="No"
                ),

                # gene list selection
                radioButtons(
                  inputId="DE_data_selection",
                  label=tags$h4("DE Data Selection"),
                  choices=c("Use DE Results", "Upload DE Data"),
                ),

                conditionalPanel(
                  condition="input.DE_data_selection == 'Upload DE Data'", 
                  # upload file
                  fileInput(
                    inputId="DE_file", 
                    label="Choose DE Data File",
                    accept=c(".csv", ".tsv", ".txt")
                  ),
                  # div(style="margin-top: -25px"),
                  # select delimiter (default is nothing until file is selected and handled in server side)
                  radioButtons(
                    inputId='sepDEButton', 
                    label='Delimiter Selector', 
                    choices=c(Default=''), 
                    selected=''
                  ),
                ),
                
                # generate subnet button
                actionButton("generate_subnet_DE", "Generate Subnetwork",),
      
                # side panel characteristics
                style="jelly", icon="NETWORK OPTIONS",
                status="primary", width="335px", size="sm",

              ),

              div(
                style="display:inline-block; float:right",

                dropdown(
                  right=TRUE,
                  tags$h4("Download Options"),
                  inputId="ASSESS_DE_download_dropdown",
                  style="minimal", icon="DOWNLOAD OPTIONS",
                  status="primary", width="300px", size="sm",
                  
                  
                  selectInput(
                    inputId="ASSESS_DE_download_format",
                    label= tags$h6("Choose Plot Download Format"),
                    choices=c(".png", ".pdf"),
                    selected=".png",
                    width="600px",
                  ),
                  selectInput(
                    inputId="ASSESS_DE_download_table_format",
                    label= tags$h6("Choose Table Download Format"),
                    choices=c(".csv", ".tsv", ".txt"),
                    selected=".png",
                    width="600px",
                  ),

                ),
              ),
              
              br(),


              navlistPanel(
                id="assessDE_navList",
                widths=c(3, 9), well=FALSE, 
                
                ########################### VIEW FILES ###########################

                tabPanel(title="View Files",
                  tabsetPanel(
                    id="subnetwork_file_tabset_DE",
                    
                    tabPanel(
                      title="DE Data", 
                      conditionalPanel(condition="input.DE_data_selection == 'Use DE Results'", 
                      dataTableOutput("AssessDE_table"),
                      
                      ),
                      conditionalPanel(condition="input.DE_data_selection == 'Upload DE Data'",
                        dataTableOutput("UIDE_loaded_Content"),
                      ),
                      
                    ),

                    # subnetwork table of assess de
                    tabPanel(
                      title="Subnetwork", 
                      
                      downloadLink("DE_subnet_download", label="Download", class="download_style"),
                      tableOutput("subnetwork_DE")
                    ),
                  ),
                ),

                ########################### CLUSTER GENES ###########################
                
                tabPanel(title="Cluster Genes",
                  
                  mainPanel(
                    h3("Cluster Genes"),
                    
                    dropdown(
                      inputId="CG_dropdown_DE",
                      style="minimal", icon="OPTIONS",
                      status="primary", width="335px", size="sm",

                      awesomeCheckboxGroup(
                        inputId="clusterPlotOptions_upreg", 
                        label=tags$h4("Upregulated"),
                        choices=c("Network", "Heatmap", "Binarized Heatmap"),
                        status=""
                      ), 

                      awesomeCheckboxGroup(
                        inputId="clusterPlotOptions_downreg", 
                        label=tags$h4("Downregulated"), 
                        choices=c("Network", "Heatmap", "Binarized Heatmap"), 
                        status=""
                      ),

                      # run button
                      actionButton(inputId="runCGDE", label="Run",),
                    ),  
                    
                    # error message
                    p(id="CG_error_DE", "Please generate a subnetwork in NETWORK OPTIONS.", style="color:red"),

                    br(),
                    
                    tabsetPanel(

                      # plots tab
                      tabPanel(title="Plots",
                      
                        conditionalPanel(
                          condition="$.inArray('Network', input.clusterPlotOptions_upreg) > -1", 
                          h5(id="CGupreg_network_text", "Network of Clustered, Upregulated Genes"), 
                          br(), 
                          plotOutput(outputId="upregNetwork", height="500px"), 
                          div(style="margin-left: 400px;", downloadLink("CG_up_network_download", label="Download", class="download_style")),
                        ),

                        conditionalPanel(
                          condition="$.inArray('Heatmap', input.clusterPlotOptions_upreg) > -1", 
                          h5(id="CGupreg_heatmap_text","Heatmap of Clustered, Upregulated Genes"), 
                          br(), 
                          plotOutput(outputId="upregHeatmap", height="500px"),
                          div(style="margin-left: 400px;", downloadLink("CG_up_heatmap_download", label="Download", class="download_style")),
                        ),

                        conditionalPanel(
                          condition="$.inArray('Binarized Heatmap', input.clusterPlotOptions_upreg) > -1", 
                          h5(id="CGupreg_bheatmap_text","Binarized Hatmap of Clustered, Upregulated Genes"), 
                          br(),
                          plotOutput(outputId="upregbinHeatmap", height="500px"), 
                          div(style="margin-left: 400px;", downloadLink("CG_up_bheatmap_download", label="Download", class="download_style")),
                        ), 

                        conditionalPanel(
                          condition="$.inArray('Network', input.clusterPlotOptions_downreg) > -1",
                          h5(id="CGdownreg_network_text", "Network of Clustered, Downregulated Genes"), 
                          br(), 
                          plotOutput(outputId="downregNetwork", height="500px"),
                          div(style="margin-left: 400px;", downloadLink("CG_down_network_download", label="Download", class="download_style")),
                        ), 

                        conditionalPanel(
                          condition="$.inArray('Heatmap', input.clusterPlotOptions_downreg) > -1", 
                          h5(id="CGdownreg_heatmap_text", "Heatmap of Clustered, Downregulated Genes"), 
                          br(), 
                          plotOutput(outputId="downregHeatmap", height="500px"),
                          div(style="margin-left: 400px;", downloadLink("CG_down_heatmap_download", label="Download", class="download_style")),
                        ), 

                        conditionalPanel(
                          condition="$.inArray('Binarized Heatmap', input.clusterPlotOptions_downreg) > -1", 
                          h5(id="CGdownreg_bheatmap_text", "Binarized Heatmap of Clustered, Downregulated Genes"), 
                          br(), 
                          plotOutput(outputId="downregbinHeatmap", height="500px"),    
                          div(style="margin-left: 400px;", downloadLink("CG_down_bheatmap_download", label="Download", class="download_style")),
                        ),
                      ),

                      tabPanel(title="Tables",

                        # clustering genes
                        h5(id="CG_table_text_upreg", "Upregulated Clustering Genes"), 
                        br(),
                        fluidRow(
                          column(11,
                                  dataTableOutput("CG_table_upreg"),
                                  downloadLink("CG_table_upreg_download", label="Download", class="download_style"),
                          )
                        ),

                        h5(id="CG_table_text_downreg", "Downregulated Clustering Genes"), 
                        br(),
                        fluidRow(
                          column(11,
                                  dataTableOutput("CG_table_downreg"),
                                  downloadLink("CG_table_downreg_download", label="Download", class="download_style"),
                          )
                        ),
                        
                      ),
                      
                    ),
                  ),
                ),
                
                ########################### GENE CONNECTIVITY ###########################
                
                tabPanel(title="Gene Connectivity",

                  mainPanel(
                    h3("Gene Connectivity"),
                    
                    # options dropdown
                    dropdown(
                      inputId="GC_dropdown_DE",
                      style="minimal", icon="OPTIONS",
                      status="primary", width="335px", size="sm",

                      # select plots
                    
                      awesomeCheckboxGroup(
                        inputId="GCPlotOptions_upreg", 
                        label=tags$h4("Upregulated"),
                        choices=c("Density", "Histogram", "Clustered Density", "Clustered Histogram"), 
                        status=""
                      ), 

                      awesomeCheckboxGroup(
                        inputId="GCPlotOptions_downreg", 
                        label=tags$h4("Downregulated"),
                        choices= c("Density", "Histogram", "Clustered Density", "Clustered Histogram"),
                        status=""

                      ),
                      
                      # histogram breaks slider
                      conditionalPanel(
                        condition="$.inArray('Histogram', input.GCPlotOptions_upreg) > -1 || $.inArray('Clustered Histogram', input.GCPlotOptions_upreg) > -1 || $.inArray('Histogram', input.GCPlotOptions_downreg) > -1 || $.inArray('Clustered Histogram', input.GCPlotOptions_downreg) > -1" ,
                        sliderInput(
                          inputId="xybreaks_DE", 
                          label="Number of breaks for histogram:",
                          min=10, max=100, value=100, step=10,
                        ),
                      ),


                      # run button
                      actionButton(inputId="runGCDE", label="Run", ),

                    ),

                    # error message
                    p(id="GC_error_DE", "Please generate a subnetwork in NETWORK OPTIONS.", style="color:red"),

                    # density - upreg 
                    conditionalPanel(
                      condition="$.inArray('Density', input.GCPlotOptions_upreg) > -1", 
                      h5(id="GCdensityG_upreg_text", "Density Plot of Upregulated Gene Connectivity"), 
                      plotOutput(outputId="GCdensityGupreg", height="500px",),
                      div(style="margin-left: 400px;", downloadLink("GC_up_density_download", label="Download", class="download_style")),
                    ),

                    # histogram - upreg 
                    conditionalPanel(
                      condition="$.inArray('Histogram', input.GCPlotOptions_upreg) > -1", 
                      h5(id="GChistogramG_upreg_text", "Histogram of Upregulated Gene Connectivity"),
                      plotOutput(outputId="GChistogramGupreg", height="500px",),
                      div(style="margin-left: 400px;", downloadLink("GC_up_histogram_download", label="Download", class="download_style")),
                    ),

                    # density (subset by clusters) - upreg 
                    conditionalPanel(
                      condition="$.inArray('Clustered Density', input.GCPlotOptions_upreg) > -1", 
                      h5(id="GCdensitySubsetG_upreg_text", "Density plot of Upregulated Gene Connectivity subset by their clusters"), 
                      plotOutput(outputId="GCdensitySubsetGupreg", height="500px",),
                      div(style="margin-left: 400px;", downloadLink("GC_up_densitySubset_download", label="Download", class="download_style")),
                    ),

                    # histogram (subset by clusters) - upreg 
                    conditionalPanel(
                      condition="$.inArray('Clustered Histogram', input.GCPlotOptions_upreg) > -1", 
                      h5(id="GChistogramSubsetG_upreg_text", "Histogram of Upregulated Gene Connectivity subset by their clusters"), 
                      plotOutput(outputId="GChistogramSubsetGupreg", height="500px",),
                      div(style="margin-left: 400px;", downloadLink("GC_up_histSubset_download", label="Download", class="download_style")),
                    ),

                    # density - downreg 
                    conditionalPanel(
                      condition="$.inArray('Density', input.GCPlotOptions_downreg) > -1", 
                      h5(id="GCdensityG_downreg_text", "Density Plot of Downregulated Gene Connectivity"), 
                      plotOutput(outputId="GCdensityGdownreg", height="500px",),
                      div(style="margin-left: 400px;", downloadLink("GC_down_density_download", label="Download", class="download_style")),
                    ),

                    # histogram - downreg
                    conditionalPanel(
                      condition="$.inArray('Histogram', input.GCPlotOptions_downreg) > -1", 
                      h5(id="GChistogramG_downreg_text", "Histogram of Downregulated Gene Connectivity"),
                      plotOutput(outputId="GChistogramGdownreg", height="500px",),
                      div(style="margin-left: 400px;", downloadLink("GC_down_hist_download", label="Download", class="download_style")),
                    ),

                    # density (subset by clusters) - downreg
                    conditionalPanel(
                      condition="$.inArray('Clustered Density', input.GCPlotOptions_downreg) > -1", 
                      h5(id="GCdensitySubsetG_downreg_text", "Density plot of Downregulated Gene Connectivity subset by their clusters"), 
                      plotOutput(outputId="GCdensitySubsetGdownreg", height="500px",),
                      div(style="margin-left: 400px;", downloadLink("GC_down_densitySubset_download", label="Download", class="download_style")),
                    ),

                    # histogram (subset by clusters) - downreg
                    conditionalPanel(
                      condition="$.inArray('Clustered Histogram', input.GCPlotOptions_downreg) > -1", 
                      h5(id="GChistogramSubsetG_downreg_text", "Histogram of Dowregulated Gene Connectivity subset by their clusters"), 
                      plotOutput(outputId="GChistogramSubsetGdownreg", height="500px",),
                      div(style="margin-left: 400px;", downloadLink("GC_down_histSubset_download", label="Download", class="download_style")),
                    ),


                  )
                ),
                
                ########################### FUNCTIONAL OUTLIERS ###########################
                
                tabPanel(title="Functional Outliers",
                  
                  mainPanel(
                    h3("Functional Outliers"),

                    # options dropdown
                    dropdown(
                      inputId="FO_dropdown_DE",
                      style="minimal", icon="OPTIONS",
                      status="primary", width="335px", size="sm",

                      awesomeCheckboxGroup(
                        inputId="FO_upreg_options", 
                        label=tags$h4("Upregulated"),
                        choices=c("Network", "Heatmap", "Genes in Module", "Functional Outliers"), 
                        status=""
                      ),

                      # select tables
                      awesomeCheckboxGroup(
                          inputId="FO_downreg_options",
                          label=tags$h4("Downregulated"),
                          choices=c("Network", "Heatmap", "Genes in Module", "Functional Outliers"),
                          status=""
                      ),

                      # other options
                      tags$h4("Other"),
                      
                      # filt_min slider
                      sliderInput("filtmin_DE", label="Number of Genes to form Module",
                          min=0, max=20, value=6, step=1
                      ),

                      br(),
                      
                      # run button
                      actionButton(inputId="runFODE", label="Run", ),

                    ),
                    
                    # error message
                    p(id="FO_error_DE", "Please generate a subnetwork in NETWORK OPTIONS.", style="color:red"),

                  ),
                  
                  br(),
                  
                  tabsetPanel(

                    # plots tab
                    tabPanel(title="Plots",                      

                      conditionalPanel(
                        condition="$.inArray('Network', input.FO_upreg_options) > -1", 
                        h5(id="FOnetwork_upreg_text", "Upregulated Network"), 
                        plotOutput(outputId="FOnetwork_upreg", height="500px"), 
                        div(style="margin-left: 400px;", downloadLink("FO_up_network_download", label="Download", class="download_style"))
                      ), 

                      conditionalPanel(
                        condition="$.inArray('Heatmap', input.FO_upreg_options) > -1",
                        h5(id="FOheatmap_upreg_text", "Upregulated Heatmap"), 
                        plotOutput(outputId="FOheatmap_upreg", height="500px"),
                        div(style="margin-left: 400px;", downloadLink("FO_up_heatmap_download", label="Download", class="download_style"))
                      ), 

                      conditionalPanel(
                        condition="$.inArray('Network', input.FO_downreg_options) > -1", 
                        h5(id="FOnetwork_downreg_text", "Downregulated Network"), 
                        plotOutput(outputId="FOnetwork_downreg", height="500px"),
                        div(style="margin-left: 400px;", downloadLink("FO_down_network_download", label="Download", class="download_style"))
                      ), 

                      conditionalPanel(
                        condition="$.inArray('Heatmap', input.FO_downreg_options) > -1", 
                        h5(id="FOheatmap_downreg_text", "Downregulated Heatmap"), 
                        plotOutput(outputId="FOheatmap_downreg", height="500px"),
                        div(style="margin-left: 400px;", downloadLink("FO_down_heatmap_download", label="Download", class="download_style"))
                      ),

                    ),

                    # tables tab
                    tabPanel(title="Tables", 

                      # selected genes table output
                      conditionalPanel(
                        condition="$.inArray('Genes in Module', input.FO_upreg_options) > -1", 
                        h5(id="genes_not_keep_table_text_upreg", "Upregulated Genes in Module"), 
                        br(),
                        fluidRow(
                          column(11,
                                  dataTableOutput("genes_not_keep_table_upreg"),
                                  downloadLink("genes_not_keep_table_upreg_download", label="Download", class="download_style"),
                          )
                        ),
                      ),

                      br(),

                      # # unselected genes table output
                      conditionalPanel(
                        condition="$.inArray('Functional Outliers', input.FO_upreg_options) > -1", 
                        h5(id="genes_keep_table_text_upreg", "Upregulated Functional Outliers"), 
                        br(),
                        fluidRow(
                          column( 11,
                                  dataTableOutput("genes_keep_table_upreg"),
                                  downloadLink("genes_keep_table_upreg_download", label="Download", class="download_style"),
                          )
                        ),
                      ),

                      br(),

                      conditionalPanel(
                        condition="$.inArray('Genes in Module', input.FO_downreg_options) > -1", 
                        h5(id="genes_not_keep_table_text_downreg", "Downregulated Genes in Module"), 
                        br(),
                        fluidRow(
                          column( 11,
                                  dataTableOutput("genes_not_keep_table_downreg"),
                                  downloadLink("genes_not_keep_table_downreg_download", label="Download", class="download_style"),
                          )
                        ),
                      ),

                      br(),

                      # # unselected genes table output
                      conditionalPanel(
                        condition="$.inArray('Functional Outliers', input.FO_downreg_options) > -1", 
                        h5(id="genes_keep_table_text_downreg", "Downregulated Functional Outliers"), 
                        br(),
                        fluidRow(
                          column( 11,
                                  dataTableOutput("genes_keep_table_downreg"),
                                  downloadLink("genes_keep_table_downreg_download", label="Download", class="download_style"),

                          )
                        ),
                      ),

                    )
                  ),
                ),

                ########################### GSEA ###########################

                tabPanel(title="Gene Set Enrichment Analysis",
                  mainPanel(
                    h3("Gene Set Enrichment Analysis"),
                    
                    # options dropdown
                    dropdown(
                      inputId="DE_GSEA_dropdown",
                      style="minimal", icon="OPTIONS",
                      status="primary", width="335px", size="sm",
                      
                      # select GSEA type
                      awesomeCheckboxGroup(
                        inputId="GSEA_type",
                        label=tags$h4("GSEA Type"),
                        choices=c("Standard GSEA", "AUCs GSEA"),
                        selected="",
                        status="",
                      ),
                      
                      # standard GSEA options
                      conditionalPanel(
                        condition="input.GSEA_type.includes('Standard GSEA')", 
                        awesomeCheckboxGroup(
                          inputId="GSEA_std_PlotOptions",
                          label=tags$h4("Standard GSEA"), 
                          choices=c("Upregulated P-value Heatmap", "Downregulated P-value Heatmap"),
                          status=""
                        ),
                      ),
                      
                      br(),

                      # run button
                      actionButton(
                        inputId="DE_GSEA_run",
                        label="Run",
                        style="color: #fff; background-color: #3E3F3A; border-color: #20201F"
                      ),

                    ),
                    

                    # error message
                    p(id="DE_GSEA_error", "Please generate a subnetwork in NETWORK OPTIONS.", style="color:red"),

                  ),
                  br(),
                  
                  tabsetPanel(

                    # Standard GSEA tab
                    tabPanel(
                      title="Standard",
                      mainPanel(

                        # upregulated heatmap
                        conditionalPanel(
                          condition="$.inArray('Upregulated P-value Heatmap', input.GSEA_std_PlotOptions) > -1", 
                          h5(id="GSEA_up_heatmap_text", "Upregulated P-value Heatmap"), 
                          br(),
                          plotOutput(outputId="GSEA_up_heatmap", height="500px"),
                          div(style="margin-left: 400px;", downloadLink("GSEA_up_heatmap_download", label="Download", class="download_style")),
                        ),

                        # downregulated heatmap
                        conditionalPanel(
                          condition="$.inArray('Downregulated P-value Heatmap', input.GSEA_std_PlotOptions) > -1", 
                          h5(id="GSEA_down_heatmap_text", "Downregulated P-value Heatmap"), 
                          br(),
                          plotOutput(outputId="GSEA_down_heatmap", height="500px"),
                          div(style="margin-left: 400px;", downloadLink("GSEA_down_heatmap_download", label="Download", class="download_style")),
                        ),

                      )
                    ),

                    # AUCs GSEA tab
                    tabPanel(
                      title="AUC",
                      mainPanel(
                        # AUROC graph
                        conditionalPanel(
                          condition="$.inArray('AUCs GSEA', input.GSEA_type) > -1", 
                          plotOutput(outputId="GSEA_auc", height="500px"),
                          div(style="margin-left: 400px;", downloadLink("GSEAauc_download", label="Download", class="download_style")),
                        ),
                        br(),

                      )
                    ),
                  )

                ),

              ),
            ),

            ######################################################################
            #                                                                    #
            #                       ASSESS GENE LIST TAB                         #
            #                                                                    #
            ######################################################################
             
            tabPanel(title="Assess Gene List",

              ########################### OPTIONS ###########################

              dropdown(

                # network selection
                tags$h4("Network Selection"),
                uiOutput("select.folder_gene_list"),

                # occr network selection
                radioButtons(
                  inputId="is_occr_gene_list",
                  label="Use occr network?",
                  choices=c("Yes", "No"),
                  selected="No"
                ),

                # gene list selection
                radioButtons(
                  inputId="gene_list_selection",
                  label=tags$h4("Gene List Selection"),
                  choices=c("Upload Gene List", "Generate Gene List"),
                  selected=""
                ),

                # generate gene list
                conditionalPanel(
                  condition="input.gene_list_selection == 'Generate Gene List'", 
                  textInput(
                    inputId='chooseChrome', 
                    label='Choose Chromosome' , 
                    placeholder="chrX"
                  ),
                  textInput(
                    inputId='chooseGeneNo', 
                    label='Choose Number of Genes',
                    placeholder="100"
                  ),
                ),

                # gene list upload
                conditionalPanel(
                  condition="input.gene_list_selection == 'Upload Gene List'", 
                  # upload file
                  fileInput(
                    inputId="DEFile", 
                    label="Choose Gene List File",
                    accept=c(".csv", ".tsv", ".txt")
                  ),

                  radioButtons( # select gene list type
                    inputId='GL_gene_list_type', 
                    label='Gene List Type', 
                    choices=c("Gene Names", "Entrez Id"), 
                    selected=''
                  ),
                ),
                
                # generate subnet button
                actionButton("generate_subnet", "Generate Subnetwork",),
      
                # side panel characteristics
                style="jelly", icon="NETWORK OPTIONS",
                status="primary", width="335px", size="sm",

              ),

              div(style="display:inline-block; float:right", dropdown(right=TRUE,
                  tags$h4("Download Options"),
                  inputId="ASSESS_GENE_LIST_download_dropdown",
                  style="minimal", icon="DOWNLOAD OPTIONS",
                  status="primary", width="300px", size="sm",
                  
                  
                  selectInput(
                    inputId="ASSESS_GENE_LIST_download_format",
                    label= tags$h6("Choose Plot Download Format"),
                    choices=c(".png", ".pdf"),
                    selected=".png",
                    width="600px",
                  ),
                  selectInput(
                    inputId="ASSESS_GENE_LIST_download_table_format",
                    label= tags$h6("Choose Table Download Format"),
                    choices=c(".csv", ".tsv", ".txt"),
                    selected=".png",
                    width="600px",
                  ),
                  
                  # run button
                  
                  
                ), 
              ),


              br(),

              navlistPanel(
                id="assessGL_navList",
                widths=c(3, 9), well=FALSE,

                ########################### VIEW FILES ###########################

                tabPanel(title="View Files",

                  tabsetPanel(
                    id="subnetwork_file_tabset",
                    
                    # view file tab
                    tabPanel(title="File",
                      uiOutput("UIDEContent"),
                    ),

                    # view subnetwork tab
                    tabPanel(title="Subnetwork", 
                      # subnetwork table of assess gene list
                      downloadLink("GL_subnet_download", label="Download", class="download_style"),
                      tableOutput("subnetwork")
                      

                    ),
                  ),
                ),

                ########################### CLUSTER GENES ###########################

                tabPanel(title="Cluster Genes",

                  mainPanel(
                    h3("Cluster Genes"),
                    dropdown(
                      inputId="CG_dropdown",
                      style="minimal", icon="OPTIONS",
                      status="primary", width="335px", size="sm",

                      # select plots
                      awesomeCheckboxGroup(
                        inputId="clusterPlotOptions_genelist",
                        label=tags$h4("Select Plots"), 
                        choices=c("Network", "Heatmap", "Binarized Heatmap"),
                        status=""
                      ),

                      # run button
                      actionButton(inputId="run", label="Run",),

                    
                    ),  

                    # error message
                    p(id="CG_error", "Please generate a subnetwork in NETWORK OPTIONS.", style="color:red"),

                    br(),

                    tabsetPanel(

                      
                      # plots tab
                      tabPanel(
                        br(),
                        title="Plots",
                        br(),

                        # network
                        conditionalPanel(
                          condition="$.inArray('Network', input.clusterPlotOptions_genelist) > -1", 
                          h5(id="CG_network_text", "Network of Clustered Genes"), 
                          br(),
                          plotOutput(outputId="network", height="500px"),

                          div(style="margin-left: 400px;", downloadLink("CG_network_download", label="Download", class="download_style")),
                          
                        ),
                        
                        # heatmap
                        conditionalPanel(
                          condition="$.inArray('Heatmap', input.clusterPlotOptions_genelist) > -1", 
                          h5(id="CG_heatmap_text", "Heatmap of Clustered Genes"),
                          br(),
                          plotOutput(outputId="heatmap", height="500px"),
                          div(style="margin-left: 400px;", downloadLink("CG_heatmap_download", label="Download", class="download_style")),
                        ),

                        # binarized heatmap
                        conditionalPanel(
                          condition="$.inArray('Binarized Heatmap', input.clusterPlotOptions_genelist) > -1", 
                          h5(id="CG_bheatmap_text", "Binarized Heatmap of Clustered Genes"), 
                          br(),
                          plotOutput(outputId="Bheatmap", height="500px"), 
                          div(style="margin-left: 400px;", downloadLink("CG_bheatmap_download", label="Download", class="download_style")),
                        ),


                      
                      ),
                      
                      # tables tab
                      tabPanel(
                        br(),
                        title="Tables",
                        br(),

                        # clustering genes
                        h5(id="CG_table_text", "Clustering Genes"), 
                        br(),
                        fluidRow(
                          column(11,
                                  dataTableOutput("CG_table"),
                                  downloadLink("table_genelist_download", label="Download", class="download_style"),
                          )
                        ),

                      ),
                    ), 
                  ), 


                ), 

                ########################### GENE CONNECTIVITY ###########################

                tabPanel(title="Gene Connectivity", 
                    mainPanel(
                    h3("Gene Connectivity"),
                    # options dropdown
                    dropdown(
                      inputId="GC_dropdown",
                      style="minimal", icon="OPTIONS",
                      status="primary", width="335px", size="sm",

                      # select plots
                      awesomeCheckboxGroup(
                        inputId="GCPlotOptions_genelist",
                        label=tags$h4("Select Plots"), 
                        choices=c("Density", "Histogram", "Clustered Density", "Clustered Histogram"),
                        status=""
                      ),
                      
                      # breaks slider
                      conditionalPanel(
                        condition="$.inArray('Histogram', input.GCPlotOptions_genelist) > -1 || $.inArray('Clustered Histogram', input.GCPlotOptions_genelist) > -1" ,
                        tags$h4("Other"),
                        sliderInput(
                          inputId="xybreaks", 
                          label="Number of breaks for histogram:",
                          min=10, max=100, value=100, step=10,
                        ),
                      ),
                      
                      br(),

                      # run button
                      actionButton(inputId="runGC", label="Run", ),
                    ),


                    # error message
                    p(id="GC_error", "Please generate a subnetwork in NETWORK OPTIONS.", style="color:red"),

                    # density
                    conditionalPanel(
                      br(),
                      condition="$.inArray('Density', input.GCPlotOptions_genelist) > -1", 
                      h5(id="GCdensityG_text", "Density Plot of Gene Connectivity"), 
                      br(),
                      plotOutput(outputId="GCdensityG", height="500px",),
                      div(style="margin-left: 400px;", downloadLink("GC_density_download", label="Download", class="download_style")),
                      br(),
                    ),

                    # histogram
                    conditionalPanel(
                      br(),
                      condition="$.inArray('Histogram', input.GCPlotOptions_genelist) > -1", 
                      h5(id="GChistogramG_text", "Histogram of Gene Connectivity"),
                      br(),
                      plotOutput(outputId="GChistogramG", height="500px",),
                      div(style="margin-left: 400px;", downloadLink("GC_histogram_download", label="Download", class="download_style")),
                      br(),
                    ),

                    # density (subset by clusters)
                    conditionalPanel(
                      br(),
                      condition="$.inArray('Clustered Density', input.GCPlotOptions_genelist) > -1", 
                      h5(id="GCdensitySubsetG_text", "Density plot of Gene Connectivity subset by their clusters"), 
                      br(),
                      plotOutput(outputId="GCdensitySubsetG", height="500px",),
                      div(style="margin-left: 400px;", downloadLink("GC_densitySubset_download", label="Download", class="download_style")),
                      br(),
                    ),

                    # histogram (subset by clusters)
                    conditionalPanel(
                      br(),
                      condition="$.inArray('Clustered Histogram', input.GCPlotOptions_genelist) > -1", 
                      h5(id="GChistogramSubsetG_text", "Histogram of Gene Connectivity subset by their clusters"), 
                      br(),
                      plotOutput(outputId="GChistogramSubsetG", height="500px",),
                      div(style="margin-left: 400px;", downloadLink("GC_histogramSubset_download", label="Download", class="download_style")),
                      br(),
                    ),

                   ),
                ), 

                ########################### FUNCTIONAL OUTLIERS ###########################

                tabPanel(title="Functional Outliers",
                  mainPanel(
                    h3("Functional Outliers"),

                    # options dropdown
                    dropdown(
                      inputId="FO_dropdown",
                      style="minimal", icon="OPTIONS",
                      status="primary", width="335px", size="sm",

                      # select plots
                      awesomeCheckboxGroup(
                        inputId="FOPlotOptions_genelist",
                        label=tags$h4("Select Plots"),
                        choices=c("Network", "Heatmap"),
                        status=""
                      ),

                      # select tables
                      awesomeCheckboxGroup(
                          inputId="FO_table_options",
                          label=tags$h4("Select Tables"),
                          choices=c("Functional Outliers", "Genes in Module"),
                          status=""
                      ),

                      # other options
                      tags$h4("Other"),
                      
                      # filt_min slider
                      sliderInput("filtmin", label="Number of Genes to form Module",
                          min=0, max=20, value=6, step=1
                      ),

                      br(),
                      
                      # run button
                      actionButton(inputId="runFO", label="Run", ),

                    ),
                    
                    # error message
                    p(id="FO_error", "Please generate a subnetwork in NETWORK OPTIONS.", style="color:red"),

                  ), 

                  br(), 

                  
                  tabsetPanel(

                    # plots tab
                    tabPanel(title="Plots",
               
                      # heatmap
                      conditionalPanel(
                        condition="$.inArray('Network', input.FOPlotOptions_genelist) > -1", 
                        h5(id="FO_network_text", "Network"), 
                        plotOutput(outputId="FO_network", height="500px"),
                        div(style="margin-left: 400px;", downloadLink("FO_network_download", label="Download", class="download_style")),
                      ),

                      # network
                      conditionalPanel(
                        condition="$.inArray('Heatmap', input.FOPlotOptions_genelist) > -1", 
                        h5(id="FO_heatmap_text", "Heatmap"), 
                        plotOutput(outputId="FO_heatmap", height="500px"),
                        div(style="margin-left: 400px;", downloadLink("FO_heatmap_download", label="Download", class="download_style")),
                      ),

                      

                    ),

                    # tables tab
                    tabPanel(title="Tables",

                      # selected genes table output
                      conditionalPanel(
                        condition="$.inArray('Genes in Module', input.FO_table_options) > -1", 
                        h5(id="genes_not_keep_table_text", "Genes in Module"), 
                        br(),
                        fluidRow(
                          column( 11,
                                  dataTableOutput("genes_not_keep_table"),
                                  downloadLink("genes_not_keep_table_download", label="Download", class="download_style"),
                          )
                        ),
                      ),

                      br(),

                      # unselected genes table output
                      conditionalPanel(
                        condition="$.inArray('Functional Outliers', input.FO_table_options) > -1", 
                        h5(id="genes_keep_table_text", "Outliers"), 
                        br(),
                        fluidRow(
                          column( 11,
                                  dataTableOutput("genes_keep_table"),
                                  downloadLink("genes_keep_table_table_download", label="Download", class="download_style"),
                          )
                        ),
                      ),
                    )
                  ),
                ),

                ########################### GSEA ###########################

                tabPanel(title="Gene Set Enrichment Analysis",

                  mainPanel(
                    h3("Gene Set Enrichment Analysis"),
                    
                    # options dropdown
                    dropdown(
                      inputId="GL_GSEA_dropdown",

                      # dropdown characteristics
                      style="minimal", icon="OPTIONS",
                      status="primary", width="335px", size="sm",
                      
                      # standard GSEA options
                      awesomeCheckboxGroup(
                        inputId="GL_GSEA_options_std",
                        label=tags$h4("Standard GSEA"), 
                        choices=c("P-value Heatmap", "Compare GSEA to heatmap"),
                        status=""
                      ),
                      
                      # run button
                      actionButton(inputId="GL_GSEA_run", label="Run"),

                    ),

                    # error message
                    p(id="GL_GSEA_error", "Please generate a subnetwork in NETWORK OPTIONS.", style="color:red"),
                    br(),
                  
                    # heatmap
                    conditionalPanel(
                      condition="$.inArray('P-value Heatmap', input.GL_GSEA_options_std) > -1", 
                      h5(id="GL_GSEA_heatmap_text", "P-value Heatmap"), 
                      plotOutput(outputId="GL_GSEA_heatmap_plot", height="500px"),
                      br(),
                    ),

                    conditionalPanel(
                      condition = "$.inArray('Compare GSEA to heatmap', input.GL_GSEA_options_std) > -1",
                      splitLayout(cellWidths = c("50%", "50%"),
                        fluidPage(
                          h5(id = "GL_GSEA_compare_heatmap_text", "P-value Heatmap"), 
                          plotOutput(outputId = "GL_GSEA_heatmap_compare", height = "500px")
                        ),

                        fluidPage(
                          h5(id = "GL_GSEA_compare_heatmap_text", "Heatmap of Clustered Genes"), 
                          plotOutput(outputId = "GL_GSEA_heatmap_cluster_compare", height = "500px"),
                        )
                      )

                  ),
                ),
                )

              ),
            ),
            

            
  )
)
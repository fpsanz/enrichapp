library(shinydashboard)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(EnsDb.Mmusculus.v79)
library(limma)
library(tidyverse)
library(DT)
library(purrr)
library(plotly)
library(chorddiag)
library(DESeq2)
library(fgsea)
source("utils.R")
options(shiny.maxRequestSize = 3000*1024^2)

 ### HEADER ############
header <- dashboardHeader(title = "RNAseq viewer and report App", 
                  titleWidth = 300, 
                  dropdownMenuOutput("messageMenu")
                  )
### SIDEBAR ##########
sidebar <- dashboardSidebar(sidebarMenu(
                            menuItem(
                            fileInput("deseqFile",
                                      "Choose RDS with Deseq object",
                                      placeholder = "RDS DEseq"))),
                            sidebarMenu(
                              menuItem("Preview dataset",
                                       tabName = "preview",
                                       icon = icon("eye")),
                              uiOutput("sampleGroup"),
                              menuItem(
                                "Kegg Enrichment",
                                tabName = "kegg",
                                icon = icon("chart-bar")
                              ),
                              menuItem(
                                "GO Enrichment",
                                tabName = "go",
                                icon = icon("chart-bar")
                              ),
                              menuItem(
                                "GSEA",
                                tabName = "gsea",
                                icon = icon("chart-line")
                              ),
                              sidebarMenu(
                              menuItem(
                              textInput("author", value="your name...", label = h4("Author report name") )),
                              downloadButton("report", "Generate report")
                              )
                            ))
### BODY ###############
body <- dashboardBody(
    shiny::tagList(shiny::tags$head(
        shiny::tags$link(rel = "stylesheet", type = "text/css", href = "busystyle.css"),
        shiny::tags$script(type = "text/javascript", src = "busy.js")
    )),
    div(
        class = "busy",
        p("Loading data and computing enrichment, please be patient..."),
        img(src = "dna.gif")
    ),
    tabItems(
        # preview tab
        tabItem(tabName = "preview",
                h3("Sample info (colData)"),
                fluidRow(
                    column(
                        width = 8,
                        offset = 2,
                        dataTableOutput("samples")
                    )
                ),
                hr(),
                h3("DE results"),
                fluidRow(
                    column(
                        width = 8,
                        offset = 2,
                        dataTableOutput("preview")
                    )
                ),
                hr(),
                fluidRow(column(
                        width = 8,
                        offset = 2,
                        plotOutput("pca", height = "600px")
                    ))
                ),
        # kegg tab content
        tabItem(tabName = "kegg",
                tabsetPanel(
                  tabPanel(
                    "Upregulated",
                    fluidRow(column(
                      width = 8,
                      offset = 2,
                      dataTableOutput("table")
                    )),
                    hr(),
                    fluidRow(
                      class = "text-center",
                      column(
                        align = "center",
                        offset = 2,
                        plotlyOutput("keggPlot"),
                        width = 5
                      ),
                      column(
                        align = "center",
                        offset = 0,
                        chorddiagOutput("keggChord", width = "500px", height = "500px"),
                        width = 4
                      )
                    ),
                    fluidRow(
                      class="text-center",
                      column(
                        align = "center",
                        offset= 2,
                        plotOutput("keggDotUp"),
                        width = 8
                      )
                    ),
                    hr(),
                    fluidRow(
                      class="text-center",
                      column(
                        align = "center",
                        offset= 2,
                        plotOutput("heatmapKeggUp"),
                        width = 8
                      )
                    ),
                    hr(),
                    fluidRow(
                      class="text-center",
                      column(
                        align = "center",
                        offset= 2,
                        plotOutput("cnetKeggUp"),
                        width = 8
                      )
                    )
                  ), # fin tabpanel upregulated
                  tabPanel(
                    "Downregulated",
                    fluidRow(column(
                      width = 8,
                      offset = 2,
                      dataTableOutput("tableDown")
                    )),
                    hr(),
                    fluidRow(
                      class = "text-center",
                      column(
                        align = "center",
                        offset = 2,
                        plotlyOutput("keggPlotDown"),
                        width = 5
                      ),
                      column(
                        align = "center",
                        offset = 0,
                        chorddiagOutput("keggChordDown", width = "500px", height = "500px"),
                        width = 4
                      )
                    ),
                    fluidRow(
                      class="text-center",
                      column(
                        align = "center",
                        offset= 2,
                        plotOutput("keggDotDown"),
                        width = 8
                      )
                    ),
                    hr(),
                    fluidRow(
                      class="text-center",
                      column(
                        align = "center",
                        offset= 2,
                        plotOutput("heatmapKeggDown"),
                        width = 8
                      )
                    ),
                    hr(),
                    fluidRow(
                      class="text-center",
                      column(
                        align = "center",
                        offset= 2,
                        plotOutput("cnetKeggDown"),
                        width = 8
                      )
                    )
                  ) # fin tabpanel
                )),
        tabItem( tabName = "go", # GO tab GO tab
          tabsetPanel(
            tabPanel("Upregulated",
            h3("Biological proccess"),
            fluidRow(column(
                # table BP
                width = 8,
                offset = 2,
                dataTableOutput("tableBP")
            )),
            hr(),
            fluidRow(
                #plot BP
                class = "text-center",
                column(
                    align = "center",
                    offset = 2,
                    plotlyOutput("plotBP"),
                    width = 8)),
            hr(),
            fluidRow(
              class="text-center",
              column(
                align = "center",
                offset= 2,
                plotOutput("BPDotUp"),
                width = 8)),
            h3("Molecular Functions"),
            fluidRow(column(
                # table MF
                width = 8,
                offset = 2,
                dataTableOutput("tableMF")
            )),
            hr(),
            fluidRow(# plot MF
                class = "text-center",
                column( align = "center",offset = 2,plotlyOutput("plotMF"),width = 8)
            ),
            hr(),
            fluidRow(
              class="text-center",
              column(
                align = "center",
                offset= 2,
                plotOutput("MFDotUp"),
                width = 8
              )
            ),
            h3("Cellular components"),
            fluidRow( # table CC
              column(width = 8,offset = 2,dataTableOutput("tableCC"))
              ),
            hr(),
            fluidRow(# plot CC
                class = "text-center",
                column(align = "center",offset = 2,plotlyOutput("plotCC"),width = 8)
            ),
            hr(),
            fluidRow(
              class="text-center",
              column(
                align = "center",
                offset= 2,
                plotOutput("CCDotUp"),
                width = 8
              )
            ) #fin fluidrow
        ),
        tabPanel("Downregulated",
                 h3("Biological proccess"),
                 fluidRow(column(# table BP
                   width = 8,
                   offset = 2,
                   dataTableOutput("tableBPdown")
                 )),
                 hr(),
                 fluidRow(#plot BP
                   class = "text-center",
                   column(
                     align = "center",
                     offset = 2,
                     plotlyOutput("plotBPdown"),
                     width = 8
                   )
                 ),
                 hr(),
                 fluidRow(
                   class="text-center",
                   column(
                     align = "center",
                     offset= 2,
                     plotOutput("BPDotDown"),
                     width = 8)),
                 h3("Molecular Functions"),
                 fluidRow(column(# table MF
                   width = 8,
                   offset = 2,
                   dataTableOutput("tableMFdown")
                 )),
                 hr(),
                 fluidRow(# plot MF
                   class = "text-center",
                   column( align = "center",offset = 2,plotlyOutput("plotMFdown"),width = 8)
                 ),
                 hr(),
                 fluidRow(
                   class="text-center",
                   column(
                     align = "center",
                     offset= 2,
                     plotOutput("MFDotDown"),
                     width = 8)),
                 h3("Cellular components"),
                 fluidRow( # table CC
                   column(width = 8,offset = 2,dataTableOutput("tableCCdown"))
                 ),
                 hr(),
                 fluidRow(# plot CC
                   class = "text-center",
                   column(align = "center",offset = 2,plotlyOutput("plotCCdown"),width = 8)
                 ),
                 hr(),
                 fluidRow(
                   class="text-center",
                   column(
                     align = "center",
                     offset= 2,
                     plotOutput("CCDotDown"),
                     width = 8)) #fin fluidrow
        ))), # tab GO
        tabItem(tabName = "gsea",
                h3("Gene Set Enrichment Analysis"),
                fluidRow(column(
                  width = 8,
                  offset = 2,
                  dataTableOutput("gseaTable")
                )),
                hr(),
                fluidRow(column(
                  width = 8,
                  offset = 2,
                  plotOutput("gseaPlot")
                ))
                ) #fin tab GSEA
    ) # fin tab items
)# fin dashboardbody

########################################## UI #################################################

ui <- dashboardPage(title="Rnaseq viewer and report",
                    header,
                    sidebar,
                    body
                    ) # fin del UI

########################################## SERVER #################################################
server <- function(input, output) {
    data <- reactiveValues()
    goDT <- reactiveValues()
    kgg <- reactiveValues()
    go <- reactiveValues()
    kggDT <- reactiveValues()
    datos <- reactiveValues(dds=NULL)
    gsea <- reactiveValues()
    
    observeEvent(input$deseqFile, {
        datos$dds <- readRDS(input$deseqFile$datapath)
        saveRDS(datos$dds, "tmpResources/deseq.Rds")
        data$genesUp <- getSigUpregulated(datos$dds)
        data$genesDown <- getSigDownregulated(datos$dds)
        saveRDS(data$genesUp, "tmpResources/genesUp.Rds")
        saveRDS(data$genesDown, "tmpResources/genesDown.Rds")
        kgg$up <- customKegg(data$genesUp, species = "Mm", species.KEGG = "mmu")
        saveRDS(kgg$up, "tmpResources/kggUp.Rds")
        kggDT$up <- kegg2DT(kgg$up, data$genesUp)
        saveRDS(kggDT$up, "tmpResources/kggDTup.Rds")
        go$up <- customGO(data$genesUp, species = "Mm")
        saveRDS(go$up, "tmpResources/goUp.Rds")
        goDT$up <- go2DT(enrichdf = go$up, data = data$genesUp )
        saveRDS(goDT$up, "tmpResources/goDTup.Rds")
        kgg$down <- customKegg(data$genesDown, species = "Mm", species.KEGG = "mmu")
        saveRDS(kgg$down, "tmpResources/kggDown.Rds")
        kggDT$down <- kegg2DT(kgg$down, data$genesDown)
        saveRDS(kggDT$down, "tmpResources/kggDTdown.Rds")
        go$down <- customGO(data$genesDown, species = "Mm")
        saveRDS(go$down, "tmpResources/goDown.Rds")
        goDT$down <- go2DT(enrichdf = go$down, data = data$genesDown )
        saveRDS(goDT$down, "tmpResources/goDTdown.Rds")
        
    })
  # generate reactive variable ###################
    rows <- reactive({input$table_rows_selected})
    bprows <- reactive({input$tableBP_rows_selected})
    mfrows <- reactive({input$tableMF_rows_selected})
    ccrows <- reactive({input$tableCC_rows_selected})
    rowsdown <- reactive({input$tableDown_rows_selected})
    bprowsdown <- reactive({input$tableBPdown_rows_selected})
    mfrowsdown <- reactive({input$tableMFdown_rows_selected})
    ccrowsdown <- reactive({input$tableCCdown_rows_selected})
    variables <- reactive({input$variables})
    gsearow <- reactive({input$gseaTable_rows_selected})
  # ui selector sample groups ###################
    output$sampleGroup <- renderUI({
        validate(need(datos$dds, ""))
        nvars <- colData(datos$dds) %>% 
                    as.data.frame() %>% 
                    select_if(is.factor) %>%
                    names()
        selectInput("variables", label="Select condition[s] to plot PCA",
                    choices = nvars,
                    multiple = TRUE)
    })
  # preview table ###################
    output$preview <- DT::renderDataTable(server=TRUE,{
        validate(need(datos$dds, "Load file to render table"))
        res <- results(datos$dds)
        res <- as.data.frame(res)
        datatable( round(res,4), 
                  filter = list(position="top", clear=FALSE),
                  options = list(
                  columnDefs = list(list(orderable = FALSE,
                                         className = "details-control",
                                         targets = 1),
                                    list(className = "dt-right", targets = 1:ncol(res))
                                    ),
                  dom = "Bfrtipl",
                  buttons = c("copy", "csv", "excel", "pdf", "print"),
                  list(pageLength = 10, white_space = "normal")
                  )
                  )
    })
  # preview samples ###################
    output$samples <- DT::renderDataTable(server = TRUE,{
        validate(need(datos$dds, "Load file to render table"))
        metadata <- as.data.frame(colData(datos$dds))
        metadata$sizeFactor <- round(metadata$sizeFactor,4)
        datatable( metadata, 
                  filter = list(position="top", clear=FALSE),
                  options = list(
                  columnDefs = list(list(orderable = FALSE,
                                         className = "details-control",
                                         targets = 1),
                                    list(className = "dt-right", targets = 1:ncol(metadata))
                                    ),
                  dom = "Bfrtipl",
                  buttons = c("copy", "csv", "excel", "pdf", "print"),
                  list(pageLength = 10, white_space = "normal")
                  )
                  )
    })
# view pca plot data ###################
    output$pca <- renderPlot( {
        validate(need(datos$dds, "Load file to render PCA"))
        validate(need(variables(),"" ) )
        plotPCA(rlog(datos$dds), intgroup = variables() )+
            theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))+
            theme(text = element_text(size=20))
        })
# KEGG table up#####################################
    output$table <- DT::renderDataTable(server=TRUE,{
        validate(need(kgg$up, "Load file to render table"))
        kgg <- kgg$up
        #genes <- data$genes
        predata <- kggDT$up
        #predata <- kegg2DT(kgg, genes)
        datatable2(
            predata,
            vars = c("genes"),
            filter = list(position="top", clear=FALSE),
            escape = FALSE,
            opts = list(pageLength = 10, white_space = "normal"))
        }) 
# KEGG barplot up################
    output$keggPlot <- renderPlotly ({
        validate(need(kgg$up, "Load file to render BarPlot"))
        kgg <- kgg$up
        nr <- rows()
        if(is.null(nr)){nr <- c(1:10)}
        plotKegg(enrichdf = kgg[nr,], nrows = length(nr))
    })
# KEGG chordiag plot up ###############
    output$keggChord <- renderChorddiag({
        validate(need(kgg$up, "Load file to render ChordPlot"))
        kgg <- kgg$up
        nr <- rows()
        if(is.null(nr)){nr <- c(1:10)}
        chordPlot(kgg[nr, ], nRows = length(nr), orderby = "P.DE")
    })
# KEGG dotplot UP ################### 
    output$keggDotUp <- renderPlot({
      validate(need(kgg$up, "Load file to render dotPlot"))
      kgg <- kgg$up
      nr <- rows()
      if(is.null(nr)){nr <- c(1:20)}
      dotPlotkegg(kgg[nr,], n = length(nr))
    })
# KEGG heatmap Up #################
    output$heatmapKeggUp <- renderPlot({
      validate(need(rows(), ""))
      validate(need(kggDT$up, ""))
      nr <- rows()
      heatmapKegg(kggDT$up, nr)
    })
# KEGG cnet Up #################
    output$cnetKeggUp <- renderPlot({
      validate(need(rows(), ""))
      validate(need(kgg$up, ""))
      nr <- rows()
      customCnetKegg(kgg$up, nr)
    })
    
# KEGG table down #####################################
    output$tableDown <- DT::renderDataTable(server=TRUE,{
      validate(need(kgg$down, "Load file to render table"))
      kgg <- kgg$down
      predata <- kggDT$down
      datatable2(
        predata,
        vars = c("genes"),
        filter = list(position="top", clear=FALSE),
        escape = FALSE,
        opts = list(pageLength = 10, white_space = "normal"))
    }) 
# KEGG barplot down ################
    output$keggPlotDown <- renderPlotly ({
      validate(need(kgg$down, "Load file to render BarPlot"))
      kgg <- kgg$down
      nrdown <- rowsdown()
      if(is.null(nrdown)){nrdown <- c(1:10)}
      plotKegg(enrichdf = kgg[nrdown,], nrows = length(nrdown))
    })
# KEGG chordiag plot down ###############
    output$keggChordDown <- renderChorddiag({
      validate(need(kgg$down, "Load file to render ChordPlot"))
      kgg <- kgg$down
      nrdown <- rowsdown()
      if(is.null(nrdown)){nrdown <- c(1:10)}
      chordPlot(kgg[nrdown, ], nRows = length(nrdown), orderby = "P.DE")
    })
# KEEGG dotplot Down ################### 
    output$keggDotDown <- renderPlot({
      validate(need(kgg$down, "Load file to render dotPlot"))
      kgg <- kgg$down
      nrdown <- rowsdown()
      if(is.null(nrdown)){nrdown <- c(1:20)}
      dotPlotkegg(kgg[nrdown,], n = length(nrdown))
    })
# KEGG heatmap Down #################
    output$heatmapKeggDown <- renderPlot({
      validate(need(rowsdown(), ""))
      validate(need(kggDT$down, ""))
      nrdown <- rowsdown()
      heatmapKegg(kggDT$down, nrdown)
    })
# KEGG cnet Up #################
    output$cnetKeggDown <- renderPlot({
      validate(need(rows(), ""))
      validate(need(kgg$up, ""))
      nrdown <- rowsdown()
      customCnetKegg(kgg$down, nrdown)
    })
# GO table BP UP#####################
    output$tableBP <- DT::renderDataTable(server=TRUE,{
    validate(need(goDT$up, "Load file to render table"))
    goDT <- goDT$up
    datatable2(goDT[goDT$Ont=="BP",], vars = c("genes"),
               filter = list(position="top", clear=FALSE),
               escape = FALSE,
               opts = list(pageLength = 10, white_space = "normal")
               )
    })
# GO plots BP UP #####################
    output$plotBP <- renderPlotly({
      validate(need(go$up, "Load file to render plot"))
       gos <- go$up
        bpnr <- bprows()
        if(is.null(bpnr)){bpnr <- c(1:10)}
        gosBP <- gos[gos$Ont=="BP",]
        plotGO(enrichdf = gosBP[bpnr, ], nrows = length(bpnr), ont="BP")
    })
# GO BP dotplot up ################### 
    output$BPDotUp <- renderPlot({
      validate(need(go$up, "Load file to render dotPlot"))
      gos <- go$up
      
      bpnr <- bprows()
      if(is.null(bpnr)){bpnr <- c(1:20)}
      gosBP <- gos[gos$Ont=="BP",]
      dotPlotGO(gosBP[bpnr,], n = length(bpnr))
    })
# GO table MF UP #####################
    output$tableMF <- DT::renderDataTable({
      validate(need(goDT$up, "Load file to render table"))
      goDT <- goDT$up
      datatable2(goDT[goDT$Ont=="MF",], vars = c("genes"),
                 filter = list(position="top", clear=FALSE),
                 escape = FALSE,
                 opts = list(pageLength = 10, white_space = "normal",
                             ajax = list(serverSide = TRUE, processing = TRUE))
      )
    })
    
    
    
# GO plots MF UP #####################
    output$plotMF <- renderPlotly({
        validate(need(go$up, "Load file to render plot"))
        gos <- go$up
        mfnr <- mfrows()
        if(is.null(mfnr)){mfnr <- c(1:10)}
        gosMF <- gos[gos$Ont=="MF",]
        plotGO(enrichdf = gosMF[mfnr, ], nrows = length(mfnr), ont = "MF")
    })
# GO MF dotplot up ################### 
    output$MFDotUp <- renderPlot({
      validate(need(go$up, "Load file to render dotPlot"))
      gos <- go$up
      mfnr <- mfrows()
      if(is.null(mfnr)){mfnr <- c(1:20)}
      gosMF <- gos[gos$Ont=="MF",]
      dotPlotGO(gosMF[mfnr,], n = length(mfnr))
    })
# GO table CC UP #####################
    output$tableCC <- DT::renderDataTable(server=TRUE,{
      validate(need(goDT$up, "Load file to render table"))
      goDT <- goDT$up
      datatable2(goDT[goDT$Ont=="CC",], vars = c("genes"),
                 filter = list(position="top", clear=FALSE),
                 escape = FALSE,
                 opts = list(pageLength = 10, white_space = "normal",
                             ajax = list(serverSide = TRUE, processing = TRUE))
      )
    })
# GO plots CC UP #####################
    output$plotCC <- renderPlotly({
        validate(need(go$up, "Load file to render plot"))
        gos <- go$up
        ccnr <- ccrows()
        if(is.null(ccnr)){ccnr <- c(1:10)}
        gosCC <- gos[gos$Ont=="CC",]
      plotGO(enrichdf = gosCC[ccnr,], nrows = length(ccnr), ont="CC")
    })
# GO CC dotplot up ################### 
    output$CCDotUp <- renderPlot({
      validate(need(go$up, "Load file to render dotPlot"))
      gos <- go$up
      ccnr <- ccrows()
      if(is.null(ccnr)){ccnr <- c(1:20)}
      gosCC <- gos[gos$Ont=="CC",]
      dotPlotGO(gosCC[ccnr,], n = length(ccnr))
    })
# GO table BP DOWN #####################
    output$tableBPdown <- DT::renderDataTable(server=TRUE,{
      validate(need(goDT$down, "Load file to render table"))
      goDT <- goDT$down
      datatable2(goDT[goDT$Ont=="BP",], vars = c("genes"),
                 filter = list(position="top", clear=FALSE),
                 escape = FALSE,
                 opts = list(pageLength = 10, white_space = "normal")
      )
    })
    # GO plots BP DOWN #####################
    output$plotBPdown <- renderPlotly({
      validate(need(go$down, "Load file to render plot"))
      gos <- go$down
      bpnrdown <- bprowsdown()
      if(is.null(bpnrdown)){bpnrdown <- c(1:10)}
      gosBP <- gos[gos$Ont=="BP",]
      plotGO(enrichdf = gosBP[bpnrdown, ], nrows = length(bpnrdown), ont="BP")
    })
# GO BP dotplot down ################### 
    output$BPDotDown <- renderPlot({
      validate(need(go$down, "Load file to render dotPlot"))
      gos <- go$down
      bpnrdown <- bprowsdown()
      if(is.null(bpnrdown)){bpnrdown <- c(1:20)}
      gosBP <- gos[gos$Ont=="BP",]
      dotPlotGO(gosBP[bpnrdown,], n = length(bpnrdown))
    })
# GO table MF DOWN #####################
    output$tableMFdown <- DT::renderDataTable({
      validate(need(goDT$down, "Load file to render table"))
      goDT <- goDT$down
      datatable2(goDT[goDT$Ont=="MF",], vars = c("genes"),
                 filter = list(position="top", clear=FALSE),
                 escape = FALSE,
                 opts = list(pageLength = 10, white_space = "normal",
                             ajax = list(serverSide = TRUE, processing = TRUE))
      )
    })
# GO plots MF DOWN #####################
    output$plotMFdown <- renderPlotly({
      validate(need(go$down, "Load file to render plot"))
      gos <- go$down
      mfnrdown <- mfrowsdown()
      if(is.null(mfnrdown)){mfnrdown <- c(1:10)}
      gosMF <- gos[gos$Ont=="MF",]
      plotGO(enrichdf = gosMF[mfnrdown, ], nrows = length(mfnrdown), ont = "MF")
    })
# GO MF dotplot down ################### 
    output$MFDotDown <- renderPlot({
      validate(need(go$down, "Load file to render dotPlot"))
      gos <- go$down
      mfnrdown <- mfrowsdown()
      if(is.null(mfnrdown)){mfnrdown <- c(1:20)}
      gosMF <- gos[gos$Ont=="MF",]
      dotPlotGO(gosMF[mfnrdown,], n = length(mfnrdown))
    })
    # GO table CC DOWN #####################
    output$tableCCdown <- DT::renderDataTable(server=TRUE,{
      validate(need(goDT$down, "Load file to render table"))
      goDT <- goDT$down
      datatable2(goDT[goDT$Ont=="CC",], vars = c("genes"),
                 filter = list(position="top", clear=FALSE),
                 escape = FALSE,
                 opts = list(pageLength = 10, white_space = "normal",
                             ajax = list(serverSide = TRUE, processing = TRUE))
      )
    })
    # GO plots CC DOWN #####################
    output$plotCCdown <- renderPlotly({
      validate(need(go$down, "Load file to render plot"))
      gos <- go$down
      ccnrdown <- ccrowsdown()
      if(is.null(ccnrdown)){ccnrdown <- c(1:10)}
      gosCC <- gos[gos$Ont=="CC",]
      plotGO(enrichdf = gosCC[ccnrdown,], nrows = length(ccnrdown), ont="CC")
    })
# GO CC dotplot down ################### 
    output$CCDotDown <- renderPlot({
      validate(need(go$down, "Load file to render dotPlot"))
      gos <- go$down
      ccnrdown <- ccrowsdown()
      if(is.null(ccnrdown)){ccnrdown <- c(1:20)}
      gosCC <- gos[gos$Ont=="CC",]
      dotPlotGO(gosCC[ccnrdown,], n = length(ccnrdown))
    })
# GSEA table ##########################
    output$gseaTable <- renderDataTable({
      validate(need(datos$dds, "Load file to render table"))
      gsea$gsea <- gseaKegg(datos$dds)
      mygsea <- gsea$gsea
      saveRDS(mygsea, "tmpResources/gsea.Rds")
      table <- mygsea@result[mygsea@result$p.adjust<=0.05 ,2:9] %>% 
        mutate_at(vars(3:7), ~round(., 3))
      DT::datatable( table,
                 rownames=FALSE,
                 filter = list(position="top", clear=FALSE),
                 options = list(
                   columnDefs = list(list(orderable = FALSE,
                                          className = "details-control",
                                          targets = 1)
                   ),
                   dom = "Bfrtipl",
                   buttons = c("copy", "csv", "excel", "pdf", "print"),
                   list(pageLength = 10, white_space = "normal")
                 )
      )
    })
# GSEA plot ##########################
    output$gseaPlot <- renderPlot({
      validate(need(gsea$gsea, "Load file to render table"))
      gseanr <- gsearow()
      if(is.null(gseanr)){gseanr <- c(1)}
        enrichplot::gseaplot2(gsea$gsea, geneSetID = gseanr, pvalue_table = TRUE, ES_geom = "line")
    })
# author name ######################
    author <- reactive({input$author})
# generate report #############################
    output$report <- downloadHandler(
        # For PDF output, change this to "report.pdf"
        filename = "report.html",
        content = function(file) {
            tempReport <- file.path(tempdir(), "report.Rmd")
            file.copy("report.Rmd", tempReport, overwrite = TRUE)
            file.copy("utils.R", file.path(tempdir(),"utils.R"), overwrite = TRUE)
            file.copy("tmpResources/", tempdir(), overwrite = TRUE, recursive = TRUE)
            do.call(file.remove, list(list.files("tmpResources/", full.names = TRUE)))
            
            # file.copy("genesUp.Rds", file.path(tempdir(), "genesUp.Rds"), overwrite = TRUE)
            # file.remove("genesUp.Rds")
            # file.copy("genesDown.Rds", file.path(tempdir(), "genesDown.Rds"), overwrite = TRUE)
            # file.remove("genesDown.Rds")
            # file.copy("kggUp.Rds", file.path(tempdir(), "kggUp.Rds"), overwrite = TRUE)
            # file.remove("kggUp.Rds")
            # file.copy("kggDown.Rds", file.path(tempdir(), "kggDown.Rds"), overwrite = TRUE)
            # file.remove("kggDown.Rds")
            # file.copy("goDown.Rds", file.path(tempdir(), "goDown.Rds"), overwrite = TRUE)
            # file.remove("goDown.Rds")
            # file.copy("goUp.Rds", file.path(tempdir(), "goUp.Rds"), overwrite = TRUE)
            # file.remove("goUp.Rds")
            # file.copy("kggDTup.Rds", file.path(tempdir(), "kggDTup.Rds"), overwrite = TRUE)
            # file.remove("kggDTup.Rds")
            # file.copy("kggDTdown.Rds", file.path(tempdir(), "kggDTdown.Rds"), overwrite = TRUE)
            # file.remove("kggDTdown.Rds")
            # file.copy("goDTup.Rds", file.path(tempdir(), "goDTup.Rds"), overwrite = TRUE)
            # file.remove("goDTup.Rds")
            # file.copy("goDTdown.Rds", file.path(tempdir(), "goDTdown.Rds"), overwrite = TRUE)
            # file.remove("goDTdown.Rds")
            # file.copy("deseq.Rds", file.path(tempdir(), "deseq.Rds"), overwrite = TRUE)
            # file.remove("deseq.Rds")
            # file.copy("gsea.Rds", file.path(tempdir(), "gsea.Rds"), overwrite = TRUE)
            # file.remove("gsea.Rds")
            

            nr <- rows()
            nrdown <- rowsdown()
            bpnr <- bprows()
            mfnr <- mfrows()
            ccnr <- ccrows()
            bpnrdown <- bprowsdown()
            mfnrdown <- mfrowsdown()
            ccnrdown <- ccrowsdown()
            variablepca <- variables()
            gseanr <- gsearow()
            if(is.null(gseanr)){gseanr <- c(1)}
            if(is.null(nr)){nr <- c(1:10)}
            if(is.null(ccnr)){ccnr <- c(1:10)}
            if(is.null(mfnr)){mfnr <- c(1:10)}
            if(is.null(bpnr)){bpnr <- c(1:10)}
            if(is.null(nrdown)){nrdown <- c(1:10)}
            if(is.null(ccnrdown)){ccnrdown <- c(1:10)}
            if(is.null(mfnrdown)){mfnrdown <- c(1:10)}
            if(is.null(bpnrdown)){bpnrdown <- c(1:10)}
            if(is.null(variablepca)){variablepca=NULL}
            params <- list(nr=nr, nrdown=nrdown, bpnr=bpnr, bpnrdown=bpnrdown,
                           mfnr=mfnr, mfnrdown=mfnrdown, ccnr=ccnr, ccnrdown=ccnrdown,
                           variablepca=variablepca, tempdir =tempdir(), gseanr=gseanr, author=author() )
            rmarkdown::render(
                tempReport,
                output_file = file,
                params = params,
                envir = new.env(parent = globalenv())
            )
        } )
}


shinyApp(ui, server)


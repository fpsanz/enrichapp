library(shinydashboard)
library(limma)
library(tidyverse)
library(DT)
library(purrr)
library(plotly)
library(chorddiag)
library(DESeq2)
source("utils.R")


#genes <- read.csv("genes.csv", header = T, stringsAsFactors = F, row.names = 1)
#genes <- genes %>% filter(!is.na(ENTREZID)& !ENTREZID=="NULL" )
#save(genes, file = "genes.Rdata")
# load("genes.Rdata")
# gos <- readRDS("gos.Rds")
# goDT <- go2DT(enrichdf = gos, data = genes)
#kgg <- customKegg(genes, species = "Mm", species.KEGG = "mmu")
#saveRDS(kgg, "kgg.Rds")
#kgg <- readRDS("kgg.Rds")
#gos <- customGO(genes,species = "Mm")
#saveRDS(gos,"gos.Rds")

  
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
                                      placeholder = "RDS file"))),
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
        tabItem(
            tabName = "kegg",
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
            )
        ),
        tabItem(
            # GO tab GO tab
            tabName = "go",
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
                    width = 8
                )
            ),
            h3("Molecular Functions"),
            fluidRow(column(
                # table MF
                width = 8,
                offset = 2,
                dataTableOutput("tableMF")
            )),
            hr(),
            fluidRow(
                # plot MF
                class = "text-center",
                column(
                    align = "center",
                    offset = 2,
                    plotlyOutput("plotMF"),
                    width = 8
                )
            ),
            h3("Cellular components"),
            fluidRow(column(
                # table CC
                width = 8,
                offset = 2,
                dataTableOutput("tableCC")
            )),
            hr(),
            fluidRow(
                # plot CC
                class = "text-center",
                column(
                    align = "center",
                    offset = 2,
                    plotlyOutput("plotCC"),
                    width = 8
                )
            ) #fin fluidrow
        ) # tab GO
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
    data <- reactiveValues(genes=NULL)
    goDT <- reactiveValues(dt=NULL)
    kgg <- reactiveValues(k=NULL)
    go <- reactiveValues(g=NULL)
    kggDT <- reactiveValues(predata=NULL)
    datos <- reactiveValues(dds=NULL)
    
    observeEvent(input$deseqFile, {
        datos$dds <- readRDS(input$deseqFile$datapath)
        #source("aux.R", local = TRUE)
        
        # data$genes <- loadGenes(input$deseqFile$datapath)
        # go$g <- customGO(data$genes, species = "Mm")
        # kgg$k <- customKegg(data$genes, species = "Mm", species.KEGG = "mmu")
        # goDT$dt <- go2DT(enrichdf = go$g, data = data$genes )
        # kggDT$predata <- kegg2DT(kgg$k, data$genes)
    })
  # generate reactive variable ###################
    rows <- reactive({input$table_rows_selected})
    bprows <- reactive({input$tableBP_rows_selected})
    mfrows <- reactive({input$tableMF_rows_selected})
    ccrows <- reactive({input$tableCC_rows_selected})
    variables <- reactive({input$variables})
  # ui selector sample groups ###################
    output$sampleGroup <- renderUI({
        validate(need(datos$dds, ""))
        nvars <- colData(datos$dds) %>% 
                    as.data.frame() %>% 
                    select_if(is.factor) %>%
                    names()
        selectInput("variables", label="Select condition[s] to plot",
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
                  dom = "B",
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
  # view kegg table #####################################
    output$table <- DT::renderDataTable(server=TRUE,{
        validate(need(data$genes, "Load file to render table"))
        kgg <- kgg$k
        #genes <- data$genes
        predata <- kggDT$predata
        #predata <- kegg2DT(kgg, genes)
        datatable2(
            predata,
            vars = c("genes"),
            filter = list(position="top", clear=FALSE),
            escape = FALSE,
            opts = list(pageLength = 10, white_space = "normal"))
    }) 
# KEGG barplot ################
    output$keggPlot <- renderPlotly ({
        validate(need(data$genes, "Load file to render BarPlot"))
        kgg <- kgg$k
        nr <- rows()
        if(is.null(nr)){nr <- c(1:10)}
        plotKegg(enrichdf = kgg[nr,], nrows = length(nr))
    })
# KEGG chordiag plot
    output$keggChord <- renderChorddiag({
        validate(need(data$genes, "Load file to render ChordPlot"))
        kgg <- kgg$k
        nr <- rows()
        if(is.null(nr)){nr <- c(1:10)}
        chordPlot(kgg[nr, ], nRows = length(nr), orderby = "P.DE")
    })
# GO table BP #####################
    output$tableBP <- DT::renderDataTable(server=TRUE,{
    validate(need(goDT$dt, "Load file to render table"))
    goDT <- goDT$dt
    datatable2(goDT[goDT$Ont=="BP",], vars = c("genes"),
               filter = list(position="top", clear=FALSE),
               escape = FALSE,
               opts = list(pageLength = 10, white_space = "normal")
               )
    })
# GO plots BP #####################
    output$plotBP <- renderPlotly({
      validate(need(go$g, "Load file to render plot"))
       gos <- go$g
        bpnr <- bprows()
        if(is.null(bpnr)){bpnr <- c(1:10)}
        gosBP <- gos[gos$Ont=="BP",]
        plotGO(enrichdf = gosBP[bpnr, ], nrows = length(bpnr), ont="BP")
    })
# GO table MF #####################
    output$tableMF <- DT::renderDataTable({
      validate(need(goDT$dt, "Load file to render table"))
      goDT <- goDT$dt
      datatable2(goDT[goDT$Ont=="MF",], vars = c("genes"),
                 filter = list(position="top", clear=FALSE),
                 escape = FALSE,
                 opts = list(pageLength = 10, white_space = "normal",
                             ajax = list(serverSide = TRUE, processing = TRUE))
      )
    })
# view go plots MF #####################
    output$plotMF <- renderPlotly({
        validate(need(go$g, "Load file to render plot"))
        gos <- go$g
        mfnr <- mfrows()
        if(is.null(mfnr)){mfnr <- c(1:10)}
        gosMF <- gos[gos$Ont=="MF",]
        plotGO(enrichdf = gosMF[mfnr, ], nrows = length(mfnr), ont = "MF")
    })
# view Go table CC #####################
    output$tableCC <- DT::renderDataTable(server=TRUE,{
      validate(need(goDT$dt, "Load file to render table"))
      goDT <- goDT$dt
      datatable2(goDT[goDT$Ont=="CC",], vars = c("genes"),
                 filter = list(position="top", clear=FALSE),
                 escape = FALSE,
                 opts = list(pageLength = 10, white_space = "normal",
                             ajax = list(serverSide = TRUE, processing = TRUE))
      )
    })
# view go plots CC #####################
    output$plotCC <- renderPlotly({
        validate(need(go$g, "Load file to render plot"))
        gos <- go$g
        ccnr <- ccrows()
        if(is.null(ccnr)){ccnr <- c(1:10)}
        gosCC <- gos[gos$Ont=="CC",]
      plotGO(enrichdf = gosCC[ccnr,], nrows = length(ccnr), ont="CC")
    })
    
        
# generate report #############################
    output$report <- downloadHandler(
        # For PDF output, change this to "report.pdf"
        filename = "report.html",
        content = function(file) {
            tempReport <- file.path(tempdir(), "report.Rmd")
            file.copy("report.Rmd", tempReport, overwrite = TRUE)
            nr <- rows()
            if(is.null(nr)){nr <- c(1:10)}
            ccnr <- ccrows()
            if(is.null(ccnr)){ccnr <- c(1:10)}
            mfnr <- mfrows()
            if(is.null(mfnr)){mfnr <- c(1:10)}
            bpnr <- bprows()
            if(is.null(bpnr)){bpnr <- c(1:10)}
            params <- list(n = nr, bp = bpnr, mf=mfnr, cc=ccnr)
            rmarkdown::render(
                tempReport,
                output_file = file,
                params = params,
                envir = new.env(parent = globalenv())
            )
        } )
}


shinyApp(ui, server)


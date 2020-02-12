library(shinydashboard)
library(limma)
library(tidyverse)
library(DT)
library(purrr)
library(plotly)
library(chorddiag)
source("utils.R")
source("global.R")
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
sidebar <- dashboardSidebar(fileInput("deseqFile",
                                      "Choose RDS with Deseq object",
                                      placeholder = "RDS file"),
                            sidebarMenu(
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
                              downloadButton("report", "Generate report")
                            ))
### BODY ###############
body <- dashboardBody(
  tabItems(
    # First tab content
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
    tabItem( # second tab GO tab
        tabName = "go",
        h3("Biological proccess"),
        fluidRow(column( # table BP
            width = 8,
            offset = 2,
            dataTableOutput("tableBP")
        )),
        hr(),
        fluidRow( #plot BP
            class = "text-center",
            column(
                align = "center",
                offset = 2,
                plotlyOutput("plotBP"),
                width = 8
            )
        ),
        h3("Molecular Functions"),
        fluidRow(column( # table MF
            width = 8,
            offset = 2,
            dataTableOutput("tableMF")
        )),
        hr(),
        fluidRow( # plot MF
            class = "text-center",
            column(
                align = "center",
                offset = 2,
                plotlyOutput("plotMF"),
                width = 8
            )),
        h3("Cellular components"),
        fluidRow(column( # table CC
            width = 8,
            offset = 2,
            dataTableOutput("tableCC")
        )),
        hr(),
        fluidRow( # plot CC
            class = "text-center",
            column(
                align = "center",
                offset = 2,
                plotlyOutput("plotCC"),
                width = 8
            ))
    )
))

########################################## UI #################################################

ui <- dashboardPage(title="Rnaseq viewer and report",
                    header,
                    sidebar,
                    body
                    ) # fin del UI

########################################## SERVER #################################################
server <- function(input, output) {
    data <- reactiveValues(genes=NULL)
    dat <- reactiveValues(kgg=NULL)
    file <- reactiveValues(f=NULL)
    observeEvent(input$deseqFile, {
        file$f <- input$deseqFile$datapath
    })
    kk <- observe({
        fileg <- file$f
        if(!is.null(fileg)){
            filegenes <- file()
            source("aux.R", local=FALSE)
            dat$kgg <- auxkgg
            data$genes <-auxgenes
        }
    })
    # generate reactive variable ###################
    rows <- reactive({input$table_rows_selected})
    bprows <- reactive({input$tableBP_rows_selected})
    mfrows <- reactive({input$tableMF_rows_selected})
    ccrows <- reactive({input$tableCC_rows_selected})
# view kegg table #####################################
    output$table <- DT::renderDataTable(server=TRUE,{
        validate(need(data$genes, "Load file to render table"))
        kgg <- dat$kgg
        genes <- data$genes
        predata <- kegg2DT(kgg, genes)
        datatable2(
            predata,
            vars = c("genes"),
            filter = list(position="top", clear=FALSE),
            escape = FALSE,
            opts = list(pageLength = 10, white_space = "normal"))
    }) 
# view kegg plot ################
    output$keggPlot <- renderPlotly ({
        validate(need(dat$kgg, "Load file to render BarPlot"))
        kgg <- dat$kgg
        nr <- rows()
        if(is.null(nr)){nr <- c(1:10)}
        plotKegg(enrichdf = kgg[nr,], nrows = length(nr))
    })
# generate view kegg chordiag plot
    output$keggChord <- renderChorddiag({
        validate(need(dat$kgg, "Load file to render ChordPlot"))
        kgg <- dat$kgg
        nr <- rows()
        if(is.null(nr)){nr <- c(1:10)}
        chordPlot(kgg[nr, ], nRows = length(nr), orderby = "P.DE")
    })
# view Go table BP #####################
    output$tableBP <- DT::renderDataTable(server=TRUE,{
    datatable2(goDT[goDT$Ont=="BP",], vars = c("genes"),
               filter = list(position="top", clear=FALSE),
               escape = FALSE,
               opts = list(pageLength = 10, white_space = "normal")
               )
    })
# view go plots BP #####################
    output$plotBP <- renderPlotly({
        bpnr <- bprows()
        if(is.null(bpnr)){bpnr <- c(1:10)}
        gosBP <- gos[gos$Ont=="BP",]
        plotGO(enrichdf = gosBP[bpnr, ], nrows = length(bpnr), ont="BP")
    })
# view Go table MF #####################
    output$tableMF <- DT::renderDataTable({
      datatable2(goDT[goDT$Ont=="MF",], vars = c("genes"),
                 filter = list(position="top", clear=FALSE),
                 escape = FALSE,
                 opts = list(pageLength = 10, white_space = "normal",
                             ajax = list(serverSide = TRUE, processing = TRUE))
      )
    })
# view go plots MF #####################
    output$plotMF <- renderPlotly({
        mfnr <- mfrows()
        if(is.null(mfnr)){mfnr <- c(1:10)}
        gosMF <- gos[gos$Ont=="MF",]
        plotGO(enrichdf = gosMF[mfnr, ], nrows = length(mfnr), ont = "MF")
    })
# view Go table CC #####################
    output$tableCC <- DT::renderDataTable(server=TRUE,{
      datatable2(goDT[goDT$Ont=="CC",], vars = c("genes"),
                 filter = list(position="top", clear=FALSE),
                 escape = FALSE,
                 opts = list(pageLength = 10, white_space = "normal",
                             ajax = list(serverSide = TRUE, processing = TRUE))
      )
    })
# view go plots CC #####################
    output$plotCC <- renderPlotly({
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


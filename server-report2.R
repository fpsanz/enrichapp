
  
output$report2 <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "report.html",
    content = function(file) {
      tempReport <- file.path(tempdir(), "flexReport.Rmd")
      file.copy("flexReport.Rmd", tempReport, overwrite = TRUE)
      file.copy("utils.R", file.path(tempdir(),"utils.R"), overwrite = TRUE)
      file.copy("tmpResources/", tempdir(), overwrite = TRUE, recursive = TRUE)
      file.copy("resources/dna-svg-small-13.gif",
                file.path(tempdir(), "tmpResources/dna-svg-small-13.gif"), overwrite = TRUE)
      #do.call(file.remove, list(list.files("tmpResources/", full.names = TRUE)))
      nrall <- rowsAll()
      nrup <- rowsUp()
      nrdown <- rowsdown()
      bpnrup <- bprowsup()
      mfnrup <- mfrowsup()
      ccnrup <- ccrowsup()
      bpnrdown <- bprowsdown()
      mfnrdown <- mfrowsdown()
      ccnrdown <- ccrowsdown()
      variablepca <- variables()
      gseanr <- gsearow()
      bpnrall <- bprowsall()
      mfnrall <- mfrowsall()
      ccnrall <- ccrowsall()
      if(is.null(gseanr)){gseanr <- c(1)}
      if(is.null(nrup)){ nrup <- ( if( dim(kggDT$up)[1]<10) 1:dim(kggDT$up)[1] else c(1:10) ) } 
      if(is.null(nrall)){ nrall <-  ( if( dim(kggDT$all)[1]<10) 1:dim(kggDT$all)[1] else c(1:10) ) }
      if(is.null(nrdown)){ nrdown <- ( if( dim(kggDT$down)[1]<10) 1:dim(kggDT$down)[1] else c(1:10) ) }
      
      if(is.null(ccnrup)){ ccnrup <- ( if (dim(goDT$up[goDT$up$Ont=="CC", ])[1]<10)
                                             1:dim(goDT$up[goDT$up$Ont=="CC", ])[1] else c(1:10) ) }
      if(is.null(mfnrup)){ mfnrup <- ( if (dim(goDT$up[goDT$up$Ont=="MF", ])[1]<10)
                                             1:dim(goDT$up[goDT$up$Ont=="MF", ])[1] else c(1:10) ) }
      if(is.null(bpnrup)){ bpnrup <- ( if (dim(goDT$up[goDT$up$Ont=="BP", ])[1]<10)
                                             1:dim(goDT$up[goDT$up$Ont=="BP", ])[1] else c(1:10) )}
      
      if(is.null(ccnrdown)){ccnrdown <- ( if (dim(goDT$down[goDT$down$Ont=="CC", ])[1]<10)
                                              1:dim(goDT$down[goDT$down$Ont=="CC", ])[1] else c(1:10) ) }
      if(is.null(mfnrdown)){mfnrdown <- ( if (dim(goDT$down[goDT$down$Ont=="MF", ])[1]<10)
                                              1:dim(goDT$down[goDT$down$Ont=="MF", ])[1] else c(1:10) ) }
      if(is.null(bpnrdown)){bpnrdown <- ( if (dim(goDT$down[goDT$down$Ont=="BP", ])[1]<10)
                                              1:dim(goDT$down[goDT$down$Ont=="BP", ])[1] else c(1:10) )}
      
      if(is.null(bpnrall)){bpnrall <- ( if (dim(goDT$all[goDT$all$Ont=="BP", ])[1]<10)
                                            1:dim(goDT$all[goDT$all$Ont=="BP", ])[1] else c(1:10)) }
      if(is.null(mfnrall)){mfnrall <- ( if (dim(goDT$all[goDT$all$Ont=="MF", ])[1]<10)
                                            1:dim(goDT$all[goDT$all$Ont=="MF", ])[1] else c(1:10)) }
      if(is.null(ccnrall)){ccnrall <- ( if (dim(goDT$all[goDT$all$Ont=="CC", ])[1]<10)
                                            1:dim(goDT$all[goDT$all$Ont=="CC", ])[1] else c(1:10))}
      
      if(is.null(variablepca)){variablepca=NULL}
      params <- list(nrup=nrup, nrdown=nrdown, bpnrup=bpnrup, bpnrdown=bpnrdown,
                     mfnrup=mfnrup, mfnrdown=mfnrdown, ccnrup=ccnrup, ccnrdown=ccnrdown,
                     variablepca=variablepca, tempdir =tempdir(),
                     gseanr=gseanr, author=author(), nrall = nrall,
                     bpnrall=bpnrall, mfnrall=mfnrall, ccnrall=ccnrall,
                     explainPreview=explainPreview(), biologicalText=biologicalText(),
                     keggAllText = keggAllText(), deseq = datos$dds, 
                     kggAll = kgg$all, kggUp = kgg$up, kggDown = kgg$down,
                     kggDTall = kggDT$all, kggDTup = kggDT$up, kggDTdown = kggDT$down,
                     goAll = go$all, goDTall = goDT$all, goUp=go$up, goDTup = goDT$up,
                     goDown = go$down, goDTdown = goDT$down, gsea = gsea$gsea)
      rmarkdown::render(
        tempReport,
        output_file = file,
        params = params,
        envir = new.env(parent = globalenv( ))
      )
    } )
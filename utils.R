# Chorplot ##########################################
chordPlot <- function(enrichdf, nRows = 10, ont=NULL,  orderby=NULL) {
  if(! "dplyr" %in% .packages()) require("dplyr")
  if(! "tidyr" %in% .packages()) require("tidyr")
  if(! "chordiag" %in% .packages()) require("chorddiag")

  name <- match.arg(names(enrichdf)[1], c("Pathway", "Term"))

  if(dim(enrichdf)[2]==7){
    if(!is.null(ont)){
      enrichdf = enrichdf[enrichdf$Ont==ont, ]
    } else{
      stop("Ontology must be provided if enrichdf is customGO object")
    }
  }

  if(!is.null(orderby)){
    orderby = match.arg(orderby, c("DE", "P.DE", "N", name))
    if(orderby=="P.DE" | orderby == name){
      enrichdf <- enrichdf %>% arrange(get(orderby))
    } else{ enrichdf <- enrichdf %>% arrange(desc(get(orderby)))}
  }

  kgg <- enrichdf
  n <- nRows

  tmp <- kgg[1:n, c(name, "genes") ]
  tmp2 <- tmp %>% separate_rows(genes, convert = TRUE)
  kk <- tidyr::pivot_wider(tmp2, names_from = name, values_from = genes)

  ns <- NULL
  nd <- NULL
  for (i in seq(1, dim(kk)[2])) {
    pt <- unlist(kk[1, i][[1]])
    for (j in seq(1, dim(kk)[2])) {
      ns <-
        append(ns, length(intersect(unlist(kk[1, i][[1]]), unlist(kk[1, j][[1]]))))
      if (i != j) {
        pt <- pt[!(pt %in% unlist(kk[1, j][[1]]))]
      }
    }
    nd <- append(nd, length(pt))
  }
  m <- matrix(ns, nrow = dim(kk)[2])
  diag(m) <- nd
  dimnames(m) <- list(have = names(kk), prefer = names(kk))

  p <- chorddiag(
    m,
    groupnamePadding = 20,
    margin = 100,
    groupColors = colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(n),
    showGroupnames = FALSE,
    tooltipFontsize = 12,
    showTicks = FALSE
  )
  return(p)
}

# plot para exportar cnet a cytoscape ###############
customCnet2Cytoscape <- function(kgg, category=NULL, nPath=NULL, byDE=FALSE){
    if(! "ggraph" %in% .packages()) require("ggraph")
    if(! "igraph" %in% .packages()) require("igraph")
    if(! "dplyr" %in% .packages()) require("dplyr")
    if(! "tidyr" %in% .packages()) require("tidyr")
    ## check if cytoscape is running
    a <- tryCatch( expr= { cytoscapePing() },
                   error = function(e){print("Error!! Cytoscape must be running")}
                   )
    ## define palette
    color_palette <- function(colors) colorRampPalette(colors)(n = 299)
    if(byDE){
        kgg <- kgg %>% arrange(-DE)
    }
    ## check parameters
    if(!is.null(category)){
        if(is.numeric(category)){
            tmp <- kgg[category, c(1,6)]
        } else{
            tmp <- kgg[kgg$Pathway %in% category, c(1,6)]
        }
    } else{
        tmp <- kgg[,c(1,6)]
    }
    if(!is.null(nPath)){
        if(is.numeric(nPath)){
            tmp <- tmp[1:nPath, ]
        } else{ stop("nPath must be numeric or NULL")}
    }
    # prepare data
    tmp2 <- tmp %>%  separate_rows( genes, convert=TRUE)
    g <- igraph::graph.data.frame(tmp2, directed=FALSE)
    size <- tmp2 %>% dplyr::select(Pathway) %>% group_by(Pathway) %>%  summarise(n=n())
    size <- left_join(tmp, size, by=c("Pathway"))
    size <- size$n
    V(g)$size <- min(size)/2
    n <- dim(tmp)[1]
    V(g)$size[1:n] <- size
    V(g)$pval <- NA
    V(g)$pval[1:n] <- kgg$P.DE[1:n]
    edge_layer <- geom_edge_link(alpha=.8, colour='darkgrey')
    fc <- V(g)$pval
    V(g)$color <- fc
    createNetworkFromIgraph(g, "customIgraph")
}

# Plot para plotear cnet para kegg ###############
customCnetKegg <- function(kgg, category=NULL, nPath=NULL, byDE=FALSE, nr){
    if(! "ggraph" %in% .packages()) require("ggraph")
    if(! "igraph" %in% .packages()) require("igraph")
    if(! "dplyr" %in% .packages()) require("dplyr")
    if(! "tidyr" %in% .packages()) require("tidyr")
    kgg <- kgg[nr,]
    color_palette <- function(colors) colorRampPalette(colors)(n = 299)
    if(byDE){
        kgg <- kgg %>% arrange(-DE)
    }
    if(!is.null(category)){
        if(is.numeric(category)){
            tmp <- kgg[category,]
        } else{
            tmp <- kgg[kgg$Pathway %in% category,]
        }
    } else{
        tmp <- kgg
    }
    if(!is.null(nPath)){
        if(is.numeric(nPath)){
            tmp <- tmp[1:nPath, ]
        } else{ stop("nPath must be numeric or NULL")}
    }
    pval <- tmp[,c(1,4)]
    tmp <- tmp[,c(1,6)]
    tmp2 <- tmp %>%  separate_rows( genes, convert=TRUE)
    g <- igraph::graph.data.frame(tmp2, directed=FALSE)
    size <- tmp2 %>% dplyr::select(Pathway) %>% group_by(Pathway) %>%  summarise(n=n())
    size <- left_join(tmp, size, by=c("Pathway"))
    pval <- left_join(tmp, pval, by=c("Pathway"))
    size <- size$n
    pval <- pval$P.DE
    V(g)$size <- min(size)/2
    n <- dim(tmp)[1]
    V(g)$size[1:n] <- size
    V(g)$pval <- NA
    V(g)$pval[1:n] <- pval
    edge_layer <- geom_edge_link(alpha=.8, colour='darkgrey')
    fc <- V(g)$pval
    V(g)$color <- fc
    palette <- color_palette(c("red", "blue"))
    p <- ggraph(g, layout="stress", circular=FALSE) +
        edge_layer +
        geom_node_point( aes_(color=~pval, size=~size) ) +
        scale_size(range=c(3, 10), breaks=unique(round(seq(min(size), max(size), length.out=4)))) +
        theme_void()
    p <- p + geom_node_text(aes_(label=~name), data = p$data[1:n,]) +
        scale_color_gradientn(name = "pval", colors=palette, na.value = "#E5C494")
    return(p)
}

# Plot para plotear cnet para GO ###############
customCnetGo <- function(gos, category=NULL, nTerm=NULL, byDE=FALSE, ont="BP"){
    if(! "ggraph" %in% .packages()) require("ggraph")
    if(! "igraph" %in% .packages()) require("igraph")
    if(! "dplyr" %in% .packages()) require("dplyr")
    if(! "tidyr" %in% .packages()) require("tidyr")
    gos <- gos[gos$Ont == ont, ]
    color_palette <- function(colors) colorRampPalette(colors)(n = 299)
    if(byDE){
        gos <- gos %>% arrange(-DE)
    }
    if(!is.null(category)){
        if(is.numeric(category)){
            tmp <- gos[category,]
        } else{
            tmp <- gos[gos$Pathway %in% category,]
        }
    } else{
        tmp <- gos
    }
    if(!is.null(nTerm)){
        if(is.numeric(nTerm)){
            tmp <- tmp[1:nTerm, ]
        } else{ stop("nTerm must be numeric or NULL")}
    }
    pval <- tmp[,c(1,5)]
    tmp <- tmp[,c(1,7)]
    tmp2 <- tmp %>%  separate_rows( genes, convert=TRUE)
    g <- igraph::graph.data.frame(tmp2, directed=FALSE)
    size <- tmp2 %>% dplyr::select(Term) %>% group_by(Term) %>%  summarise(n=n())
    size <- left_join(tmp, size, by=c("Term"))
    pval <- left_join(tmp, pval, by=c("Term"))
    size <- size$n
    pval <- pval$P.DE
    V(g)$size <- min(size)/2
    n <- dim(tmp)[1]
    V(g)$size[1:n] <- size
    V(g)$pval <- NA
    V(g)$pval[1:n] <- pval
    edge_layer <- geom_edge_link(alpha=.8, colour='darkgrey')
    fc <- V(g)$pval
    V(g)$color <- fc
    palette <- color_palette(c("red", "blue"))
    p <- ggraph(g, layout="stress", circular=FALSE) +
        edge_layer +
        geom_node_point( aes_(color=~pval, size=~size) ) +
        scale_size(range=c(3, 10), breaks=unique(round(seq(min(size), max(size), length.out=4)))) +
        theme_void()
    p <- p + geom_node_text(aes_(label=~name), data = p$data[1:n,]) +
        scale_color_gradientn(name = "pval", colors=palette, na.value = "#E5C494")
    return(p)
}

# Función para hacer enrich GO ################
customGO <- function(data, universe = NULL, species = "Hs", prior.prob = NULL,
    covariate = NULL, plot = FALSE, coef = 1, FDR = 0.05, golevelFile="resources/GOlevel.Rds") {
    if (!is.data.frame(data)) {
        stop("de should be a data.frame with firt column as symbol
        Id and second column as entrez Id.")
    }
    if( dim(data)[2] != 2){
        stop("de should be a data.frame with firt column as symbol
        Id and second column as entrez Id.")
    }
    names(data) <- c("SYMBOL","ENTREZID")
    de <- data
    de <- as.character(de$ENTREZID)
    suppressPackageStartupMessages(OK <- requireNamespace("GO.db",quietly = TRUE))
    if (!OK) stop("GO.db package required but not installed (or can't be loaded)")
    suppressPackageStartupMessages(OK <- requireNamespace("AnnotationDbi",quietly = TRUE))
    if (!OK) stop("AnnotationDbi package required but not installed (or can't be loaded)")
    orgPkg <- paste0("org.", species, ".eg.db")
    suppressPackageStartupMessages(OK <- requireNamespace(orgPkg, quietly = TRUE))
    if (!OK)
        stop(orgPkg, " package required but not not installed (or can't be loaded)")
    obj <- paste0("org.", species, ".egGO2ALLEGS")
    egGO2ALLEGS <- tryCatch(getFromNamespace(obj, orgPkg), error = function(e) FALSE)
    if (is.logical(egGO2ALLEGS))
        stop("Can't find gene ontology mappings in package ", orgPkg)
    if (is.null(universe)) {
        #GeneID.PathID <- AnnotationDbi::toTable(egGO2ALLEGS)[, c("gene_id",
         #  "go_id", "Ontology")]
        file1 <- paste0("resources/",species,"GOlinks.Rds")
        GeneID.PathID <- readRDS(file1)
        i <- !duplicated(GeneID.PathID[, c("gene_id", "go_id")])
        GeneID.PathID <- GeneID.PathID[i, ]
        universe <- unique(GeneID.PathID[, 1])
        prior.prob <- covariate <- NULL
    } else {
        universe <- as.character(universe)
        lu <- length(universe)
        if (!is.null(prior.prob) && length(prior.prob) != lu)
            stop("universe and prior.prob must have same length")
        if (!is.null(covariate) && length(covariate) != lu)
            stop("universe and covariate must have same length")
        if (anyDuplicated(universe)) {
            i <- !duplicated(universe)
            if (!is.null(covariate))
                covariate <- covariate[i]
            if (!is.null(prior.prob))
                prior.prob <- prior.prob[i]
            universe <- universe[i]
        }
        i <- (universe %in% AnnotationDbi::Lkeys(egGO2ALLEGS))
        universe <- universe[i]
        if (!is.null(covariate))
            covariate <- covariate[i]
        if (!is.null(prior.prob))
            prior.prob <- prior.prob[i]
        AnnotationDbi::Lkeys(egGO2ALLEGS) <- universe
        GeneID.PathID <- AnnotationDbi::toTable(egGO2ALLEGS)[, c("gene_id",
            "go_id", "Ontology")]
        d <- duplicated(GeneID.PathID[, c("gene_id", "go_id")])
        GeneID.PathID <- GeneID.PathID[!d, ]
    }
    if (is.list(de)) {
        if (is.data.frame(de))
            stop("de should be a list of character vectors. It should not be a data.frame.")
    } else {
        de <- list(DE = de)
    }
    nsets <- length(de)
    if (!all(vapply(de, is.vector, TRUE)))
        stop("components of de should be vectors")
    for (s in 1:nsets) de[[s]] <- unique(as.character(de[[s]]))
    names(de) <- trimWhiteSpace(names(de))
    NAME <- names(de)
    i <- which(NAME == "" | is.na(NAME))
    NAME[i] <- paste0("DE", i)
    names(de) <- makeUnique(NAME)
    NGenes <- length(universe)
    if (NGenes < 1L)
        stop("No annotated genes found in universe")
    for (s in 1:nsets) de[[s]] <- de[[s]][de[[s]] %in% universe]
    i <- GeneID.PathID[, 1] %in% universe
    if (sum(i) == 0L)
        stop("Pathways do not overlap with universe")
    GeneID.PathID <- GeneID.PathID[i, ]
    if (!is.null(covariate)) {
        if (!is.null(prior.prob))
            message("prior.prob being recomputed from covariate")
        covariate <- as.numeric(covariate)
        isDE <- (universe %in% unlist(de))
        o <- order(covariate)
        prior.prob <- covariate
        span <- approx(x = c(20, 200), y = c(1, 0.5), xout = sum(isDE),
            rule = 2)$y
        prior.prob[o] <- tricubeMovingAverage(isDE[o], span = span)
        if (plot)
            barcodeplot(covariate, index = isDE, worm = TRUE, span.worm = span,
                main = "DE status vs covariate")
    }
    if (is.null(prior.prob)) {
        X <- matrix(1, nrow(GeneID.PathID), nsets + 1)
        colnames(X) <- c("N", names(de))
    } else {
        names(prior.prob) <- universe
        X <- matrix(1, nrow(GeneID.PathID), nsets + 2)
        X[, nsets + 2] <- prior.prob[GeneID.PathID[, 1]]
        colnames(X) <- c("N", names(de), "PP")
    }
    for (s in 1:nsets) X[, s + 1] <- (GeneID.PathID[, 1] %in% de[[s]])
    S <- rowsum(X, group = GeneID.PathID[, 2], reorder = FALSE)
    PValue <- matrix(0, nrow = nrow(S), ncol = nsets)
    colnames(PValue) <- paste("P", names(de), sep = ".")
    nde <- lengths(de, use.names = FALSE)
    if (!is.null(prior.prob)) {
        SumPP <- sum(prior.prob)
        M2 <- NGenes - S[, "N"]
        Odds <- S[, "PP"]/(SumPP - S[, "PP"]) * M2/S[, "N"]
        if (!requireNamespace("BiasedUrn", quietly = TRUE))
            stop("BiasedUrn package required but is not installed (or can't be loaded)")
        for (j in seq_len(nsets)) for (i in seq_len(nrow(S))) PValue[i,
            j] <- BiasedUrn::pWNCHypergeo(S[i, 1L + j], S[i, "N"],
            M2[i], nde[j], Odds[i], lower.tail = FALSE) + BiasedUrn::dWNCHypergeo(S[i,
            1L + j], S[i, "N"], M2[i], nde[j], Odds[i])
        S <- S[, -ncol(S)]
    } else {
        for (j in seq_len(nsets)) PValue[, j] <- phyper(S[, 1L + j] -
            0.5, nde[j], NGenes - nde[j], S[, "N"], lower.tail = FALSE)
    }
    GOID <- rownames(S)
    TERM <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = GOID,
        columns = "TERM"))
    m <- match(GOID, GeneID.PathID[, 2])
    Ont <- GeneID.PathID[m, 3]
    ## aportacion
    kk <- GeneID.PathID  # [ which(GeneID.PathID$gene_id %in% de[[1]]), ]
    kkb <- left_join(data, kk, by = c("ENTREZID" = "gene_id"))
    kkb <- kkb[which(!is.na(kkb$ENTREZID)), ]
    kkk <- kkb %>% group_by(go_id) %>%
        summarize(genes = paste(ENTREZID, collapse = ", "))
    ##
    Results <- data.frame(Term = TERM[, 2], Ont = Ont, S, PValue, stringsAsFactors = FALSE)
    Results$go_id <- rownames(Results)
    resultado <- left_join(Results, kkk, by = c(go_id = "go_id"))
    resultado <- resultado[which(resultado$P.DE < 0.05), ]
    GOlevel <- readRDS(golevelFile)
    GOlevel$id <- as.character(GOlevel$id)
    resultado <- left_join(resultado, GOlevel, by = c("go_id"="id"))
    return(resultado)
}

# Función para hacer enrich kegg ################
customKegg <- function(data, universe = NULL, restrict.universe = FALSE,
    species = "Hs", species.KEGG = "hsa", convert = FALSE, gene.pathway = NULL,
    pathway.names = NULL, prior.prob = NULL, covariate = NULL,
    plot = FALSE, ...) {
    if( !is.data.frame(data)){
        stop("de should be a data.frame with firt column as symbol Id and second
             column as entrez Id.")
    }
    if( dim(data)[2] != 2){
        stop("de should be a data.frame with firt column as symbol
        Id and second column as entrez Id.")
    }
    names(data) <- c("SYMBOL","ENTREZID")
    de <- data
    de <- as.vector(de[,2])
    de <- list(DE = de)
    nsets <- length(de)
    if (!all(vapply(de, is.vector, TRUE))) {
        stop("components of de should be vectors")
    }
    for (s in 1:nsets) de[[s]] <- unique(as.character(de[[s]]))
    names(de) <- trimWhiteSpace(names(de))
    NAME <- names(de)
    i <- which(NAME == "" | is.na(NAME))
    NAME[i] <- paste0("DE", i)
    names(de) <- makeUnique(NAME)
    if (is.null(species.KEGG)) {
        species <- match.arg(species, c("Ag", "At", "Bt", "Ce",
            "Dm", "Dr", "EcK12", "EcSakai", "Gg", "Hs", "Mm",
            "Mmu", "Pf", "Pt", "Rn", "Ss", "Xl"))
        species.KEGG <- switch(species, Ag = "aga", At = "ath",
            Bt = "bta", Ce = "cel", Cf = "cfa", Dm = "dme", Dr = "dre",
            EcK12 = "eco", EcSakai = "ecs", Gg = "gga", Hs = "hsa",
            Mm = "mmu", Mmu = "mcc", Pf = "pfa", Pt = "ptr",
            Rn = "rno", Ss = "ssc", Xl = "xla")
    }
    if (is.null(gene.pathway)) {
        #GeneID.PathID <- getGeneKEGGLinks(species.KEGG, convert = convert)
         file1 <- paste0("resources/",species.KEGG,"KeggLinks.Rds")
         GeneID.PathID <- readRDS(file1)
    } else {
        GeneID.PathID <- gene.pathway
        d <- dim(GeneID.PathID)
        if (is.null(d)) {
            stop("gene.pathway must be data.frame or matrix")
        }
        if (d[2] < 2) {
            stop("gene.pathway must have at least 2 columns")
        }
        isna <- rowSums(is.na(GeneID.PathID[, 1:2])) > 0.5
        GeneID.PathID <- GeneID.PathID[!isna, ]
        ID.ID <- paste(GeneID.PathID[, 1], GeneID.PathID[, 2],
            sep = ".")
        if (anyDuplicated(ID.ID)) {
            d <- duplicated(ID.ID)
            GeneID.PathID <- GeneID.PathID[!d, ]
        }
    }
    if (is.null(pathway.names)) {
        #PathID.PathName <- getKEGGPathwayNames(species.KEGG,
         #   remove.qualifier = TRUE)
        file2 <- paste0("resources/",species.KEGG,"KeggPathwayNames.Rds")
        PathID.PathName <- readRDS(file2)
    } else {
        PathID.PathName <- pathway.names
        d <- dim(PathID.PathName)
        if (is.null(d)) {
            stop("pathway.names must be data.frame or matrix")
        }
        if (d[2] < 2) {
            stop("pathway.names must have at least 2 columns")
        }
        isna <- rowSums(is.na(PathID.PathName[, 1:2])) > 0.5
        PathID.PathName <- PathID.PathName[!isna, ]
    }
    if (is.null(universe)) {
        universe <- unique(GeneID.PathID[, 1])
        prior.prob <- covariate <- NULL
    } else {
        universe <- as.character(universe)
        lu <- length(universe)
        if (!lu) {
            stop("No genes in universe")
        }
        if (!is.null(prior.prob) && length(prior.prob) != lu) {
            stop("universe and prior.prob must have same length")
        }
        if (!is.null(covariate) && length(covariate) != lu) {
            stop("universe and covariate must have same length")
        }
        if (restrict.universe) {
            i <- universe %in% GeneID.PathID[, 1]
            universe <- universe[i]
            if (!is.null(prior.prob)) {
                prior.prob <- prior.prob[i]
            }
            if (!is.null(covariate)) {
                covariate <- covariate[i]
            }
        }
    }
    if (anyDuplicated(universe)) {
        d <- duplicated(universe)
        if (!is.null(covariate)) {
            covariate <- rowsum(covariate, group = universe,
                reorder = FALSE)
            n <- rowsum(rep_len(1L, length(universe)), group = universe,
                reorder = FALSE)
            covariate <- covariate/n
        }
        if (!is.null(prior.prob)) {
            prior.prob <- rowsum(prior.prob, group = universe,
                reorder = FALSE)
            n <- rowsum(rep_len(1L, length(universe)), group = universe,
                reorder = FALSE)
            prior.prob <- prior.prob/n
        }
        universe <- universe[!d]
    }
    NGenes <- length(universe)
    if (NGenes < 1L) {
        stop("No annotated genes found in universe")
    }
    for (s in 1:nsets) de[[s]] <- de[[s]][de[[s]] %in% universe]
    i <- GeneID.PathID[, 1] %in% universe
    if (sum(i) == 0L) {
        stop("Pathways do not overlap with universe")
    }
    GeneID.PathID <- GeneID.PathID[i, ]
    if (!is.null(covariate)) {
        if (!is.null(prior.prob)) {
            message("prior.prob being recomputed from covariate")
        }
        covariate <- as.numeric(covariate)
        isDE <- (universe %in% unlist(de))
        o <- order(covariate)
        prior.prob <- covariate
        span <- approx(x = c(20, 200), y = c(1, 0.5), xout = sum(isDE),
            rule = 2)$y
        prior.prob[o] <- tricubeMovingAverage(isDE[o], span = span)
        if (plot) {
            barcodeplot(covariate, index = isDE, worm = TRUE,
                span.worm = span, main = "DE status vs covariate")
        }
    }
    if (is.null(prior.prob)) {
        X <- matrix(1, nrow(GeneID.PathID), nsets + 1)
        colnames(X) <- c("N", names(de))
    } else {
        names(prior.prob) <- universe
        X <- matrix(1, nrow(GeneID.PathID), nsets + 2)
        X[, nsets + 2] <- prior.prob[GeneID.PathID[, 1]]
        colnames(X) <- c("N", names(de), "PP")
    }
    for (s in 1:nsets) {
        X[, s + 1] <- (GeneID.PathID[, 1] %in% de[[s]])
    }
    S <- rowsum(X, group = GeneID.PathID[, 2], reorder = FALSE)
    PValue <- matrix(0, nrow = nrow(S), ncol = nsets)
    colnames(PValue) <- paste("P", names(de), sep = ".")
    nde <- lengths(de, use.names = FALSE)
    if (!is.null(prior.prob)) {
        SumPP <- sum(prior.prob)
        M2 <- NGenes - S[, "N"]
        Odds <- S[, "PP"]/(SumPP - S[, "PP"]) * M2/S[, "N"]
        if (!requireNamespace("BiasedUrn", quietly = TRUE)) {
            stop("BiasedUrn package required but is not installed (or can't be loaded)")
        }
        for (j in seq_len(nsets)) {
            for (i in seq_len(nrow(S))) {
                PValue[i, j] <- BiasedUrn::pWNCHypergeo(S[i,
                  1L + j], S[i, "N"], M2[i], nde[j], Odds[i],
                  lower.tail = FALSE) + BiasedUrn::dWNCHypergeo(S[i,
                  1L + j], S[i, "N"], M2[i], nde[j], Odds[i])
            }
        }
        S <- S[, -ncol(S)]
    } else {
        for (j in seq_len(nsets)) {
            PValue[, j] <- phyper(S[, 1L + j] - 0.5, nde[j],
                NGenes - nde[j], S[, "N"], lower.tail = FALSE)
        }
    }

    g <- rownames(S)
    m <- match(g, PathID.PathName[, 1])
    path <- PathID.PathName[m, 2]
    ##
    kk <- GeneID.PathID
    data <- data.frame(ENTREZID = data[,2])
    data$ENTREZID <- as.character(data$ENTREZID)
    kkb <- inner_join(kk, data, by = c(GeneID = "ENTREZID"))
    kkb <- kkb[!duplicated(kkb), ]
    # kkb <- kkb[ which( !is.na( kkb$SYMBOL ) ), ]
    kkk <- kkb %>% group_by(PathwayID) %>% summarize(genes = paste(GeneID,
        collapse = ","))
    ##

    Results <- data.frame(Pathway = path, S, PValue, stringsAsFactors = FALSE)
    rownames(Results) <- g
    Results$pathID <- rownames(Results)
    resultado <- left_join(Results, kkk, by = c(pathID = "PathwayID"))
    resultado <- resultado[which(resultado$P.DE < 0.05), ]
    return(resultado)
}

# Función para crear tablas  con desplegable de genes ############
datatable2 <- function(x, vars = NULL, opts = NULL, ...) {
  names_x <- names(x)
  if (is.null(vars))
    stop("'vars' must be specified!")
  pos <- match(vars, names_x)
  if (any(map_chr(x[, pos], typeof) == "list"))
    stop("list columns are not supported in datatable2()")
  pos <- pos[pos <= ncol(x)] + 1
  rownames(x) <- NULL
  if (nrow(x) > 0)
    x <- cbind(` ` = "&oplus;", x)
  if(dim(x)[2]==7){
      js <- c("function(row, data) {",
                "if (data[5]>1000 | data[5]<1){",
                "$('td:eq('+(4)+')', row).html(data[5].toExponential(3));",
                "}",
                "}")
  }else if(dim(x)[2]==9){
      js <- c("function(row, data) {",
                "if (data[6]>1000 | data[6]<1){",
                "$('td:eq('+(5)+')', row).html(data[6].toExponential(3));",
                "}",
                "}")
  }
  opts <- c(opts,
            list(
              columnDefs = list(list(visible = FALSE, targets = c(0,pos)),
                                list(orderable = FALSE,
                                     className = "details-control", targets = 1),
                                list(className = "dt-left", targets = 1:2),
                                list(className = "dt-right",targets = 3:ncol(x) 
                                     )
                                ),
              rowCallback = JS(js),
              dom = "Bfrtipl",
              buttons = c("copy", "csv", "excel", "pdf", "print") ) )
  datatable(x, ..., extensions = "Buttons", options = opts, 
            callback = JS(.callback2(x = x, pos = c(0, pos) ) ) )
}

# funcion auxiliar de datatable2 ##############################
.callback2 <- function(x, pos = NULL) {
    part1 <- "table.column(1).nodes().to$().css({cursor: 'pointer'});"
    part2 <- .child_row_table2(x, pos = pos)
    part3 <- "
   table.on('click', 'td.details-control', function() {
    var td = $(this), row = table.row(td.closest('tr'));
    if (row.child.isShown()) {
      row.child.hide();
      td.html('&oplus;');
    } else {
      row.child(format(row.data())).show();
      td.html('&ominus;');
    }
  });"
    paste(part1, part2, part3)
}

# funcion auxiliar de datatable2 ##############################
.child_row_table2 <- function(x, pos = NULL) {
    names_x <- paste0(names(x), ":")
    text <- "var format = function(d) {
    text = '<div><table >' +"
    for (i in seq_along(pos)) {
        text <- paste(text, glue::glue("'<tr>' +
          '<td>' + '{names_x[pos[i]]}' + '</td>' +
          '<td>' + d[{pos[i]}] + '</td>' +
        '</tr>' + "))
    }
    paste0(text, "'</table></div>'
      return text;};")
}

# funcion que preparar los datos de enrich go para pasárlos a datatable2 ###############
go2DT <- function(enrichdf, data, orderby = NULL, nrows = NULL) {
    if(!is.data.frame(enrichdf) | !is.data.frame(data)){
        stop("enrichdf and data should be data.frame")
    }
    names(data) <- c("SYMBOL", "ENTREZID")
    if(!is.null(orderby)){
        orderby = match.arg(orderby, c("DE", "P.DE", "N", "Term"))
        if(orderby=="P.DE" | orderby =="Term"){
            enrichdf <- enrichdf %>% arrange(get(orderby))
        } else{ enrichdf <- enrichdf %>% arrange(desc(get(orderby)))}
    }
    enrichdf2 <- enrichdf %>%
        mutate(url = paste0("<a href='http://amigo.geneontology.org/amigo/search/bioentity?q=",
                            go_id,"' target='_blank'>",go_id,"</a>"))
    CAup <- enrichdf2[, c(1, 2, 3, 4, 5, 7, 8, 9)]
    CAup$genes <- gsub(",", ", ", CAup$genes)
    for (i in seq(1:length(CAup$genes))) {
       mg <- as.numeric(unlist(strsplit(CAup$genes[i], ", ")))
       mg2 <- match(mg, data$ENTREZID)
       CAup$genes[i] <- paste0(data[mg2, 1], collapse = ", ")
    }
    splitGenes <- strsplit(CAup$genes, ", ")
    CAup$genes <- lapply(
        splitGenes, function(x){
            paste(sort(x), collapse = ", ")
            }
        )
    if(!is.null(nrows) & is.numeric(nrows)){
        CAup <- CAup[1:nrows, ]
    }
    #CAup <- CAup %>% mutate(P.DE = format(P.DE, scientific = T, digits = 4))
    return(CAup)
}

# Recupera todos los ids de GO y el nivel al que pertenecen
# Ejemplo de uso:
# GOlevel = getGOlevel()
# Guardar lo que genera en resources/GOlevel.Rds #############
getGOlevel <- function(){
    bp <- "GO:0008150"
    mf <- "GO:0003674"
    cc <- "GO:0005575"
    getAllBPChildren <- function(goids){
    ans <- unique(unlist(mget(goids, GOBPCHILDREN), use.names=FALSE))
    ans <- ans[!is.na(ans)]
    return(ans)
    }
    getAllCCChildren <- function(goids){
    ans <- unique(unlist(mget(goids, GOCCCHILDREN), use.names=FALSE))
    ans <- ans[!is.na(ans)]
    return(ans)
    }
    getAllMFChildren <- function(goids){
    ans <- unique(unlist(mget(goids, GOMFCHILDREN), use.names=FALSE))
    ans <- ans[!is.na(ans)]
    return(ans)
    }
    reskk <- list()
    dfBP <- data.frame()
    dfCC <- data.frame()
    dfMF <- data.frame()
    kk <- bp
    for (i in seq_len(100)) {
        reskk[[i]] <- getAllBPChildren(kk)
        if (is.logical(reskk[[i]])) {
            reskk <- reskk[-i]
            break()
        }
        dfBP <- rbind(dfBP, data.frame(id = reskk[[i]], level = i))
        kk <- reskk[[i]]
    }
    kk <- cc
    reskk <- list()
    for (i in seq_len(100)) {
        reskk[[i]] <- getAllCCChildren(kk)
        if (is.logical(reskk[[i]])) {
            reskk <- reskk[-i]
            break()
        }
        dfCC <- rbind(dfCC, data.frame(id = reskk[[i]], level = i))
        kk <- reskk[[i]]
    }
    reskk <- list()
    kk <- mf
    for (i in seq_len(100)) {
        reskk[[i]] <- getAllMFChildren(kk)
        if (is.logical(reskk[[i]])) {
            reskk <- reskk[-i]
            break()
        }
        dfMF <- rbind(dfMF, data.frame(id = reskk[[i]], level = i))
        kk <- reskk[[i]]
    }
    GOlevel <- data.frame(rbind(dfBP,dfCC, dfMF))
    GOlevel <- rbind(
        GOlevel,
        data.frame(
            id = c("GO:0008150","GO:0003674","GO:0005575"),
            level= 0))
    GOlevel <- aggregate(level~id, GOlevel, function(x)x[which.min(abs(x))])
    return(GOlevel)
}

# funcion que preparar los datos de enrich kegg para pasárlos a datatable2 ###############
kegg2DT <- function(enrichdf, data, orderby = NULL, nrows = NULL) {
    if(!is.data.frame(enrichdf) | !is.data.frame(data)){
        stop("enrichdf and data should be data.frame")
    }
    names(data) <- c("SYMBOL", "ENTREZID")
    if(!is.null(orderby)){
        orderby = match.arg(orderby, c("DE", "P.DE", "N", "Pathway"))
        if(orderby=="P.DE" | orderby =="Pathway"){
            enrichdf <- enrichdf %>% arrange(get(orderby))
        } else{ enrichdf <- enrichdf %>% arrange(desc(get(orderby)))}
    }
    enrichdf2 <- enrichdf %>%
        mutate(pathway = limma::strsplit2(enrichdf$pathID, ":")[, 2],
               genesURL = gsub(",", "+", genes)) %>%
        mutate(url = paste0("<a href='https://www.kegg.jp/kegg-bin/show_pathway?",
        pathway, "/", genesURL,"' target='_blank'>", pathway, "</a>"))
    CAup <- enrichdf2[, c(1, 2, 3, 4, 6, 9)]
    CAup$genes <- gsub(",", ", ", CAup$genes)
    for (i in seq(1:length(CAup$genes))) {
        mg <- as.numeric(unlist(strsplit(CAup$genes[i], ", ")))
        mg2 <- match(mg, data$ENTREZID)
        CAup$genes[i] <- paste0(data[mg2, 1], collapse = ", ")
    }
    splitGenes <- strsplit(CAup$genes, ", ")
    CAup$genes <- lapply(
        splitGenes, function(x){
            paste(sort(x), collapse = ", ")
            }
        )
    if(!is.null(nrows) & is.numeric(nrows)){
        CAup <- CAup[1:nrows, ]
    }
    #CAup <- CAup %>% mutate(P.DE = format(P.DE, scientific = T, digits = 4))
    return(CAup)
}

# Plot barrar de GO ####################
plotGO <- function(enrichdf, nrows = 30, orderby=NULL, ont){
    require(plotly)
    if(!is.data.frame(enrichdf)){
        stop("enrichdf should be data.frame")
    }
    if(!exists("ont") | !(ont %in% c("BP","MF","CC")) ){
        stop("A valid value should be provided for 'ont'")
    }
    dataTitle <- list(BP=c("Biological Process", 'rgb(227,74,51)'),
                      MF=c("Molecular Function",'rgb(31,119,180)'),
                      CC=c("Cellular Component", 'rgb(49,163,84)'))
    enrichdf <- enrichdf[enrichdf$Ont == ont, ]
    if(!is.null(orderby)){
        orderby = match.arg(orderby, c("DE", "P.DE", "N", "Term"))
        if(orderby=="P.DE" | orderby =="Term"){
            enrichdf <- enrichdf %>% arrange(get(orderby))
        } else{ enrichdf <- enrichdf %>% arrange(desc(get(orderby)))}
    }
    p <- enrichdf[1:nrows,] %>%
        plot_ly(x=~DE, y=~go_id, text=~Term, type = "bar",
                marker = list(color=dataTitle[[ont]][2]),
                orientation = "v") %>%
        layout(margin = list(l=100), yaxis = list(title=""),
               title=dataTitle[[ont]][1])
    return(p)
}

# Plot barrar de Kegg #############################333
plotKegg <- function(enrichdf, nrows = 30, orderby=NULL){
    require(plotly)
    if(!is.data.frame(enrichdf)){
        stop("enrichdf should be data.frame")
    }
    if(!is.null(orderby)){
        orderby = match.arg(orderby, c("DE", "P.DE", "N", "Pathway"))
        if(orderby=="P.DE" | orderby =="Pathway"){
            enrichdf <- enrichdf %>% arrange(get(orderby))
        } else{ enrichdf <- enrichdf %>% arrange(desc(get(orderby)))}
    }
    p <- enrichdf[1:nrows,] %>%
        plot_ly(x=~DE, y=~pathID, text=~Pathway, type = "bar",
                marker = list(color='rgb(117,107,177)'),
                orientation = "v") %>%
        layout(margin = list(l=100), yaxis = list(title=""),
               title="Kegg pathways")
    return(p)
}

# Función sin uso actualmente -creo- #############
loadGenes <- function(filegenes){
  load(filegenes)
  auxgenes <- genes
}

    # PCA de un objeto DESeq #####################

plotPCA = function(object, intgroup = "condition", ntop = 500, returnData = TRUE){
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  # select the ntop genes by variance
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select, ]))
  # the contribution to the total variance for each component
  percentVar <- pca$sdev ^ 2 / sum(pca$sdev ^ 2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <-
    as.data.frame(colData(object)[, intgroup, drop = FALSE])
  # add the intgroup factors together to create a new grouping factor
  # group <- if (length(intgroup) > 1) {
  #     factor(apply(intgroup.df, 1, paste, collapse = ":"))
  # } else {
  #     colData(object)[[intgroup]]
  # }
  if(length(intgroup)>1){
    colgroup <- factor(intgroup.df[ ,intgroup[1] ] )
    shapegroup <- factor(intgroup.df[ ,intgroup[2] ] )
    d <-
      data.frame(
        PC1 = pca$x[, 1],
        PC2 = pca$x[, 2],
        group = colgroup,
        shape = shapegroup,
        intgroup.df,
        name = colnames(object)
      )
  } else{
    colgroup <- factor(intgroup.df[ ,intgroup[1] ] )
    d <-
      data.frame(
        PC1 = pca$x[, 1],
        PC2 = pca$x[, 2],
        group = colgroup,
        intgroup.df,
        name = colnames(object)
      )
  }
  # assembly the data for the plot
  
  getPalette <- colorRampPalette(c("#008000","#800080"))
  colours <- getPalette(length(levels(d$group)))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    #return(d)
  }
  
  if(length(intgroup)>1){
    p <- ggplot(data = d,
                aes_string(x = "PC1", y = "PC2", color = "group", shape = "shape")) +
      geom_point(size = 3) +
      ggtitle("PCA for VST data transformation") +
      xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
      ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
      scale_color_manual(values = colours, name = intgroup[1]) +
      scale_shape_manual(values = seq_len(length(d$shape)), name=intgroup[2] )+
      coord_fixed() +
      ggrepel::geom_text_repel(aes(label = paste("",d$name, sep = ""))) +
      theme(text = element_text(size=20))}
  else{
    p <- ggplot(data = d,
                aes_string(x = "PC1", y = "PC2", color = "group")) +
      geom_point(size = 3) +
      ggtitle("PCA for transformed data according to variance") +
      xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
      ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
      scale_color_manual(values = colours, name = intgroup[1]) +
      coord_fixed() +
      ggrepel::geom_text_repel(aes(label = paste("",d$name, sep = ""))) +
      theme(text = element_text(size=20))
  }
  return(p)
}




# Función para recuperar los genes up de un objeto DEseq #############
# actualmente para p-val <= 0.05 fijo.
getSigUpregulated <- function(dds, pval=0.05, logfc=0){
  res.sh <- lfcShrink(dds, coef=2, type="apeglm", res = results(dds))
  rk <- as.data.frame(res.sh)
  rk <- rk[rk$log2FoldChange >logfc & rk$padj<=pval,]
  rk <- rk[ order(rk$pvalue, decreasing = TRUE), ]
  annot <- geneIdConverter(rownames(rk))
  return(data.frame(SYMBOL = annot$consensus, ENTREZID = annot$ENTREZID, stringsAsFactors = F) )
}

# Función para recuperar los genes down de un objeto DEseq #############
# actualmente para p-val <= 0.05 fijo.
getSigDownregulated <- function(dds, pval=0.05, logfc=0){
  res.sh <- lfcShrink(dds, coef=2, type="apeglm", res = results(dds))
  rk <- as.data.frame(res.sh)
  rk <- rk[rk$log2FoldChange <logfc & rk$padj<=pval,]
  rk <- rk[ order(rk$pvalue, decreasing = TRUE), ]
  annot <- geneIdConverter(rownames(rk))
  return(data.frame(SYMBOL = annot$consensus, ENTREZID = annot$ENTREZID, stringsAsFactors = F) )
}

# Convertidor de nombres de genes ###################
# Se le pasa un vector de ensembl y devuelve un df con varios nombres
geneIdConverter <- function(genes){ # genes = vector of ensembl gene ids (sólo para Mm por ahora)
  require("EnsDb.Mmusculus.v79")
  require("org.Mm.eg.db")
  annot <- NULL
  annot$ENSEMBL <- genes
  annot$SYMBOL <-  mapIds(EnsDb.Mmusculus.v79, keys=genes, column="SYMBOL",keytype="GENEID")
  annot$SYMBOL1 <- mapIds(org.Mm.eg.db, keys = genes, column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first') 
  annot$description <- mapIds(org.Mm.eg.db, keys = genes, column = 'GENENAME', keytype = 'ENSEMBL', multiVals = 'first')
  annot <- as.data.frame(annot)
  consensus <- data.frame('Symbol'= ifelse(!is.na(annot$SYMBOL), as.vector(annot$SYMBOL),
                                           ifelse(!is.na(annot$SYMBOL1),as.vector(annot$SYMBOL1),
                                                  as.vector(annot$ENSEMBL))), stringsAsFactors = F)
  annot$consensus <- consensus$Symbol
  entrez1 <- mapIds(org.Mm.eg.db, keys = annot$consensus, column = "ENTREZID", keytype = "SYMBOL")
  entrez2 <- mapIds(org.Mm.eg.db, keys = as.character(annot$ENSEMBL),
                    column = "ENTREZID", keytype = "ENSEMBL")
  annot$entrez1 <- entrez1
  annot$entrez2 <- entrez2
  ENTREZID <- ifelse(!is.na(annot$entrez1), annot$entrez1, annot$entrez2)
  annot$ENTREZID <- ENTREZID
  return(annot)
}

# Dotplot de objeto enrich kegg ##########################
dotPlotkegg <- function(data, n = 20){
  data$ratio <- data$DE/data$N
  data <- data[order(data$ratio, decreasing = F), ]
  data <- data[seq_len(n),]
  data$Pathway <- factor(data$Pathway, levels = data$Pathway)
  p <- ggplot(data, aes(y=Pathway, x=ratio, color=P.DE))+
    geom_point(aes(size=DE), stat="identity")+
    theme_bw()+
    labs(x = "ratio (DE/N)") +
    scale_color_continuous(low = "red", high = "blue", 
                           guide = guide_colorbar(reverse = TRUE))+
    theme(text = element_text(size=20))
  return(p)
}

# Dotplot de objeto enrich GO ##########################
dotPlotGO <- function(data, n = 20){
  data$ratio <- data$DE/data$N
  data <- data[order(data$ratio, decreasing = F), ]
  data <- data[seq_len(n),]
  data$Term <- factor(data$Term, levels = data$Term)
  p <- ggplot(data, aes(y=Term, x=ratio, color=P.DE))+
    geom_point(aes(size=DE), stat="identity")+
    theme_bw()+
    labs(x = "ratio (DE/N)") +
    scale_color_continuous(low = "red", high = "blue", 
                           guide = guide_colorbar(reverse = TRUE))+
    theme(text = element_text(size=20))
  return(p)
}

# Heatmap de objeto enrich kegg ##########################
heatmapKegg <- function(kdt, nr){
  kdt <- kdt[nr, ]
  colourCount <- length(unique(kdt$DE)) # number of levels
  getPalette <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))
  kdt %>% dplyr::select(Pathway, genes,DE) %>% 
    separate_rows(genes) %>%
    mutate(Pathway = fct_inorder(Pathway)) %>% 
    mutate(Pathway = fct_rev(Pathway)) %>% 
    mutate(genes = fct_infreq(genes)) %>% 
    mutate(DE = factor(DE)) %>% 
    ggplot(aes_(~genes, ~Pathway)) + 
    geom_tile(aes_(fill = ~DE), color = 'black', size =0.2) +
    xlab(NULL) + ylab(NULL) +
    theme_minimal() +
    theme(panel.grid.major = element_line(colour = "gray88", size = 0.8),
          axis.text.x=element_text(angle = 55, hjust = 1, vjust = 0.5, size=10))+
    # scale_fill_continuous(low="blue", high="red", name = "N")
    # scale_fill_brewer(palette = "YlOrRd")
    scale_fill_manual(values = getPalette(colourCount))+
    theme(text = element_text(size=16, angle=0))
  
}

# Función para crear dataset para hacer GESA pathway ##################
buildKeggDataset <- function(specie="mmu"){
  GeneID.PathID <- getGeneKEGGLinks(specie, convert = FALSE)
  PathName <- getKEGGPathwayNames(specie,remove.qualifier = TRUE)
  PathName$Id <- paste(PathName$PathwayID,PathName$Description,sep="_")
  dataSet <- left_join(GeneID.PathID, PathName, by = c("PathwayID"="PathwayID"))
  dataSet$Id <- gsub("path:","",dataSet$Id)
  dataSet <- dataSet[,c(4,1)]
  saveRDS(dataSet,"resources/keggDataGSEA.Rds")
  }

# Función para hacer GSEA pathway #################################
gseaKegg <- function(dds){
  pathwayDataSet <- readRDS("resources/keggDataGSEA.Rds")
  
  res.sh <- as.data.frame(lfcShrink(dds, coef=2, type="apeglm", res = results(dds)))
  res.sh <- res.sh[order(res.sh$log2FoldChange, decreasing = TRUE), ]
  res.sh$ENSEMBL <- rownames(res.sh)
  geneRank <- geneIdConverter( res.sh$ENSEMBL)
  resRank <- left_join(res.sh, geneRank, by=c("ENSEMBL"="ENSEMBL"))
  resRank <- resRank[!is.na(resRank$ENTREZID), c("ENTREZID","log2FoldChange") ]
  vectRank <- resRank$log2FoldChange
  attr(vectRank, "names") <- as.character(resRank$ENTREZID)
  mygsea <- clusterProfiler::GSEA(vectRank, 
                                  TERM2GENE = pathwayDataSet, 
                                  by="fgsea", pvalueCutoff = 0.1)
  mygsea <- DOSE::setReadable(mygsea, "org.Mm.eg.db", "ENTREZID")
  return(mygsea)
}

# Función para actualizar las bases de datos de kegg y GO #############
# esto mejora la velocidad de los enrich en unos 10 segs
updateDatabases <- function(species){
    species.KEGG <- NULL
        if (is.null(species.KEGG)) {
        species <- match.arg(species, c("Ag", "At", "Bt", "Ce",
            "Dm", "Dr", "EcK12", "EcSakai", "Gg", "Hs", "Mm",
            "Mmu", "Pf", "Pt", "Rn", "Ss", "Xl"))
        species.KEGG <- switch(species, Ag = "aga", At = "ath",
            Bt = "bta", Ce = "cel", Cf = "cfa", Dm = "dme", Dr = "dre",
            EcK12 = "eco", EcSakai = "ecs", Gg = "gga", Hs = "hsa",
            Mm = "mmu", Mmu = "mcc", Pf = "pfa", Pt = "ptr",
            Rn = "rno", Ss = "ssc", Xl = "xla")
        }
    # Kegg
    GeneID.PathID <- getGeneKEGGLinks(species.KEGG, convert = FALSE)
    filename <- paste0(species.KEGG,"KeggLinks.Rds")
    saveRDS(GeneID.PathID, filename)
    
    PathID.PathName <- getKEGGPathwayNames(species.KEGG,
           remove.qualifier = TRUE)
    filename <- paste0(species.KEGG,"KeggPathwayNames.Rds")
    saveRDS(PathID.PathName, filename)
    # GO
    orgPkg <- paste0("org.", species, ".eg.db")
    obj <- paste0("org.", species, ".egGO2ALLEGS")
    egGO2ALLEGS <- tryCatch(getFromNamespace(obj, orgPkg), error = function(e) FALSE)
    GeneID.PathID <- AnnotationDbi::toTable(egGO2ALLEGS)[, c("gene_id",
            "go_id", "Ontology")]
    filename <- paste0(species,"GOlinks.Rds")
    saveRDS(GeneID.PathID, filename)
}

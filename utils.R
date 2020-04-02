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
    species = "Hs", species.KEGG = NULL, convert = FALSE, gene.pathway = NULL,
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

# Plot barras de GO ####################
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
# Plot barras de GOAll ####################
plotGOAll <- function(enrichdf, nrows = 30, orderby=NULL,
                      ont, genesUp = NULL, genesDown = NULL){
    require(plotly)
    require(ggplot2)
    if(!is.data.frame(enrichdf)){
        stop("enrichdf should be data.frame")
    }
    if(!exists("ont") | !(ont %in% c("BP","MF","CC")) ){
        stop("A valid value should be provided for 'ont'")
    }
    dataTitle <- list(BP=c("Biological Process", '#3e90bd', "#eb3f3f"),
                      MF=c("Molecular Function",'#3e90bd',"#eb3f3f"),
                      CC=c("Cellular Component",'#3e90bd',"#eb3f3f" ))
    enrichdf <- enrichdf[enrichdf$Ont == ont, ]
    if(!is.null(orderby)){
        orderby = match.arg(orderby, c("DE", "P.DE", "N", "Term"))
        if(orderby=="P.DE" | orderby =="Term"){
            enrichdf <- enrichdf %>% arrange(get(orderby))
        } else{ enrichdf <- enrichdf %>% arrange(desc(get(orderby)))}
    }
    enrichAll <- enrichdf
    enrichAll <- enrichAll %>% rowwise() %>%
        mutate(numUp = length(which(strsplit(genes,", ")[[1]] %in% genesUp$ENTREZID ))) %>% 
        mutate(numDown = length(which(strsplit(genes,", ")[[1]] %in% genesDown$ENTREZID ))) %>% 
        mutate(numDownNeg = -length(which(strsplit(genes,", ")[[1]] %in% genesDown$ENTREZID )))
    enrichAll <- enrichAll %>% dplyr::select(go_id, numUp, numDown, numDownNeg)
    df <- data.frame(
        Regulation = rep(c("Up", "Down"), each = nrows),
        goId = enrichAll$go_id,
        x = rep(1:nrows, 2),
        DE = c(enrichAll$numUp, enrichAll$numDown),
        DENeg = c(enrichAll$numUp, enrichAll$numDownNeg)
    )
    colorfill <- c(dataTitle[[ont]][2:3])
    r <- ggplot(df, aes(x = goId, y = DENeg, fill = Regulation)) +
        geom_bar(stat = "identity", position = "identity") + coord_flip() +
        theme(axis.text.y = element_text(angle = 0, hjust = 1)) + theme_bw() +
        scale_fill_manual(values = colorfill) +
        theme(panel.grid.major.y  = element_blank(),
              axis.title.y = element_blank())
    r <- r %>% plotly::ggplotly()
    p <- ggplot(df, aes(fill = Regulation, y = DE, x = goId)) +
        geom_bar(position = "dodge", stat = "identity") + coord_flip() +
        theme(axis.text.y = element_text(angle = 0, hjust = 1)) + theme_bw() +
        scale_fill_manual(values = colorfill) +
        theme(panel.grid.major.y  = element_blank(),
              axis.title.y = element_blank())
    p <- p %>% ggplotly()
    q <- ggplot(df, aes(fill = Regulation, y = DE, x = goId)) +
        geom_bar(position = "stack", stat = "identity") + coord_flip() +
        theme(axis.text.y = element_text(angle = 0, hjust = 1)) + theme_bw() +
        scale_fill_manual(values = colorfill) +
        theme(panel.grid.major.y  = element_blank(),
              axis.title.y = element_blank())
    q <- q %>% ggplotly()
    return(list(p,q,r) ) 
}

# Plot barras de Kegg ###########################
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

# Plot barras de KeggALL ###################
plotKeggAll <- function(enrichdf, nrows = 10, orderby = NULL, 
                        genesUp = NULL, genesDown = NULL){
    require(plotly)
    require(ggplot2)
        if(!is.data.frame(enrichdf)){
        stop("enrichdf should be data.frame")
    }
    if(!is.null(orderby)){
        orderby = match.arg(orderby, c("DE", "P.DE", "N", "Pathway"))
        if(orderby=="P.DE" | orderby =="Pathway"){
            enrichdf <- enrichdf %>% arrange(get(orderby))
        } else{ enrichdf <- enrichdf %>% arrange(desc(get(orderby)))}
    }
    enrichAll <- enrichdf
    enrichAll <- enrichAll %>% rowwise() %>%
        mutate(numUp = length(which(strsplit(genes,",")[[1]] %in% genesUp$ENTREZID ))) %>% 
        mutate(numDown = length(which(strsplit(genes,",")[[1]] %in% genesDown$ENTREZID ))) %>% 
        mutate(numDownNeg = -length(which(strsplit(genes,",")[[1]] %in% genesDown$ENTREZID )))
    enrichAll <- enrichAll %>% dplyr::select(pathID, numUp, numDown, numDownNeg)
    enrichAll$pathID <- sub(  "path:", "", enrichAll$pathID )
    df <- data.frame(
        Regulation = rep(c("Up", "Down"), each = nrows),
        pathId = enrichAll$pathID,
        x = rep(1:nrows, 2),
        DE = c(enrichAll$numUp, enrichAll$numDown),
        DENeg = c(enrichAll$numUp, enrichAll$numDownNeg)
    )
    colorfill <- c("#3e90bd","#eb3f3f")
    r <- ggplot(df, aes(x = pathId, y = DENeg, fill = Regulation)) +
        geom_bar(stat = "identity", position = "identity") + coord_flip() +
        theme(axis.text.y = element_text(angle = 0, hjust = 1)) + theme_bw() +
        scale_fill_manual(values = colorfill) +
        theme(panel.grid.major.y  = element_blank(),
              axis.title.y = element_blank())
    r <- r %>% plotly::ggplotly()
    p <- ggplot(df, aes(fill = Regulation, y = DE, x = pathId)) +
        geom_bar(position = "dodge", stat = "identity") + coord_flip() +
        theme(axis.text.y = element_text(angle = 0, hjust = 1)) + theme_bw() +
        scale_fill_manual(values = colorfill) +
        theme(panel.grid.major.y  = element_blank(),
              axis.title.y = element_blank())
    p <- p %>% ggplotly()
    q <- ggplot(df, aes(fill = Regulation, y = DE, x = pathId)) +
        geom_bar(position = "stack", stat = "identity") + coord_flip() +
        theme(axis.text.y = element_text(angle = 0, hjust = 1)) + theme_bw() +
        scale_fill_manual(values = colorfill) +
        theme(panel.grid.major.y  = element_blank(),
              axis.title.y = element_blank())
    q <- q %>% ggplotly()

    return(list(p, q, r))
}

# Función sin uso actualmente -creo- #############
loadGenes <- function(filegenes){
  load(filegenes)
  auxgenes <- genes
}

# PCA de un objeto DESeq #####################

plotPCA = function(object, intgroup = "condition", ntop = 500,
                   returnData = TRUE, labels = NULL){
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  # select the ntop genes by variance
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select, ]))
  # the contribution to   the total variance for each component
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
        name = colnames(object),
        labels = colData(object)[[labels]]
      )
  } else{
    colgroup <- factor(intgroup.df[ ,intgroup[1] ] )
    d <-
      data.frame(
        PC1 = pca$x[, 1],
        PC2 = pca$x[, 2],
        group = colgroup,
        intgroup.df,
        name = colnames(object),
        labels = colData(object)[[labels]]
      )
  }
  # assembly the data for the plot
  
  getPalette <- colorRampPalette(c("#f7837b","#1cc3c8"))
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
      #coord_fixed() +
      ggrepel::geom_text_repel(aes(label = labels, #paste("",d$name, sep = ""),
                                   size = "tam"),
                               show.legend = FALSE, size=3, nudge_y = 0.1) +
      # scale_size_manual("tam", c(1)) +
      theme(text = element_text(size=20))}
  else{
    p <- ggplot(data = d,
                aes_string(x = "PC1", y = "PC2", color = "group")) +
      geom_point(size = 3) +
      ggtitle("PCA for VST data transformation") +
      xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
      ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
      scale_color_manual(values = colours, name = intgroup[1]) +
      #coord_fixed() +
      ggrepel::geom_text_repel(aes(label = labels,# paste("",d$name, sep = ""),
                                   size = "tam"),
                               show.legend = FALSE, nudge_y = 0.1) +
      # scale_size_manual(labels = c("tam"), values = c(1)) +
      theme(text = element_text(size=20))
  }

  return(p)
}




# Función para recuperar los genes up de un objeto DEseq #############
# actualmente para p-val <= 0.05 fijo.
getSigUpregulated <- function(dds, pval=0.05, logfc=0, specie="Mm"){
  #res.sh <- lfcShrink(dds, coef=2, type="apeglm", res = results(dds))
  rk <- as.data.frame(dds)
  rk <- rk[rk$log2FoldChange >logfc & rk$padj<=pval,]
  rk <- rk[ order(rk$padj, decreasing = TRUE), ]
  annot <- geneIdConverter(rownames(rk), specie)
  return(data.frame(SYMBOL = annot$consensus, ENTREZID = annot$ENTREZID, stringsAsFactors = F) )
}

# Función para recuperar los genes down de un objeto DEseq #############
# actualmente para p-val <= 0.05 fijo.
getSigDownregulated <- function(dds, pval=0.05, logfc=0, specie="Mm"){
  #res.sh <- lfcShrink(dds, coef=2, type="apeglm", res = results(dds))
  rk <- as.data.frame(dds)
  rk <- rk[rk$log2FoldChange <logfc & rk$padj<=pval,]
  rk <- rk[ order(rk$padj, decreasing = TRUE), ]
  annot <- geneIdConverter(rownames(rk), specie)
  return(data.frame(SYMBOL = annot$consensus, ENTREZID = annot$ENTREZID, stringsAsFactors = F) )
}

# Convertidor de nombres de genes ###################
# Se le pasa un vector de ensembl y devuelve un df con varios nombres
geneIdConverter <- function(genes, specie="Mm"){ # genes = vector of ensembl gene ids (sólo para Mm por ahora)
  require("EnsDb.Mmusculus.v79")
  require("org.Mm.eg.db")
  require("EnsDb.Hsapiens.v86")
  require("org.Hs.eg.db")
  if(specie=="Mm"){
      ensdb <- EnsDb.Mmusculus.v79
      orgdb <- org.Mm.eg.db
  }
    else{
        ensdb <- EnsDb.Hsapiens.v86
        orgdb <- org.Hs.eg.db
    }
  annot <- NULL
  annot$ENSEMBL <- genes
  annot$SYMBOL <-  mapIds(ensdb, keys=genes, column="SYMBOL",keytype="GENEID")
  annot$SYMBOL1 <- mapIds(orgdb, keys = genes, column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first') 
  annot$description <- mapIds(orgdb, keys = genes, column = 'GENENAME', keytype = 'ENSEMBL', multiVals = 'first')
  annot <- as.data.frame(annot)
  consensus <- data.frame('Symbol'= ifelse(!is.na(annot$SYMBOL), as.vector(annot$SYMBOL),
                                           ifelse(!is.na(annot$SYMBOL1),as.vector(annot$SYMBOL1),
                                                  as.vector(annot$ENSEMBL))), stringsAsFactors = F)
  annot$consensus <- consensus$Symbol
  entrez1 <- mapIds(orgdb, keys = annot$consensus, column = "ENTREZID", keytype = "SYMBOL")
  entrez2 <- mapIds(orgdb, keys = as.character(annot$ENSEMBL),
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
          axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5, size=10))+
    # scale_fill_continuous(low="blue", high="red", name = "N")
    # scale_fill_brewer(palette = "YlOrRd")
    scale_fill_manual(values = getPalette(colourCount))+
    theme(text = element_text(size=16, angle=0), plot.margin = unit(c(15,25,15,15), "pt"))
  
}

# Función para crear dataset para hacer GESA pathway ##################
buildKeggDataset <- function(specie="mmu"){
  GeneID.PathID <- getGeneKEGGLinks(specie, convert = FALSE)
  PathName <- getKEGGPathwayNames(specie,remove.qualifier = TRUE)
  PathName$Id <- paste(PathName$PathwayID,PathName$Description,sep="_")
  dataSet <- left_join(GeneID.PathID, PathName, by = c("PathwayID"="PathwayID"))
  dataSet$Id <- gsub("path:","",dataSet$Id)
  dataSet <- dataSet[,c(4,1)]
  saveRDS(dataSet,paste0("resources/",specie,"keggDataGSEA.Rds"))
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

# Customized Volcano Plot ###############
CustomVolcano <- function (toptable, lab, x, y, selectLab = NULL, xlim = c(min(toptable[[x]], 
                           na.rm = TRUE), max(toptable[[x]], na.rm = TRUE)), 
                           ylim = c(0, max(-log10(toptable[[y]]), na.rm = TRUE) + 5), xlab = bquote(~Log[2] ~ "fold change"), 
                           ylab = bquote(~-Log[10] ~ italic(P)), axisLabSize = 18, 
                           title = "Volcano plot", subtitle = "", caption = paste0("Total = ", 
                           nrow(toptable), " variables"), titleLabSize = 18, subtitleLabSize = 14, 
                           captionLabSize = 14, pCutoff = 1e-05, pLabellingCutoff = pCutoff, 
                           FCcutoffDOWN = -1, FCcutoffUP = 1 , cutoffLineType = "longdash", cutoffLineCol = "black", 
                           cutoffLineWidth = 0.4, transcriptPointSize = 0.8, transcriptLabSize = 3, 
                           transcriptLabCol = "black", transcriptLabFace = "plain", 
                           transcriptLabhjust = 0, transcriptLabvjust = 1.5, pointSize = 2, 
                           labSize = 3, labCol = "black", labFace = "plain", labhjust = 0, 
                           labvjust = 1.5, boxedlabels = FALSE, boxedLabels = FALSE, 
                           shape = 19, shapeCustom = NULL, col = c("grey30", "forestgreen", 
                           "royalblue", "red2"), colCustom = NULL, colAlpha = 1/2, 
                           legend = c("NS", "Log2 FC", "P", "P & Log2 FC"), legendLabels = c("NS", 
                           expression(Log[2]~FC), "p-adjusted", expression(p-adjusted ~ and ~ log[2]~FC)), 
                           legendPosition = "top", legendLabSize = 14, 
                           legendIconSize = 4, legendVisible = TRUE, shade = NULL, 
                           shadeLabel = NULL, shadeAlpha = 1/2, shadeFill = "grey", 
                           shadeSize = 0.01, shadeBins = 2, drawconnectors = FALSE, 
                           drawConnectors = FALSE, widthConnectors = 0.5, typeConnectors = "closed", 
                           endsConnectors = "first", lengthConnectors = unit(0.01, "npc"), colConnectors = "grey10", 
                           hline = NULL, hlineType = "longdash", 
                           hlineCol = "black", hlineWidth = 0.4, vline = NULL, vlineType = "longdash", 
                           vlineCol = "black", vlineWidth = 0.4, gridlines.major = TRUE, 
                           gridlines.minor = TRUE, border = "partial", borderWidth = 0.8, 
                           borderColour = "black") 
{
  
  if (!is.numeric(toptable[[x]])) {
    stop(paste(x, " is not numeric!", sep = ""))
  }
  if (!is.numeric(toptable[[y]])) {
    stop(paste(y, " is not numeric!", sep = ""))
  }
  i <- xvals <- yvals <- Sig <- NULL
  if (!missing("transcriptPointSize")) {
    warning(paste0("transcriptPointSize argument deprecated in v1.4", 
                   " - please use pointSize"))
    pointSize <- transcriptPointSize
  }
  if (!missing("transcriptLabSize")) {
    warning(paste0("transcriptLabSize argument deprecated in v1.4", 
                   " - please use labSize"))
    labSize <- transcriptLabSize
  }
  if (!missing("transcriptLabCol")) {
    warning(paste0("transcriptLabCol argument deprecated in v1.4", 
                   " - please use labCol"))
    labCol <- transcriptLabCol
  }
  if (!missing("transcriptLabFace")) {
    warning(paste0("transcriptLabFace argument deprecated in v1.4", 
                   " - please use labFace"))
    labFace <- transcriptLabFace
  }
  if (!missing("transcriptLabhjust")) {
    warning(paste0("transcriptLabhjust argument deprecated in v1.4", 
                   " - please use labhjust"))
    labhjust <- transcriptLabhjust
  }
  if (!missing("transcriptLabvjust")) {
    warning(paste0("transcriptLabvjust argument deprecated in v1.4", 
                   " - please use labvjust"))
    labvjust <- transcriptLabvjust
  }
  if (!missing("boxedlabels")) {
    warning(paste0("boxedlabels argument deprecated in v1.4", 
                   " - please use boxedLabels"))
    boxedLabels <- boxedlabels
  }
  if (!missing("drawconnectors")) {
    warning(paste0("drawconnectors argument deprecated since v1.2", 
                   " - please use drawConnectors"))
    drawConnectors <- drawconnectors
  }
  toptable <- as.data.frame(toptable)
  toptable$Sig <- "NS"
  toptable$Sig[toptable[[x]] >= FCcutoffUP] <- "FC"
  toptable$Sig[toptable[[x]] <= FCcutoffDOWN] <- "FC"
  toptable$Sig[(toptable[[y]] < pCutoff)] <- "P"
  toptable$Sig[(toptable[[y]] < pCutoff) & (toptable[[x]] >= FCcutoffUP)] <- "FC_P"
  toptable$Sig[(toptable[[y]] < pCutoff) & (toptable[[x]] <= FCcutoffDOWN)] <- "FC_P"
  toptable$Sig <- factor(toptable$Sig, levels = c("NS", "FC", "P", "FC_P"))
  
  if (min(toptable[[y]], na.rm = TRUE) == 0) {
    warning(paste("One or more p-values is 0.", "Converting to 10^-1 * current", 
                  "lowest non-zero p-value..."), call. = FALSE)
    toptable[which(toptable[[y]] == 0), y] <- min(toptable[which(toptable[[y]] != 0), y], na.rm = TRUE) * 10^-1
  }
  toptable$lab <- lab
  toptable$xvals <- toptable[[x]]
  toptable$yvals <- toptable[[y]]
  if (!is.null(selectLab)) {
    names.new <- rep(NA, length(toptable$lab))
    indices <- which(toptable$lab %in% selectLab)
    names.new[indices] <- toptable$lab[indices]
    toptable$lab <- names.new
  }
  th <- theme_bw(base_size = 24) + theme(legend.background = element_rect(), 
          plot.title = element_text(angle = 0, size = titleLabSize, 
          face = "bold", vjust = 1), plot.subtitle = element_text(angle = 0, 
          size = subtitleLabSize, face = "plain", vjust = 1), 
          plot.caption = element_text(angle = 0, size = captionLabSize, 
          face = "plain", vjust = 1), axis.text.x = element_text(angle = 0, size = axisLabSize, vjust = 1), 
          axis.text.y = element_text(angle = 0, 
          size = axisLabSize, vjust = 1), axis.title = element_text(size = axisLabSize), 
          legend.position = legendPosition, legend.key = element_blank(), 
          legend.key.size = unit(0.5, "cm"), legend.text = element_text(size = legendLabSize), 
          title = element_text(size = legendLabSize), legend.title = element_blank())
  if (!is.null(colCustom) & !is.null(shapeCustom)) {
    plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
      th + guides(colour = guide_legend(order = 1, override.aes = list(size = legendIconSize)), 
                  shape = guide_legend(order = 2, override.aes = list(size = legendIconSize))) + 
      geom_point(aes(color = factor(names(colCustom)), 
                     shape = factor(names(shapeCustom))), alpha = colAlpha, 
                 size = pointSize, na.rm = TRUE) + scale_color_manual(values = colCustom) + 
      scale_shape_manual(values = shapeCustom)
  }
  else if (!is.null(colCustom) & is.null(shapeCustom) & length(shape) == 
           1) {
    plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
      th + guides(colour = guide_legend(order = 1, override.aes = list(size = legendIconSize)), 
                  shape = guide_legend(order = 2, override.aes = list(size = legendIconSize))) + 
      geom_point(aes(color = factor(names(colCustom))), 
                 alpha = colAlpha, shape = shape, size = pointSize, 
                 na.rm = TRUE) + scale_color_manual(values = colCustom) + 
      scale_shape_manual(guide = TRUE)
  }
  else if (!is.null(colCustom) & is.null(shapeCustom) & length(shape) == 
           4) {
    plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
      th + guides(colour = guide_legend(order = 1, override.aes = list(size = legendIconSize)), 
                  shape = guide_legend(order = 2, override.aes = list(size = legendIconSize))) + 
      geom_point(aes(color = factor(names(colCustom)), 
                     shape = factor(Sig)), alpha = colAlpha, size = pointSize, 
                 na.rm = TRUE) + scale_color_manual(values = colCustom) + 
      scale_shape_manual(values = c(NS = shape[1], FC = shape[2], 
                                    P = shape[3], FC_P = shape[4]), labels = c(NS = legendLabels[1], 
                                    FC = legendLabels[2], P = legendLabels[3], FC_P = legendLabels[4]), 
                                    guide = TRUE)
  }
  else if (is.null(colCustom) & !is.null(shapeCustom)) {
    plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
      th + guides(colour = guide_legend(order = 1, override.aes = list(size = legendIconSize)), 
                  shape = guide_legend(order = 2, override.aes = list(size = legendIconSize))) + 
      geom_point(aes(color = factor(Sig), shape = factor(names(shapeCustom))), 
                 alpha = colAlpha, size = pointSize, na.rm = TRUE) + 
      scale_color_manual(values = c(NS = col[1], FC = col[2], 
                                    P = col[3], FC_P = col[4]), labels = c(NS = legendLabels[1], 
                                                                           FC = legendLabels[2], P = legendLabels[3], FC_P = legendLabels[4])) + 
      scale_shape_manual(values = shapeCustom)
  }
  else if (is.null(colCustom) & is.null(shapeCustom) & length(shape) == 
           1) {
    plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
      th + guides(colour = guide_legend(order = 1, override.aes = list(shape = shape, 
               size = legendIconSize))) + geom_point(aes(color = factor(Sig)), 
               alpha = colAlpha, shape = shape, size = pointSize, 
               na.rm = TRUE, show.legend = legendVisible) + scale_color_manual(values = c(NS = col[1], 
               FC = col[2], P = col[3], FC_P = col[4]), labels = c(NS = legendLabels[1], 
               FC = legendLabels[2], P = legendLabels[3], FC_P = legendLabels[4]))
  }
  else if (is.null(colCustom) & is.null(shapeCustom) & length(shape) == 4) {
    plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
      th + guides(colour = guide_legend(order = 1, override.aes = list(shape = c(NS = shape[1], 
      FC = shape[2], P = shape[3], FC_P = shape[4]), size = legendIconSize))) + 
      geom_point(aes(color = factor(Sig), shape = factor(Sig)), 
                 alpha = colAlpha, size = pointSize, na.rm = TRUE, 
                 show.legend = legendVisible) + scale_color_manual(values = c(NS = col[1], 
                 FC = col[2], P = col[3], FC_P = col[4]), labels = c(NS = legendLabels[1], 
                 FC = legendLabels[2], P = legendLabels[3], FC_P = legendLabels[4])) + 
      scale_shape_manual(values = c(NS = shape[1], FC = shape[2], 
                                    P = shape[3], FC_P = shape[4]), guide = FALSE)
  }
  plot <- plot + xlab(xlab) + ylab(ylab) + xlim(xlim[1], xlim[2]) + 
    ylim(ylim[1], ylim[2]) + geom_vline(xintercept = c(FCcutoffDOWN,FCcutoffUP), linetype = cutoffLineType, colour = cutoffLineCol, 
                                        size = cutoffLineWidth) + geom_hline(yintercept = -log10(as.numeric(pCutoff)), 
                                        linetype = cutoffLineType, colour = cutoffLineCol, size = cutoffLineWidth)
  plot <- plot + labs(title = title, subtitle = subtitle, caption = caption)
  if (!is.null(vline)) {
    plot <- plot + geom_vline(xintercept = vline, linetype = vlineType, 
                              colour = vlineCol, size = vlineWidth)
  }
  if (!is.null(hline)) {
    plot <- plot + geom_hline(yintercept = -log10(hline), 
                              linetype = hlineType, colour = hlineCol, size = hlineWidth)
  }
  if (border == "full") {
    plot <- plot + theme(panel.border = element_rect(colour = borderColour, 
                                                     fill = NA, size = borderWidth))
  }
  else if (border == "partial") {
    plot <- plot + theme(axis.line = element_line(size = borderWidth, 
                                                  colour = borderColour), panel.border = element_blank(), 
                         panel.background = element_blank())
  }
  else {
    stop("Unrecognised value passed to 'border'. Must be 'full' or 'partial'")
  }
  if (gridlines.major == TRUE) {
    plot <- plot + theme(panel.grid.major = element_line())
  }
  else {
    plot <- plot + theme(panel.grid.major = element_blank())
  }
  if (gridlines.minor == TRUE) {
    plot <- plot + theme(panel.grid.minor = element_line())
  }
  else {
    plot <- plot + theme(panel.grid.minor = element_blank())
  }
  if (boxedLabels == FALSE) {
    if (drawConnectors == TRUE && is.null(selectLab)) {
      plot <- plot + geom_text_repel(data = subset(toptable, 
                  toptable[[y]] < pLabellingCutoff & ((toptable[[x]] > FCcutoffUP) | (toptable[[x]] < FCcutoffDOWN))),
                  aes(label = subset(toptable, toptable[[y]] < pLabellingCutoff & ((toptable[[x]] > FCcutoffUP) | (toptable[[x]] < FCcutoffDOWN)))[["lab"]]), 
                  size = labSize, segment.color = colConnectors, 
                  segment.size = widthConnectors, arrow = arrow(length = lengthConnectors, 
                                     type = typeConnectors, ends = endsConnectors), 
                                     hjust = labhjust, vjust = labvjust, colour = labCol, 
                                     fontface = labFace, na.rm = TRUE)
    }
    else if (drawConnectors == TRUE && !is.null(selectLab)) {
      plot <- plot + geom_text_repel(data = subset(toptable, 
                                     !is.na(toptable[["lab"]])), aes(label = subset(toptable, 
                                     !is.na(toptable[["lab"]]))[["lab"]]), size = labSize, 
                                     segment.color = colConnectors, segment.size = widthConnectors, 
                                     arrow = arrow(length = lengthConnectors, type = typeConnectors, 
                                                   ends = endsConnectors), hjust = labhjust, 
                                     vjust = labvjust, colour = labCol, fontface = labFace, 
                                     na.rm = TRUE)
    }
    else if (drawConnectors == FALSE && !is.null(selectLab)) {
      plot <- plot + geom_text(data = subset(toptable, 
                               !is.na(toptable[["lab"]])), aes(label = subset(toptable, 
                               !is.na(toptable[["lab"]]))[["lab"]]), size = labSize, 
                               check_overlap = TRUE, hjust = labhjust, vjust = labvjust, 
                               colour = labCol, fontface = labFace, na.rm = TRUE)
    }
    else if (drawConnectors == FALSE && is.null(selectLab)) {
      plot <- plot + geom_text(data = subset(toptable, 
                                             toptable[[y]] < pLabellingCutoff & ((toptable[[x]] > FCcutoffUP) | (toptable[[x]] < FCcutoffDOWN))), 
                                             aes(label = subset(toptable, toptable[[y]] < 
                                             pLabellingCutoff & ((toptable[[x]] > FCcutoffUP) | (toptable[[x]] < FCcutoffDOWN)))[["lab"]]), 
                               size = labSize, check_overlap = TRUE, hjust = labhjust, 
                               vjust = labvjust, colour = labCol, fontface = labFace, 
                               na.rm = TRUE)
    }
  }
  else {
    if (drawConnectors == TRUE && is.null(selectLab)) {
      plot <- plot + geom_label_repel(data = subset(toptable, 
                                                    toptable[[y]] < pLabellingCutoff & (toptable[[x]] > FCcutoffUP) | (toptable[[x]] < FCcutoffDOWN)), 
                                                    aes(label = subset(toptable, 
                                                    toptable[[y]] < pLabellingCutoff & ((toptable[[x]] > FCcutoffUP) | (toptable[[x]] < FCcutoffDOWN)))[["lab"]]), 
                                      size = labSize, segment.color = colConnectors, 
                                      segment.size = widthConnectors, arrow = arrow(length = lengthConnectors, 
                                                                                    type = typeConnectors, ends = endsConnectors), 
                                      hjust = labhjust, vjust = labvjust, colour = labCol, 
                                      fontface = labFace, na.rm = TRUE)
    }
    else if (drawConnectors == TRUE && !is.null(selectLab)) {
      plot <- plot + geom_label_repel(data = subset(toptable, 
                                      !is.na(toptable[["lab"]])), aes(label = subset(toptable, 
                                      !is.na(toptable[["lab"]]))[["lab"]]), size = labSize, 
                                      segment.color = colConnectors, segment.size = widthConnectors, 
                                      arrow = arrow(length = lengthConnectors, type = typeConnectors, 
                                                    ends = endsConnectors), hjust = labhjust, 
                                      vjust = labvjust, colour = labCol, fontface = labFace, 
                                      na.rm = TRUE)
    }
    else if (drawConnectors == FALSE && !is.null(selectLab)) {
      plot <- plot + geom_label(data = subset(toptable, 
                                !is.na(toptable[["lab"]])), aes(label = subset(toptable, 
                                !is.na(toptable[["lab"]]))[["lab"]]), size = labSize, 
                                hjust = labhjust, vjust = labvjust, colour = labCol, 
                                fontface = labFace, na.rm = TRUE)
    }
    else if (drawConnectors == FALSE && is.null(selectLab)) {
      plot <- plot + geom_label(data = subset(toptable, 
                                toptable[[y]] < pLabellingCutoff & ((toptable[[x]] > FCcutoffUP) | (toptable[[x]] < FCcutoffDOWN))),
                                aes(label = subset(toptable, toptable[[y]] < pLabellingCutoff & ((toptable[[x]] > FCcutoffUP) | (toptable[[x]] < FCcutoffDOWN)))[["lab"]]), 
                                size = labSize, hjust = labhjust, vjust = labvjust, 
                                colour = labCol, fontface = labFace, na.rm = TRUE)
    }
  }
  if (!is.null(shade)) {
    plot <- plot + stat_density2d(data = subset(toptable, 
                                                rownames(toptable) %in% shade), fill = shadeFill, 
                                  alpha = shadeAlpha, geom = "polygon", contour = TRUE, 
                                  size = shadeSize, bins = shadeBins, show.legend = FALSE, 
                                  na.rm = TRUE) + scale_fill_identity(name = shadeLabel, 
                                                                      labels = shadeLabel, guide = "legend")
  }
  return(plot)
}


# MA plot ##################
.levels <- function (x) 
    {
    if (!is.factor(x)) 
        x <- as.factor(x)
    levels(x)
    }

MA <- function (data, fdr = 0.05, fcDOWN = -1, fcUP = 1, genenames = NULL, detection_call = NULL, 
          size = NULL, font.label = c(12, "plain", "black"), label.rectangle = FALSE, 
          palette = c("#f7837b", "#1cc3c8", "darkgray"), top = 15, 
          select.top.method = c("padj", "fc"), main = NULL, xlab = "Log2 mean expression", 
          ylab = "Log2 fold change", ggtheme = theme_classic(), ...) 
{
  if (!base::inherits(data, c("matrix", "data.frame", "DataFrame", 
                              "DE_Results", "DESeqResults"))) 
    stop("data must be an object of class matrix, data.frame, DataFrame, DE_Results or DESeqResults")
  if (!is.null(detection_call)) {
    if (nrow(data) != length(detection_call)) 
      stop("detection_call must be a numeric vector of length = nrow(data)")
  }
  else if ("detection_call" %in% colnames(data)) {
    detection_call <- as.vector(data$detection_call)
  }
  else detection_call = rep(1, nrow(data))
  if (is.null(list(...)$legend)) 
    legend <- c(0.12, 0.9)
  ss <- base::setdiff(c("baseMean", "log2FoldChange", "padj"), 
                      colnames(data))
  if (length(ss) > 0) 
    stop("The colnames of data must contain: ", paste(ss, collapse = ", "))
  if (is.null(genenames)) 
    genenames <- rownames(data)
  else if (length(genenames) != nrow(data)) 
    stop("genenames should be of length nrow(data).")
  sig <- rep(3, nrow(data))
  sig[which(data$padj <= fdr & data$log2FoldChange < 0 & data$log2FoldChange <= (as.numeric(fcDOWN)) & detection_call == 1)] = 2
  sig[which(data$padj <= fdr & data$log2FoldChange > 0 & data$log2FoldChange >= (as.numeric(fcUP)) & detection_call == 1)] = 1
  data <- data.frame(name = genenames, mean = data$baseMean, lfc = data$log2FoldChange, padj = data$padj, sig = sig)
  . <- NULL
  data$sig <- as.factor(data$sig)
  .lev <- .levels(data$sig) %>% as.numeric()
  palette <- palette[.lev]
  new.levels <- c(paste0("Up: ", sum(sig == 1)), paste0("Down: ", sum(sig == 2)), "NS") %>% .[.lev]
  data$sig <- factor(data$sig, labels = new.levels)
  select.top.method <- match.arg(select.top.method)
  if (select.top.method == "padj") 
    data <- data[order(data$padj), ]
  else if (select.top.method == "fc") 
    data <- data[order(abs(data$lfc), decreasing = TRUE), 
                 ]
  labs_data <- stats::na.omit(data)
  labs_data <- subset(labs_data, padj <= fdr & name != "" & 
                        (lfc >= fcUP | lfc <=fcDOWN) )
  labs_data <- utils::head(labs_data, top)
  font.label <- ggpubr:::.parse_font(font.label)
  font.label$size <- ifelse(is.null(font.label$size), 12, 
                            font.label$size)
  font.label$color <- ifelse(is.null(font.label$color), "black", 
                             font.label$color)
  font.label$face <- ifelse(is.null(font.label$face), "plain", 
                            font.label$face)
  set.seed(42)
  mean <- lfc <- sig <- name <- padj <- NULL
  p <- ggplot(data, aes(x = log2(mean + 1), y = lfc)) + geom_point(aes(color = sig), size = size)
  if (label.rectangle) {
    p <- p + ggrepel::geom_label_repel(data = labs_data, 
                                       mapping = aes(label = name), box.padding = unit(0.35,"lines"), point.padding = unit(0.3, "lines"), 
                                       force = 1, fontface = font.label$face, size = font.label$size/3, 
                                       color = font.label$color)
  }
  else {
    p <- p + ggrepel::geom_text_repel(data = labs_data, 
                                      mapping = aes(label = name), box.padding = unit(0.35,  "lines"), point.padding = unit(0.3, "lines"), 
                                      force = 1, fontface = font.label$face, size = font.label$size/3, 
                                      color = font.label$color)
  }
  p <- p + scale_x_continuous(breaks = seq(0, max(log2(data$mean + 
                                                         1)), 2)) + labs(x = xlab, y = ylab, title = main, color = "") + 
    geom_hline(yintercept = c(0, (fcDOWN), (fcUP)), linetype = c(1, 2, 2), color = c("black", "black", "black"))
  p <- ggpar(p, palette = palette, ggtheme = ggtheme, ...)
  return(p)
}


# VST ###############  SIN USO CREO !!!!!
VST <- function (object, blind = TRUE, nsub = 1000, fitType = "parametric") 
{
  if (nrow(object) < nsub) {
    stop("less than 'nsub' rows,\n  it is recommended to use varianceStabilizingTransformation directly")
  }
  if (is.null(colnames(object))) {
    colnames(object) <- seq_len(ncol(object))
  }
  if (is.matrix(object)) {
    matrixIn <- TRUE
    object <- DESeqDataSetFromMatrix(object, DataFrame(row.names = colnames(object)), 
                                     ~1)
  }
  else {
    if (blind) {
      design(object) <- ~1
    }
    matrixIn <- FALSE
  }
  if (is.null(sizeFactors(object)) & is.null(normalizationFactors(object))) {
    object <- estimateSizeFactors(object)
  }
  baseMean <- rowMeans(counts(object, normalized = TRUE))
  if (sum(baseMean > 5) < nsub) {
    stop("less than 'nsub' rows with mean normalized count > 5, \n  it is recommended to use varianceStabilizingTransformation directly")
  }
  object.sub <- object[baseMean > 5, ]
  baseMean <- baseMean[baseMean > 5]
  o <- order(baseMean)
  idx <- o[round(seq(from = 1, to = length(o), length = nsub))]
  object.sub <- object.sub[idx, ]
  object.sub <- estimateDispersionsGeneEst(object.sub, quiet = TRUE)
  object.sub <- estimateDispersionsFit(object.sub, fitType = fitType, 
                                       quiet = TRUE)
  suppressMessages({
    dispersionFunction(object) <- dispersionFunction(object.sub)
  })
  vsd <- varianceStabilizingTransformation(object, blind = FALSE)
  if (matrixIn) {
    return(assay(vsd))
  }
  else {
    return(vsd)
  }
}

# Heatmap #############

heat <- function (vsd, n = 40, intgroup = "condition", sampleName = "condition", specie="Mm") 
{
  require("EnsDb.Mmusculus.v79")
  require("org.Mm.eg.db")
  require("EnsDb.Hsapiens.v86")
  require("org.Hs.eg.db")
  if(specie=="Mm"){
    ensdb <- EnsDb.Mmusculus.v79
    orgdb <- org.Mm.eg.db
  }
  else{
    ensdb <- EnsDb.Hsapiens.v86
    orgdb <- org.Hs.eg.db
  }
#vsd <- vst(data)
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), n)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
if (!all(intgroup %in% names(colData(vsd)))) {
  stop("the argument 'intgroup' should specify columns of colData(dds)")
}
df <- as.data.frame(colData(vsd)[, intgroup, drop = FALSE])

annot <- NULL
annot$ENSEMBL <- rownames(mat)
annot$SYMBOL <-  mapIds(ensdb, keys=rownames(mat), column="SYMBOL",keytype="GENEID")
annot$SYMBOL1 <- mapIds(orgdb, keys = rownames(mat), column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first') 
annot$description <- mapIds(orgdb, keys = rownames(mat), column = 'GENENAME', keytype = 'ENSEMBL', multiVals = 'first')
annot <- as.data.frame(annot)
consensus <- data.frame('Symbol'= ifelse(!is.na(annot$SYMBOL), as.vector(annot$SYMBOL),
                                         ifelse(!is.na(annot$SYMBOL1),as.vector(annot$SYMBOL1),
                                                as.vector(annot$ENSEMBL))), stringsAsFactors = F)

pheatmap(mat, cluster_rows=TRUE, cluster_cols=TRUE,
         show_colnames=TRUE, show_rownames = TRUE, annotation_col = df,
         labels_col = as.character(vsd[[sampleName]]),
         labels_row = as.character(consensus$Symbol),
         main = "Heatmap top genes")
}


# cluster #############

cluster <- function(vsd, intgroup = "condition")
  {
  #vsd <- vst(data)
  sampleDists_vsd <- dist(t(assay(vsd)))
  sampleDistMatrix_vsd <- as.matrix( sampleDists_vsd )
  rownames(sampleDistMatrix_vsd) <- vsd[[intgroup]]
  colnames(sampleDistMatrix_vsd) <- vsd[[intgroup]]
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix_vsd,
           clustering_distance_rows = sampleDists_vsd,
           clustering_distance_cols = sampleDists_vsd,
           col = colors, main = 'Heatmap clustering')
}






#############  TOP6 genes #########################


plotCountsSymbol <- function (dds, gene, res, intgroup = "condition", normalized = TRUE,
                              transform = TRUE, main, xlab = "group", returnData = FALSE,
                              replaced = FALSE, pc, ...){
  stopifnot(length(gene) == 1 & (is.character(gene) | (is.numeric(gene) &
                                                         (gene >= 1 & gene <= nrow(dds)))))
  if (!all(intgroup %in% names(colData(dds))))
    stop("all variables in 'intgroup' must be columns of colData")
  if (!returnData) {
    if (!all(sapply(intgroup, function(v) is(colData(dds)[[v]],
                                             "factor")))) {
      stop("all variables in 'intgroup' should be factors, or choose returnData=TRUE and plot manually")
    }
  }
  if (missing(pc)) {
    pc <- if (transform)
      0.5
    else 0
  }
  if (is.null(sizeFactors(dds)) & is.null(normalizationFactors(dds))) {
    dds <- estimateSizeFactors(dds)
  }
  cnts <- counts(dds, normalized = normalized, replaced = replaced)[gene,]
  group <- if (length(intgroup) == 1) {
    colData(dds)[[intgroup]]
  }
  else if (length(intgroup) == 2) {
    lvls <- as.vector(t(outer(levels(colData(dds)[[intgroup[1]]]),
                              levels(colData(dds)[[intgroup[2]]]), function(x,
                                                                            y) paste(x, y, sep = ":"))))
    droplevels(factor(apply(as.data.frame(colData(dds)[,
                                                       intgroup, drop = FALSE]), 1, paste, collapse = ":"),
                      levels = lvls))
  }
  else {
    factor(apply(as.data.frame(colData(dds)[, intgroup,
                                            drop = FALSE]), 1, paste, collapse = ":"))
  }
  #rownames(cnts) <- res$GeneName_Symbol
  data <- data.frame(count = cnts + pc, group = as.integer(group))
  logxy <- if (transform)
    "y"
  else ""
  ylab <- ifelse(normalized, "normalized counts")
  #colors = c("#008000","#800080")
  if (returnData)
    return(data.frame(count = data$count, colData(dds)[intgroup]))
  plot(data$group + runif(ncol(dds), -0.05, 0.05), data$count, #col=colors[(dds)[intergroup]],
       xlim = c(0.5, max(data$group) + 0.5), log = logxy, xaxt = "n",
       xlab = xlab, ylab = ylab, main = "Top 6 significant gene", ...)
  axis(1, at = seq_along(levels(group)), levels(group))
  #text(data$group + runif(ncol(dds), -0.05, 0.05), data$count, labels=colnames(dds))
}






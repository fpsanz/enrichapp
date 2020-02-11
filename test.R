library("limma") #dependencia de customkegg y customgo
library("purrr") #dependencia de datatable2
library("DT") #dependencia de datatable2
library("tidyverse")
library("GO.db")
genes <- read.csv("genes.csv", row.names = 1, stringsAsFactors = FALSE)

kgg <- customKegg(genes, species="Mm",species.KEGG = "mmu")
kggDT <- kegg2DT(enrichdf = kgg, data = genes, nrows = 10)
datatable2(kggDT, vars = c("genes"),
           escape = FALSE,
           filter = list(position="top", clear=FALSE),
           opts = list(pageLength = 10, white_space = "normal"))
plotKegg(enrichdf = kgg, nrows = 10)

gos <- customGO(genes, species = "Mm")
goDT <- go2DT(enrichdf = gos, data = genes, nrows = 1000)
datatable2(goDT, vars = c("genes"),
           escape = FALSE,
           filter = list(position="top", clear=FALSE),
           opts = list(pageLength = 10, white_space = "normal"))

plotGO(enrichdf = gos, nrows = 10, ont="BP")

chordKgg <- chordPlot(kgg, nRows = 10, orderby = "P.DE")
chordGO <- chordPlot(gos, nRows = 10, ont="MF", orderby = "P.DE")



kggDT %>%
    filter(Pathway %in% rutas) %>%
    dplyr::select(Pathway, genes,DE) %>%
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
          axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5, size=5))+
    # scale_fill_continuous(low="blue", high="red", name = "N")
    # scale_fill_brewer(palette = "YlOrRd")
    scale_fill_manual(values = getPalette(colourCount))



gos

library(GO.db)
getAllBPChildren <- function(goids)
{
    ans <- unique(unlist(mget(goids, GOBPCHILDREN), use.names=FALSE))
    ans <- ans[!is.na(ans)]
    return(ans)
}



saveRDS(GOlevel, "GOlevel.Rds")
gos_level <- left_join(gos, GOlevel, by = c("go_id"="id"))


gokk <- gos %>% filter(go_id %in% kk)





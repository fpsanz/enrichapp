load(filegenes)
auxgenes <- genes
gos <- readRDS("gos.Rds")
auxkgg <- readRDS("kgg.Rds")
goDT <- go2DT(enrichdf = gos, data = auxgenes)

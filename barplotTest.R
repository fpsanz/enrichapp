dds <- readRDS("deseq.Rds")
library(DESeq2)
library(dplyr)
res <- as.data.frame(lfcShrink(dds, coef=2, type="apeglm", res = results(dds)))
        conversion <- geneIdConverter(rownames(res), "Mm" )
        res$baseMean <- round(res$baseMean,4)
        res$lfcSE <- round(res$lfcSE,4)
        res$log2FoldChange <- round(res$log2FoldChange,4)
        res <- cbind(`Description`=conversion$description, res)
        res <- cbind(`GeneName_Symbol`=conversion$consensus, res)
        res <-  res %>% dplyr::select(-c(pvalue))
        spc = "Mus_musculus"
        
        links = paste0("<a href='http://www.ensembl.org/",spc,",/Gene/Summary?db=core;g=",
                       rownames(res),"' target='_blank'>",rownames(res),"</a>")
        res <- cbind(`GeneName_Ensembl`= links, res)
genesUp <- getSigUpregulated(res, 0.05, 0.5, "Mm" ) 
    genesDown <- getSigDownregulated(res, 0.05, -0.05, "Mm" ) 
    genesall <- rbind(genesUp, genesDown)
    library(limma)
go <- customGO(genesall, species = "Mm")

##nuevo
go2 <- go %>% rowwise() %>% mutate(numUp = length(which(strsplit(genes,", ")[[1]] %in% genesUp$ENTREZID ))) %>% 
    mutate(numDown = -length(which(strsplit(genes,", ")[[1]] %in% genesDown$ENTREZID )))
go2 <- go2[1:10, c(1,5,6,9,10)]
library(ggplot2)
library(plotly)
#df <- data.frame(group = rep(c("Above", "Below"), each=10), x = rep(1:10, 2), y = c(runif(10, 0, 1), runif(10, -1, 0)))

df <- data.frame(Regulation = rep(c("Up", "Down"), each=10), 
                 pval=go2$P.DE, goid = go2$go_id,
                 x = rep(1:10, 2), 
                 DE = c(go2$numUp,go2$numDown))
colorfill <- c("#e27d60", "#41B3A3")
p <- ggplot(df, aes(x=goid, y=DE, fill=Regulation ) ) + 
	    geom_bar(stat="identity", position="identity") + coord_flip()+ 
            theme(axis.text.y=element_text(angle=0, hjust=1))+
        theme(axis.title.y = element_blank()) + theme(legend.position = "none") +
    scale_fill_manual(values=colorfill)
print(p) 
p %>% plotly::ggplotly() %>% style( hoverlabel = list(bgcolor="#000"))



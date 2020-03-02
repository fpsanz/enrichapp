# Variables enrichApp

## Pestaña preview

| Elemento         | Variables   |
| ---------------- | ----------- |
| tabla samples    | **samples** |
| tabla resultados | **preview** |
| gráfico PCA      | **pca**     |

## Pestaña kegg

| Elemento  | All_together       | upRegulated       | downregulated       |
| --------- | ------------------ | ----------------- | ------------------- |
| Tabla     | **tableAll**       | **table**         | **tableDown**       |
| Barplot   | **keggPlotAll**    | **keggPlot**      | **keggPlotDown**    |
| Chordplot | **keggChordAll**   | **keggChord**     | **keggChordDown**   |
| Dotplot   | **keggDotAll**     | **keggDotUp**     | **keggDotDown**     |
| HeatMap   | **heatmapKeggAll** | **heatmapKeggUp** | **heatmapKeggDown** |
| NetPlot   | **cnetKeggAll**    | **cnetKeggUp**    | **cnetKeggDown**    |
|           |                    |                   |                     |

## Pestaña GO BP

| ELEMENTO | all_together | upregulated | downregulated   |
| -------- | ------------ | ----------- | --------------- |
| Tabla    |              | **tableBP** | **tableBPdown** |
| Barplot  |              | **plotBP**  | **plotBPdown**  |
| Dotplot  |              | **BPDotUp** | **BPDotdown**   |

## Pestaña GO MF

| ELEMENTO | all_together | upregulated | downregulated   |
| -------- | ------------ | ----------- | --------------- |
| Tabla    |              | **tableMF** | **tableMFdown** |
| Barplot  |              | **plotMF**  | **plotMFdown**  |
| Dotplot  |              | **MFDotUp** | **MFDotdown**   |

## Pestaña GO CC

| ELEMENTO | all_together | upregulated | downregulated   |
| -------- | ------------ | ----------- | --------------- |
| Tabla    |              | **tableCC** | **tableCCdown** |
| Barplot  |              | **plotCC**  | **plotCCdown**  |
| Dotplot  |              | **CCDotUp** | **CCDotdown**   |



## Pestaña GSEA

| ELEMENTO | All genes deseq |
| -------- | --------------- |
| Tabla    | **gseaTable**   |
| Barplot  | **gseaPlot**    |



## Variables observeEvent

| variable | valores                         | contenido              |
| -------- | ------------------------------- | ---------------------- |
| data     | $genesUp, $genesDown, $genesall | listados de genes      |
| goDT     | $up, $down                      | pretabla GO            |
| go       | $up, $down                      | tabla enrich GO        |
| kgg      | $up, $down, $all                | tabla enrich Kegg      |
| kggDT    | $up, $down. $all                | pretabla enrich Kegg   |
| datos    | $dds                            | objeto DESeq importado |
| gsea     | $gsea                           | objeto GSEA            |

## Variable reactive

| Variable   | contenido                          |
| ---------- | ---------------------------------- |
| rowsAll    | filas seleccionadas de tableAll    |
| rows       | filas seleccionadas de table       |
| rowdown    | filas seleccionadas de tableDown   |
| bprows     | filas seleccionadas de tableBP     |
| mfrows     | filas seleccionadas de tableMF     |
| ccrows     | filas seleccionadas de tableCC     |
| bprowndown | filas seleccionadas de tableBPdown |
| mfrowsdown | filas seleccionadas de tableMFdown |
| ccrowsdown | filas seleccionadas de tableCCdown |
| variables  | Variables seleccionadas para PCA   |
| gsearow    | filas seleccionadas para GSEA      |



## Diario de acciones

01/03/2020: FPS: 

* Arreglado tabla samples para que no salgan las dos últimas columnas
* Arreglado tabla results para que salga con notación científica p-val ordenable
* Añadida y funcional, pestaña kegg: All DE genes
* Centrado boton download report

## TODOs

* Hacer all together en GO
* Poner en marcha report
* Sugerencias??¿??¿?
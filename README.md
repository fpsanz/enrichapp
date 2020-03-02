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
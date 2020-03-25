
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
                  toptable[[y]] < pLabellingCutoff & ((toptable[[x]] > FCcutoffUP) & (toptable[[x]] < FCcutoffDOWN))),
                  aes(label = subset(toptable, toptable[[y]] < pLabellingCutoff & ((toptable[[x]] > FCcutoffUP) & (toptable[[x]] < FCcutoffDOWN)))[["lab"]]), 
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
                                                    toptable[[y]] < pLabellingCutoff & (toptable[[x]] > FCcutoffUP) & (toptable[[x]] < FCcutoffDOWN)), 
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
                                toptable[[y]] < pLabellingCutoff & ((toptable[[x]] > FCcutoffUP) & (toptable[[x]] < FCcutoffDOWN))),
                                aes(label = subset(toptable, toptable[[y]] < pLabellingCutoff & ((toptable[[x]] > FCcutoffUP) & (toptable[[x]] < FCcutoffDOWN)))[["lab"]]), 
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

print("hkk")
p <- CustomVolcano(res.sh, x = 'log2FoldChange', y = "padj", lab = res.sh$GeneName_Symbol)
print(p)


  if (!requireNamespace('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

  BiocManager::install('EnhancedVolcano')


devtools::install_github('kevinblighe/EnhancedVolcano')


library(airway)
library(magrittr)
library(EnhancedVolcano)
library(DESeq2)
library(pasilla)
library(gridExtra)
library(grid)

data('airway')
  airway$dex %<>% relevel('untrt')



  dds <- DESeqDataSet(airway, design = ~ cell + dex)
  dds <- DESeq(dds, betaPrior=FALSE)
  res1 <- results(dds,
    contrast = c('dex','trt','untrt'))
  res1 <- lfcShrink(dds,
    contrast = c('dex','trt','untrt'), res=res1, type = 'normal')
  res2 <- results(dds,
    contrast = c('cell', 'N061011', 'N61311'))
  res2 <- lfcShrink(dds,
    contrast = c('cell', 'N061011', 'N61311'), res=res2, type = 'normal')


EnhancedVolcano(res1,
    lab = rownames(res1),
    x = 'log2FoldChange',
    y = 'pvalue',
    xlim = c(-5, 8))

#Modify cut-offs for log2FC and P value; specify title; adjust point and label size

 EnhancedVolcano(res2,
    lab = rownames(res2),
    x = 'log2FoldChange',
    y = 'pvalue',
    xlim = c(-8, 8),
    title = 'N061011 versus N61311',
    pCutoff = 10e-16,
    FCcutoff = 1.5,
    pointSize = 3.0,
    labSize = 3.0)

#Adjust colour and alpha for point shading

  EnhancedVolcano(res2,
    lab = rownames(res2),
    x = 'log2FoldChange',
    y = 'pvalue',
    xlim = c(-8, 8),
    title = 'N061011 versus N61311',
    pCutoff = 10e-16,
    FCcutoff = 1.5,
    pointSize = 3.0,
    labSize = 3.0,
    col=c('black', 'black', 'black', 'red3'),
    colAlpha = 1)


#Adjust shape of plotted points

 EnhancedVolcano(res2,
    lab = rownames(res2),
    x = 'log2FoldChange',
    y = 'pvalue',
    xlim = c(-8, 8),
    title = 'N061011 versus N61311',
    pCutoff = 10e-16,
    FCcutoff = 1.5,
    pointSize = 4.0,
    labSize = 3.0,
    shape = 8,
    colAlpha = 1)




#or

  EnhancedVolcano(res2,
    lab = rownames(res2),
    x = 'log2FoldChange',
    y = 'pvalue',
    xlim = c(-8, 8),
    title = 'N061011 versus N61311',
    pCutoff = 10e-16,
    FCcutoff = 1.5,
    pointSize = 3.0,
    labSize = 3.0,
    shape = c(1, 4, 23, 25),
    colAlpha = 1)

#Adjust cut-off lines and add extra threshold lines


EnhancedVolcano(res2,
    lab = rownames(res2),
    x = 'log2FoldChange',
    y = 'pvalue',
    xlim = c(-6, 6),
    title = 'N061011 versus N61311',
    pCutoff = 10e-12,
    FCcutoff = 1.5,
    pointSize = 3.0,
    labSize = 3.0,
    colAlpha = 1,
    cutoffLineType = 'blank',
    cutoffLineCol = 'black',
    cutoffLineWidth = 0.8,
    hline = c(10e-12, 10e-36, 10e-60, 10e-84),
    hlineCol = c('grey0', 'grey25','grey50','grey75'),
    hlineType = 'longdash',
    hlineWidth = 0.8,
    gridlines.major = FALSE,
    gridlines.minor = FALSE)


#Adjust legend position, size, and text

EnhancedVolcano(res2,
    lab = rownames(res2),
    x = 'log2FoldChange',
    y = 'pvalue',
    xlim = c(-6, 6),
    pCutoff = 10e-12,
    FCcutoff = 1.5,
    cutoffLineType = 'twodash',
    cutoffLineWidth = 0.8,
    pointSize = 4.0,
    labSize = 4.0,
    colAlpha = 1,
    legendLabels=c('Not sig.','Log (base 2) FC','p-value',
      'p-value & Log (base 2) FC'),
    legendPosition = 'right',
    legendLabSize = 16,
    legendIconSize = 5.0)



#Fit more labels by adding connectors


EnhancedVolcano(res2,
    lab = rownames(res2),
    x = 'log2FoldChange',
    y = 'pvalue',
    xlim = c(-6,6),
    xlab = bquote(~Log[2]~ 'fold change'),
    pCutoff = 10e-14,
    FCcutoff = 2.0,
    pointSize = 4.0,
    labSize = 4.0,
    colAlpha = 1,
    legendPosition = 'right',
    legendLabSize = 12,
    legendIconSize = 4.0,
    drawConnectors = TRUE,
    widthConnectors = 0.2,
    colConnectors = 'grey30')


#Only label key variables

EnhancedVolcano(res2,
    lab = rownames(res2),
    x = 'log2FoldChange',
    y = 'pvalue',
    selectLab = c('ENSG00000106565','ENSG00000187758'),
    xlim = c(-6,7),
    xlab = bquote(~Log[2]~ 'fold change'),
    pCutoff = 10e-14,
    FCcutoff = 2.0,
    pointSize = 4.0,
    labSize = 5.0,
    shape = c(4, 35, 17, 18),
    colAlpha = 1,
    legendPosition = 'right',
    legendLabSize = 14,
    legendIconSize = 5.0)

#Draw labels in boxes
EnhancedVolcano(res2,
    lab = rownames(res2),
    x = 'log2FoldChange',
    y = 'pvalue',
    selectLab = c('ENSG00000106565','ENSG00000187758',
      'ENSG00000230795', 'ENSG00000164530',
      'ENSG00000143153'),
    xlim = c(-5.5,8),
    xlab = bquote(~Log[2]~ 'fold change'),
    pCutoff = 10e-14,
    FCcutoff = 2.0,
    pointSize = 4.0,
    labSize = 5.0,
    labCol = 'black',
    labFace = 'bold',
    boxedLabels = TRUE,
    colAlpha = 4/5,
    legendPosition = 'right',
    legendLabSize = 14,
    legendIconSize = 4.0,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black')


#Over-ride colouring scheme with custom key-value pairs

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
  # this can be achieved with nested ifelse statements
  keyvals <- ifelse(
    res2$log2FoldChange < -2.5, 'royalblue',
      ifelse(res2$log2FoldChange > 2.5, 'gold',
        'black'))
  keyvals[is.na(keyvals)] <- 'black'
  names(keyvals)[keyvals == 'gold'] <- 'high'
  names(keyvals)[keyvals == 'black'] <- 'mid'
  names(keyvals)[keyvals == 'royalblue'] <- 'low'


 p1 <- EnhancedVolcano(res2,
    lab = rownames(res2),
    x = 'log2FoldChange',
    y = 'pvalue',
    selectLab = rownames(res2)[which(names(keyvals) %in% c('high', 'low'))],
    xlim = c(-6.5,6.5),
    xlab = bquote(~Log[2]~ 'fold change'),
    title = 'Custom colour over-ride',
    pCutoff = 10e-14,
    FCcutoff = 1.0,
    pointSize = 4.5,
    labSize = 4.5,
    shape = c(6, 4, 2, 11),
    colCustom = keyvals,
    colAlpha = 1,
    legendPosition = 'left',
    legendLabSize = 15,
    legendIconSize = 5.0,
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    colConnectors = 'grey50',
    gridlines.major = TRUE,
    gridlines.minor = FALSE,
    border = 'partial',
    borderWidth = 1.5,
    borderColour = 'black')

  p2 <- EnhancedVolcano(res2,
    lab = rownames(res2),
    x = 'log2FoldChange',
    y = 'pvalue',
    selectLab = rownames(res2)[which(names(keyvals) %in% c('high', 'low'))],
    xlim = c(-6.5,6.5),
    xlab = bquote(~Log[2]~ 'fold change'),
    title = 'No custom colour over-ride',
    pCutoff = 10e-14,
    FCcutoff = 1.0,
    pointSize = 4.5,
    labSize = 4.5,
    colCustom = NULL,
    colAlpha = 1,
    legendPosition = 'right',
    legendLabSize = 15,
    legendIconSize = 5.0,
    drawConnectors = FALSE,
    widthConnectors = 0.5,
    colConnectors = 'grey50',
    gridlines.major = TRUE,
    gridlines.minor = FALSE,
    border = 'full',
    borderWidth = 1.0,
    borderColour = 'black')

  library(gridExtra)
  library(grid)
  grid.arrange(p1, p2,
    ncol=2,
    top = textGrob('EnhancedVolcano',
      just = c('center'),
      gp = gpar(fontsize = 32)))


#Over-ride colour and/or shape scheme with custom key-value pairs

# define different cell-types that will be shaded
  celltype1 <- c('ENSG00000106565', 'ENSG00000002933',
    'ENSG00000165246', 'ENSG00000224114')
  celltype2 <- c('ENSG00000230795', 'ENSG00000164530',
    'ENSG00000143153', 'ENSG00000169851',
    'ENSG00000231924', 'ENSG00000145681')

  # create custom key-value pairs for different cell-types
  # this can be achieved with nested ifelse statements
  keyvals.shape <- ifelse(
    rownames(res2) %in% celltype1, 17,
      ifelse(rownames(res2) %in% celltype2, 64,
        3))
  keyvals.shape[is.na(keyvals.shape)] <- 3
  names(keyvals.shape)[keyvals.shape == 3] <- 'PBMC'
  names(keyvals.shape)[keyvals.shape == 17] <- 'Cell-type 1'
  names(keyvals.shape)[keyvals.shape == 64] <- 'Cell-type 2'


  p1 <- EnhancedVolcano(res2,
    lab = rownames(res2),
    x = 'log2FoldChange',
    y = 'pvalue',
    selectLab = rownames(res2)[which(names(keyvals) %in% c('high', 'low'))],
    xlim = c(-6.5,6.5),
    xlab = bquote(~Log[2]~ 'fold change'),
    title = 'Custom shape over-ride',
    pCutoff = 10e-14,
    FCcutoff = 1.0,
    pointSize = 4.5,
    labSize = 4.5,
    shapeCustom = keyvals.shape,
    colCustom = NULL,
    colAlpha = 1,
    legendLabSize = 15,
    legendPosition = 'left',
    legendIconSize = 5.0,
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    colConnectors = 'grey50',
    gridlines.major = TRUE,
    gridlines.minor = FALSE,
    border = 'partial',
    borderWidth = 1.5,
    borderColour = 'black')

  # create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
  # this can be achieved with nested ifelse statements
  keyvals.colour <- ifelse(
    res2$log2FoldChange < -2.5, 'royalblue',
      ifelse(res2$log2FoldChange > 2.5, 'gold',
        'black'))
  keyvals.colour[is.na(keyvals.colour)] <- 'black'
  names(keyvals.colour)[keyvals.colour == 'gold'] <- 'high'
  names(keyvals.colour)[keyvals.colour == 'black'] <- 'mid'
  names(keyvals.colour)[keyvals.colour == 'royalblue'] <- 'low'

  p2 <- EnhancedVolcano(res2,
    lab = rownames(res2),
    x = 'log2FoldChange',
    y = 'pvalue',
    selectLab = rownames(res2)[which(names(keyvals) %in% c('High', 'Low'))],
    xlim = c(-6.5,6.5),
    xlab = bquote(~Log[2]~ 'fold change'),
    title = 'Custom shape & colour over-ride',
    pCutoff = 10e-14,
    FCcutoff = 1.0,
    pointSize = 5.5,
    labSize = 0.0,
    shapeCustom = keyvals.shape,
    colCustom = keyvals.colour,
    colAlpha = 1,
    legendPosition = 'right',
    legendLabSize = 15,
    legendIconSize = 5.0,
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    colConnectors = 'grey50',
    gridlines.major = TRUE,
    gridlines.minor = FALSE,
    border = 'full',
    borderWidth = 1.0,
    borderColour = 'black')

  library(gridExtra)
  library(grid)
  grid.arrange(p1, p2,
    ncol=2,
    top = textGrob('EnhancedVolcano',
      just = c('center'),
      gp = gpar(fontsize = 32)))





#Highlighting key variables via custom point sizes

  library("pasilla")
  pasCts <- system.file("extdata", "pasilla_gene_counts.tsv",
    package="pasilla", mustWork=TRUE)
  pasAnno <- system.file("extdata", "pasilla_sample_annotation.csv",
    package="pasilla", mustWork=TRUE)
  cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
  coldata <- read.csv(pasAnno, row.names=1)
  coldata <- coldata[,c("condition","type")]
  rownames(coldata) <- sub("fb", "", rownames(coldata))
  cts <- cts[, rownames(coldata)]
  library("DESeq2")
  dds <- DESeqDataSetFromMatrix(countData = cts,
    colData = coldata,
    design = ~ condition)

  featureData <- data.frame(gene=rownames(cts))
  mcols(dds) <- DataFrame(mcols(dds), featureData)
  dds <- DESeq(dds)
  res <- results(dds)

  p1 <- EnhancedVolcano(res,
    lab = rownames(res),
    x = "log2FoldChange",
    y = "pvalue",
    pCutoff = 10e-4,
    FCcutoff = 2,
    xlim = c(-5.5, 5.5),
    ylim = c(0, -log10(10e-12)),
    pointSize = c(ifelse(res$log2FoldChange>2, 8, 1)),
    labSize = 4.0,
    shape = c(6, 6, 19, 16),
    title = "DESeq2 results",
    subtitle = "Differential expression",
    caption = bquote(~Log[2]~ "fold change cutoff, 2; p-value cutoff, 10e-4"),
    legendPosition = "right",
    legendLabSize = 14,
    col = c("grey30", "forestgreen", "royalblue", "red2"),
    colAlpha = 0.9,
    drawConnectors = TRUE,
    hline = c(10e-8),
    widthConnectors = 0.5)

  p1


#Change to continuous colour scheme

  p1 <- EnhancedVolcano(res,
    lab = rownames(res),
    x = "log2FoldChange",
    y = "pvalue",
    pCutoff = 10e-4,
    FCcutoff = 2,
    xlim = c(-5.5, 5.5),
    ylim = c(0, -log10(10e-12)),
    pointSize = c(ifelse(res$log2FoldChange>2, 8, 1)),
    labSize = 4.0,
    shape = c(6, 6, 19, 16),
    title = "DESeq2 results",
    subtitle = "Differential expression",
    caption = bquote(~Log[2]~ "fold change cutoff, 2; p-value cutoff, 10e-4"),
    legendPosition = "right",
    legendLabSize = 14,
    colAlpha = 0.9,
    colGradient = c('red3', 'royalblue'),
    drawConnectors = TRUE,
    hline = c(10e-8),
    widthConnectors = 0.5)

  p1


#Custom axis tick marks
  p1 +
    ggplot2::coord_cartesian(xlim=c(-6, 6)) +
    ggplot2::scale_x_continuous(
      breaks=seq(-6,6, 1))



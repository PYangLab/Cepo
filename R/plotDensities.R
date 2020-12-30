#' Plot densities
#' @return A \code{\link{ggplot}} object
#' with cell-type specific densities for a gene.
#' @param x a \code{\linkS4class{SummarizedExperiment}} or a \code{\linkS4class{SingleCellExperiment}} object.
#' @param cepoOutput an output from Cepo or doLimma/doVoom/doTtest/doWilcoxon functions
#' @param assay a character ('logcounts' by default),
#' indicating the name of the assays(x) element which stores the expression data (i.e., assays(x)$name_assays_expression).
#' We strongly encourage using normalized data, such as counts per million (CPM) or log-CPM.
#' @param celltype a character, indicating the name of the cell type to plot. Default is NULL which selects all celltypes in the cepoOutput.
#' @param celltypeColumn a character, indicating the name of the name of the cell type column in the colData(x).
#' @param genes a character vector, indicating the name of the genes to plot. Default to NULL, so that 2 top genes from each celltype will be plotted.
#' @param nGenes number of top genes from each celltype to plot. Default to 2.
#' @param color a named color vector. The names should correspond to the `celltype` argument above
#' @param plotType Either 'histogram' or 'density'
#' @import ggplot2
#' @importFrom methods is
#' @importFrom SingleCellExperiment colData
#' @importFrom SummarizedExperiment assays
#' @importFrom grDevices rainbow
#' @importFrom patchwork plot_layout wrap_plots
#' @importFrom rlang .data
#' @importFrom reshape2 melt
#' @return A \code{\link{ggplot}} object.
#' @export
#' @examples
#' library(SingleCellExperiment)
#' data('cellbench', package = 'Cepo')
#' cellbench
#' cepoOutput <- Cepo(logcounts(cellbench), cellbench$celltype)
#'
#' plotDensities(
#'   x = cellbench,
#'   cepoOutput = cepoOutput,
#'   assay = 'logcounts',
#'   plotType = 'histogram',
#'   celltypeColumn = 'celltype'
#' )
#'
#' plotDensities(
#'   x = cellbench,
#'   cepoOutput = cepoOutput,
#'   genes = c('PLTP', 'CPT1C', 'MEG3', 'SYCE1', 'MICOS10P3', 'HOXB7'),
#'   assay = 'logcounts',
#'   plotType = 'histogram',
#'   celltypeColumn = 'celltype'
#' )
plotDensities <- function(x, cepoOutput, nGenes = 2, assay = "logcounts", celltypeColumn, 
    celltype = NULL, genes = NULL, plotType = c("histogram", "density"), color = NULL) {
    if (!is.null(genes)) {
        nGenes <- NULL
    }
    if (is.null(celltype)) {
        celltype <- colnames(cepoOutput$stats)
    }
    
    plotType <- match.arg(plotType)
    
    stopifnot((is(x, "SummarizedExperiment") | is(x, "SingleCellExperiment")), is.character(assay), 
        length(assay) == 1L, is.character(celltypeColumn), length(celltypeColumn) == 
            1L, length(plotType) == 1L)
    
    # exprsMatrix:
    if (assay %notin% SummarizedExperiment::assayNames(x)) {
        message("'assay' not found in assayNames(x)")
        return(NULL)
    }
    exprsMatrix <- SummarizedExperiment::assays(x)[[assay]]
    
    ## extract cell type labels
    if (celltypeColumn %notin% colnames(colData(x))) {
        message("'cell-type label' not found in colnames(colData(x))")
        return(NULL)
    }
    celltype_label <- colData(x)[[celltypeColumn]]
    if (any(celltype_label %notin% celltype)) {
        # if idx_celltype FALSE:
        message("Some supplied celltypes are not found in the data")
        return(NULL)
    }
    
    ## If the genes argument is supplied (user-defined plotting genes)
    if (!is.null(genes)) {
        if (length(genes) > 10) {
            message("Maximum number of genes to plot is 10!")
            return(NULL)
        }
        
        if (any(genes %notin% rownames(x))) {
            # if some idx_gene FALSE:
            message("one or more genes are not found in `rownames(x)`")
            return(NULL)
        }
    }
    
    ## If the nGenes argument is supplied (so we will pick out the top genes in the
    ## computed object)
    if (!is.null(nGenes)) {
        if (nGenes > 10) {
            message("Maximum number of genes to plot is 10!")
            return(NULL)
        }
        
        genes = unlist(topGenes(object = cepoOutput, n = nGenes, returnValues = FALSE))
        message(paste(genes, collapse = ", "), " will be plotted")
        
    }
    
    # idx_gene = rownames(x) %in% genes
    
    if (!is.null(color)) {
        ## If user defined some colors
        if (length(color) != length(unique(celltype))) {
            message("`color` is of length ", length(color), ", but there are ", length(celltype), 
                " cell types")
            return(NULL)
        }
        if (length(names(color)) == 0L) {
            warning("`color` should be a named vector with names corresponding to the cell types. \n Colors sequential to the ordering of cell types will be used.")
            names(color) <- celltype
        }
    } else {
        color <- rainbow(length(unique(celltype)))
        names(color) <- unique(celltype)
    }
    
    plotmat <- exprsMatrix[genes, celltype_label %in% celltype, drop = FALSE]
    
    annotateDf <- data.frame(
      cellNames = colnames(plotmat),
      celltypeAnnotate = celltype_label[celltype_label %in% celltype]
    )
    
    plotdf <- merge(
      x = reshape2::melt(plotmat, varnames = c("geneNames", "cellNames")), 
      y = annotateDf, 
      by = "cellNames")
    
    plotdfSplit <- split.data.frame(x = plotdf, f = plotdf$geneNames)
    
    ggList <- lapply(
      X = plotdfSplit, 
      FUN = function(thisCellTypePlotDf){
        ggBase <- ggplot2::ggplot(data = thisCellTypePlotDf, ggplot2::aes(x = .data$value)) +
          ggplot2::facet_wrap(~ celltypeAnnotate, scales = "free_y", ncol = 1) + 
          ggplot2::theme_classic() + 
          ggplot2::labs(title = unique(thisCellTypePlotDf$geneNames),
                        fill = "Cell type") +
          ggplot2::scale_fill_manual(values = color[celltype]) +
          ggplot2::theme(axis.title.x = element_blank(), 
                         axis.title.y = element_blank(), 
                         strip.background = element_blank(),
                         strip.text.x = element_blank(),
                         legend.position = "bottom",
                         plot.title = element_text(hjust = 0.5))
        
        if(plotType == "histogram"){
          result <- ggBase +
            ggplot2::geom_histogram(aes(y = .data$..density..), 
                                    colour = "black", fill = "white", bins = 15) + 
            ggplot2::geom_density(aes(fill = .data$celltypeAnnotate), alpha = 0.2)
        } else if (plotType == "density"){
          result <- ggBase +
            ggplot2::geom_density(aes(fill = .data$celltypeAnnotate), alpha = 0.2)
        }
        return(result)
      })
    
    finalPlot <- patchwork::wrap_plots(ggList, nrow = 1) + 
      patchwork::plot_layout(guides = "collect") & 
      ggplot2::theme(legend.position = 'bottom')
    
    return(finalPlot)
}

`%notin%` <- Negate(`%in%`)

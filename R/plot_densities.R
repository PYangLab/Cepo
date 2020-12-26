#' Plot densities 
#' @return A \code{\link{ggplot}} object
#' with cell-type specific densities for a gene.
#' @param x a \code{\linkS4class{SummarizedExperiment}} or a \code{\linkS4class{SingleCellExperiment}} object.
#' @param dx_list an output from Cepo or doLimma/doVoom/doTtest/doWilcoxon functions
#' @param assay a character ("logcounts" by default), 
#' indicating the name of the assays(x) element which stores the expression data (i.e., assays(x)$name_assays_expression).
#' We strongly encourage using normalized data, such as counts per million (CPM) or log-CPM.
#' @param celltype a character, indicating the name of the cell type to plot. Default is NULL which selects all celltypes in the dx_list. 
#' @param celltype_col a character, indicating the name of the name of the cell type column in the colData(x).
#' @param genes a character vector, indicating the name of the genes to plot. Default to NULL, so that 2 top genes from each celltype will be plotted. 
#' @param n_genes number of top genes from each celltype to plot. Default to 2. 
#' @param color a named color vector. The names should correspond to the `celltype` argument above
#' @param plotType Either "histogram" or "density"
#' @import ggplot2
#' @importFrom ggpubr ggarrange
#' @importFrom methods is
#' @importFrom SingleCellExperiment colData
#' @importFrom SummarizedExperiment assays
#' @importFrom grDevices rainbow
#' @importFrom rlang .data
#' @return A \code{\link{ggplot}} object.
#' @export
#' @examples
#' library(SingleCellExperiment)
#' data("cellbench", package = "Cepo")
#' cellbench
#' ds_res <- Cepo(logcounts(cellbench), cellbench$celltype)
#' 
#' plot_densities(x = cellbench,
#'                dx_list = ds_res, 
#'                assay = "logcounts", 
#'                celltype_col = "celltype")
#' 
#' plot_densities(x = cellbench,
#'                dx_list = ds_res,
#'                genes = c("PLTP", "CPT1C", "MEG3", "SYCE1", "MICOS10P3", "HOXB7"),
#'                assay = "logcounts", 
#'                celltype_col = "celltype")
plot_densities = function(x, 
                          dx_list,
                          n_genes = 2,
                          assay = "logcounts",
                          celltype_col,
                          celltype = NULL,
                          genes = NULL,
                          plotType = c("histogram", "density"),
                          color = NULL){
  # browser()
  if(!is.null(genes)){n_genes = NULL}
  if(is.null(celltype)) {celltype = colnames(dx_list$stats)}
  
  plotType <- match.arg(plotType)

  stopifnot(
    (is(x, "SummarizedExperiment") | is(x, "SingleCellExperiment")),
    is.character(assay), length(assay) == 1L,
    is.character(celltype_col), length(celltype_col) == 1L,
    length(plotType) == 1L)
  
  # exprsMatrix:
  if(assay %notin% SummarizedExperiment::assayNames(x)){
    message("'assay' not found in assayNames(x)")
    return(NULL)
  }
  exprsMatrix = SummarizedExperiment::assays(x)[[assay]]
  
  ## extract cell type labels
  if(celltype_col %notin% colnames(colData(x))){
    message("'cell-type label' not found in colnames(colData(x))")
    return(NULL)
  }
  celltype_label = colData(x)[[celltype_col]]
  if(any(celltype_label %notin% celltype)) { # if idx_celltype FALSE:
    message("Some supplied celltypes are not found in the data")
    return(NULL)
  }
  
  ## If the genes argument is supplied (user-defined plotting genes)
  if(!is.null(genes)) {
    
    if (length(genes) > 10) {
      message("Maximum number of genes to plot is 10!")
      return(NULL)
    }
    
    if(any(genes %notin% rownames(x))){ # if some idx_gene FALSE:
      message("one or more genes are not found in `rownames(x)`")
      return(NULL)
    }
  }
  
  ## If the n_genes argument is supplied (so we will pick out the top genes in the computed object)
  if(!is.null(n_genes)) {
    if (n_genes > 10) {
      message("Maximum number of genes to plot is 10!")
      return(NULL)
    }
    
    genes <- unlist(topGenes(object = dx_list, n = n_genes, returnValues = FALSE))
    message(paste(genes, collapse = ", "), " will be plotted")
  }
  
  # idx_gene = rownames(x) %in% genes
  
  if(!is.null(color)) {
    ## If user defined some colors
    if(length(color) != length(unique(celltype))){
      message("`color` is of length ", length(color), ", but there are ", length(celltype), " cell types")
      return(NULL)
    }
    if(length(names(color)) == 0L){
      warning("`color` should be a named vector with names corresponding to the cell types. \n Colors sequential to the ordering of cell types will be used.")
      names(color) = celltype
    }
    
  } else {
    color = rainbow(length(unique(celltype)))
    names(color) <- unique(celltype)
  }
  

    # Density plot:
    gg.density.all <- lapply(celltype, function(cty) {
      Dx = genes
      n = length(Dx)
      df <- as.data.frame(t(exprsMatrix[Dx, celltype_label %in% cty, drop = FALSE]))
      df$celltype <- as.character(cty)
      
      my_col <- color
      my_col <- my_col[[cty]]
      
      listGGDensity <- lapply(1:n, function(i) {
        
        dftoplot <- df[,c(i, ncol(df)), drop = FALSE]
        genename <- colnames(dftoplot)[[1]]
        colnames(dftoplot) <- c("gene", "celltype")
        
      
        ggBase <- ggplot2::ggplot(dftoplot, aes(x = .data$gene)) + 
            ggplot2::facet_wrap(~celltype, scales = "free_y", dir  ="v") + 
            ggplot2::ggtitle(genename) + 
            ggplot2::theme_classic() + 
            ggplot2::theme(axis.title.x = element_blank(),
                           axis.title.y = element_blank(), 
                           strip.background = element_blank(), 
                           strip.text.x = element_blank(),
                           legend.position = "bottom")
        
        if(plotType == "histogram"){
          ggDensity = ggBase +
            ggplot2::geom_histogram(aes(y = .data$..density..), colour = "black", fill = "white", bins = 15) +
            ggplot2::geom_density(alpha = 0.2, fill = my_col)
        } else {
          ggDensity = ggBase +
            ggplot2::geom_density(alpha = 0.2, fill = my_col) 
        }

        return(ggDensity)
      })
      return(listGGDensity)
    })
    n = unique(sapply(gg.density.all, length))
    gg.list <- lapply(gg.density.all, function(x) do.call(ggpubr::ggarrange, c(x, ncol=n)))
    g <- do.call(ggpubr::ggarrange, c(gg.list, nrow = length(celltype)))
    
    g

}

`%notin%` <- Negate(`%in%`)

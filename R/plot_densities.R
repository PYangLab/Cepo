#' Plot densities 
#'
#' @author Hani Jieun Kim \email{hani.kim@sydney.edu.au}
#' 
#' \code{plot_density} returns a \code{\link{ggplot}} object
#' with cell-type specific densities for a gene.
#' 
#' @param x a \code{\linkS4class{SummarizedExperiment}} or a \code{\linkS4class{SingleCellExperiment}} object.
#' @param dx_list a list containing the ordered output of Dx runs.
#' @param assay a character ("logcounts" by default), 
#' indicating the name of the assays(x) element which stores the expression data (i.e., assays(x)$name_assays_expression).
#' We strongly encourage using normalized data, such as counts per million (CPM) or log-CPM.
#' @param celltype a character, indicating the name of the cluster to plot. Default is all clusters.
#' @param celltype_col X.
#' @param gene a character, indicating the name of the gene to plot. Default is NULL. If is.null(gene), this argument overrides the `dx_list`.
#' @param hist a logical, indicating whether to plot histogram (if TRUE) or not (if FALSE).
#' @param color a named color vector.
#' @import ggplot2
#' @importFrom ggpubr ggarrange
#' @importFrom methods is
#' @importFrom SingleCellExperiment colData
#' 
#' @return A \code{\link{ggplot}} object.
#' @examples
#' data("cellbench", package = "Cepo")
#' cellbench
#' ds_res <- Cepo(logcounts(cellbench), cellbench$celltype)
#' 
#' plot_densities(x = cellbench,
#'                dx_list = ds_res, 
#'                n_genes = 2,
#'                assay = "logcounts", 
#'                celltype_col = "celltype",
#'                celltype = names(ds_res))
#' 
#' plot_densities(x = cellbench,
#'                dx_list=ds_res,
#'                genes = unlist(lapply(ds_res, function(x) names(x)[10:11])),
#'                assay = "logcounts", 
#'                celltype_col = "celltype",
#'                celltype = names(ds_res))
#' 
#' 
#' @export
plot_densities = function(x, 
                          dx_list,
                          n_genes = NULL,
                          assay = "logcounts",
                          celltype_col,
                          celltype,
                          genes = NULL,
                          hist_plot = TRUE,
                          color=NULL){
  
  if (is.null(n_genes)) { n_genes = 2 } else { n_genes }
  
  stopifnot(
    ( is(x, "SummarizedExperiment") | is(x, "SingleCellExperiment") ),
    is.numeric(n_genes), length(n_genes) == 1L,
    is.character(assay), length(assay) == 1L,
    is.character(celltype_col), length(celltype_col) == 1L,
    is.character(celltype), length(celltype) > 0L,
    is.logical(hist_plot), length(hist_plot) == 1L)
  
  # exprsMatrix:
  sel = which(names(assays(x)) == assay)
  if( length(sel) == 0 ){
    message("'assay' not found in names(assays(x))")
    return(NULL)
  }
  if( length(sel) > 1 ){
    message("multiple 'assay' found in names(assays(x))")
    return(NULL)
  }
  exprsMatrix = assays(x)[[sel]]
  
  # extract cell type labels
  sel = which(names(colData(x)) == celltype_col)
  if( length(sel) == 0 ){
    message("'cell-type label' not found in names(colData(x))")
    return(NULL)
  }
  if( length(sel) > 1 ){
    message("more than one 'cell-type label' found in names(colData(x))")
    return(NULL)
  }
  celltype_label = colData(x)[[sel]]
  
  # extract genes
  if(!is.null(genes)) {
    
    if (length(genes) > 10) {
      message("number of `genes` too many to plot")
      return(NULL)
    } 
    
    idx_gene = genes %in% rownames(x)
    idx_gene[is.na(idx_gene)] = FALSE
    
    if( sum(!idx_gene) > 0 ){ # if some idx_gene FALSE:
      message("one or more 'genes' not found in `rownames(x)`")
      return(NULL)
    }
  }
    
  if (!is.null(n_genes)) {

    if (n_genes > 10) {
      message("number of `genes` too many to plot")
      return(NULL)
    } 
    
  }
  
  # genes to plot
  if (!is.null(genes)) {
    
    idx_gene = rownames(x) %in% genes

  } else {
    
    idx_gene_list <- lapply(dx_list, function(x) names(x)[1:n_genes])
    idx_gene_list <- idx_gene_list[celltype]

  }
  
  # cell types
  idx_celltype = celltype_label %in% celltype
  idx_celltype[is.na(idx_celltype)] = FALSE
  
  if( sum(!idx_celltype) > 0) { # if idx_celltype FALSE:
    message("'celltype' not found in `colData(x)[[celltype_col]]`")
    return(NULL)
  }
  
  if( !missing(color) ) { 
    
    # colors
    idx_color = names(color) %in% celltype
    idx_color[is.na(idx_color)] = FALSE
    
    if (length(idx_color) == length(unique(celltype))) {
      message("'color' does not match number of `celltype` to plot")
      return(NULL)
    } else if (sum(!idx_color) > 0) {
      message("names of 'color' does not match `celltype` labels")
      return(NULL)
    }
    
  } else {
    
    color = rainbow(length(unique(celltype)))
    names(color) <- unique(celltype)

  }
  

    # Density plot:
    gg.density.all <- lapply(celltype, function(cty) {
      
      if(is.null(genes)){
      
      Dx <- unlist(idx_gene_list) 
      
      } else {
        
      Dx <- genes
      
      }
      n = length(Dx)
      df <- as.data.frame(t(exprsMatrix[Dx,celltype_label %in% cty]))
      df$celltype <- as.character(cty)
      
      my_col <- color
      my_col <- my_col[[cty]]
      
      gg.density <- lapply(1:n, function(i) {
        
        dftoplot <- df[,c(i,ncol(df))]
        genename <- colnames(dftoplot)[[1]]
        colnames(dftoplot) <- c("gene", "celltype")
        
        if (hist_plot == TRUE) {
          gg.den <- ggplot2::ggplot(dftoplot, aes(x=gene)) + 
            ggplot2::geom_histogram(aes(y=..density..), colour="black", fill="white", bins=15) +
            ggplot2::facet_wrap(~celltype, scales="free_y", dir="v") + 
            ggplot2::ggtitle(genename) + 
            ggplot2::geom_density(alpha=.2, fill=my_col) + 
            ggplot2::theme_classic() + ggplot2::theme(axis.title.x=element_blank(), axis.title.y=element_blank(), 
                                    strip.background = element_blank(), strip.text.x = element_blank())
          return(gg.den)
        } else {
          gg.den <- ggplot2::ggplot(dftoplot, aes(x=gene)) + 
            ggplot2::facet_wrap(~celltype, scales="free_y", dir="v") + 
            ggplot2::ggtitle(genename) + 
            gggplot2::eom_density(alpha=.2, fill=my_col) + 
            ggplot2::theme_classic() + ggplot2::theme(axis.title.x=element_blank(), axis.title.y=element_blank(), 
                                    strip.background = element_blank(), strip.text.x = element_blank())
          return(gg.den)
        }
      })
      return(gg.density)
    })
    
    n = unique(sapply(gg.density.all, length))
    gg.list <- lapply(gg.density.all, function(x) do.call(ggpubr::ggarrange, c(x, ncol=n)))
    g <- do.call(ggpubr::ggarrange, c(gg.list, nrow=length(celltype)))
    
    g

}



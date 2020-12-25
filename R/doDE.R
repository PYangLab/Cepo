#' @title Computing (trend-)Limma cell identity genes
#' @param exprsMat expression matrix where columns denote cells and rows denote genes
#' @param cellTypes vector of cell type labels
#' @importFrom DelayedMatrixStats rowSums2 rowMeans2
#' @importFrom DelayedArray cbind
#' @importFrom limma lmFit eBayes topTable
#' @importFrom methods new
#' @return Returns a list, each element corresponds to statistics (and p-values) associated with a cell-type
#' @description exprsMat accepts various matrix objects, including DelayedArray and HDF5Array for 
#' out-of-memory computations. 
#' @rdname DE-cell-identity-methods
#' @export
#' @examples 
#' set.seed(1234)
#' n = 50 ## genes, rows
#' p = 100 ## cells, cols
#' exprsMat = matrix(rpois(n*p, lambda = 5), nrow = n)
#' rownames(exprsMat) = paste0("gene", 1:n)
#' colnames(exprsMat) = paste0("cell", 1:p)
#' cellTypes = sample(letters[1:3], size = p, replace = TRUE)
#' 
#' doLimma(exprsMat = exprsMat, cellTypes = cellTypes)
#' doVoom(exprsMat = exprsMat, cellTypes = cellTypes)
#' doTtest(exprsMat = exprsMat, cellTypes = cellTypes)
#' doWilcoxon(exprsMat = exprsMat, cellTypes = cellTypes)
doLimma <- function(exprsMat, cellTypes){
  cellTypes <- droplevels(as.factor(cellTypes))
  tt <- list()
  for (i in 1:nlevels(cellTypes)) {
    tmp_celltype <- (ifelse(cellTypes == levels(cellTypes)[i], 1, 0))
    design <- stats::model.matrix(~tmp_celltype)
    meanExprs <- do.call(DelayedArray::cbind, lapply(c(0,1), function(i){
      DelayedMatrixStats::rowMeans2(exprsMat[, tmp_celltype == i, drop = FALSE])
    }))
    meanPct <- do.call(DelayedArray::cbind, lapply(c(0,1), function(i){
      DelayedMatrixStats::rowSums2(exprsMat[, tmp_celltype == i, drop = FALSE] > 0)/sum(tmp_celltype == i)
    }))
    rownames(meanExprs) = rownames(meanPct) = rownames(exprsMat)
    y <- methods::new("EList")
    y$E <- exprsMat
    fit <- limma::lmFit(y, design = design)
    fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
    tt[[i]] <- limma::topTable(fit, n = Inf, adjust.method = "BH", coef = 2)
    if (!is.null(tt[[i]]$ID)) {
      tt[[i]] <- tt[[i]][!duplicated(tt[[i]]$ID),]
      rownames(tt[[i]]) <- tt[[i]]$ID
    }
    tt[[i]]$meanExprs.1 <- meanExprs[rownames(tt[[i]]), 1]
    tt[[i]]$meanExprs.2 <- meanExprs[rownames(tt[[i]]), 2]
    tt[[i]]$meanPct.1 <- meanPct[rownames(tt[[i]]), 1]
    tt[[i]]$meanPct.2 <- meanPct[rownames(tt[[i]]), 2]
  }
  names(tt) <- levels(cellTypes)
  
  statsResult = lapply(tt, function(x){
      stats <- x$t
      names(stats) <- rownames(x)
      stats <- sort(stats, decreasing = TRUE, na.last = TRUE)
      return(stats)
    })
  
  pvalsResult = lapply(tt, function(x) {
    pvals <- x$P.Value
    names(pvals) <- rownames(x)
    return(pvals)
  })
  
  result = list(
    stats = S4Vectors::DataFrame(sortList(statsResult)),
    pvalues = S4Vectors::DataFrame(sortList(pvalsResult))
  )
  return(result)
}

#' @title Computing Limma-Voom cell identity genes
#' @importFrom DelayedMatrixStats rowSums2 rowMeans2
#' @importFrom DelayedArray cbind
#' @importFrom limma lmFit eBayes topTable voom
#' @rdname DE-cell-identity-methods
#' @export
doVoom <- function(exprsMat, cellTypes) {
  # input must be normalised, log-transformed data
  cty <- droplevels(as.factor(cellTypes))
  names(cty) <- colnames(exprsMat)
  
  tt <- list()
  for (i in 1:nlevels(cty)) {
    tmp_celltype <- (ifelse(cty == levels(cty)[i], 1, 0))
    
    design <- stats::model.matrix(~tmp_celltype)
    
    y <- methods::new("EList")
    y$E <- exprsMat
    vm <- limma::voom(y, design = design)
    fit <- limma::lmFit(vm, design = design)
    fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
    tt[[i]] <- limma::topTable(fit, n = Inf, adjust.method = "BH", coef = 2)
    if (!is.null(tt[[i]]$ID)) {
      tt[[i]] <- tt[[i]][!duplicated(tt[[i]]$ID),]
      rownames(tt[[i]]) <- tt[[i]]$ID
    }
    
  }
  names(tt) <- levels(cty)
  statsResult = lapply(tt, function(x){
    stats <- x$t
    names(stats) <- rownames(x)
    stats <- sort(stats, decreasing = TRUE, na.last = TRUE)
    return(stats)
  })
  
  pvalsResult = lapply(tt, function(x) {
    pvals <- x$P.Value
    names(pvals) <- rownames(x)
    return(pvals)
  })
  
  result = list(
    stats = S4Vectors::DataFrame(sortList(statsResult)),
    pvalues = S4Vectors::DataFrame(sortList(pvalsResult))
  )
  return(result)
}

#' @title Computing t-test identity genes
#' @importFrom DelayedArray apply
#' @rdname DE-cell-identity-methods
#' @export
doTtest <- function(exprsMat, cellTypes) {
  # input must be normalised, log-transformed data
  cty <- droplevels(as.factor(cellTypes))
  names(cty) <- colnames(exprsMat)
  
  tt <- list()
  for (i in 1:nlevels(cty)) {
    tmp_celltype <- (ifelse(cty == levels(cty)[i], 1, 0))
    
    tt[[i]] <- t(DelayedArray::apply(exprsMat, 1, function(x) {
      x1 <- x[tmp_celltype == 0]
      x2 <- x[tmp_celltype == 1]
      
      res <- stats::t.test(x2, y=x1)
      return(c(stats=res$statistic,
               pvalue=res$p.value))
    }))
    tt[[i]] <- as.data.frame(tt[[i]])
    tt[[i]]$adj.pvalue <- stats::p.adjust(tt[[i]]$pvalue, method = "BH")
  }
  names(tt) <- levels(cty)
  statsResult = lapply(tt, function(x){
    stats <- x$stats.t
    names(stats) <- rownames(x)
    stats <- sort(stats, decreasing = TRUE, na.last = TRUE)
    return(stats)
  })
  
  pvalsResult = lapply(tt, function(x) {
    pvals <- x$pvalue
    names(pvals) <- rownames(x)
    return(pvals)
  })
  
  result = list(
    stats = S4Vectors::DataFrame(sortList(statsResult)),
    pvalues = S4Vectors::DataFrame(sortList(pvalsResult))
  )
  return(result)
}

#' @title Computing Wilcoxon identity genes
#' @importFrom DelayedArray apply
#' @rdname DE-cell-identity-methods
#' @export
doWilcoxon <- function(exprsMat, cellTypes) {
  if(!is.matrix(exprsMat)){
    stop("Wilcoxon method only supports matrix input at the moment")
  }
  # input must be normalised, log-transformed data
  cty <- droplevels(as.factor(cellTypes))
  names(cty) <- colnames(exprsMat)
  
  tt <- list()
  for (i in 1:nlevels(cty)) {
    tmp_celltype <- (ifelse(cty == levels(cty)[i], 1, 0))
    
    tt[[i]] <- t(DelayedArray::apply(exprsMat, 1, function(x) {
      res <- stats::wilcox.test(x ~ tmp_celltype)
      c(stats=res$statistic,
        pvalue=res$p.value)
    }))
    tt[[i]] <- as.data.frame(tt[[i]])
    tt[[i]]$adj.pvalue <- stats::p.adjust(tt[[i]]$pvalue, method = "BH")
  }
  names(tt) <- levels(cty)
  
  statsResult = lapply(tt, function(x){
    stats <- x$stats.W
    names(stats) <- rownames(x)
    stats <- sort(stats, decreasing = TRUE, na.last = TRUE)
    return(stats)
  })
  
  pvalsResult = lapply(tt, function(x) {
    pvals <- x$pvalue
    names(pvals) <- rownames(x)
    return(pvals)
  })
  
  result = list(
    stats = S4Vectors::DataFrame(sortList(statsResult)),
    pvalues = S4Vectors::DataFrame(sortList(pvalsResult))
  )
  
  return(result)
}
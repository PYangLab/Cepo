#' @title Extract differential statistics p-values from DX methods
#' @param object a computed result from Cepo()
#' @param returnDF Whether return a data.frame, default to TRUE.
#' @importFrom S4Vectors DataFrame
#' @return If returnDF is set to TRUE, then returns a data frame. Otherwise returns a list.
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
#' cepo_output = Cepo(exprsMat = exprsMat, cellTypes = cellTypes)
#' limma_output = doLimma(exprsMat = exprsMat, cellTypes = cellTypes)
#' voom_output = doVoom(exprsMat = exprsMat, cellTypes = cellTypes)
#' ttest_output = doTtest(exprsMat = exprsMat, cellTypes = cellTypes)
#' wilcoxon_output = doWilcoxon(exprsMat = exprsMat, cellTypes = cellTypes)
#' 
#' getPvalues(cepo_output)
#' getPvalues(limma_output)
#' getPvalues(voom_output)
#' getPvalues(ttest_output)
#' getPvalues(wilcoxon_output)
getPvalues <- function(object, returnDF = TRUE){
  method = attr(x = object, "differential_method")
  
  listResult = switch(
    method, 
    # Cepo = getPvalues_cepo(object),
    limma = getPvalues_limma(object = object),
    voom = getPvalues_limma(object = object),
    ttest = getPvalues_tTest(object = object),
    wilcoxon = getPvalues_wilcoxon(object = object))
  
  if(returnDF){ ## if TRUE, turn results to a data.frame
    return(listUnsorted2DF(listResult))
  } else{ ## else, just return the lists
    return(listResult)
  }
}

getPvalues_limma = function(object){
  listResult <- lapply(object, function(x) {
    pvals <- x$P.Value
    names(pvals) <- rownames(x)
    pvals <- sort(pvals, decreasing = TRUE, na.last = TRUE)
    return(pvals)
  })
  return(listResult)
}

getPvalues_tTest = function(object){
  listResult <- lapply(object, function(x) {
    pvals <- x$pvalue
    names(pvals) <- rownames(x)
    pvals <- sort(pvals, decreasing = TRUE, na.last = TRUE)
    return(pvals)
  })
  return(listResult)
}

getPvalues_wilcoxon = function(object){
  listResult <- lapply(object, function(x) {
    pvals <- -x$pvalue
    names(pvals) <- rownames(x)
    pvals <- sort(pvals, decreasing = TRUE, na.last = TRUE)
    return(pvals)
  })
  return(listResult)
}

# getPvalues_cepo = function(object, exprsMat, cellTypes, exprs_pct, times = 100){
#   list_cepo_outputs = lapply(
#     X = seq_len(times),
#     FUN = function(i){
#       Cepo(exprsMat = exprsMat, cellTypes = sample(cellTypes), exprs_pct = exprs_pct)}
#     )
#   
#   list_pvals = vector("list", length = length(object))
#   names(list_pvals) = names(object)
#   
#   for(i in names(object)){
#     list_binary = lapply(list_cepo_outputs, function(this_run){
#       object[[i]] >= this_run[[i]][names(object[[i]])]
#     })
#     
#     list_pvals[[i]] = DelayedMatrixStats::colMeans2(do.call(rbind, list_binary))
#     names(list_pvals[[i]]) = names(object[[i]])
#   }
#   
#   return(list_pvals)
# }
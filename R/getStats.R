#' @title Extract differential statistics from DX methods
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
#' 
#' getStats(cepo_output)
#' getStats(limma_output)
#' getStats(voom_output)
getStats <- function(object, returnDF = TRUE){
  method = attr(x = object, "differential_method")
  
  listResult = switch(
    method, 
    Cepo = object,
    limma = getStats_limma(object = object),
    voom = getStats_limma(object = object)
  )
  
  if(returnDF){ ## if TRUE, turn results to a data.frame
    return(listUnsorted2DF(listResult))
  } else{ ## else, just return the lists
    return(listResult)
  }
}

getStats_limma = function(object){
  listResult <- lapply(object, function(x) {
    stats <- x$t
    names(stats) <- rownames(x)
    stats <- sort(stats, decreasing = TRUE, na.last = TRUE)
    return(stats)
  })
  return(listResult)
}

listUnsorted2DF = function(object){
  geneNames = sort(names(object[[1]]))
  objectSorted = lapply(object, function(thisCellTypeResult){thisCellTypeResult[geneNames]})
  result = data.frame(do.call(cbind, objectSorted))
  result = S4Vectors::DataFrame(result)
  # result = cbind(geneName = rownames(result),
  #                result)
  return(result)
}
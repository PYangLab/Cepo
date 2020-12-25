#' @title Extract the top genes from the Cepo output
#' @param object Output from the Cepo function
#' @param n Number of top genes to extract
#' @param returnValues Whether to return the numeric value associated with the top selected genes
#' @return Returns a list of key genes. 
#' @export
#' @examples
#' set.seed(1234)
#' n = 50 ## genes, rows
#' p = 100 ## cells, cols
#' exprsMat = matrix(rpois(n*p, lambda = 5), nrow = n)
#' rownames(exprsMat) = paste0("gene", 1:n)
#' colnames(exprsMat) = paste0("cell", 1:p)
#' cellTypes = sample(letters[1:3], size = p, replace = TRUE)
#' cepo_output = Cepo(exprsMat = exprsMat, cellTypes = cellTypes)
#' cepo_output
#' topGenes(cepo_output, n = 10)
#' topGenes(cepo_output, n = 10, returnValues = TRUE)
topGenes <- function(object, n = min(10, nrow(object$stats)),
                     returnValues = FALSE){
  geneNames <- rownames(object$stats)
  
  stopifnot(n <= nrow(object$stats))
  
  if(returnValues){
    result <- lapply(object$stats, function(this_column){
      this_column[order(this_column, decreasing = TRUE) <= n]
    })
  } else {
    result <- lapply(object$stats, function(this_column){
      geneNames[order(this_column, decreasing = TRUE) <= n]
    })
  }
  
  return(result)
}
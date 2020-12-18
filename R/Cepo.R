#' @title Computing Cepo cell identity genes
#' @param exprsMat expression matrix where columns denote cells and rows denote genes
#' @param cellTypes vector of cell type labels
#' @param exprs_pct Percentage of lowly expressed genes to remove. Default to NULL to not remove any genes. 
#' @importFrom DelayedMatrixStats rowSums2 rowMeans2 rowSds
#' @importFrom DelayedArray cbind 
#' @return Returns a list of key genes. 
#' @description exprsMat accepts various matrix objects, including DelayedArray and HDF5Array for 
#' out-of-memory computations. 
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
#' Cepo(exprsMat = exprsMat, cellTypes = cellTypes)
Cepo <- function(exprsMat, cellTypes, exprs_pct = NULL) {
    cts <- names(table(cellTypes))
    
    if (!is.null(exprs_pct)) {
        
        meanPct.list <- list()
        for(i in 1:length(cts)){
            idx <- which(cellTypes == cts[i])
            meanPct.list[[i]] <- (rowSums_withnames(exprsMat[, idx, drop = FALSE] > 0)/sum(cellTypes == cts[i])) > exprs_pct 
        }
        names(meanPct.list) <- cts
        keep = rowSums_withnames(do.call(DelayedArray::cbind, meanPct.list)) == length(cts) 
        exprsMat <- exprsMat[keep,]
        
    }
    segIdx.list <- list()
    for(i in 1:length(cts)){
        idx <- which(cellTypes == cts[i])
        segIdx.list[[i]] <- segIndex(exprsMat[,idx])
    }
    names(segIdx.list) <- cts
    
    segMat <- segIdxList2Mat(segIdx.list)
    segGenes <- consensusSegIdx(segMat)
    names(segGenes) <- colnames(segMat)
    result = segGenes
    class(result) = c("Cepo", class(result))
    attr(result, "differential_method") = "Cepo"
    return(result)
}

print.Cepo = function(x, ...){
 cat("Top 6 Cepo genes in each celltype: \n")
 for(j in names(x)){
     cat(j, "\n")
     print(utils::head(x[[j]]))
 }
}

segIndex <- function(mat){
    nz <- rowMeans_withnames(mat != 0)
    ms <- rowMeans_withnames(mat)
    sds <- rowSds_withnames(mat)
    cvs <- sds/ms
    
    x1 <- rank(nz)/(length(nz)+1)
    x2 <- 1 - rank(cvs)/(length(cvs)+1)
    
    segIdx <- rowMeans_withnames(DelayedArray::cbind(x1, x2))
    return(segIdx)
}

rowMeans_withnames = function(mat){
    result = DelayedMatrixStats::rowMeans2(mat)
    names(result) = rownames(mat)
    return(result)
}

rowSds_withnames = function(mat){
    result = DelayedMatrixStats::rowSds(mat)
    names(result) = rownames(mat)
    return(result)
}

rowSums_withnames = function(mat){
    result = DelayedMatrixStats::rowMeans2(mat)
    names(result) = rownames(mat)
    return(result)
}

segIdxList2Mat <- function(segIdx.list) {
    allGenes <- unique(unlist(lapply(segIdx.list, names)))
    segMat <- matrix(0, nrow=length(allGenes), ncol=length(segIdx.list))
    rownames(segMat) <- allGenes
    colnames(segMat) <- names(segIdx.list)
    
    for(i in 1:length(segIdx.list)) {
        si <- segIdx.list[[i]]
        segMat[names(si),i] <- si
    }
    return(segMat)
}

consensusSegIdx <- function(mat) {
    tt <- mat
    
    CIGs <- list()
    for (i in 1:ncol(tt)) {
        avgRank <- c()
        for(j in 1:ncol(tt)) {
            if(i == j){next}
            avgRank <- cbind(avgRank, tt[,i] - tt[,j])
        }
        
        meanAvgRank = DelayedMatrixStats::rowMeans2(avgRank)
        names(meanAvgRank) <- rownames(avgRank)
        CIGs[[i]] <- sort(meanAvgRank, decreasing = TRUE)
    }
    return(CIGs)
}

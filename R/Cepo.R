#' @title Computing Cepo cell identity genes
#' @param exprsMat expression matrix where columns denote cells and rows denote genes
#' @param cellTypes vector of cell type labels
#' @param exprs_pct Percentage of lowly expressed genes to remove. Default to NULL to not remove any genes. 
#' @param BPPARAM A BiocParallelParam object specifying how parallelization should be performed.
#' Default to BiocParallel::SerialParam()
#' @importFrom DelayedMatrixStats rowSums2 rowMeans2 rowSds
#' @importFrom DelayedArray cbind blockApply 
#' @importFrom BiocParallel SerialParam
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
Cepo <- function(exprsMat, cellTypes, exprs_pct = NULL, BPPARAM = BiocParallel::SerialParam()) {
    cts <- names(table(cellTypes))
    
    if (!is.null(exprs_pct)) {
        
        meanPct.list <- list()
        for(i in 1:length(cts)){
            idx <- which(cellTypes == cts[i])
            meanPct.list[[i]] <- (block_rowSums(exprsMat[, idx, drop = FALSE] > 0)/sum(cellTypes == cts[i])) > exprs_pct 
        }
        names(meanPct.list) <- cts
        keep = block_rowSums(do.call(DelayedArray::cbind, meanPct.list)) == length(cts) 
        exprsMat <- exprsMat[keep,]
        
    }
    segIdx.list <- list()
    for(i in 1:length(cts)){
        idx <- which(cellTypes == cts[i])
        segIdx.list[[i]] <- segIndex(exprsMat[,idx], BPPARAM = BPPARAM)
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

segIndex <- function(mat, BPPARAM){
    nz <- block_rowMeans(mat != 0, BPPARAM = BPPARAM)
    ms <- block_rowMeans(mat, BPPARAM = BPPARAM)
    sds <- block_rowSds(mat, BPPARAM = BPPARAM)
    cvs <- sds/ms
    names(nz) = names(sds) = names(cvs) = names(ms) = rownames(mat)
    
    x1 <- rank(nz)/(length(nz)+1)
    x2 <- 1 - rank(cvs)/(length(cvs)+1)
    
    segIdx <- block_rowMeans(DelayedArray::cbind(x1, x2), BPPARAM = BPPARAM)
    names(segIdx) = rownames(mat)
    return(segIdx)
}

block_rowMeans = function(mat, BPPARAM = BiocParallel::SerialParam()){
    result = DelayedArray::blockApply(x = mat, FUN = DelayedMatrixStats::rowMeans2,
                                      BPPARAM = BPPARAM)[[1]]
    names(result) = rownames(mat)
    return(result)
}

block_rowSds = function(mat, BPPARAM = BiocParallel::SerialParam()){
    result = DelayedArray::blockApply(x = mat, FUN = DelayedMatrixStats::rowSds,
                                      BPPARAM = BPPARAM)[[1]]
    names(result) = rownames(mat)
    return(result)
}

block_rowSums = function(mat, BPPARAM = BiocParallel::SerialParam()){
    result = DelayedArray::blockApply(x = mat, FUN = DelayedMatrixStats::rowSums2,
                                      BPPARAM = BPPARAM)[[1]]
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

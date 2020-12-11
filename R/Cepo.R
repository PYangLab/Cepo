#' @title Functions to generate differentially expressed genes
#' @param exprsMat expression matrix where columns denote cells and rows denote genes
#' @param cellTypes vector of cell type labels
#' @examples 
#' 
#' ####### an example where `sce`` is a SingleCellExperiment object
#' mat <- logcounts(sce)
#' labels <- sce$lineage
#' 
#' ds_res <- Cepo(mat, labels, exprs_pct=0.05, filter=FALSE)
#' 
#' ####### A quick simulation. Note that we need to have colnames and rownames on matrix ##########
#' set.seed(1234)
#' n = 1000 ## genes, rows
#' p = 10000 ## cells, cols
#' exprsMat = matrix(rpois(n*p, lambda = 5), nrow = n)
#' rownames(exprsMat) = paste0("gene", 1:n)
#' colnames(exprsMat) = paste0("cell", 1:p)
#' cellTypes = sample(letters[1:3], size = p, replace = TRUE)
#' mat_output = Cepo(exprsMat = exprsMat, cellTypes = cellTypes)
#' 
#' library(DelayedArray)
#' library(DelayedMatrixStats)
#' library(HDF5Array)
#' DelayedArray:::set_verbose_block_processing(TRUE)
#' ## DelayedArray::setAutoBlockSize(size = 1e3)
#' ## DelayedArray::getAutoBlockSize()
#' exprsMat_da = DelayedArray::DelayedArray(realize(exprsMat, "HDF5Array"))
#' class(seed(exprsMat_da))
#' da_output = Cepo(exprsMat = exprsMat_da, cellTypes = cellTypes)
#' identical(mat_output, da_output)
#' 
#' library(microbenchmark)
#' microbenchmark::microbenchmark(
#' mat_output = Cepo(exprsMat = exprsMat, cellTypes = cellTypes),
#' da_output = Cepo(exprsMat = exprsMat_da, cellTypes = cellTypes), 
#' times = 10)
Cepo <- function(exprsMat, cellTypes, exprs_pct=0.05, filter=FALSE) {
    cts <- names(table(cellTypes))
    
    if (filter==TRUE) {
        
        meanPct.list <- list()
        for(i in 1:length(cts)){
            idx <- which(cellTypes == cts[i])
            meanPct.list[[i]] <- (DelayedMatrixStats::rowSums2(exprsMat[, idx, drop = FALSE] > 0)/sum(cellTypes == cts[i])) > exprs_pct 
        }
        names(meanPct.list) <- cts
        keep = rowSums(do.call(cbind, meanPct.list)) == length(cts) 
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
    
    return(segGenes)
}

segIndex <- function(mat){
    # nz <- (rowSums(mat != 0) / ncol(mat))
    # ms <- apply(mat, 1, function(x){mean(x)})
    # cvs <- apply(mat, 1, function(x){
    #     return(sd(x)/mean(x))
    # })
    
    # nz <- Matrix::rowMeans(mat != 0)
    # ms <- Matrix::rowMeans(mat)
    # sds <- matrixStats::rowSds(mat)
    nz <- DelayedMatrixStats::rowMeans2(mat != 0)
    ms <- DelayedMatrixStats::rowMeans2(mat)
    sds <- DelayedMatrixStats::rowSds(mat)
    cvs <- sds/ms
    names(nz) = names(sds) = names(cvs) = names(ms) = rownames(mat)
    
    x1 <- rank(nz)/(length(nz)+1)
    x2 <- 1 - rank(cvs)/(length(cvs)+1)
    
    # segIdx <- apply(cbind(x1, x2), 1, mean)
    segIdx <- DelayedMatrixStats::rowMeans2(DelayedArray::cbind(x1, x2))
    names(segIdx) = rownames(mat)
    return(segIdx)
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
        
        CIGs[[i]] <- sort(rowMeans(avgRank), decreasing = TRUE)
    }
    return(CIGs)
}
###############




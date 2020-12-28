#' @title Computing Cepo cell identity genes
#' @param exprsMat expression matrix where columns denote cells and rows denote genes
#' @param cellTypes vector of cell type labels
#' @param exprsPct Percentage of lowly expressed genes to remove. Default to NULL to not remove any genes.
#' @param computePvalue Whether to compute p-values using bootstrap test. Default to NULL to not make computations.
#' Set this to an integer to set the number of bootstraps needed (recommend to be at least 100).
#' @param workers Number of cores to use. Default to 1, which invokes `BiocParallel::SerialParam`.
#' For workers greater than 1, see the `workers` argument in `BiocParallel::MulticoreParam` and `BiocParallel::SnowParam`. 
#' @param ... Additional arguments passed to `BiocParallel::MulticoreParam` and `BiocParallel::SnowParam`. 
#' @importFrom DelayedMatrixStats rowSums2 rowMeans2 rowSds
#' @importFrom DelayedArray cbind
#' @importFrom HDF5Array HDF5Array
#' @importFrom BiocParallel SerialParam MulticoreParam SnowParam bplapply
#' @return Returns a list of key genes.
#' @description exprsMat accepts various matrix objects, including DelayedArray and HDF5Array for
#' out-of-memory computations. See vignette.
#' @export
#' @examples
#' set.seed(1234)
#' n <- 50 ## genes, rows
#' p <- 100 ## cells, cols
#' exprsMat <- matrix(rpois(n * p, lambda = 5), nrow = n)
#' rownames(exprsMat) <- paste0('gene', 1:n)
#' colnames(exprsMat) <- paste0('cell', 1:p)
#' cellTypes <- sample(letters[1:3], size = p, replace = TRUE)
#'
#' Cepo(exprsMat = exprsMat, cellTypes = cellTypes)
#' Cepo(exprsMat = exprsMat, cellTypes = cellTypes, computePvalue = 100)
Cepo <- function(exprsMat, cellTypes, exprsPct = NULL, computePvalue = NULL,
                 workers = 1L, ...) {
    if (is.null(rownames(exprsMat))) {
        ## Add rownames if missing
        nGenes <- nrow(exprsMat)
        rownames(exprsMat) <- sprintf(paste0("gene%0", log10(nGenes) + 1, "d"), seq_len(nGenes))
    }
    cellTypes <- as.character(cellTypes)
    
    ## A single run of oneCepo gives a list of output
    singleResult <- oneCepo(exprsMat = exprsMat, cellTypes = cellTypes, exprsPct = exprsPct,
                            workers = workers, ...)
    ## Export Cepo outputs as a DataFrame
    singleStatsResult <- S4Vectors::DataFrame(sortList(singleResult))
    
    if (is.null(computePvalue)) {
        ## If no need to compute p-values, then that is it.
        result <- list(stats = singleStatsResult, pvalues = NULL)
    } else {
        ## Coerce the input to an integer
        times <- as.integer(computePvalue)
        ## Pass onto a boot function for computation
        listPvals <- bootCepo(exprsMat = exprsMat, cellTypes = cellTypes, exprsPct = exprsPct, 
            singleResult = singleResult, times = times, workers = workers, ...)
        ## The output has two components, one DataFrame of stats and another for p-values
        result <- list(stats = singleStatsResult, pvalues = S4Vectors::DataFrame(sortList(listPvals)))
    }  ## End else
    class(result) <- c("Cepo", class(result))
    return(result)
}

print.Cepo <- function(x) {
    cat("Computed statistics: \n \n")
    print(x$stats)
    cat("Computed p-values: \n")
    if (is.null(x$pvalues)) {
        cat("Note: a valid value for `computePvalue` argument is needed to get p-values when running the `Cepo` function")
    } else {
        print(x$pvalues)
    }
}

bootCepo <- function(exprsMat, cellTypes, exprsPct, singleResult, times, workers = 1L, ...) {
    ## Running multiple runs of Cepo based on bootstrap
    listCepoOutputs <- BiocParallel::bplapply(X = seq_len(times), FUN = function(i) {
        oneCepo(exprsMat = exprsMat, cellTypes = sample(cellTypes), exprsPct = exprsPct, workers = 1L)
    }, BPPARAM = setCepoBPPARAM(workers = workers, ...))
    
    ## Initialise p-value calculations
    listPvals <- vector("list", length = length(singleResult))
    names(listPvals) <- names(singleResult)
    for (i in names(singleResult)) {
        ## For each celltype and each gene, calculate the proportion of times that the
        ## gene exceeds the statistics value under bootstrap runs.
        listBinary <- lapply(listCepoOutputs, function(this_run) {
            singleResult[[i]] >= this_run[[i]][names(singleResult[[i]])]
        })
        listPvals[[i]] <- DelayedMatrixStats::colMeans2(do.call(rbind, listBinary))
        names(listPvals[[i]]) <- names(singleResult[[i]])
    }
    return(listPvals)
}

oneCepo <- function(exprsMat, cellTypes, exprsPct = NULL, workers = 1L, ...) {
    cts <- names(table(cellTypes))
    
    if (!is.null(exprsPct)) {
        meanPct.list <- list()
        for (i in seq_along(cts)) {
            idx <- which(cellTypes == cts[i])
            meanPct.list[[i]] <- (rowSums_withnames(exprsMat[, idx, drop = FALSE] > 
                0)/sum(cellTypes == cts[i])) > exprsPct
        }
        names(meanPct.list) <- cts
        keep <- rowSums_withnames(do.call(DelayedArray::cbind, meanPct.list)) == 
            length(cts)
        exprsMat <- exprsMat[keep, , drop = FALSE]
    }
    
    listIndexCelltypes <- lapply(cts, function(thisCelltypeName){which(cellTypes == thisCelltypeName)})
    
    segIdx.list <- BiocParallel::bplapply(
        X = listIndexCelltypes,
        FUN = function(thisCelltypeIndices){
            if(length(thisCelltypeIndices) < 20){ 
                ## If there is insufficient number of cells, we return NA's
                message("Less than 20 cells found in some cell types, returning NA")
                return(NA)
            } else {
                return(segIndex(exprsMat[, thisCelltypeIndices, drop = FALSE]))
            }}, BPPARAM = setCepoBPPARAM(workers = workers, ...))
    names(segIdx.list) <- cts
    
    segMat <- segIdxList2Mat(segIdx.list)
    segGenes <- consensusSegIdx(segMat)
    names(segGenes) <- colnames(segMat)
    result <- segGenes
    return(result)
}

segIndex <- function(mat) {
    nz <- rowMeans_withnames((mat != 0) + 0L)
    ms <- rowMeans_withnames(mat)
    sds <- rowSds_withnames(mat)
    cvs <- sds/ms
    
    x1 <- rank(nz)/(length(nz) + 1)
    x2 <- 1 - rank(cvs)/(length(cvs) + 1)
    
    segIdx <- rowMeans_withnames(DelayedArray::cbind(x1, x2))
    return(segIdx)
}

rowMeans_withnames <- function(mat) {
    result <- DelayedMatrixStats::rowMeans2(mat)
    names(result) <- rownames(mat)
    return(result)
}

rowSds_withnames <- function(mat) {
    result <- DelayedMatrixStats::rowSds(mat)
    names(result) <- rownames(mat)
    return(result)
}

rowSums_withnames <- function(mat) {
    result <- DelayedMatrixStats::rowMeans2(mat)
    names(result) <- rownames(mat)
    return(result)
}

segIdxList2Mat <- function(segIdx.list) {
    allGenes <- unique(unlist(lapply(segIdx.list, names)))
    segMat <- matrix(0, nrow = length(allGenes), ncol = length(segIdx.list))
    rownames(segMat) <- allGenes
    colnames(segMat) <- names(segIdx.list)
    
    for (i in seq_along(segIdx.list)) {
        si <- segIdx.list[[i]]
        segMat[names(si), i] <- si
    }
    return(segMat)
}

consensusSegIdx <- function(mat) {
    tt <- mat
    
    CIGs <- list()
    for (i in seq_len(ncol(tt))) {
        avgRank <- c()
        for (j in seq_len(ncol(tt))) {
            if (i == j) {
                next
            }
            avgRank <- cbind(avgRank, tt[, i] - tt[, j])
        }
        
        meanAvgRank <- DelayedMatrixStats::rowMeans2(avgRank)
        names(meanAvgRank) <- rownames(avgRank)
        CIGs[[i]] <- sort(meanAvgRank, decreasing = TRUE)
    }
    return(CIGs)
}

## Sorts every element of the list (assumed each element is a vector) by the names
## of the first element.
sortList <- function(listResult) {
    result <- lapply(listResult, function(thisElement) {
        thisElement[names(listResult[[1]])]
    })
    return(result)
}

#' @title Setting parallel params based on operating platform
#' @param workers Number of cores to use. Default to 1, which invokes `BiocParallel::SerialParam`.
#' For workers greater than 1, see the `workers` argument in `BiocParallel::MulticoreParam` and `BiocParallel::SnowParam`. 
#' @param ... Additional arguments passed to `BiocParallel::MulticoreParam` and `BiocParallel::SnowParam`. 
#' @examples 
#' # system.time(BiocParallel::bplapply(1:3, FUN = function(i){Sys.sleep(i)}, 
#' # BPPARAM = setCepoBPPARAM(workers = 1)))
#' # system.time(BiocParallel::bplapply(1:3, FUN = function(i){Sys.sleep(i)}, 
#' # BPPARAM = setCepoBPPARAM(workers = 3)))
setCepoBPPARAM = function(workers = 1L, ...){
    if(workers == 1){
        return(BiocParallel::SerialParam())
    } else if (.Platform$OS.type == "windows") {
        return(BiocParallel::SnowParam(workers = workers, ...))
    } else {
        return(BiocParallel::MulticoreParam(workers = workers, ...))
    } 
}
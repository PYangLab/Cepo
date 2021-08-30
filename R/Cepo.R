#' @title Computing Cepo cell identity genes
#' @param exprsMat expression matrix where columns denote cells and rows denote genes
#' @param cellTypes vector of cell type labels
#' @param exprsPct Percentage of lowly expressed genes to remove. Default to NULL to not remove any genes.
#' @param computePvalue Whether to compute p-values using bootstrap test. Default to NULL to not make computations.
#' Set this to an integer to set the number of bootstraps needed (recommend to be at least 100).
#' @param variability a character indicating the stability measure (CV, IQR, MAD, SD). Default is set to CV.
#' @param workers Number of cores to use. Default to 1, which invokes `BiocParallel::SerialParam`.
#' For workers greater than 1, see the `workers` argument in `BiocParallel::MulticoreParam` and `BiocParallel::SnowParam`. 
#' @param block vector of batch labels
#' @param minCelltype integer indicating the minimum number of cell types required in each batch
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
#' library(SingleCellExperiment)
#' data('cellbench', package = 'Cepo')
#' cellbench
#' cepoOutput <- Cepo(logcounts(cellbench), cellbench$celltype)
#' cepoOutput
Cepo <- function(exprsMat, cellTypes, exprsPct = NULL, computePvalue = NULL, variability = "CV",
                 workers = 1L, block = NULL, minCelltype = 3, ...) {
    stopifnot(ncol(exprsMat) == length(cellTypes))
    
    variability <- match.arg(variability, c("CV", "SD", "MAD", "IQR"),
                               several.ok = FALSE)
    
    if (!is.null(block)) {
        stopifnot(ncol(exprsMat) == length(block))
    }
        
    if(is.null(rownames(exprsMat))) {
        ## Add rownames if missing
        message("Gene names are missing in the input matrix, automatically adding `rownames`.")
        nGenes <- nrow(exprsMat)
        rownames(exprsMat) <- sprintf(paste0("gene%0", ceiling(log10(nGenes)) + 1L, "d"), seq_len(nGenes))
    }
    
    if(is.null(colnames(exprsMat))) {
        ## Add colnames if missing
        message("Cell names are missing in the input matrix, automatically adding `colnames`.")
        nCells <- ncol(exprsMat)
        colnames(exprsMat) <- sprintf(paste0("cell%0", ceiling(log10(nCells)) + 1, "d"), seq_len(nCells))
    }
    cellTypes <- as.character(cellTypes)
    
    if (!is.null(block)) {
        
        ## Select only batches with more than `minCelltype` number of cell types
        batches = names(which(rowSums(table(block, cellTypes) > 10) >= minCelltype))
        
        ## Run Cepo by batch
        batch_result <- lapply(batches, function(batch) {
            
            exprsMat_batch <- exprsMat[, block == batch]
            cellTypes_batch <- cellTypes[block == batch]
            
            ## A single run of oneCepo gives a list of output
            singleResult <- oneCepo(exprsMat = exprsMat_batch, cellTypes = cellTypes_batch, exprsPct = exprsPct, variability = variability,
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
                listPvals <- bootCepo(exprsMat = exprsMat_batch, cellTypes = cellTypes_batch, exprsPct = exprsPct, variability = variability,
                                      singleResult = singleResult, times = times, workers = workers, ...)
                ## The output has two components, one DataFrame of stats and another for p-values
                result <- list(stats = singleStatsResult, pvalues = S4Vectors::DataFrame(sortList(listPvals)))
            }  ## End else
            class(result) <- c("Cepo", class(result))
            return(result)
        })
        names(batch_result) = batches
        
        types = unique(unlist(lapply(batch_result, function(x) names(x$stats@listData))))
        idx = Reduce(intersect, lapply(batch_result, function(x) x$stats@rownames))
        
        averageCepo <- lapply(types, function(celltype) {
            mat <- do.call(cbind, lapply(batch_result, function(x) {
                x$stats@listData[[celltype]][idx]
            }))
            return(rowMeans(mat))
        })
        names(averageCepo) = types
        averageStatsResult <- S4Vectors::DataFrame(sortList(averageCepo))
        
        
        if (is.null(computePvalue)) {
        
            averageResult <- list(stats = averageStatsResult, pvalues = NULL)
            
        } else {
            
            averageCepoPvals <- lapply(types, function(celltype) {
                mat <- do.call(cbind, lapply(batch_result, function(x) {
                    x$pvalues@listData[[celltype]][idx]
                }))
                return(mat)
            })
            names(averageCepoPvals) = types
            averagePvalResult <- S4Vectors::DataFrame(sortList(averageCepoPvals))
            averageResult <- list(stats = averageStatsResult, pvalues = averagePvalResult)
        }
        
        batch_result$average <- averageResult
        
        return(batch_result)
        
    } else {
        
        ## A single run of oneCepo gives a list of output
        singleResult <- oneCepo(exprsMat = exprsMat, cellTypes = cellTypes, exprsPct = exprsPct, variability = variability,
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
            listPvals <- bootCepo(exprsMat = exprsMat, cellTypes = cellTypes, exprsPct = exprsPct, variability = variability,
                                  singleResult = singleResult, times = times, workers = workers, ...)
            ## The output has two components, one DataFrame of stats and another for p-values
            result <- list(stats = singleStatsResult, pvalues = S4Vectors::DataFrame(sortList(listPvals)))
        }  ## End else
        class(result) <- c("Cepo", class(result))
        return(result)
        
    }
    
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

bootCepo <- function(exprsMat, cellTypes, exprsPct, variability, singleResult, times, workers = 1L, ...) {
    ## Running multiple runs of Cepo based on bootstrap
    listCepoOutputs <- BiocParallel::bplapply(X = seq_len(times), FUN = function(i) {
        oneCepo(exprsMat = exprsMat, variability = variability, cellTypes = sample(cellTypes, replace = TRUE), exprsPct = exprsPct, workers = 1L)
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

oneCepo <- function(exprsMat, cellTypes, variability = variability, exprsPct = NULL, workers = 1L, ...) {
    cts <- names(table(cellTypes))
    
    if (!is.null(exprsPct)) {
        meanPct.list <- list()
        for (i in seq_along(cts)) {
            idx <- which(cellTypes == cts[i])
            meanPct.list[[i]] <- (rowSums_withnames(exprsMat[, idx, drop = FALSE] > 
                0)/sum(cellTypes == cts[i])) > exprsPct
        }
        names(meanPct.list) <- cts
        keep <- rowSums_withnames(do.call(DelayedArray::cbind, meanPct.list)) > 0
        exprsMat <- exprsMat[keep, , drop = FALSE]
    }
    
    listIndexCelltypes <- lapply(cts, function(thisCelltypeName){which(cellTypes == thisCelltypeName)})

    segIdx.list <- BiocParallel::bplapply(
        X = listIndexCelltypes,
        FUN = function(thisCelltypeIndices){
            if(length(thisCelltypeIndices) < 10){ 
                ## If there is insufficient number of cells, we return NA's
                message("Less than 10 cells found in some cell types, returning NA")
                return(NA)
            } else {
                return(segIndex(exprsMat[, thisCelltypeIndices, drop = FALSE], variability))
            }}, BPPARAM = setCepoBPPARAM(workers = workers, ...))
    names(segIdx.list) <- cts
    
    segMat <- segIdxList2Mat(segIdx.list)
    segGenes <- consensusSegIdx(segMat)
    result <- segGenes
    return(result)
}

segIndex <- function(mat, stability) {
    nz <- rowMeans_withnames((mat != 0) + 0L)
    ms <- rowMeans_withnames(mat)
    sds <- rowSds_withnames(mat)
    
    if (stability == "CV") {
        cvs <- sds/ms
        s <- cvs
    } else if (stability == "MAD") {
        s <- rowMads_withnames(mat)
    } else if (stability == "IQR") {
        s <- rowIqrs_withnames(mat)
    } else if (stability == "SD") {
        s <- sds
    }
    
    x1 <- rank(nz)/(length(nz) + 1)
    x2 <- 1 - rank(s)/(length(s) + 1)
    
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

rowMads_withnames <- function(mat) {
    result <- DelayedMatrixStats::rowMads(mat)
    names(result) <- rownames(mat)
    return(result)
}

rowIqrs_withnames <- function(mat) {
    result <- DelayedMatrixStats::rowIQRs(mat)
    names(result) <- rownames(mat)
    return(result)
}

rowSums_withnames <- function(mat) {
    result <- DelayedMatrixStats::rowSums2(mat)
    names(result) <- rownames(mat)
    return(result)
}

segIdxList2Mat <- function(segIdx.list) {
    allGenes <- unique(unlist(lapply(segIdx.list, names)))
    
    segMat <- matrix(NA, nrow = length(allGenes), ncol = length(segIdx.list))
    rownames(segMat) <- allGenes
    colnames(segMat) <- names(segIdx.list)

    for (i in seq_along(segIdx.list)) {
        si <- segIdx.list[[i]]
        segMat[names(si), i] <- si
    }
    return(segMat)
}

consensusSegIdx <- function(mat) {
    CIGs2 <- lapply(seq_len(ncol(mat)),
                   function(i){
                       meanAvgRank <- DelayedMatrixStats::rowMeans2(-(mat[, -i, drop = FALSE] - mat[, i, drop = TRUE]))
                       names(meanAvgRank) <- rownames(mat)
                       return(sort(meanAvgRank, decreasing = TRUE))
                   })
    names(CIGs2) <- colnames(mat)
    return(CIGs2)
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


#' @title Computing Cepo cell identity genes
#' @param exprsMat Expression matrix where columns denote cells and rows 
#' denote genes
#' @param cellTypes Vector of cell type labels
#' @param exprsPct Percentage of lowly expressed genes to remove. 
#' Default to NULL to not remove any genes.
#' @param prefilter_sd Numeric value indicating threshold relating to standard 
#' deviation of genes. Used with prefilter_zeros. 
#' @param prefilter_pzeros Numeric value indicating threshold relating to the 
#' percentage of zero expression of genes. Used with prefilter_sd. 
#' @param logfc Numeric value indicating the threshold of log fold-change 
#' to use to filter genes.
#' @param computePvalue Whether to compute p-values using bootstrap test. 
#' Default to NULL to not make computations.
#' Set this to an integer to set the number of bootstraps needed 
#' (recommend to be at least 100).
#' @param computeFastPvalue Logical vector indicating whether to perform a 
#' faster version of p-value calculation. Set to TRUE by default.
#' @param variability A character indicating the stability measure 
#' (CV, IQR, MAD, SD). Default is set to CV.
#' @param workers Number of cores to use. Default to 1, which invokes 
#' `BiocParallel::SerialParam`.
#' For workers greater than 1, see the `workers` argument in 
#' `BiocParallel::MulticoreParam` and `BiocParallel::SnowParam`. 
#' @param method Character indicating the method for integration 
#' the two stability measures. By default this is set to
#' 'weightedMean' with equal weights.
#' @param weight Vector of two values indicating the weights for each stability 
#' measure. By default this value is c(0.5, 0.5).
#' @param block Vector of batch labels
#' @param minCelltype Integer indicating the minimum number of cell types 
#' required in each batch
#' @param minCells Integer indicating the minimum number of cells required 
#' within a cell type
#' @param ... Additional arguments passed to `BiocParallel::MulticoreParam` 
#' and `BiocParallel::SnowParam`. 
#' @importFrom DelayedMatrixStats rowSums2 rowMeans2 rowSds
#' @importFrom DelayedArray cbind
#' @importFrom HDF5Array HDF5Array
#' @importFrom BiocParallel SerialParam MulticoreParam SnowParam bplapply
#' @importFrom stats pnorm pchisq
#' @import dplyr
#' @return Returns a list of key genes.
#' @description ExprsMat accepts various matrix objects, 
#' including DelayedArray and HDF5Array for
#' out-of-memory computations. See vignette.
#' @export
#' @examples
#' library(SingleCellExperiment)
#' data('cellbench', package = 'Cepo')
#' cellbench
#' cepoOutput <- Cepo(logcounts(cellbench), cellbench$celltype)
#' cepoOutput
#' 
Cepo <- function(exprsMat, cellTypes, minCells = 20, minCelltype = 3, exprsPct = NULL, 
                 prefilter_sd = NULL, prefilter_pzero = NULL,
                 logfc = NULL, computePvalue = NULL, computeFastPvalue = TRUE,
                 variability = "CV", method = "weightedMean",
                 weight = c(0.5, 0.5), workers = 1L, block = NULL, ...) {
    
    stopifnot(ncol(exprsMat) == length(cellTypes))
    
    variability <- match.arg(variability, c("CV", "SD", "MAD", "IQR"),
                             several.ok = FALSE)
    
    method <- match.arg(method, c("weightedMean", "Stouffer", "OSP", "Fisher",
                                  "maxP"), several.ok = FALSE)
    
    if (!is.null(block)) {
        stopifnot(ncol(exprsMat) == length(block))
    }
    
    if (is.null(rownames(exprsMat))) {
        ## Add rownames if missing
        message("Gene names are missing in the input matrix, 
                automatically adding `rownames`.")
        nGenes <- nrow(exprsMat)
        rownames(exprsMat) <- sprintf(paste0("gene%0", ceiling(log10(nGenes)) +
                                                 1L, "d"), seq_len(nGenes))
    }
    
    if (is.null(colnames(exprsMat))) {
        ## Add colnames if missing
        message("Cell names are missing in the input matrix, 
                automatically adding `colnames`.")
        nCells <- ncol(exprsMat)
        colnames(exprsMat) <- sprintf(paste0("cell%0", ceiling(log10(nCells)) +
                                                 1, "d"), seq_len(nCells))
    }
    cellTypes <- as.character(cellTypes)
    
    if (!(is.null(prefilter_sd) && is.null(prefilter_pzero))) {
        
        pzeros = rowSums_withnames(exprsMat == 0)/ncol(exprsMat)
        
        rowSdsRes <- lapply(names(table(cellTypes)), function(i) {
            rowSds_withnames(exprsMat[, cellTypes == i] )
        })
        rowSdsRes <- do.call(cbind, rowSdsRes)
        rowSdsResNorm <- t(apply(rowSdsRes, 1, function(i) {
            mean((i - mean(i))*10^16)
        }))[1,]
        
        if (is.null(prefilter_sd)) { prefilter_sd <- max(rowSdsResNorm)}
        if (is.null(prefilter_pzero)) { prefilter_pzero <- max(pzeros)}
        
        ## Remove genes where there is HIGH expression and LOW difference in sd
        genesFilter = names(pzeros[which(abs(rowSdsResNorm) < prefilter_sd 
                                         & pzeros < prefilter_pzero)])
        
        exprsMat = exprsMat[!rownames(exprsMat) %in% genesFilter,]
        message(paste0("Prefiltering ", length(genesFilter), " genes...."))
        
    }
    
    if (!is.null(block)) {
        
        block <- as.character(block)
        
        ## Select only batches with more than `minCelltype` number of
        ## cell types
        batches <- names(which(rowSums(table(block, cellTypes) > minCells) >=
                                   minCelltype))
        
        ## Run Cepo by batch
        batch_result <- lapply(batches, function(batch) {
            
            exprsMat_batch <- exprsMat[, block == batch]
            cellTypes_batch <- cellTypes[block == batch]
            
            singleBatch <- singleBatchCepo(exprsMat = exprsMat_batch, cellTypes = cellTypes_batch,
                                           minCells = minCells, exprsPct = exprsPct, logfc = logfc,
                                           method = method, weight = weight, variability = variability,
                                           computePvalue = computePvalue, 
                                           computeFastPvalue = computeFastPvalue,
                                           workers = workers, ...)
            
            return(singleBatch)
            
        })
        names(batch_result) <- batches
        
        types <- unique(unlist(lapply(batch_result, function(x) {
            colnames(x$stats)
        })))
        idx <- Reduce(intersect, lapply(batch_result, function(x) {
            rownames(x$stats)
        }))
        
        averageCepo <- lapply(types, function(celltype) {
            mat <- do.call(cbind, lapply(batch_result, function(x) {
                x$stats[idx ,celltype]
            }))
            return(rowMeans(mat))
        })
        names(averageCepo) <- types
        averageStatsResult <- S4Vectors::DataFrame(sortList(averageCepo))
        
        
        if (is.null(computePvalue)) {
            
            averageResult <- list(stats = averageStatsResult, pvalues = NULL)
            
        } else {
            
            averageCepoPvals <- lapply(types, function(celltype) {
                mat <- do.call(cbind, lapply(batch_result, function(x) {
                    x$pvalues[idx, celltype]
                }))
                return(mat)
            })
            names(averageCepoPvals) <- types
            averagePvalResult <- S4Vectors::DataFrame(sortList(averageCepoPvals))
            averageResult <- list(stats = averageStatsResult, pvalues = averagePvalResult)
        }
        
        batch_result$average <- averageResult
        
        return(batch_result)
        
    } else {
        
        singleBatch <- singleBatchCepo(exprsMat = exprsMat, cellTypes = cellTypes,
                                       minCells = minCells, exprsPct = exprsPct, logfc = logfc, method = method,
                                       weight = weight, variability = variability, 
                                       computePvalue = computePvalue,
                                       computeFastPvalue = computeFastPvalue,
                                       workers = workers, ...)
        
        return(singleBatch)
    }
    
}

print.Cepo <- function(x) {
    cat("Computed statistics: \n \n")
    print(x$stats)
    cat("Computed p-values: \n")
    if (is.null(x$pvalues)) {
        cat("Note: a valid value for `computePvalue` argument is needed to 
            get p-values when running the `Cepo` function")
    } else {
        print(x$pvalues)
    }
}

bootCepo <- function(exprsMat, cellTypes, minCells, exprsPct, logfc, variability,
                     method, weight, singleResult, times, workers = 1L, ...) {
    ## Running multiple runs of Cepo based on bootstrap
    listCepoOutputs <- BiocParallel::bplapply(X = seq_len(times), FUN = function(i) {
        oneCepo(exprsMat = exprsMat, minCells = minCells, variability = variability,
                logfc = logfc, method = method, weight = weight, cellTypes = sample(cellTypes,
                                                                                    replace = TRUE), exprsPct = exprsPct, workers = 1L)
    }, BPPARAM = setCepoBPPARAM(workers = workers, ...))
    
    ## Initialise p-value calculations
    listPvals <- vector("list", length = length(singleResult))
    names(listPvals) <- names(singleResult)
    for (i in names(singleResult)) {
        ## For each celltype and each gene, calculate the proportion
        ## of times that the gene exceeds the statistics value under
        ## bootstrap runs.
        listBinary <- lapply(listCepoOutputs, function(this_run) {
            singleResult[[i]] <= this_run[[i]][names(singleResult[[i]])]
        })
        listPvals[[i]] <- DelayedMatrixStats::colMeans2(do.call(rbind,
                                                                listBinary))
        names(listPvals[[i]]) <- names(singleResult[[i]])
    }
    return(listPvals)
}

bootFastCepo <- function(exprsMat, cellTypes, minCells, exprsPct, logfc, variability,
                         method, weight, singleResult, times, workers = 1L, ...) {
    
    ## 
    cts <- names(table(cellTypes))
    
    geneNames <- sort(names(singleResult[[1]]))
    singleResult <- lapply(singleResult, function(x) x[geneNames])
    
    ## Running multiple runs of Cepo based on bootstrap
    sampled_cepo_stats <- BiocParallel::bplapply(X = seq_len(times), FUN = function(i) {
        sampleResult <- oneCepo(exprsMat = exprsMat, minCells = minCells, variability = variability,
                                logfc = logfc, method = method, weight = weight, cellTypes = sample(cellTypes,
                                                                                                    replace = TRUE), exprsPct = exprsPct, workers = 1L)
        
        sampleResult <- S4Vectors::DataFrame(sortList(sampleResult))
        return(sampleResult)
    }, BPPARAM = setCepoBPPARAM(workers = workers, ...))
    
    sampled_cepo_stats = sampled_cepo_stats %>% 
        purrr::map(.f = ~ .x[geneNames,])
    
    celltype_norm_stats = purrr::map(
        .x = cts, 
        .f = function(this_celltype){
            this_celltype_stats_matrix = purrr::map(.x = sampled_cepo_stats, .f = function(this_sampled_stats){
                this_sampled_stats[,this_celltype]
            }) %>% do.call(cbind, .)
            this_celltype_stats_mean = Cepo:::rowMeans_withnames(this_celltype_stats_matrix)
            this_celltype_stats_sd = Cepo:::rowSds_withnames(this_celltype_stats_matrix)
            return(data.frame(mean = this_celltype_stats_mean, sd = this_celltype_stats_sd))
        })
    names(celltype_norm_stats) = cts
    
    celltype_norm_pvalues = purrr::map2(
        .x = celltype_norm_stats, 
        .y = cts,
        .f = function(this_celltype_stats, this_celltype_name){
            pnorm(q = singleResult[[this_celltype_name]], 
                  mean = this_celltype_stats$mean,
                  sd = this_celltype_stats$sd, lower.tail = FALSE)
        }) %>% 
        do.call(cbind, .)
    colnames(celltype_norm_pvalues) = cts
    rownames(celltype_norm_pvalues) = geneNames
    
    return(list(approx_pvalues = celltype_norm_pvalues, 
                mean_sd_stats = celltype_norm_stats))
}

oneCepo <- function(exprsMat, cellTypes, minCells = minCells, variability = variability,
                    exprsPct = NULL, logfc = NULL, method = method, weight = weight, workers = 1L,
                    ...) {
    cts <- names(table(cellTypes))
    
    if (!is.null(exprsPct)) {
        meanPct.list <- list()
        for (i in seq_along(cts)) {
            idx <- which(cellTypes == cts[i])
            meanPct.list[[i]] <- (rowSums_withnames(exprsMat[, idx, drop = FALSE] >
                                                        0)/sum(cellTypes == cts[i])) > exprsPct
        }
        names(meanPct.list) <- cts
        keep <- rowSums_withnames(do.call(DelayedArray::cbind, meanPct.list)) >
            0
        exprsMat <- exprsMat[keep, , drop = FALSE]
    }
    
    if (!is.null(logfc)) {
        logfc.list <- list()
        for (i in seq_along(cts)) {
            idx <- which(cellTypes == cts[i])
            logfc.list[[i]] <- abs(log(x = rowMeans_withnames(mat = expm1(x = exprsMat[,
                                                                                       idx, drop = FALSE])) + 1, base = 2) - log(x = rowMeans_withnames(mat = expm1(x = exprsMat[,
                                                                                                                                                                                 -idx, drop = FALSE])) + 1, base = 2)) > logfc
        }
        names(logfc.list) <- cts
        keep <- rowSums_withnames(do.call(DelayedArray::cbind, logfc.list)) >
            0
        exprsMat <- exprsMat[keep, , drop = FALSE]
    }
    
    listIndexCelltypes <- lapply(cts, function(thisCelltypeName) {
        which(cellTypes == thisCelltypeName)
    })
    
    segIdx.list <- BiocParallel::bplapply(X = listIndexCelltypes, FUN = function(thisCelltypeIndices) {
        if (length(thisCelltypeIndices) < minCells) {
            ## If there is insufficient number of cells, we return
            ## NA's
            message("Less than ", minCells, " cells found in some cell types, returning NA")
            return(NA)
        } else {
            return(segIndex(exprsMat[, thisCelltypeIndices, drop = FALSE],
                            stability = variability, method = method, weight = weight))
        }
    }, BPPARAM = setCepoBPPARAM(workers = workers, ...))
    names(segIdx.list) <- cts
    
    segMat <- segIdxList2Mat(segIdx.list)
    segGenes <- consensusSegIdx(segMat)
    result <- segGenes
    return(result)
}

singleBatchCepo <- function(exprsMat, cellTypes, minCells = minCells, variability = variability,
                            exprsPct = NULL, logfc = NULL, method = method, weight = weight, 
                            computePvalue = NULL, computeFastPvalue = computeFastPvalue,
                            workers = 1L, ...) {
    
    singleResult <- oneCepo(exprsMat = exprsMat, cellTypes = cellTypes,
                            minCells = minCells, exprsPct = exprsPct, logfc = logfc, method = method,
                            weight = weight, variability = variability, workers = workers,
                            ...)
    ## Export Cepo outputs as a DataFrame
    singleStatsResult <- S4Vectors::DataFrame(sortList(singleResult))
    
    if (is.null(computePvalue)) {
        ## If no need to compute p-values, then that is it.
        result <- list(stats = singleStatsResult, pvalues = NULL)
    } else if (computeFastPvalue == TRUE) {
        ## Coerce the input to an integer
        times <- as.integer(computePvalue)
        ## Pass onto a boot function for computation
        listPvals <- bootFastCepo(exprsMat = exprsMat, cellTypes = cellTypes,
                                  minCells = minCells, exprsPct = exprsPct, weight = weight,
                                  method = method, logfc = logfc, variability = variability,
                                  singleResult = singleResult, times = times, workers = workers,
                                  ...)
        listPvals$mean_sd_stats <- lapply(listPvals$mean_sd_stats, function(x) {
            S4Vectors::DataFrame(x[rownames(singleStatsResult), ])
        })
        ## The output has two components, one DataFrame of stats and
        ## another for p-values
        result <- list(stats = singleStatsResult, 
                       pvalues = S4Vectors::DataFrame(listPvals$approx_pvalues[rownames(singleStatsResult),]),
                       geneStatistics = listPvals$mean_sd_stats)
    } else {
        ## Coerce the input to an integer
        times <- as.integer(computePvalue)
        ## Pass onto a boot function for computation
        listPvals <- bootCepo(exprsMat = exprsMat, cellTypes = cellTypes,
                              minCells = minCells, exprsPct = exprsPct, weight = weight,
                              method = method, logfc = logfc, variability = variability,
                              singleResult = singleResult, times = times, workers = workers,
                              ...)
        ## The output has two components, one DataFrame of stats and
        ## another for p-values
        result <- list(stats = singleStatsResult, pvalues = S4Vectors::DataFrame(sortList(listPvals)))
    } ## End else
    class(result) <- c("Cepo", class(result))
    return(result)
    
}

segIndex <- function(mat, stability, method, weight) {
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
    
    # segIdx <- rowMeans_withnames(DelayedArray::cbind(x1, x2))
    
    if (method == "weightedMean") {
        segIdx <- rowMeansWeighted_withnames(DelayedArray::cbind(x1, x2),
                                             weight)
    } else if (method == "Stouffer") {
        segIdx <- apply(-DelayedArray::cbind(x1, x2), 1, geneStats, method = "Stouffer")
    } else if (method == "OSP") {
        segIdx <- apply(-DelayedArray::cbind(x1, x2), 1, geneStats, method = "OSP")
    } else if (method == "Fisher") {
        segIdx <- apply(-DelayedArray::cbind(x1, x2), 1, geneStats, method = "Fisher")
    } else if (method == "maxP") {
        segIdx <- apply(-DelayedArray::cbind(x1, x2), 1, geneStats, method = "maxP")
    }
    
    return(segIdx)
}

rowMeans_withnames <- function(mat) {
    result <- DelayedMatrixStats::rowMeans2(mat)
    names(result) <- rownames(mat)
    return(result)
}

rowMeansWeighted_withnames <- function(mat, weight) {
    result <- DelayedMatrixStats::rowWeightedMeans(mat, w = weight)
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
    CIGs2 <- lapply(seq_len(ncol(mat)), function(i) {
        meanAvgRank <- DelayedMatrixStats::rowMeans2(-(mat[, -i, drop = FALSE] -
                                                           mat[, i, drop = TRUE]))
        names(meanAvgRank) <- rownames(mat)
        return(sort(meanAvgRank, decreasing = TRUE))
    })
    names(CIGs2) <- colnames(mat)
    return(CIGs2)
}

## Sorts every element of the list (assumed each element is a vector)
## by the names of the first element.
sortList <- function(listResult) {
    result <- lapply(listResult, function(thisElement) {
        thisElement[names(listResult[[1]])]
    })
    return(result)
}

geneStats <- function(Tstat, method = "OSP") {
    pvalue <- 0
    if (method == "Stouffer") {
        pvalue <- stats::pnorm(sum(Tstat), 0, sqrt(length(Tstat)), lower.tail = FALSE)
    } else if (method == "OSP") {
        p <- stats::pnorm(Tstat, lower.tail = TRUE)
        pvalue <- stats::pchisq(-2 * sum(log(p)), 2 * length(p), lower.tail = TRUE)
    } else if (method == "Fisher") {
        p <- stats::pnorm(Tstat, lower.tail = FALSE)
        pvalue <- stats::pchisq(-2 * sum(log(p)), 2 * length(p), lower.tail = FALSE)
    } else if (method == "maxP") {
        pvalue <- stats::pnorm(max(Tstat), lower.tail = FALSE)
    }
    return(pvalue)
}


#' @title Setting parallel params based on operating platform
#' @param workers Number of cores to use. Default to 1, which invokes 
#' `BiocParallel::SerialParam`.
#' For workers greater than 1, see the `workers` argument in 
#' `BiocParallel::MulticoreParam` and `BiocParallel::SnowParam`. 
#' @param ... Additional arguments passed to `BiocParallel::MulticoreParam` 
#' and `BiocParallel::SnowParam`. 
#' @return Parameters for parallel computing depending on OS
#' @examples 
#' # system.time(BiocParallel::bplapply(1:3, FUN = function(i){Sys.sleep(i)}, 
#' # BPPARAM = setCepoBPPARAM(workers = 1)))
#' # system.time(BiocParallel::bplapply(1:3, FUN = function(i){Sys.sleep(i)}, 
#' # BPPARAM = setCepoBPPARAM(workers = 3)))
setCepoBPPARAM <- function(workers = 1L, ...) {
    if (workers == 1) {
        return(BiocParallel::SerialParam())
    } else if (.Platform$OS.type == "windows") {
        return(BiocParallel::SnowParam(workers = workers, ...))
    } else {
        return(BiocParallel::MulticoreParam(workers = workers, ...))
    }
}



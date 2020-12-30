# library(DelayedArray)
# library(DelayedMatrixStats)
# library(HDF5Array)
# ## DelayedArray:::set_verbose_block_processing(TRUE)
# ## DelayedArray::setAutoBlockSize(size = 1e8)
# ########### Simulation ##############
# set.seed(1234)
# n <- 40000 ## genes, rows
# p <- 100 ## cells, cols
# exprsMat <- matrix(rpois(n * p, lambda = 5), nrow = n)
# rownames(exprsMat) <- paste0("gene", 1:n)
# colnames(exprsMat) <- paste0("cell", 1:p)
# cellTypes <- sample(letters[1:2], size = p, replace = TRUE)
# onecore_output <- Cepo(exprsMat = exprsMat, cellTypes = cellTypes, computePvalue = 100)
# twocore_output <- Cepo(exprsMat = exprsMat, cellTypes = cellTypes, computePvalue = 100, workers = 2)

# system.time(onecore_output <- Cepo(exprsMat = exprsMat, cellTypes = cellTypes, computePvalue = 100))
# system.time(twocore_output <- Cepo(exprsMat = exprsMat, cellTypes = cellTypes, computePvalue = 100, workers = 4, progressbar = TRUE))
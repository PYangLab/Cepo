library(DelayedArray)
library(DelayedMatrixStats)
library(HDF5Array)
## DelayedArray:::set_verbose_block_processing(TRUE)
## DelayedArray::setAutoBlockSize(size = 1e8)
########### Simulation ##############
set.seed(1234)
n <- 50 ## genes, rows
p <- 100 ## cells, cols
exprsMat <- matrix(rpois(n * p, lambda = 5), nrow = n)
rownames(exprsMat) <- paste0("gene", 1:n)
colnames(exprsMat) <- paste0("cell", 1:p)
cellTypes <- sample(letters[1:3], size = p, replace = TRUE)
da_exprsMat <- DelayedArray::DelayedArray(DelayedArray::realize(exprsMat, "HDF5Array"))
########### Normal matrix outputs ####################
cepo_output <- Cepo(exprsMat = exprsMat, cellTypes = cellTypes)
########### End normal matrix outputs ####################

########### Delayed matrix outputs ####################
da_cepo_output <- Cepo(exprsMat = da_exprsMat, cellTypes = cellTypes)
########### End delayed matrix outputs ####################

########### Compare outputs ####################
testthat::expect_identical(cepo_output, da_cepo_output)

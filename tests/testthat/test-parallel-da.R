library(DelayedArray)
library(DelayedMatrixStats)
library(HDF5Array)
library(BiocParallel)
BPPARAM <- if (.Platform$OS.type == "windows") {
  BiocParallel::SnowParam(workers = 2)
} else {
  BiocParallel::MulticoreParam(workers = 2)
}
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
########### Comparing matrix outputs ####################
cepo_1core_output <- Cepo(exprsMat = exprsMat, cellTypes = cellTypes)
da_cepo_1core_output <- Cepo(exprsMat = da_exprsMat, cellTypes = cellTypes)
DelayedArray::setAutoBPPARAM(BPPARAM = BPPARAM) ## Setting two cores for computation
da_cepo_2core_output <- Cepo(exprsMat = da_exprsMat, cellTypes = cellTypes)
DelayedArray::setAutoBPPARAM(BPPARAM = SerialParam()) ## Revert back to only one core
testthat::expect_identical(cepo_1core_output, da_cepo_1core_output)
testthat::expect_identical(cepo_1core_output, da_cepo_2core_output)

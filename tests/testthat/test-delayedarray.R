set.seed(1234)
n = 100 ## genes, rows
p = 1000 ## cells, cols
exprsMat = matrix(rpois(n*p, lambda = 5), nrow = n)
rownames(exprsMat) = paste0("gene", 1:n)
colnames(exprsMat) = paste0("cell", 1:p)
cellTypes = sample(letters[1:3], size = p, replace = TRUE)
mat_output = Cepo(exprsMat = exprsMat, cellTypes = cellTypes)

library(DelayedArray)
library(DelayedMatrixStats)
library(HDF5Array)
## DelayedArray:::set_verbose_block_processing(TRUE)
## DelayedArray::setAutoBlockSize(size = 1e8)
exprsMat_da = DelayedArray::DelayedArray(DelayedArray::realize(exprsMat, "HDF5Array"))
da_output = Cepo(exprsMat = exprsMat_da, cellTypes = cellTypes)
testthat::expect_identical(mat_output, da_output)

library(DelayedArray)
library(DelayedMatrixStats)
library(HDF5Array)
## DelayedArray:::set_verbose_block_processing(TRUE)
## DelayedArray::setAutoBlockSize(size = 1e8)
########### Simulation ##############
set.seed(1234)
n = 50 ## genes, rows
p = 100 ## cells, cols
exprsMat = matrix(rpois(n*p, lambda = 5), nrow = n)
rownames(exprsMat) = paste0("gene", 1:n)
colnames(exprsMat) = paste0("cell", 1:p)
cellTypes = sample(letters[1:3], size = p, replace = TRUE)
da_exprsMat = DelayedArray::DelayedArray(DelayedArray::realize(exprsMat, "HDF5Array"))
########### Normal matrix outputs ####################
cepo_output = Cepo(exprsMat = exprsMat, cellTypes = cellTypes)
limma_output = doLimma(exprsMat = exprsMat, cellTypes = cellTypes)
voom_output = doVoom(exprsMat = exprsMat, cellTypes = cellTypes)
ttest_output = doTtest(exprsMat = exprsMat, cellTypes = cellTypes)
wilcoxon_output = doWilcoxon(exprsMat = exprsMat, cellTypes = cellTypes)
########### End normal matrix outputs ####################

########### Delayed matrix outputs ####################
da_cepo_output = Cepo(exprsMat = da_exprsMat, cellTypes = cellTypes)
da_limma_output = doLimma(exprsMat = da_exprsMat, cellTypes = cellTypes)
da_voom_output = doVoom(exprsMat = da_exprsMat, cellTypes = cellTypes)
da_ttest_output = doTtest(exprsMat = da_exprsMat, cellTypes = cellTypes)

## Expecting error due to non-support of DA
da_wilcoxon_output = expect_error(doWilcoxon(exprsMat = da_exprsMat, cellTypes = cellTypes))
########### End delayed matrix outputs ####################

########### Compare outputs ####################
testthat::expect_identical(cepo_output, da_cepo_output)
testthat::expect_identical(limma_output, da_limma_output)
testthat::expect_identical(voom_output, da_voom_output)
testthat::expect_identical(ttest_output, da_ttest_output)

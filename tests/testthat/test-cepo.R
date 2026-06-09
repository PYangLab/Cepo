set.seed(1234)
n <- 50 ## genes, rows
p <- 100 ## cells, cols
exprsMat <- matrix(rpois(n * p, lambda = 5), nrow = n)
rownames(exprsMat) <- paste0('gene', 1:n)
colnames(exprsMat) <- paste0('cell', 1:p)
cellTypes <- sample(letters[1:3], size = p, replace = TRUE)

cepoOutput_usual <- Cepo(exprsMat = exprsMat, cellTypes = cellTypes)

## Without computePvalues, the pvalues should be NULL
expect_null(cepoOutput_usual$pvalues)

## Input without rownames and colnames
expect_message(
  Cepo(exprsMat = unname(exprsMat), cellTypes = cellTypes)
)

## Input using factors
cepoOutput_fct <- Cepo(exprsMat = exprsMat, cellTypes = cellTypes)
expect_identical(cepoOutput_usual, cepoOutput_fct)

## Testing exprsPct argument
cepoOutput_perc <- Cepo(exprsMat = exprsMat, cellTypes = cellTypes,
                        exprsPct = 0.5)

## Bootstrap ouput
cepoOutput_boot <- Cepo(exprsMat = exprsMat, cellTypes = cellTypes, computePvalue = 10)
## Testing print method
print(cepoOutput_boot)

## P-value resampling should preserve cell-type sizes so classes just above
## minCells do not become all-NA during null runs.
set.seed(1234)
imbalancedCellTypes <- c(rep("a", 21), rep("b", 30), rep("c", 44))
imbalancedExprsMat <- matrix(rpois(n * length(imbalancedCellTypes), lambda = 5), nrow = n)
rownames(imbalancedExprsMat) <- paste0("gene", 1:n)
colnames(imbalancedExprsMat) <- paste0("cell", seq_along(imbalancedCellTypes))
cepoOutput_imbalanced <- Cepo(
  exprsMat = imbalancedExprsMat,
  cellTypes = imbalancedCellTypes,
  computePvalue = 10
)
expect_false(all(is.na(as.matrix(cepoOutput_imbalanced$pvalues))))
expect_false(any(vapply(
  cepoOutput_imbalanced$geneStatistics,
  function(x) all(is.na(x$mean)) || all(is.na(x$sd)),
  logical(1)
)))

## Testing insufficient number of cells 
# cellTypes_insufficient = c(rep("a", 10), rep("b", p - 10))
# expect_message(
#   cepoOutput_insufficient <- Cepo(exprsMat = exprsMat, cellTypes = cellTypes_insufficient)
# )

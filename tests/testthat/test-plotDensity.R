library(SingleCellExperiment)
data('cellbench', package = 'Cepo')
cellbench
ds_res <- Cepo(logcounts(cellbench), cellbench$celltype)

expect_message(plotDensities(
  x = cellbench,
  cepoOutput = ds_res,
  assay = 'logcounts',
  celltypeColumn = 'celltype'
))

plotDensities(
  x = cellbench,
  cepoOutput = ds_res,
  genes = c('PLTP', 'CPT1C', 'MEG3', 'SYCE1', 'MICOS10P3', 'HOXB7'),
  assay = 'logcounts',
  celltypeColumn = 'celltype'
)
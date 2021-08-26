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

plotDensities(
  x = cellbench,
  cepoOutput = ds_res,
  genes = c('PLTP', 'CPT1C', 'MEG3', 'SYCE1', 'MICOS10P3', 'HOXB7'),
  assay = 'logcounts',
  plotType = "density",
  celltypeColumn = 'celltype'
)

expect_null(
  plotDensities(
    x = cellbench,
    cepoOutput = ds_res,
    assay = 'not_exist',
    celltypeColumn = 'celltype'
  ))

expect_null(
  plotDensities(
    x = cellbench,
    cepoOutput = ds_res,
    assay = 'logcounts',
    celltypeColumn = 'not_exist'
  ))

## Too many genes
expect_null(
  plotDensities(
    x = cellbench,
    cepoOutput = ds_res,
    genes = rownames(ds_res$stats)[1:11],
    assay = 'logcounts',
    celltypeColumn = 'celltype'
  ))


## genes missing in data
expect_null(
  plotDensities(
    x = cellbench,
    cepoOutput = ds_res,
    genes = c("not_exist"),
    assay = 'logcounts',
    celltypeColumn = 'celltype'
  ))

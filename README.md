# CCIGSVA

CCIGSVA calculates sample-level pathway enrichment of cell-cell interactions
inferred from single-cell transcriptome data by applying GSVA to ligand-receptor
interaction scores for each source-target cell-type pair.

The CCI-GSVA method was introduced in
[Jeong et al. (2023)](https://doi.org/10.1038/s41598-023-46350-2).

## Installation

Requires R >= 4.3 and GSVA >= 1.50.

```r
install.packages(c("BiocManager", "remotes"))
BiocManager::install(c("GSVA", "SummarizedExperiment", "S4Vectors", "BiocParallel"))
remotes::install_github("Khreat0205/CCIGSVA", upgrade = "never")
```

## Quick start

```r
library(CCIGSVA)
ex <- example_ccigsva()
result <- cci_gsva(ex$cci, ex$lr_sets)
SummarizedExperiment::assay(result, "ccigsva")
SummarizedExperiment::rowData(result)
```

The assay has source-target-pathway features in rows and samples in columns.
The example uses synthetic interaction scores.

## Use CellChat results

Start with one CSV per sample exported using `CellChat::subsetCommunication()`
and a GMT file mapping ligand-receptor interactions to signaling pathways.

```r
samples <- c("sample1", "sample2", "sample3") # complete cohort, including empty exports
tables <- setNames(lapply(samples, function(s) {
  read.csv(file.path("Res_CellChat", paste0("cellchat_", s, ".csv")))
}), samples)
cci <- from_cellchat(tables)
lr_sets <- read_lr_gmt("lr_pathways.gmt")

# Cell counts must come from the cell metadata, not significant CCI rows.
cell_counts <- table(obj$Anno1, obj$Individual)
coverage <- filter_cell_types(cell_counts[, samples, drop = FALSE],
                              min_cells = 10, min_coverage = 0.8)
result <- cci_gsva(cci, lr_sets, samples = samples,
                   cell_types = coverage$cell_type[coverage$keep], min_size = 2)
```

Here, `obj` is a Seurat object with cell-type labels in `Anno1` and sample IDs
in `Individual`. Replace these names with the corresponding metadata columns.

- `min_cells = 10`: minimum cells for a cell type to be present in a sample.
- `min_coverage = 0.8`: retain cell types present in at least 80% of samples.
  Use `0.7` for 70% or `1` for all samples.
- `min_size = 2`: minimum matched, nonconstant LR interactions per pathway.

Alternatively, construct sets directly from a CellChat database table:

```r
lr_sets <- prepare_lr_sets(CellChatDB.human$interaction)
```

CellChat is not required to install CCIGSVA. `from_cellchat()` accepts exported
communication tables.
Other tools can supply a long table through `prepare_cci()`, provided scores
and LR identifiers have a compatible interpretation.

## Processing rules

- Union LR interactions across samples within each observed source-target pair.
- Fill absent entries with zero; reject explicit NA/Inf scores and duplicate keys.
- Preserve every explicitly supplied sample, including samples with no CCI rows.
- Apply no score normalization, transformation or scaling.
- Remove constant LR rows and apply pathway minimum size after matching.
- Use Gaussian GSVA with `tau = 1`, `maxDiff = TRUE`, `absRanking = FALSE`.
- Exclude sets covering the full variable LR universe (no background), with an audit report.

Zero filling is an input convention. A missing cell type, an untested edge and
a nonsignificant edge are not necessarily biologically equivalent. Choose the
coverage threshold upstream and retain the complete sample list. GSVA scores
depend on the input cohort; rerunning on a different cohort can change scores.

```r
S4Vectors::metadata(result)$pair_report
S4Vectors::metadata(result)$effective_lr_sets
S4Vectors::metadata(result)$versions
```

## Reference and license

Jeong K, Kim Y, Jeon J, Kim K. Subtyping of COVID-19 samples based on
cell-cell interaction in single cell transcriptomes. *Scientific Reports*
13, 19629 (2023). DOI: [10.1038/s41598-023-46350-2](https://doi.org/10.1038/s41598-023-46350-2).

GPL-3.

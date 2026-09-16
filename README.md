# CCIGSVA

CCIGSVA calculates sample-level pathway enrichment from ligand-receptor (LR)
communication scores, separately for each source-target cell-type pair.

**Development version 0.1.0.** Based on the CCI-GSVA approach in
[Jeong et al. (2023)](https://doi.org/10.1038/s41598-023-46350-2) and an
author-supplied analysis script. Exact reproduction of the published results
has **not** been established: the original GMT, reference inputs/results and
historical software versions are still needed. See [implementation decisions](METHODS.md).

## Installation

Requires R >= 4.3 and GSVA >= 1.50 (parameter-object API).

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
The example is synthetic and contains no patient data.

## Use existing CellChat exports and GMT

```r
samples <- c("sample1", "sample2", "sample3") # complete cohort, including empty exports
tables <- setNames(lapply(samples, function(s) {
  read.csv(file.path("Res_CellChat", paste0("cellchat_", s, ".csv")))
}), samples)
cci <- from_cellchat(tables)
lr_sets <- read_lr_gmt("typeB_interactionSet.gmt")

# Cell counts must come from the cell metadata, not significant CCI rows.
cell_counts <- table(obj$Anno1, obj$Individual)
coverage <- filter_cell_types(cell_counts[, samples, drop = FALSE],
                              min_cells = 10, min_coverage = 0.8)
result <- cci_gsva(cci, lr_sets, samples = samples,
                   cell_types = coverage$cell_type[coverage$keep], min_size = 2)
```

Set `min_coverage = 0.7` for 70%, or `1` for presence in every sample.
The historical script used `rowSums(cell_counts >= 10) == 16`; this is only
equivalent to 100% coverage if the cohort contains exactly 16 samples.
The script passed `min.sz = 1`; use `min_size = 1` for that setting.
The new package defaults to 2 at the author's request. GMT construction may
have imposed a separate minimum, which is not yet known.

Alternatively, construct sets directly from a CellChat database table:

```r
lr_sets <- prepare_lr_sets(CellChatDB.human$interaction)
```

This alternative is not assumed equivalent to the historical GMT.
CellChat is not required to install CCIGSVA. Run CCI inference upstream;
`from_cellchat()` accepts exported tables, not CellChat S4 objects.
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
- Keep source, target and pathway labels separate; never parse them using dots.

Zero filling is an input convention. A missing cell type, an untested edge and
a nonsignificant edge are not necessarily biologically equivalent. Choose the
coverage threshold upstream and retain the complete sample list. GSVA scores
depend on the input cohort; rerunning on a different cohort can change scores.

```r
S4Vectors::metadata(result)$pair_report
S4Vectors::metadata(result)$effective_lr_sets
S4Vectors::metadata(result)$versions
```

COVID-19 severity feature selection, clustering, IGP and subtype labeling are
application analyses and are outside this package. No historical regression
test is claimed until genuine reference inputs and outputs are available.

## Reference and license

Jeong K, Kim Y, Jeon J, Kim K. Subtyping of COVID-19 samples based on
cell-cell interaction in single cell transcriptomes. *Scientific Reports*
13, 19629 (2023). DOI: [10.1038/s41598-023-46350-2](https://doi.org/10.1038/s41598-023-46350-2).

GPL-3. This repository contains newly organized package code and synthetic
examples; it does not redistribute CellChatDB, the historical GMT or cohort data.

# Source audit and implementation decisions

Audit date: 2026-09-16. Source: an R script supplied directly by the author.
The supplied script uses Tonsil-named inputs; it is evidence for implementation
details, not by itself proof of the exact COVID-19 analysis environment.

## Confirmed in the supplied script

- CellChat runs per `Individual`, with `filterCommunication(min.cells = 10)`.
- Input rows are exported by `subsetCommunication(cellchat)`, so zero filling
  operates on exported interactions, not necessarily every raw CellChat probability.
- Cell counts come from `table(obj$Anno1, obj$Individual)`.
- Cell-type eligibility is `rowSums(ct_ind >= 10) == 16` (exactly 16 samples).
- Directed source-target combinations include self interactions.
- Interaction identifiers are `interaction_name`; scores are `prob`.
- `dcast(..., value.var = "prob", fill = 0)` forms the observed union and fills gaps.
- Per-pair matrices have LR interactions in rows and samples in columns.
- There is no additional score transformation before GSVA.
- The sets are read from `typeB_interactionSet.gmt`, whose contents and creation
  code were not supplied.
- The call explicitly uses `method = "gsva"`, `min.sz = 1`, `max.sz = Inf`,
  `abs.ranking = FALSE`, `parallel.sz = 1`, `verbose = FALSE`.
- `ssgsea.norm = TRUE` belongs to the ssGSEA method, not the selected GSVA method;
  it is not carried into the new GSVA parameter object.
- Gaussian kernel, tau 1 and mx.diff TRUE are consistent with historical defaults;
  the installed historical GSVA version is not available. The paper explicitly
  reports a Gaussian kernel.

## Deliberate package choices

- Default `min_size = 2` implements the author's stated desired threshold.
  `min_size = 1` implements the script's argument; neither establishes the GMT's
  original construction threshold. Size means unique matched, nonconstant LRs.
- Coverage is a separate count-based helper, default 80%, configurable to 70%
  or 100%. The scorer does not silently infer cell presence from CCI edges.
- An explicit sample list retains empty samples that `dcast` could otherwise drop.
  This can change scores relative to the old script if those samples were absent.
- Exact equality on source and target replaces substring matching with `grep`.
  Opaque feature IDs avoid greedy regex parsing and collisions in dotted labels.
- Duplicate interaction keys fail instead of allowing reshape aggregation/counts.
- Explicit missing/non-finite scores fail; absent rows alone are zero-filled.
- Constant features are removed before size filtering, following dense GSVA.
- Sets containing the entire variable LR universe are skipped and reported.
  This is an explicit guard, not a claim of historical behavior.
- Empty pairs are skipped. A run with no eligible pathways fails with a message.
- GSVA uses serial execution and dense matrices; package and GSVA versions,
  parameters, effective sets and per-pair counts are retained in metadata.

## Script-only problems not copied

- `df.net.target_cts` is referenced before assignment in an unused check.
- The first sample is processed once before, then again inside, the loop.
- Fixed-substring matching can select unintended source-target features.
- Dot concatenation and greedy `gsub` can corrupt labels containing dots.
- `gsva_res_dir` is assigned but an undefined `share_dir` is created/used.
- The last `list_res_gsva` aggregation discards original set membership provenance.

## Needed before claiming historical reproduction

1. The actual COVID-19 analysis script and its sample/cell-type selection.
2. The original GMT and its construction rule (including any LR-count filter).
3. A small shareable input and saved output from the original implementation.
4. Historical GSVA, CellChat and database versions plus relevant session information.
5. Numerical comparison with a justified tolerance across the intended GSVA versions.

Synthetic unit/integration tests check package behavior and agreement with a
direct call to the installed GSVA. They are not a 2023 paper regression test.

Sources: [paper](https://doi.org/10.1038/s41598-023-46350-2),
[GSVA manual](https://bioconductor.org/packages/release/bioc/manuals/GSVA/man/GSVA.pdf),
[GSVA vignette](https://bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.html).

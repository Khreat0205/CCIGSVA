test_that("scores agree with a direct modern GSVA call", {
    ex <- example_ccigsva()
    result <- cci_gsva(ex$cci, ex$lr_sets)
    input <- build_lr_matrices(ex$cci)
    p <- GSVA::gsvaParam(input$matrices[[1]], ex$lr_sets, minSize = 2,
                         kcdf = "Gaussian", tau = 1, maxDiff = TRUE, absRanking = FALSE)
    reference <- GSVA::gsva(p, verbose = FALSE, BPPARAM = BiocParallel::SerialParam())
    rd <- as.data.frame(SummarizedExperiment::rowData(result))
    observed <- SummarizedExperiment::assay(result)[rd$source == input$pairs$source[1], , drop = FALSE]
    expect_equal(unname(observed), unname(reference), tolerance = 1e-12)
    expect_s4_class(result, "SummarizedExperiment")
    expect_equal(dim(result), c(4L, 6L))
    expect_equal(rd$n_lr_pairs, rep(3L, 4))
    expect_equal(S4Vectors::metadata(result)$parameters$min_size, 2)
})

test_that("constant rows are removed before pathway size filtering", {
    ex <- example_ccigsva()
    ex$cci$score[ex$cci$interaction == "LR1"] <- 0
    result <- cci_gsva(ex$cci, ex$lr_sets)
    rd <- as.data.frame(SummarizedExperiment::rowData(result))
    expect_equal(rd$n_lr_pairs[rd$pathway == "Pathway_A"], c(2L, 2L))
    expect_error(cci_gsva(ex$cci, list(single = "LR2")), "No pathways")
    # One pair emits one expected GSVA singleton warning.
    single <- ex$cci[ex$cci$source == "B", ]
    expect_warning(one <- cci_gsva(single, list(single = "LR2"), min_size = 1),
                   "Some gene sets have size one")
    expect_equal(nrow(one), 1L)
    expect_error(cci_gsva(ex$cci, list(all = paste0("LR", 1:8))), "No pathways")
})

test_that("sample order and label boundaries are preserved", {
    ex <- example_ccigsva()
    ex$cci$source[ex$cci$source == "B"] <- "B.with.dot"
    result <- cci_gsva(ex$cci, ex$lr_sets, samples = paste0("S", 6:1))
    expect_equal(colnames(result), paste0("S", 6:1))
    expect_true("B.with.dot" %in% SummarizedExperiment::rowData(result)$source)
    expect_error(cci_gsva(ex$cci, ex$lr_sets, samples = "S1"), "every observed")
    expect_error(cci_gsva(ex$cci, ex$lr_sets, min_size = 0), "positive integer")
})

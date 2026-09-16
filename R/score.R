#' Build interaction by sample matrices for each cell-type pair
#' @param cci Standardized table from [prepare_cci()] or [from_cellchat()].
#' @param samples Complete sample identifiers in the desired order. Defaults to
#'   the samples attribute, then the observed samples. Supply explicitly to retain
#'   samples without communication rows.
#' @param cell_types Optional cell types retained for both source and target.
#' @return A list containing matrices, a pairs data.frame, and samples. Matrix
#'   names are opaque pair IDs; source and target labels remain separate columns.
#' @details The LR universe is the union across samples within each observed
#'   pair. Absent rows are zero-filled, including a sample missing the entire
#'   pair. This is a numerical convention, not evidence of biological absence.
#'   Completely unobserved pairs produce no matrix. No normalization is applied.
#' @export
build_lr_matrices <- function(cci, samples = NULL, cell_types = NULL) {
    if (is.null(samples)) samples <- attr(cci, "samples")
    cci <- prepare_cci(cci, interaction = "interaction", score = "score")
    if (is.null(samples)) samples <- unique(cci$sample)
    samples <- .ids(samples, "samples", TRUE)
    if (!all(cci$sample %in% samples))
        stop("samples must include every observed sample.", call. = FALSE)
    if (!is.null(cell_types)) {
        cell_types <- .ids(cell_types, "cell_types", TRUE)
        cci <- cci[cci$source %in% cell_types & cci$target %in% cell_types, , drop = FALSE]
    }
    if (!nrow(cci)) stop("No communication rows remain after filtering.", call. = FALSE)
    pairs <- unique(cci[c("source", "target")])
    pairs <- pairs[order(pairs$source, pairs$target), , drop = FALSE]
    rownames(pairs) <- NULL
    pairs$pair_id <- sprintf("pair_%04d", seq_len(nrow(pairs)))
    matrices <- lapply(seq_len(nrow(pairs)), function(i) {
        d <- cci[cci$source == pairs$source[i] & cci$target == pairs$target[i], , drop = FALSE]
        lr <- sort(unique(d$interaction))
        m <- matrix(0, length(lr), length(samples), dimnames = list(lr, samples))
        m[cbind(match(d$interaction, lr), match(d$sample, samples))] <- d$score
        m
    })
    names(matrices) <- pairs$pair_id
    list(matrices = matrices, pairs = pairs, samples = samples)
}

#' Calculate pathway enrichment of communication profiles
#' @param cci Standardized long communication table.
#' @param lr_sets Named list of pathway interaction identifiers.
#' @param samples,cell_types See [build_lr_matrices()].
#' @param min_size Minimum number of matched nonconstant LR interactions per
#'   pathway. Defaults to 2 for the new package; the supplied historical script
#'   used min.sz = 1. Set this explicitly to 1 when comparing to that script.
#' @param kcdf GSVA kernel: Gaussian (default), Poisson, none or auto.
#' @param tau GSVA random-walk weight; defaults to 1.
#' @param max_diff,abs_ranking GSVA options; defaults TRUE and FALSE.
#' @return A SummarizedExperiment with an assay named ccigsva, rowData containing
#'   source, target, pathway and n_lr_pairs, and colData containing sample.
#'   Metadata records parameters, package versions, effective LR sets and a
#'   per-pair processing report, including pairs with no scored pathways.
#' @details Uses the GSVA parameter-object API and serial execution. Constant
#'   LR rows are removed before pathway size filtering, matching dense GSVA
#'   preprocessing. Pathways covering the entire remaining LR universe are
#'   skipped because the enrichment statistic has no out-of-set background.
#'   Scores depend on the cohort supplied, and are not a frozen prediction model.
#' @export
cci_gsva <- function(cci, lr_sets, samples = NULL, cell_types = NULL,
                     min_size = 2L, kcdf = "Gaussian", tau = 1,
                     max_diff = TRUE, abs_ranking = FALSE) {
    if (!is.list(lr_sets) || is.data.frame(lr_sets) || !length(lr_sets))
        stop("lr_sets must be a nonempty named list.", call. = FALSE)
    .ids(names(lr_sets), "Pathway names", TRUE)
    lr_sets <- lapply(lr_sets, function(x) unique(.ids(x, "LR set")))
    if (length(min_size) != 1L || !is.numeric(min_size) || !is.finite(min_size) ||
        min_size < 1 || min_size != floor(min_size))
        stop("min_size must be a positive integer.", call. = FALSE)
    if (length(tau) != 1L || !is.numeric(tau) || !is.finite(tau) || tau <= 0)
        stop("tau must be a positive finite number.", call. = FALSE)
    for (v in list(max_diff, abs_ranking))
        if (!is.logical(v) || length(v) != 1L || is.na(v))
            stop("max_diff and abs_ranking must be single logical values.", call. = FALSE)
    kcdf <- match.arg(kcdf, c("Gaussian", "Poisson", "none", "auto"))
    input <- build_lr_matrices(cci, samples, cell_types)
    if (length(input$samples) < 2L) stop("GSVA requires at least two samples.", call. = FALSE)
    results <- rows <- reports <- used_sets <- list()
    for (i in seq_along(input$matrices)) {
        m <- input$matrices[[i]]
        pair <- input$pairs[i, ]
        n_input <- nrow(m)
        varying <- apply(m, 1L, function(x) any(x != x[1L]))
        m <- m[varying, , drop = FALSE]
        sets <- lapply(lr_sets, intersect, y = rownames(m))
        size_ok <- lengths(sets) >= min_size
        full <- size_ok & lengths(sets) == nrow(m)
        sets <- sets[size_ok & !full]
        report <- data.frame(pair, n_lr_input = n_input, n_lr_variable = nrow(m),
                             n_pathways = length(sets), n_full_universe = sum(full))
        reports[[i]] <- report
        used_sets[[pair$pair_id]] <- sets
        if (nrow(m) < 2L || !length(sets)) next
        param <- GSVA::gsvaParam(exprData = m, geneSets = sets, minSize = min_size,
                                 maxSize = Inf, kcdf = kcdf, tau = tau,
                                 maxDiff = max_diff, absRanking = abs_ranking)
        score <- GSVA::gsva(param, verbose = FALSE, BPPARAM = BiocParallel::SerialParam())
        if (any(!is.finite(score)))
            stop("GSVA returned non-finite scores for ", pair$pair_id, call. = FALSE)
        score <- score[, input$samples, drop = FALSE]
        rows[[length(rows) + 1L]] <- data.frame(source = pair$source, target = pair$target,
            pathway = rownames(score), n_lr_pairs = unname(lengths(sets[rownames(score)])))
        results[[length(results) + 1L]] <- score
    }
    if (!length(results))
        stop("No pathways remain: check varying LR scores, set overlap, min_size and background.", call. = FALSE)
    score <- do.call(rbind, results)
    row_data <- do.call(rbind, rows)
    rownames(score) <- rownames(row_data) <- sprintf("feature_%06d", seq_len(nrow(score)))
    report <- do.call(rbind, reports)
    if (any(report$n_pathways == 0L) || any(report$n_full_universe > 0L))
        warning("Some pairs/pathways were skipped; see metadata(result)$pair_report.", call. = FALSE)
    SummarizedExperiment::SummarizedExperiment(
        assays = list(ccigsva = score),
        rowData = S4Vectors::DataFrame(row_data),
        colData = S4Vectors::DataFrame(sample = input$samples, row.names = input$samples),
        metadata = list(parameters = list(min_size = min_size, kcdf = kcdf, tau = tau,
                         max_diff = max_diff, abs_ranking = abs_ranking,
                         missing_rows = "zero", normalization = "none"),
                         pair_report = report, effective_lr_sets = used_sets,
                         versions = c(CCIGSVA = as.character(utils::packageVersion("CCIGSVA")),
                                      GSVA = as.character(utils::packageVersion("GSVA")))))
}

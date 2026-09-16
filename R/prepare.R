.ids <- function(x, label, unique_required = FALSE) {
    if (!(is.character(x) || is.factor(x)))
        stop(label, " must contain character identifiers.", call. = FALSE)
    x <- as.character(x)
    if (!length(x) || anyNA(x) || any(!nzchar(trimws(x))))
        stop(label, " must be nonempty and cannot contain missing/blank identifiers.", call. = FALSE)
    if (unique_required && anyDuplicated(x))
        stop(label, " must be unique.", call. = FALSE)
    x
}

.columns <- function(x, cols) {
    if (!is.data.frame(x) || anyDuplicated(names(x)) ||
        !all(cols %in% names(x)))
        stop("Input must be a data.frame with columns: ", paste(cols, collapse = ", "), call. = FALSE)
    as.data.frame(x)
}

#' Standardize a long ligand-receptor communication table
#'
#' @param data A data.frame with one row per sample, source, target and interaction.
#' @param sample,source,target,interaction,score Input column names.
#' @return A data.frame with columns sample, source, target, interaction and score.
#' @details Duplicate keys and non-finite scores are rejected. Absent rows are
#'   filled with zero by [build_lr_matrices()], but explicit NA scores are errors.
#' @export
prepare_cci <- function(data, sample = "sample", source = "source",
                        target = "target", interaction = "interaction_name",
                        score = "prob") {
    cols <- c(sample, source, target, interaction, score)
    if (length(cols) != 5L || anyDuplicated(cols))
        stop("Specify five distinct column names.", call. = FALSE)
    data <- .columns(data, cols)
    out <- data[cols]
    names(out) <- c("sample", "source", "target", "interaction", "score")
    for (nm in names(out)[1:4]) out[[nm]] <- .ids(out[[nm]], nm)
    if (!is.numeric(out$score) || any(!is.finite(out$score)))
        stop("Scores must be finite numeric values.", call. = FALSE)
    if (anyDuplicated(out[1:4]))
        stop("Duplicate sample/source/target/interaction keys; resolve upstream.", call. = FALSE)
    rownames(out) <- NULL
    out
}

#' Construct pathway interaction sets
#' @param data A mapping data.frame, for example CellChatDB.human$interaction.
#' @param interaction,pathway Mapping column names.
#' @return A named list of unique interaction identifiers per pathway.
#' @details No minimum size is applied here. Size filtering is performed on
#'   the matched, nonconstant LR universe within each source-target pair.
#' @export
prepare_lr_sets <- function(data, interaction = "interaction_name",
                            pathway = "pathway_name") {
    data <- .columns(data, c(interaction, pathway))
    ids <- .ids(data[[interaction]], "interaction")
    pathways <- .ids(data[[pathway]], "pathway")
    lapply(split(ids, pathways), unique)
}

#' Read ligand-receptor pathway sets in GMT format
#' @param file Path to a tab-delimited GMT file.
#' @return A named list of interaction identifiers. The second GMT field is ignored.
#' @export
read_lr_gmt <- function(file) {
    lines <- readLines(file, warn = FALSE)
    lines <- lines[nzchar(trimws(lines))]
    fields <- strsplit(lines, "\t", fixed = TRUE)
    if (!length(fields) || any(lengths(fields) < 3L))
        stop("Each nonempty GMT line needs a name, description and interaction.", call. = FALSE)
    sets <- lapply(fields, function(x) unique(.ids(x[-c(1, 2)], "GMT interactions")))
    names(sets) <- .ids(vapply(fields, `[`, character(1), 1L), "GMT pathway names", TRUE)
    sets
}

#' Filter cell types by sample coverage
#' @param cell_counts Named numeric matrix with cell types in rows and all samples
#'   in columns; entries are nonnegative integer cell counts.
#' @param min_cells Minimum cells for a cell type to be present in a sample.
#' @param min_coverage Required fraction of all samples, between zero and one.
#' @return A data.frame with cell_type, n_samples, coverage and keep columns.
#' @details Uses counts, not the presence of significant communication edges.
#'   To require presence in every sample, set min_coverage = 1.
#' @export
filter_cell_types <- function(cell_counts, min_cells = 10L, min_coverage = 0.8) {
    if (!is.matrix(cell_counts) || !is.numeric(cell_counts) ||
        !all(dim(cell_counts) > 0L) || any(!is.finite(cell_counts)) ||
        any(cell_counts < 0 | cell_counts != floor(cell_counts)))
        stop("cell_counts must be a nonnegative integer-valued numeric matrix.", call. = FALSE)
    ct <- .ids(rownames(cell_counts), "Cell-type row names", TRUE)
    .ids(colnames(cell_counts), "Sample column names", TRUE)
    if (length(min_cells) != 1L || !is.numeric(min_cells) || !is.finite(min_cells) ||
        min_cells < 1 || min_cells != floor(min_cells))
        stop("min_cells must be a positive integer.", call. = FALSE)
    if (length(min_coverage) != 1L || !is.numeric(min_coverage) ||
        !is.finite(min_coverage) || min_coverage <= 0 || min_coverage > 1)
        stop("min_coverage must be in (0, 1].", call. = FALSE)
    n <- rowSums(cell_counts >= min_cells)
    data.frame(cell_type = ct, n_samples = unname(n), coverage = unname(n / ncol(cell_counts)),
               keep = unname(n / ncol(cell_counts) >= min_coverage))
}

#' Combine exported CellChat communication tables
#' @param tables Named list of data.frames returned by CellChat::subsetCommunication().
#' @return A standardized communication table. The samples attribute retains
#'   names of samples with empty tables.
#' @details CellChat itself is not a dependency. Empty tables are permitted,
#'   provided at least one sample has communication rows.
#' @export
from_cellchat <- function(tables) {
    if (!is.list(tables) || is.data.frame(tables))
        stop("tables must be a named list of data.frames.", call. = FALSE)
    samples <- .ids(names(tables), "Table sample names", TRUE)
    rows <- lapply(seq_along(tables), function(i) {
        d <- .columns(tables[[i]], c("source", "target", "interaction_name", "prob"))
        if (!nrow(d)) return(NULL)
        d$sample <- samples[i]
        prepare_cci(d)
    })
    out <- do.call(rbind, rows)
    if (is.null(out)) stop("No communication rows in any sample.", call. = FALSE)
    attr(out, "samples") <- samples
    out
}

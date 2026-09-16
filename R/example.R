#' Small synthetic communication example
#' @return A list with cci, lr_sets and cell_counts. No patient data are included.
#' @examples
#' ex <- example_ccigsva()
#' result <- cci_gsva(ex$cci, ex$lr_sets)
#' result
#' @export
example_ccigsva <- function() {
    d <- expand.grid(sample = paste0("S", seq_len(6)),
                     source = c("B", "CD4T"), target = "Mono",
                     interaction_name = paste0("LR", seq_len(8)),
                     stringsAsFactors = FALSE)
    i <- seq_len(nrow(d))
    d$prob <- ((i * 17L + (i %% 5L) * 13L) %% 101L) / 1000
    d <- d[i %% 7L != 0L, ]
    counts <- matrix(20, 3L, 6L,
                     dimnames = list(c("B", "CD4T", "Mono"), paste0("S", seq_len(6))))
    list(cci = prepare_cci(d), lr_sets = list(Pathway_A = c("LR1", "LR2", "LR3"),
         Pathway_B = c("LR4", "LR5", "LR6")), cell_counts = counts)
}

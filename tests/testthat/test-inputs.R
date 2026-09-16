test_that("union, zeros and empty samples have explicit semantics", {
    d <- data.frame(sample = c("s1", "s2"), source = "A.B", target = "C",
                    interaction_name = c("l1", "l2"), prob = c(0.2, 0.3))
    x <- build_lr_matrices(prepare_cci(d), samples = c("s2", "s1", "empty"))
    expect_equal(x$matrices[[1]], matrix(c(0, .3, .2, 0, 0, 0), 2,
        dimnames = list(c("l1", "l2"), c("s2", "s1", "empty"))))
    expect_equal(x$pairs$source, "A.B")
    d2 <- d; d2$source <- "A"; d2$target <- "B.C"
    expect_length(build_lr_matrices(prepare_cci(rbind(d, d2)))$matrices, 2)
    expect_error(prepare_cci(rbind(d, d)), "Duplicate")
    d$prob[1] <- NA_real_
    expect_error(prepare_cci(d), "finite")
})

test_that("coverage counts cells over the complete sample denominator", {
    counts <- matrix(c(10, 9, 10, 10, 10, 10, 0, 10), 2,
                      dimnames = list(c("A", "B"), paste0("s", 1:4)))
    x <- filter_cell_types(counts, min_coverage = .7)
    expect_equal(x$coverage, c(.75, .75))
    expect_true(all(x$keep))
    expect_false(any(filter_cell_types(counts, min_coverage = .8)$keep))
    expect_error(filter_cell_types(counts, min_coverage = 1.1), "min_coverage")
})

test_that("adapters retain empty samples and GMT membership", {
    d <- data.frame(source = "A", target = "B", interaction_name = "l1", prob = .1)
    cci <- from_cellchat(list(s1 = d, s2 = d[FALSE, ]))
    expect_equal(colnames(build_lr_matrices(cci)$matrices[[1]]), c("s1", "s2"))
    f <- tempfile(); on.exit(unlink(f))
    writeLines(c("p1\tdesc\tl1\tl2\tl1", "p2\tdesc\tl3"), f)
    expect_equal(read_lr_gmt(f), list(p1 = c("l1", "l2"), p2 = "l3"))
    mapping <- data.frame(interaction_name = c("l1", "l2", "l1"), pathway_name = "p1")
    expect_equal(prepare_lr_sets(mapping), list(p1 = c("l1", "l2")))
})

test_that("QC detects incompatible pigeonhole/crossing pairs", {
  data(cluster_data_ind2, package = "CloneTree")
  qc <- clone_qc(cluster_data_ind2, subject = "IGC-11-1130")
  expect_s3_class(qc, "clonetree_qc")
  expect_equal(qc_solutions(qc)$excluded, "C11")
  bad_pairs <- qc$pairwise[qc$pairwise$incompatible, c("cluster_a", "cluster_b")]
  expect_true(any(bad_pairs$cluster_a == "C5" & bad_pairs$cluster_b == "C11"))
  expect_true(any(bad_pairs$cluster_a == "C8" & bad_pairs$cluster_b == "C11"))
})

test_that("QC keeps multiple best solutions when exclusions tie", {
  df <- data.frame(
    subject = rep("TumorA", 6),
    sampleID = rep(c("S1", "S2"), each = 3),
    cluster.no = rep(c("1", "2", "3"), times = 2),
    Clone = rep(c("Y", "N", "N"), times = 2),
    no.of.mutations = rep(c(100, 10, 10), times = 2),
    CCF = c(1, 0.7, 0.4, 1, 0.2, 0.9)
  )
  qc <- clone_qc(df, subject = "TumorA", objective = "max_clusters")
  expect_equal(length(qc$solutions), 2)
  excluded <- sort(vapply(qc$solutions, function(x) x$excluded, character(1)))
  expect_equal(excluded, c("C2", "C3"))
})

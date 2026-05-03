test_that("tree inference nests compatible clusters", {
  data(cluster_data_ind2, package = "CloneTree")
  qc <- clone_qc(cluster_data_ind2, subject = "IGC-11-1130")
  tree <- infer_clone_tree(qc)
  expect_s3_class(tree, "clonetree")
  expect_equal(unname(tree$parent["C10"]), "C8")
  expect_equal(unname(tree$parent["C7"]), "C5")
  expect_true("C11" %in% tree$excluded)
})

test_that("plotting functions return ggplot objects", {
  data(cluster_data_ind2, package = "CloneTree")
  qc <- clone_qc(cluster_data_ind2, subject = "IGC-11-1130")
  tree <- infer_clone_tree(qc)
  expect_s3_class(plot_clone_qc(qc), "ggplot")
  expect_s3_class(plot_clone_ccf(qc), "ggplot")
  expect_s3_class(plot_clone_tree(tree), "ggplot")
  expect_s3_class(plot_qc_rules_guide(), "ggplot")
})

test_that("layout separates sibling branches from the same parent", {
  data(cluster_data_ind2, package = "CloneTree")
  qc <- clone_qc(cluster_data_ind2, subject = "IGC-11-1130")
  tree <- infer_clone_tree(qc)
  layout <- layout_clone_tree(tree)
  edges_by_parent <- split(layout$edges, layout$edges$parent)

  for (edge_group in edges_by_parent[lengths(edges_by_parent) > 1]) {
    angles <- round(
      atan2(edge_group$xend - edge_group$x, edge_group$y - edge_group$yend),
      4
    )
    expect_equal(length(unique(angles)), nrow(edge_group))
  }
})

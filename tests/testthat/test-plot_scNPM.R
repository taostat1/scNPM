# Test 1: plot_group_network functionality
test_that("plot_group_network works with built-in top50 results", {
  data(top50_Data_result_with_protein, package = "scNPM")
  gene_names <- gene_names
  treat_edge_mat <- edge_treat_est[, , 1]
  ctrl_edge_mat <-edge_ctrl_est[, , 1]

  p_treat <- plot_group_network(
    edge_mat = treat_edge_mat,
    gene_names = gene_names,
    group_type = "treatment"
  )
  expect_true("ggplot" %in% class(p_treat))

  p_ctrl <- plot_group_network(
    edge_mat = ctrl_edge_mat,
    gene_names = gene_names,
    group_type = "control"
  )
  expect_true("ggplot" %in% class(p_treat))
})

# Test 2: plot_diff_network functionality
test_that("plot_diff_network works with built-in top50 differential data", {
  data(top50_Data_result_with_protein, package = "scNPM")

  gene_names <- gene_names
  diff_mat <- edge_diff_networks$treat_1_ctrl_2$diff_matrix

  p_combined <- plot_diff_network(
    diff_mat = diff_mat,
    gene_names = gene_names,
    output_mode = "combined",
    solid_threshold = 0.5,
    dashed_threshold = -0.5
  )
  expect_true("ggplot" %in% class(p_combined) || is.null(p_combined))

  p_separate <- plot_diff_network(
    diff_mat = diff_mat,
    gene_names = gene_names,
    output_mode = "separate",
    solid_threshold = 0.5,
    dashed_threshold = -0.5
  )
  expect_type(p_separate, "list")
  expect_true(all(c("solid", "dashed") %in% names(p_separate)))
  if (!is.null(p_separate$solid)) expect_true("ggplot" %in% class(p_separate$solid))
  if (!is.null(p_separate$dashed)) expect_true("ggplot" %in% class(p_separate$dashed))
})

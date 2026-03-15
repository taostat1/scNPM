# ==============================================================================
# scNPM Package Demo Script
# Purpose: Demonstrate the complete scNPM analysis pipeline using built-in data
# Author: scNPM Development Team
# Date: 2026-03
# ==============================================================================

# --------------------------
# 1. Load scNPM package
# --------------------------
# Use devtools::load_all() for development mode; use library(scNPM) if installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::load_all()  # Load package in development mode
# library(scNPM)  # Load installed package (choose one)

# --------------------------
# 2. Load built-in example data
# --------------------------
# Load gene expression matrix (from package's data directory)
data("top50_goal_gene_mt", package = "scNPM")
goal_gene_mt <- as.matrix(goal_gene_mt)

# Create cell metadata matrix (cell_t)
# First column: 0 = control, 1 = treatment
# Create cell metadata dataframe (match cell order with expression matrix)
cell_t <- matrix(
  c(
    ifelse(cell$PR_or_MR == "Yes", 1, 0),
    scale(cell$scater_qc.all.total_features_by_counts)[,1],
    scale(cell$scater_qc.all.total_counts)[,1],
    scale(cell$scater_qc.endogenous.pct_counts)[,1],
    scale(cell$scater_qc.feature_control_MT.pct_counts)[,1],
    scale(cell$scater_qc.all.pct_counts_in_top_50_features)[,1]
  ),
  ncol = 6, dimnames = list(NULL, c("treatment","complete","expression","reliability","apoptosis","homogeneity"))
)

# Extract protein interaction matrix (PPI)
tsv_path <- system.file(
  "data",
  "top50_string_interactions_short.tsv",
  package = "scNPM"
)
tsv_lines <- readLines(tsv_path)
tsv_lines_processed <- gsub("^\\s*#", "", tsv_lines)
tsv_lines_processed <- tsv_lines_processed[tsv_lines_processed != ""]
string_interactions_short <- read.table(
  text = paste(tsv_lines_processed, collapse = "\n"),
  sep = "\t", header = TRUE, stringsAsFactors = FALSE,
  check.names = FALSE
)
# --------------------------
# 3. Run core scNPM analysis
# --------------------------
result_list <- run_scNPM(
  goal_gene_mt = goal_gene_mt,
  cell_t = cell_t,
  K = 4,                   # Number of cell clusters
  ds_seed = 20250102,      # Random seed for reproducibility
  num_iter = 10000,         # Total MCMC iterations
  num_save = 5000,          # Post-burn-in iterations to retain
  protein_input = string_interactions_short,  # Protein interaction data
  edge_prob_threshold = 0.6,      # Edge probability threshold
  edge_prob_adj = 0.75            # Edge probability adjustment factor
)
cat("scNPM analysis completed! Runtime:", round(result_list$runtime, 2), "minutes\n")
save_path <- "out/top50_Data_result_with_protein.RData"
# Reduce Memory
result_list <- result_list[names(result_list) != "Result"]
save(result_list, file = save_path)
# --------------------------
# 4. Visualize key results
# --------------------------
cat("\n📊 Generating visualizations...\n")
load("out/Data_result_with_protein.RData")
# Create output directory
if (!dir.exists("out")) {
  dir.create("out", recursive = TRUE)
}

# 4.1 Plot treatment group network (cluster 1)
p_treat <- plot_group_network(
  edge_mat = result_list$edge_treat_est[, , 1],
  gene_names = result_list$gene_names,
  group_type = "treatment"
)
ggplot2::ggsave("out/treatment_cluster_1_network.pdf", plot = p_treat, width = 10, height = 8, dpi = 300)

# 4.2 Plot control group network (cluster 1)
p_ctrl <- plot_group_network(
  edge_mat = result_list$edge_ctrl_est[, , 1],
  gene_names = result_list$gene_names,
  group_type = "control"
)
ggsave("output/control_cluster_1_network.png", plot = p_ctrl, width = 10, height = 8, dpi = 300)

# 4.3 Plot differential network (treatment vs control)
diff_mat <- result_list$edge_diff_networks[[1]]$diff_matrix
p_diff <- plot_diff_network(
  diff_mat = diff_mat,
  gene_names = result_list$gene_names,
  output_mode = "separate"
)
ggplot2::ggsave("out/differential_network_cluster_1.pdf", plot = p_diff$solid, width = 12, height = 10, dpi = 300)

K <- result_list$K
gene_names <- result_list$gene_names

for (cluster_id in 1:K) {
  if (!as.character(cluster_id) %in% names(result_list$match_map)) next
  # Treatment network
  p_treat <- plot_group_network(
    edge_mat = result_list$edge_treat_est[, , cluster_id],
    gene_names = gene_names,
    group_type = "treatment"
  )
  ggplot2::ggsave(paste0("out/treatment_cluster_", cluster_id, "_network.pdf"),
         plot = p_treat, width = 10, height = 6.5)

  # Control network
  ctrl_id <- result_list$match_map[as.character(cluster_id)]
  p_ctrl <- plot_group_network(
    edge_mat = result_list$edge_ctrl_est[, , ctrl_id],
    gene_names = gene_names,
    group_type = "control"
  )
  ggplot2::ggsave(paste0("out/control_cluster_", ctrl_id, "_network.pdf"),
         plot = p_ctrl, width = 10, height = 6.5)
  # Differential network
  diff_mat <- result_list$edge_diff_networks[[cluster_id]]$diff_matrix
  p_diff <- plot_diff_network(
    diff_mat = diff_mat,
    gene_names = gene_names,
    output_mode = "separate"
  )
  plot_solid <- p_diff$solid
  plot_dashed <- p_diff$dashed
  ggplot2::ggsave(paste0("out/differential_network_solid_", cluster_id, ".pdf"),
                  plot = plot_solid, width = 10, height = 6.5)
  ggplot2::ggsave(paste0("out/differential_network_dashed_", cluster_id, ".pdf"),
                  plot = plot_dashed, width = 10, height = 6.5)
}


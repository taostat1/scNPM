#' Run Single-Cell Network Propensity Matching (scNPM)
#'
#' This function encapsulates the complete scNPM analysis pipeline, including data preprocessing,
#' MCMC inference, and result extraction for single-cell RNA-seq data. Protein interaction matrix
#' is optional (auto-initialized to zero matrix if not provided). The core input is gene expression
#' matrix and cell metadata matrix (cell_t), where only the first column of cell_t (0/1 for control/treatment)
#' is necessary.
#'
#' @param goal_gene_mt Gene expression matrix (required, rows = genes, columns = cells)
#' @param cell_t Cell metadata matrix/data frame (required, first column must be binary 0/1: 0=control, 1=treatment;
#'        other columns are optional covariates with no format restrictions)
#' @param K Number of cell types/clusters for MCMC inference (required, e.g., 4)
#' @param ds_seed Random seed for reproducibility (default = 20250102)
#' @param num_iter Total MCMC iterations (default = 1000)
#' @param num_save Number of post-burn-in iterations to retain (default = 500)
#' @param protein_input Protein-protein interaction (PPI) input (default = NULL, auto-initialize zero matrix).
#'        Supports 3 types of input:
#'        1. Local TSV file path (e.g., "ppi.tsv");
#'        2. TSV text content (character vector, one line per row);
#'        3. Pre-loaded data frame of PPI (columns: node1, node2, combined_score).
#'        Automatically removes lines starting with '#' (comment lines) or leading '#' symbols.
#' @param ssp_v0 Spike-slab prior parameter v0 (default = 0.02)
#' @param ssp_v1 Spike-slab prior parameter v1 (default = 1)
#' @param ssp_l Spike-slab prior parameter l (default = 1)
#' @param edge_prob_threshold Edge probability threshold for network construction (default = 0.6)
#' @param edge_prob_adj Edge probability adjustment factor (default = 0.75)
#' @param epsilon_theta Step size for Hamiltonian Monte Carlo (HMC) sampling of theta (default = 0.2)
#' @param num_step_theta Number of steps for HMC sampling of theta (default = 20)
#' @param epsilon_lam Step size for HMC sampling of lambda (default = 0.01)
#' @param num_step_lam Number of steps for HMC sampling of lambda (default = 10)
#'
#' @return A named list containing all key analysis results (can be saved with `save()`):
#' \item{runtime}{MCMC runtime in minutes (difftime object)}
#' \item{K}{Number of cell clusters used in analysis}
#' \item{G}{Number of genes (rows in goal_gene_mt)}
#' \item{N}{Number of cells (columns in goal_gene_mt)}
#' \item{lam0_est}{Posterior estimates of lambda0 parameter (length = G)}
#' \item{lam1_est}{Posterior estimates of lambda1 parameter (length = G)}
#' \item{mu_treat_est}{Posterior mean estimates for treatment group (G × K matrix)}
#' \item{mu_ctrl_est}{Posterior mean estimates for control group (G × K matrix)}
#' \item{invcov_treat_est}{Posterior precision matrix estimates for treatment group (G × G × K array)}
#' \item{invcov_ctrl_est}{Posterior precision matrix estimates for control group (G × G × K array)}
#' \item{edge_treat_est}{Posterior edge matrix estimates for treatment group (G × G × K array)}
#' \item{edge_ctrl_est}{Posterior edge matrix estimates for control group (G × G × K array)}
#' \item{group_treat_est}{Cluster labels for treatment cells (length = number of treatment cells)}
#' \item{group_ctrl_est}{Cluster labels for control cells (length = number of control cells)}
#' \item{pi_treat_est}{Posterior cell type probability estimates for treatment group (length = K)}
#' \item{pi_ctrl_est}{Posterior cell type probability estimates for control group (length = K)}
#' \item{theta_treat_est}{Posterior expression matrix estimates for treatment group}
#' \item{theta_ctrl_est}{Posterior expression matrix estimates for control group}
#' \item{theta_est}{Posterior overall expression matrix estimates (cells × genes)}
#' \item{matching_result_df}{Cluster matching table (treatment_cluster vs matched_control_cluster)}
#' \item{smd_results}{Standardized Mean Difference (SMD) values for each treatment cluster}
#' \item{match_map}{Named vector mapping treatment clusters to matched control clusters}
#' \item{edge_diff_networks}{List of per-cluster differential networks (solid/dashed edges based on 0.5 threshold)}
#' \item{Result}{Raw output from C++ MCMC function (full iteration results)}
#'
#' @examples
#' \dontrun{
#' # Step 1: Prepare input data
#' load("top50_goal_gene_mt.RData")
#' goal_gene_mt <- as.matrix(goal_gene_mt)
#'
#' # Step 2: Create cell_t matrix (only first column is mandatory 0/1)
#' # Example: First column = treatment/control (0/1), other columns as optional covariates
#' cell_t <- matrix(
#'   c(
#'     ifelse(cell$PR_or_MR == "Yes", 1, 0),
#'     scale(cell$scater_qc.all.total_counts),
#'     scale(cell$scater_qc.all.total_features_by_counts)
#'   ),
#'   ncol = 3,
#'   dimnames = list(NULL, c("treatment", "total_counts", "total_features"))
#' )
#'
#' # Step 3: Load and process protein interaction TSV file
#' # Read TSV lines and remove leading # (keep line content if no #)
#' tsv_lines <- readLines("top50_string_interactions_short.tsv")
#' tsv_lines_processed <- gsub("^\\s*#", "", tsv_lines)  # Remove leading # only
#' tsv_lines_processed <- tsv_lines_processed[tsv_lines_processed != ""]  # Remove empty lines
#' # Read processed lines into data frame (preserve original column names)
#' top50_string_interactions_short <- read.table(
#'   text = paste(tsv_lines_processed, collapse = "\n"),
#'   sep = "\t", header = TRUE, stringsAsFactors = FALSE,
#'   check.names = FALSE
#' )
#'
#' # Step 4: Run scNPM function
#' result_list <- run_scNPM(
#'   goal_gene_mt = goal_gene_mt,
#'   cell_t = cell_t,
#'   K = 4,
#'   protein_input = top50_string_interactions_short
#' )
#'
#' # Step 4: Save results
#' save(result_list, file = "output/top50_Data_result_with_protein.RData")
#' }
#' @export
run_scNPM <- function(
    goal_gene_mt,
    cell_t,
    K,
    ds_seed = 20250102,
    num_iter = 1000,
    num_save = 500,
    protein_input = NULL,
    ssp_v0 = 0.02,
    ssp_v1 = 1,
    ssp_l = 1,
    edge_prob_threshold = 0.6,
    edge_prob_adj = 0.75,
    epsilon_theta = 0.2,
    num_step_theta = 20,
    epsilon_lam = 0.01,
    num_step_lam = 10
) {
  # --------------------------
  # 0. Load dependencies
  # --------------------------
  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    stop("Package 'Rcpp' is required but not installed.")
  }
  if (!requireNamespace("RcppArmadillo", quietly = TRUE)) {
    stop("Package 'RcppArmadillo' is required but not installed.")
  }
  if (!requireNamespace("RcppParallel", quietly = TRUE)) {
    stop("Package 'RcppParallel' is required but not installed.")
  }

  # --------------------------
  # 1. Validate input parameters
  # --------------------------
  # Validate gene expression matrix
  if (!is.matrix(goal_gene_mt) && !is.data.frame(goal_gene_mt)) {
    stop("'goal_gene_mt' must be a matrix or data frame (rows=genes, columns=cells).")
  }
  X <- as.matrix(goal_gene_mt)  # Convert to matrix for consistency
  G <- dim(X)[1]     # Number of genes (rows)
  N <- dim(X)[2]     # Number of cells (columns)
  if (G == 0 || N == 0) {
    stop("'goal_gene_mt' is empty (no genes or cells detected).")
  }

  # Set gene names (auto-generate if missing)
  gene_names <- rownames(X)
  if (is.null(gene_names)) {
    warning("'goal_gene_mt' has no row names - auto-generating gene names (Gene_1, Gene_2, ...)")
    gene_names <- paste0("Gene_", 1:G)
    rownames(X) <- gene_names
  }

  # Validate cell_t matrix (core requirement: first column = 0/1)
  if (!is.matrix(cell_t) && !is.data.frame(cell_t)) {
    stop("'cell_t' must be a matrix or data frame (rows=cells, columns=covariates).")
  }
  cell_t <- as.matrix(cell_t)  # Convert to matrix for consistency
  if (nrow(cell_t) != N) {
    stop(paste("'cell_t' has", nrow(cell_t), "rows, but 'goal_gene_mt' has", N, "columns (cells) - dimension mismatch."))
  }

  # Check first column of cell_t (must be binary 0/1)
  treatment_col <- cell_t[, 1]
  if (!all(treatment_col %in% c(0, 1))) {
    stop("First column of 'cell_t' must contain only 0 (control) and 1 (treatment) values.")
  }

  # Validate treatment/control groups (non-empty)
  treat_cell_idx <- which(treatment_col == 1)
  ctrl_cell_idx <- which(treatment_col == 0)
  if (length(treat_cell_idx) == 0) {
    stop("No treatment cells found (first column of 'cell_t' has no 1 values).")
  }
  if (length(ctrl_cell_idx) == 0) {
    stop("No control cells found (first column of 'cell_t' has no 0 values).")
  }

  # --------------------------
  # 2. Prepare cell metadata
  # --------------------------
  # Standardize optional covariates (columns 2+) if present
  if (ncol(cell_t) > 1) {
    # Standardize columns 2 to ncol(cell_t)
    for (i in 2:ncol(cell_t)) {
      if (is.numeric(cell_t[, i])) {
        cell_t[, i] <- scale(cell_t[, i])  # Z-score normalization
        cell_t[is.na(cell_t[, i]), i] <- 0  # Replace NA with 0 after scaling
      } else {
        warning(paste("Column", i, "of 'cell_t' is non-numeric - kept as original values"))
      }
    }
  }

  # --------------------------
  # 3. Initialize protein interaction matrix
  # --------------------------
  if (is.null(protein_input)) {
    # Default: initialize zero matrix with gene names
    cat("Without protein interaction data provided - initializing zero matrix\n")
    protein_mat <- matrix(0, nrow = G, ncol = G, dimnames = list(gene_names, gene_names))
  } else {
    if (is.data.frame(protein_input)) {
      cat("Using pre-loaded protein interaction data frame\n")
      string_df <- protein_input
    }
    else {
      stop("'protein_input' must be a file path, TSV text, or data frame")
    }
    required_cols <- c("node1", "node2", "combined_score")
    if (!all(required_cols %in% colnames(string_df))) {
      stop(paste(
        "Protein interaction data must contain columns:",
        paste(required_cols, collapse = ", ")
      ))
    }
    protein_mat <- matrix(0, nrow = G, ncol = G, dimnames = list(gene_names, gene_names))
    for (i in 1:nrow(string_df)) {
      g1 <- string_df$node1[i]
      g2 <- string_df$node2[i]
      score <- string_df$combined_score[i]

      # Map gene names to matrix indices
      idx1 <- which(gene_names == g1)
      idx2 <- which(gene_names == g2)

      # Only add interactions for genes present in expression matrix
      if (length(idx1) > 0 && length(idx2) > 0) {
        protein_mat[idx1, idx2] <- score
        protein_mat[idx2, idx1] <- score  # Symmetric matrix for undirected network
      }
    }
    # Binarize protein matrix (non-zero = 1, zero = 0)
    protein_mat[protein_mat != 0] <- 1
    cat(paste0("Protein interaction matrix loaded (dimensions: ", nrow(protein_mat), "×", ncol(protein_mat), ")\n"))
  }

  # --------------------------
  # 4. Independent K-means clustering
  # --------------------------
  set.seed(ds_seed)
  X_log <- log(X + 1)  # Log-transform expression matrix for better clustering

  # Clustering for treatment group (cells × genes matrix)
  X_treat_log <- t(X_log[, treat_cell_idx])
  kmeans_treat <- kmeans(X_treat_log, centers = K, nstart = 5)
  group_treat <- kmeans_treat$cluster

  # Clustering for control group (different seed to avoid correlation)
  set.seed(ds_seed + 1)
  X_ctrl_log <- t(X_log[, ctrl_cell_idx])
  kmeans_ctrl <- kmeans(X_ctrl_log, centers = K, nstart = 5)
  group_ctrl <- kmeans_ctrl$cluster

  # --------------------------
  # 5. Theta_t initialization
  # --------------------------
  ind_zero <- (X == 0)  # Zero value indicator matrix
  theta_t <- X          # Initialize theta_t with original expression values

  # Stratified zero imputation (by treatment/control and cluster)
  for (g in 1:G) {
    ind_zero_g <- ind_zero[g, ]
    X_g <- X[g, ]

    # Impute treatment group cells
    for (k in 1:K) {
      treat_k_idx <- intersect(treat_cell_idx, which(group_treat == k))
      if (length(treat_k_idx) == 0) next

      X_gk <- X_g[treat_k_idx]
      ind_zero_gk <- ind_zero_g[treat_k_idx]

      # Impute zeros with 5th percentile of non-zero values (if available)
      if (sum(!ind_zero_gk) > 0) {
        theta_gk <- X_gk
        theta_gk[ind_zero_gk] <- quantile(X_gk[!ind_zero_gk], na.rm = TRUE, probs = 0.05)
        theta_t[g, treat_k_idx] <- theta_gk
      }
    }

    # Impute control group cells
    for (k in 1:K) {
      ctrl_k_idx <- intersect(ctrl_cell_idx, which(group_ctrl == k))
      if (length(ctrl_k_idx) == 0) next

      X_gk <- X_g[ctrl_k_idx]
      ind_zero_gk <- ind_zero_g[ctrl_k_idx]

      # Impute zeros with 5th percentile of non-zero values (if available)
      if (sum(!ind_zero_gk) > 0) {
        theta_gk <- X_gk
        theta_gk[ind_zero_gk] <- quantile(X_gk[!ind_zero_gk], na.rm = TRUE, probs = 0.05)
        theta_t[g, ctrl_k_idx] <- theta_gk
      }
    }
  }

  # Log transformation and NA/Inf handling
  theta_t <- log(theta_t)
  theta_t[is.na(theta_t) | theta_t == -Inf] <- 0  # Replace NA/-Inf with 0
  cat("Theta_t initialization completed (stratified zero imputation)\n")

  # --------------------------
  # 6. Initialize mean/covariance/precision matrices
  # --------------------------
  # Treatment group matrices
  mu_treat <- matrix(0, nrow = G, ncol = K)
  cov_treat <- array(diag(1, G), dim = c(G, G, K))
  invcov_treat <- array(diag(1, G), dim = c(G, G, K))

  for (k in 1:K) {
    treat_k_idx <- intersect(treat_cell_idx, which(group_treat == k))
    if (length(treat_k_idx) == 0) next

    theta_k <- theta_t[, treat_k_idx]
    mu_treat[, k] <- rowMeans(theta_k)  # Cluster-specific mean

    # Gene-specific variance (diagonal covariance matrix)
    gene_var <- diag((theta_k - mu_treat[, k]) %*% t(theta_k - mu_treat[, k]) / (length(treat_k_idx) - 1))
    gene_var[gene_var == 0] <- 1  # Avoid division by zero
    diag(cov_treat[, , k]) <- gene_var
    diag(invcov_treat[, , k]) <- 1.0 / gene_var  # Precision = inverse variance
  }

  # Control group matrices
  mu_ctrl <- matrix(0, nrow = G, ncol = K)
  cov_ctrl <- array(diag(1, G), dim = c(G, G, K))
  invcov_ctrl <- array(diag(1, G), dim = c(G, G, K))

  for (k in 1:K) {
    ctrl_k_idx <- intersect(ctrl_cell_idx, which(group_ctrl == k))
    if (length(ctrl_k_idx) == 0) next

    theta_k <- theta_t[, ctrl_k_idx]
    mu_ctrl[, k] <- rowMeans(theta_k)  # Cluster-specific mean

    # Gene-specific variance (diagonal covariance matrix)
    gene_var <- diag((theta_k - mu_ctrl[, k]) %*% t(theta_k - mu_ctrl[, k]) / (length(ctrl_k_idx) - 1))
    gene_var[gene_var == 0] <- 1  # Avoid division by zero
    diag(cov_ctrl[, , k]) <- gene_var
    diag(invcov_ctrl[, , k]) <- 1.0 / gene_var  # Precision = inverse variance
  }

  # --------------------------
  # 7. Initialize edge matrices
  # --------------------------
  edge_treat <- array(0, dim = c(G, G, K))
  edge_ctrl <- array(0, dim = c(G, G, K))

  # Assign protein matrix to all clusters
  for (k in 1:K) {
    edge_treat[, , k] <- protein_mat
    edge_ctrl[, , k] <- protein_mat
  }
  cat("Edge matrices initialized (based on protein interaction matrix)\n")

  # --------------------------
  # 8. Final parameter preparation
  # --------------------------
  theta_t <- t(theta_t)       # Transpose to cells × genes (C++ input format)
  ind_zero <- t(ind_zero)     # Match theta_t dimensions
  gam <- rep(1, K)            # Dirichlet prior weights for cell type probabilities

  # Cell type probability initialization
  pi_treat <- update_pi_R(group_treat[treat_cell_idx], gam, K)
  pi_ctrl <- update_pi_R(group_ctrl[ctrl_cell_idx], gam, K)

  # Spike-slab prior parameters for lambda
  lambda0_t <- rnorm(G, mean = 3, sd = 0.2)
  lambda1_t <- rnorm(G, mean = -2, sd = 0.2)
  ssp_xi <- 2 / (G - 1)        # Spike-slab prior parameter xi

  # --------------------------
  # 9. Run C++ MCMC function
  # --------------------------
  start_time <- Sys.time()
  Result <- MCMC_with_difficulty_based_matching(
    num_iter = num_iter,
    num_save = num_save,
    theta_t = theta_t,
    ind_zero = ind_zero,
    mu_treat = mu_treat,
    mu_ctrl = mu_ctrl,
    invcov_treat = invcov_treat,
    cov_treat = cov_treat,
    invcov_ctrl = invcov_ctrl,
    cov_ctrl = cov_ctrl,
    edge_treat = edge_treat,
    edge_ctrl = edge_ctrl,
    group_treat = group_treat,
    group_ctrl = group_ctrl,
    lambda0_t = lambda0_t,
    lambda1_t = lambda1_t,
    pi_treat = pi_treat,
    pi_ctrl = pi_ctrl,
    gam = gam,
    cell_t = cell_t,
    G = G,
    N = N,
    K = K,
    ssp_v0 = ssp_v0,
    ssp_v1 = ssp_v1,
    ssp_l = ssp_l,
    ssp_xi = ssp_xi,
    epsilon_theta = epsilon_theta,
    num_step_theta = num_step_theta,
    eta_mu = 0,
    tau_sq_mu = 1,
    lam0_0 = 3,
    lam1_0 = -2,
    sigma2_lam0 = 0.2,
    sigma2_lam1 = 0.2,
    edge_prob_threshold = edge_prob_threshold,
    edge_prob_adj = edge_prob_adj,
    epsilon_lam = epsilon_lam,
    num_step_lam = num_step_lam
  )

  # Calculate MCMC runtime
  end_time <- Sys.time()
  runtime <- difftime(end_time, start_time, units = "mins")
  cat(paste0("\nMCMC Runtime: ", round(runtime, 2), " minutes\n"))

  # --------------------------
  # 10. Extract posterior estimates
  # --------------------------
  # Core posterior estimates
  lam0_est <- rowMeans(Result$lam0_iter)
  lam1_est <- rowMeans(Result$lam1_iter)
  mu_treat_est <- Result$mu_treat_est
  mu_ctrl_est <- Result$mu_ctrl_est
  invcov_treat_est <- Result$invcov_treat_est
  invcov_ctrl_est <- Result$invcov_ctrl_est
  edge_treat_est <- Result$edge_treat_est
  edge_ctrl_est <- Result$edge_ctrl_est
  group_treat_est <- Result$group_treat_est
  group_ctrl_est <- Result$group_ctrl_est
  group_est <- Result$group_est
  pi_treat_est <- rowMeans(Result$pi_treat_iter)
  pi_ctrl_est <- rowMeans(Result$pi_ctrl_iter)
  theta_treat_est <- Result$theta_treat_est
  theta_ctrl_est <- Result$theta_ctrl_est
  theta_est <- Result$theta_est

  # Cluster matching results
  matching_result_df <- Result$matching_results$matching_result_df
  smd_results <- Result$matching_results$smd_results

  # Create cluster matching map (treatment → control)
  match_map <- setNames(
    matching_result_df$matched_control_cluster,
    matching_result_df$treatment_cluster
  )
  match_map <- match_map[!is.na(match_map)]  # Remove unmatched clusters

  # --------------------------
  # 11. Per-cluster differential network analysis
  # --------------------------
  edge_diff_networks <- list()

  # Iterate over each matched treatment-control cluster pair
  for (treat_clust in names(match_map)) {
    treat_clust <- as.numeric(treat_clust)
    ctrl_clust <- as.numeric(match_map[as.character(treat_clust)])

    # Extract cluster-specific edge matrices
    edge_treat_clust <- edge_treat_est[, , treat_clust]
    edge_ctrl_clust <- edge_ctrl_est[, , ctrl_clust]

    # Calculate differential network (treatment - control)
    edge_diff <- edge_treat_clust - edge_ctrl_clust

    # Initialize solid/dashed edge matrices and edge lists
    solid_mat <- matrix(0, nrow = G, ncol = G)  # Solid edges (difference > 0.5)
    dashed_mat <- matrix(0, nrow = G, ncol = G) # Dashed edges (difference < -0.5)
    solid_edges <- data.frame(from = character(), to = character(), value = numeric(), stringsAsFactors = FALSE)
    dashed_edges <- data.frame(from = character(), to = character(), value = numeric(), stringsAsFactors = FALSE)

    # Identify solid/dashed edges (upper triangle only for undirected network)
    for (i in 1:nrow(edge_diff)) {
      for (j in i:ncol(edge_diff)) {
        val <- edge_diff[i, j]

        # Solid edge (positive difference > 0.5)
        if (val > 0.5) {
          solid_mat[i, j] <- 1
          solid_mat[j, i] <- 1  # Symmetric for undirected network
          solid_edges <- rbind(solid_edges, data.frame(
            from = gene_names[i],
            to = gene_names[j],
            value = val,
            stringsAsFactors = FALSE
          ))
        }
        # Dashed edge (negative difference < -0.5)
        else if (val < -0.5) {
          dashed_mat[i, j] <- 1
          dashed_mat[j, i] <- 1  # Symmetric for undirected network
          dashed_edges <- rbind(dashed_edges, data.frame(
            from = gene_names[i],
            to = gene_names[j],
            value = val,
            stringsAsFactors = FALSE
          ))
        }
      }
    }

    # Store per-cluster differential network results
    edge_diff_networks[[paste0("treat_", treat_clust, "_ctrl_", ctrl_clust)]] <- list(
      pair_info = data.frame(
        treatment_cluster = treat_clust,
        matched_control_cluster = ctrl_clust,
        smd_value = smd_results[treat_clust]
      ),
      diff_matrix = edge_diff,
      solid_network = list(
        mat = solid_mat,
        edges = solid_edges,
        count = nrow(solid_edges)
      ),
      dashed_network = list(
        mat = dashed_mat,
        edges = dashed_edges,
        count = nrow(dashed_edges)
      )
    )
  }

  # --------------------------
  # 12. Compile final results list
  # --------------------------
  result_list <- list(
    # Basic metadata
    runtime = runtime,
    K = K,
    G = G,
    N = N,
    # Core posterior parameter estimates
    lam0_est = lam0_est,
    lam1_est = lam1_est,
    gene_names = gene_names,
    mu_treat_est = mu_treat_est,
    mu_ctrl_est = mu_ctrl_est,
    invcov_treat_est = invcov_treat_est,
    invcov_ctrl_est = invcov_ctrl_est,
    edge_treat_est = edge_treat_est,
    edge_ctrl_est = edge_ctrl_est,
    group_treat_est = group_treat_est,
    group_ctrl_est = group_ctrl_est,
    group_est = group_est,
    pi_treat_est = pi_treat_est,
    pi_ctrl_est = pi_ctrl_est,
    theta_treat_est = theta_treat_est,
    theta_ctrl_est = theta_ctrl_est,
    theta_est = theta_est,
    # Cluster matching results
    matching_result_df = matching_result_df,
    smd_results = smd_results,
    match_map = match_map,
    # Per-cluster differential networks
    edge_diff_networks = edge_diff_networks,
    # Raw MCMC output (full iteration data)
    Result = Result
  )

  # Return complete results list
  return(result_list)
}

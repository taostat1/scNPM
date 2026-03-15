#include "scNPM_types.h"
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <Rmath.h>
#include <trng/yarn2.hpp>
#include <trng/uniform_dist.hpp>
#include <trng/normal_dist.hpp>
#include <trng/lcg64.hpp>
#include <trng/discrete_dist.hpp>
#include <trng/gamma_dist.hpp>
#include <omp.h>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <cmath>

using namespace arma;
using namespace Rcpp;
using namespace std;
using namespace trng;

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(rTRNG)]]
// TRNG >= 4.22 requires C++11
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends("RcppArmadillo")]]

vec mvnrnd(const vec& mu, const mat& sigma, int n) {
  int d = mu.n_elem;
  mat Z = randn(d, n);
  mat L = chol(sigma, "lower");
  mat X = repmat(mu, 1, n) + L * Z;
  return X.col(0);
}

double U_theta(const vec &theta_miss, const vec &mu, const mat &prec, const vec &lambda0_miss, const vec &lambda1_miss) {
  vec U = (theta_miss - mu).t() * prec * (theta_miss - mu) / 2.0 -
    sum(log(normcdf(lambda0_miss + lambda1_miss % theta_miss)));
  return U(0);
}

double calculate_edge_probability(double omega, double v0, double v1, double xi) {
  double term1 = omega * omega * (1.0/(2.0*v1*v1) - 1.0/(2.0*v0*v0));
  double term2 = log(v1) - log(v0) + log(1.0 - xi);
  double log_p = term1 + term2;
  double r = 1.0 / (1.0 + exp(log_p));
  return r;
}

void U_theta_grad(vec &U_grad, const vec &theta_miss, const vec &mu, const mat &prec, const vec &lambda0_miss, const vec &lambda1_miss) {
  U_grad = prec * (theta_miss - mu) -
    (lambda1_miss % normpdf(lambda0_miss + lambda1_miss % theta_miss)) / normcdf(lambda0_miss + lambda1_miss % theta_miss);
}

void HMC_theta(vec &theta_new, double &log_r, const vec &theta_miss, const vec &mu, const mat &prec, const vec &lambda0_miss,
               const vec &lambda1_miss, const vec &p_rnorm, const double &epsilon, const int &num_step) {
  //The kinetic energy has the simplest form sum(p^2)/2
  vec U_grad;

  theta_new = theta_miss;
  vec p_new = p_rnorm;

  // a half step for momentum at the beginning
  U_theta_grad(U_grad, theta_new, mu, prec, lambda0_miss, lambda1_miss);
  p_new -= epsilon *  U_grad / 2.0;

  // full steps for position and momentum
  for (int i = 0; i < (num_step-1); i++) {
    theta_new += epsilon * p_new;

    U_theta_grad(U_grad, theta_new, mu, prec, lambda0_miss, lambda1_miss);
    p_new -= epsilon * U_grad;
  }
  theta_new += epsilon * p_new;

  U_theta_grad(U_grad, theta_new, mu, prec, lambda0_miss, lambda1_miss);
  p_new -= epsilon * U_grad / 2.0;

  p_new = -p_new;
  log_r = U_theta(theta_miss, mu, prec, lambda0_miss, lambda1_miss) - U_theta(theta_new, mu, prec, lambda0_miss, lambda1_miss) +
    sum(p_rnorm % p_rnorm) / 2.0 - sum(p_new % p_new) / 2.0;
}

void update_theta(mat &theta_t, const mat &ind_zero, const mat &mu_t, const cube &invcov_t, const cube &cov_t, const vec &lambda0_t,
                  const vec &lambda1_t, const vec &group_t, const int &N, const unsigned long &seed, const double &epsilon = 0.01,
                  const int &num_step = 50) {
  int i;

#pragma omp parallel shared(theta_t) private(i) num_threads(24)
{
  trng::yarn2 rx;
  rx.seed(seed);
  int size = omp_get_num_threads();
  int rank = omp_get_thread_num();
  rx.split(size, rank);
  trng::uniform_dist<> r_unif(0,1);
  trng::normal_dist<> r_norm(0,1);
#pragma omp for schedule(auto)
  for (i = 0; i < N; i++) {

    vec ind_zero_i = ind_zero.row(i).t();
    uvec ind_0 = find(ind_zero_i == true);

    if (ind_0.n_elem > 0) {
      vec theta_i = theta_t.row(i).t();
      uvec ind_obs = find(ind_zero_i == false);

      vec theta_i_obs = theta_i(ind_obs);

      vec theta_i_0 = theta_i(ind_0);

      vec mu_i = mu_t.col(group_t(i));
      mat cov_i = cov_t.slice(group_t(i));
      mat invcov_i = invcov_t.slice(group_t(i));

      //calculate the acceptance probability
      vec mu_obs = mu_i(ind_obs);
      vec mu_0 = mu_i(ind_0);

      mat invcov_i_21 = invcov_i.submat(ind_obs, ind_0);

      mat prec_cond = invcov_i.submat(ind_0, ind_0);

      mat cov_obs_inv = invcov_i.submat(ind_obs, ind_obs) - invcov_i_21 * inv(prec_cond) * invcov_i_21.t();

      mat cov_21 = cov_i.submat(ind_0, ind_obs);
      mat cov_0 = cov_i.submat(ind_0, ind_0);

      vec mu_cond = mu_0 + cov_21 * cov_obs_inv * (theta_i_obs - mu_obs);

      vec p_rnorm(ind_0.n_elem);
      for (uword t = 0; t < ind_0.n_elem; t++) {
        p_rnorm(t) = r_norm(rx);
      }
      double tmp_unif = r_unif(rx);

      vec lambda0_miss = lambda0_t(ind_0);
      vec lambda1_miss = lambda1_t(ind_0);

      vec theta_star_0;
      double log_r;
      HMC_theta(theta_star_0, log_r, theta_i_0, mu_cond, prec_cond, lambda0_miss, lambda1_miss, p_rnorm, epsilon, num_step);

      if (tmp_unif < exp(log_r)) {
        theta_i(ind_0) = theta_star_0;
        theta_t.row(i) = theta_i.t();
      }
    }
  }
}
}

void update_mu(mat &mu_t, const mat &theta_t, const cube &invcov_t, const vec &group_t, const int &G, const int &K, const double &eta_mu = 0, const double &tau_sq_mu = 1) {
  mat I_tau_sq = eye(G,G);
  I_tau_sq /= tau_sq_mu;
  for (int k = 0; k < K; k++) {
    uvec ind_k = find(group_t == k);
    vec tmp1 = sum(theta_t.rows(ind_k), 0).t();

    mat invcov_k = invcov_t.slice(k);

    mat COV = inv(ind_k.n_elem * invcov_k + I_tau_sq);
    vec MU = COV * (invcov_k * tmp1 + eta_mu/tau_sq_mu);

    vec mu_k = mvnrnd(MU, COV, 1);
    mu_t.col(k) = mu_k;
  }
}

void update_invcov(cube &invcov_t, cube &cov_t, const cube &edge_t, const mat &theta_t, const mat &mu_t, const vec &group_t,
                   const double &ssp_v0, const double &ssp_v1, const double &ssp_l, const int &G, const int &K, const unsigned long &seed) {
  int k;
#pragma omp parallel shared(invcov_t, cov_t) private(k) num_threads(K)
{
  trng::yarn2 rx;
  rx.seed(seed);
  int size = omp_get_num_threads();
  int rank = omp_get_thread_num();
  rx.split(size, rank);
  trng::normal_dist<> r_norm(0,1);
#pragma omp for schedule(auto)
  for (k = 0; k < K; k++) {
    uvec ind_k = find(group_t == k);
    mat theta_k = theta_t.rows(ind_k).t();

    mat theta_mu = theta_k.each_col() - mu_t.col(k);
    mat S = theta_mu * theta_mu.t();

    mat edge_k = edge_t.slice(k);
    mat V = edge_k * ssp_v1 * ssp_v1;
    uvec ind_n_v1 = find(edge_k == 0);
    V(ind_n_v1).fill(ssp_v0 * ssp_v0);
    V.diag().fill(0);

    mat invcov_k = invcov_t.slice(k);
    mat cov_k = cov_t.slice(k);

    vec G_vec = regspace(0, G-1);

    for (int g = 0; g < G; g++) {
      uvec ind_ng = find(G_vec != g);
      vec v12_g = V.col(g);
      mat v12_inv = diagmat(1/v12_g(ind_ng));
      vec cov_k_g = cov_k.col(g);
      mat w11_inv = cov_k.submat(ind_ng, ind_ng) - cov_k_g(ind_ng) * cov_k_g(ind_ng).t() / cov_k_g(g);
      mat C = (S(g,g) + ssp_l) * w11_inv + v12_inv;

      mat C_chol_inv = inv(chol(C));

      vec s12_g = S.col(g);
      vec mu_w12 = -C_chol_inv * C_chol_inv.t() * s12_g(ind_ng);

      vec rnorm_vec(G-1);
      for (int t = 0; t < (G-1); t++) {
        rnorm_vec(t) = r_norm(rx);
      }
      vec w12 = mu_w12 + C_chol_inv * rnorm_vec;

      trng::gamma_dist<> r_gamma(ind_k.n_elem/2.0 + 1.0, 2.0/(S(g,g)+ssp_l));

      double w_v = r_gamma(rx);
      vec w11_inv_w12 = w11_inv * w12;
      vec w22 = w_v + w12.t() * w11_inv_w12;

      vec w_update(G);
      w_update(ind_ng) = w12;
      w_update(g) = w22(0);

      invcov_k.col(g) = w_update;
      invcov_k.row(g) = w_update.t();

      cov_k.submat(ind_ng, ind_ng) = w11_inv + w11_inv_w12 * w11_inv_w12.t() / w_v;
      vec cov_k_ng = - w11_inv_w12 / w_v;

      vec cov_update(G);
      cov_update(ind_ng) = cov_k_ng;
      cov_update(g) = 1 / w_v;

      cov_k.col(g) = cov_update;
      cov_k.row(g) = cov_update.t();
    }
    cov_t.slice(k) = cov_k;
    invcov_t.slice(k) = invcov_k;
  }
}
}

vec calculate_propensity_score(const mat& X, const vec& treatment_status) {
  Function glm("glm");
  List glm_model = glm(
    _["formula"] = Formula("y ~ ."),
    _["data"] = DataFrame::create(_["y"] = treatment_status, _["X"] = X),
    _["family"] = "binomial"
  );
  return as<vec>(glm_model["fitted.values"]);
}

vec calculate_quadratic_form(
  const mat& theta, const cube& invcov, const vec& group,
  const uvec& idx, const int& K, const int& G
) {
  vec qf(idx.n_elem, fill::zeros);
  #pragma omp parallel for schedule(dynamic, 32) num_threads(24)
  for (uword i = 0; i < idx.n_elem; i++) {
    const uword k = static_cast<uword>(clamp(static_cast<int>(group(i)), 0, K-1));
    if (k >= static_cast<uword>(K) || k >= invcov.n_slices) continue;
    const uword cell_global_idx = idx(i);
    if (cell_global_idx >= theta.n_rows) continue;
    const rowvec theta_row = theta.row(cell_global_idx);
    qf(i) = dot(theta_row, theta_row * invcov.slice(k));
    if (isnan(qf(i)) || isinf(qf(i))) qf(i) = 0.0;
  }
  return qf;
}

List network_propensity_matching_by_difficulty(
  const mat& theta_t, const mat& cell_t,
  const vec& group_treat, const vec& group_ctrl,
  const uvec& treat_idx, const uvec& ctrl_idx,
  const cube& invcov_treat, const cube& invcov_ctrl,
  const int& G, const int& K) {

  uword N = cell_t.n_rows;
  vec treatment_status = cell_t.col(0);
  uvec treatment_indices = treat_idx;
  uvec control_indices = ctrl_idx;

  vec smd_results(K, fill::zeros);
  vec qf_treat = calculate_quadratic_form(theta_t, invcov_treat, group_treat, treat_idx, K, G);
  vec qf_ctrl = calculate_quadratic_form(theta_t, invcov_ctrl, group_ctrl, ctrl_idx, K, G);
  vec qf_global(N, fill::zeros);
  qf_global.elem(treat_idx) = qf_treat;
  qf_global.elem(ctrl_idx) = qf_ctrl;
  mat X = join_horiz(ones<vec>(N), cell_t.cols(1, cell_t.n_cols-1), qf_global);
  vec base_propensity = calculate_propensity_score(X, treatment_status);
  if (base_propensity.n_elem != N) {
    Rcpp::stop("Propensity score length mismatch with cell count!");
  }
  vec prop_treat = base_propensity.elem(treatment_indices);
  vec prop_ctrl = base_propensity.elem(control_indices);

  vec cluster_mean_treat(K, fill::zeros);
  vec difficulty_treat(K, fill::zeros);
  for (uword k = 0; k < (uword)K; k++) {
    uvec ind_k = find(group_treat == k);
    if (ind_k.is_empty() || ind_k.max() >= prop_treat.n_elem) {
      cluster_mean_treat(k) = mean(prop_treat);
      difficulty_treat(k) = 0;
      continue;
    }
    cluster_mean_treat(k) = mean(prop_treat.elem(ind_k));
    difficulty_treat(k) = abs(cluster_mean_treat(k) - mean(prop_treat));
  }

  vec cluster_mean_ctrl(K, fill::zeros);
  vec cluster_var_treat(K, fill::zeros);
  vec cluster_var_ctrl(K, fill::zeros);
  for (uword c = 0; c < (uword)K; c++) {
    uvec ind_k = find(group_treat == c);
    if (!ind_k.is_empty() && ind_k.max() < prop_treat.n_elem) {
      cluster_var_treat(c) = var(prop_treat.elem(ind_k));
    } else {
      cluster_var_treat(c) = var(prop_treat);
    }
    uvec ind_c = find(group_ctrl == c);
    if (ind_c.is_empty() || ind_c.max() >= prop_ctrl.n_elem) {
      cluster_mean_ctrl(c) = mean(prop_ctrl);
      cluster_var_ctrl(c) = var(prop_ctrl);
      continue;
    }
    cluster_mean_ctrl(c) = mean(prop_ctrl.elem(ind_c));
    cluster_var_ctrl(c) = var(prop_ctrl.elem(ind_c));
  }

  uvec sorted_clusters = sort_index(difficulty_treat);
  sorted_clusters = reverse(sorted_clusters);
  sorted_clusters = sorted_clusters.elem(find(difficulty_treat.elem(sorted_clusters) > 0));
  if (sorted_clusters.is_empty()) sorted_clusters = regspace<uvec>(0, K-1);

  uvec ctrl_matched_flag(K, fill::zeros);
  uvec treat2ctrl_map(K);
  treat2ctrl_map.fill(K);

  for (uword s = 0; s < sorted_clusters.n_elem; s++) {
    int k = sorted_clusters(s);
    if (k < 0 || k >= K) continue;

    double min_smd = 1e9;
    int best_c = -1;
    for (uword c = 0; c < (uword)K; c++) {
      if (ctrl_matched_flag(c) == 1) continue;

      double mean_diff = abs(cluster_mean_treat(k) - cluster_mean_ctrl(c));
      double var_sum = cluster_var_treat(k) + cluster_var_ctrl(c);
      double smd = 0.0;
      if (var_sum < 1e-10) {
        smd = (mean_diff < 1e-10) ? 0.0 : 1e9;
      } else {
        smd = mean_diff / sqrt(var_sum / 2.0);
      }

      // Select the smallest SMD control group cluster
      if (smd < min_smd) {
        min_smd = smd;
        best_c = c;
      }
    }

    if (best_c != -1 && best_c < K) {
      treat2ctrl_map(k) = best_c;
      ctrl_matched_flag(best_c) = 1;
      smd_results(k) = min_smd;
    }
  }

  int unmatched_treat_clusters = sum(treat2ctrl_map == K);
  int matched_pairs = K - unmatched_treat_clusters;

  vec treat_cluster(K);
  vec ctrl_cluster(K);
  for (uword k = 0; k < (uword)K; k++) {
    treat_cluster(k) = k + 1;
    if (treat2ctrl_map(k) < (uword)K) {
      ctrl_cluster(k) = treat2ctrl_map(k) + 1;
    } else {
      ctrl_cluster(k) = NA_REAL;
    }
  }
  DataFrame matching_result_df = DataFrame::create(
    Named("treatment_cluster") = treat_cluster,
    Named("matched_control_cluster") = ctrl_cluster
  );

  return List::create(
    _["smd_results"] = smd_results,
    _["quadratic_form_global"] = qf_global,
    _["treat2ctrl_cluster_map"] = treat2ctrl_map,
    _["matched_pairs"] = matched_pairs,
    _["ctrl_matched_flag"] = ctrl_matched_flag,
    _["matching_result_df"] = matching_result_df
  );
}

void update_edge_with_matching_by_difficulty(
  cube &edge_treat, cube &edge_ctrl,
  const cube &invcov_treat, const cube &invcov_ctrl,
  double ssp_v0, double ssp_v1, double ssp_xi, double edge_prob_threshold, double edge_prob_adj,
  const int &G, const int &K, const uword &N,
  const unsigned long &seed,
  const mat &theta_t, const mat &cell_t,
  const vec &group_treat, const vec &group_ctrl,
  const uvec &treat_idx, const uvec &ctrl_idx) {
  double edge_prob_threshold_up = std::max(edge_prob_threshold, 1.0 - edge_prob_threshold);
  double edge_prob_threshold_down = std::min(edge_prob_threshold, 1.0 - edge_prob_threshold);
  List matching_result = network_propensity_matching_by_difficulty(
    theta_t, cell_t, group_treat, group_ctrl, treat_idx, ctrl_idx,
    invcov_treat, invcov_ctrl,
    G, K);

  uvec treat2ctrl_map = as<uvec>(matching_result["treat2ctrl_cluster_map"]);
  uvec ctrl2treat_map(K);
  ctrl2treat_map.fill(K);
  for (uword k = 0; k < (uword)K; k++) {
    if (treat2ctrl_map(k) < (uword)K) {
      ctrl2treat_map(treat2ctrl_map(k)) = k;
    }
  }

  const uword G_u = static_cast<uword>(G);
  const uword K_u = static_cast<uword>(K);
  const double ssp_xi_adj1 = ssp_xi * edge_prob_adj;
  const double ssp_xi_adj4 = ssp_xi / edge_prob_adj;

  omp_set_num_threads(24);
  #pragma omp parallel sections shared(edge_treat, edge_ctrl) num_threads(24)
  {
    #pragma omp section
    {
      omp_set_num_threads(12);
      trng::yarn2 rx;
      trng::uniform_dist<> r_unif(0, 1);

      uword g, k, i;
      #pragma omp parallel for schedule(static,64) private(g, k, i) \
        shared(edge_treat, invcov_treat, edge_ctrl, treat2ctrl_map, seed, G_u, K_u, ssp_v0, ssp_v1, ssp_xi, ssp_xi_adj1, ssp_xi_adj4) \
        num_threads(12)
      for (g = 0; g < G_u; g++) {
        rx.seed(seed + 1000 + omp_get_thread_num() * G_u + g);
        for (k = 0; k < K_u; k++) {
          if (k >= edge_treat.n_slices || k >= invcov_treat.n_slices) continue;
          const mat& invcov_k = invcov_treat.slice(k);

          int matched_ctrl_cluster = treat2ctrl_map(k);
          bool is_matched = (matched_ctrl_cluster >= 0 && matched_ctrl_cluster < K);

          for (i = g+1; i < G_u; i++) {
            double omega_treat = 0.0;
            if (g < invcov_k.n_rows && i < invcov_k.n_cols) {
              omega_treat = invcov_k(g, i);
            }

            const double original_prob = calculate_edge_probability(omega_treat, ssp_v0, ssp_v1, ssp_xi);
            const double tmp = r_unif(rx);

            if (is_matched) {
              const uword c = static_cast<uword>(matched_ctrl_cluster);
              if (original_prob > edge_prob_threshold_up || original_prob < edge_prob_threshold_down) {
                edge_treat(g, i, k) = (tmp < original_prob) ? 1 : 0;
              } else if (c < edge_ctrl.n_slices) {
                double ctrl_edge_val = 0.0;
                if (g < edge_ctrl.n_rows && i < edge_ctrl.n_cols && c < edge_ctrl.n_slices) {
                  ctrl_edge_val = edge_ctrl(g, i, c);
                }
                const double adjusted_xi = (ctrl_edge_val == 1) ? ssp_xi_adj1 : ssp_xi_adj4;
                const double adjusted_prob = calculate_edge_probability(omega_treat, ssp_v0, ssp_v1, adjusted_xi);
                edge_treat(g, i, k) = (tmp < adjusted_prob) ? 1 : 0;
              } else {
                edge_treat(g, i, k) = (tmp < original_prob) ? 1 : 0;
              }
            } else {
              edge_treat(g, i, k) = (tmp < original_prob) ? 1 : 0;
            }
            edge_treat(i, g, k) = edge_treat(g, i, k);
          }
        }
      }
    }

    #pragma omp section
    {
      omp_set_num_threads(12);
      trng::yarn2 rx;
      trng::uniform_dist<> r_unif(0, 1);

      uword g, k, i;
      #pragma omp parallel for schedule(static,64) private(g, k, i) \
        shared(edge_ctrl, invcov_ctrl, edge_treat, ctrl2treat_map, seed, G_u, K_u, ssp_v0, ssp_v1, ssp_xi, ssp_xi_adj1, ssp_xi_adj4) \
        num_threads(12)
      for (g = 0; g < G_u; g++) {
        rx.seed(seed + 2000 + omp_get_thread_num() * G_u + g);
        for (k = 0; k < K_u; k++) {
          if (k >= edge_ctrl.n_slices || k >= invcov_ctrl.n_slices) continue;
          const mat& invcov_k = invcov_ctrl.slice(k);

          int matched_treat_cluster = ctrl2treat_map(k);
          bool is_matched = (matched_treat_cluster >= 0 && matched_treat_cluster < K);

          for (i = g+1; i < G_u; i++) {
            double omega_ctrl = 0.0;
            if (g < invcov_k.n_rows && i < invcov_k.n_cols) {
              omega_ctrl = invcov_k(g, i);
            }
            const double original_prob = calculate_edge_probability(omega_ctrl, ssp_v0, ssp_v1, ssp_xi);
            const double tmp = r_unif(rx);

            if (is_matched) {
              const uword t = static_cast<uword>(matched_treat_cluster);
              if (original_prob > edge_prob_threshold_up || original_prob < edge_prob_threshold_down) {
                edge_ctrl(g, i, k) = (tmp < original_prob) ? 1 : 0;
              } else if (t < edge_treat.n_slices) {
                double treat_edge_val = 0.0;
                if (g < edge_treat.n_rows && i < edge_treat.n_cols && t < edge_treat.n_slices) {
                  treat_edge_val = edge_treat(g, i, t);
                }
                const double adjusted_xi = (treat_edge_val == 1) ? ssp_xi_adj1 : ssp_xi_adj4;
                const double adjusted_prob = calculate_edge_probability(omega_ctrl, ssp_v0, ssp_v1, adjusted_xi);
                edge_ctrl(g, i, k) = (tmp < adjusted_prob) ? 1 : 0;
              } else {
                edge_ctrl(g, i, k) = (tmp < original_prob) ? 1 : 0;
              }
            } else {
              edge_ctrl(g, i, k) = (tmp < original_prob) ? 1 : 0;
            }
            edge_ctrl(i, g, k) = edge_ctrl(g, i, k);
          }
        }
      }
    }
  }
}

double U_lam(const vec &lambda, const vec &theta_miss, const vec &theta_obs,
             const double &lam0_0, const double &lam1_0, const double &sigma2_lam0, const double &sigma2_lam1) {
  vec U = - sum(log(1.0 - normcdf(lambda(0) + lambda(1) * theta_obs))) - sum(log(normcdf(lambda(0) + lambda(1) * theta_miss))) +
    (lambda(0) - lam0_0) * (lambda(0) - lam0_0) / (2.0 * sigma2_lam0) + (lambda(1) - lam1_0) * (lambda(1) - lam1_0) / (2.0 * sigma2_lam1);
  return U(0);
}

void U_lam_grad(vec &U_grad, const vec &lambda, const vec &theta_miss, const vec &theta_obs,
                const double &lam0_0, const double &lam1_0, const double &sigma2_lam0, const double &sigma2_lam1) {
  vec part1 = normpdf(lambda(0) + lambda(1) * theta_obs) / (1.0 - normcdf(lambda(0) + lambda(1) * theta_obs));
  vec part2 = - normpdf(lambda(0) + lambda(1) * theta_miss) / normcdf(lambda(0) + lambda(1) * theta_miss);

  U_grad(0) = sum(part1) + sum(part2) + (lambda(0) - lam0_0) / sigma2_lam0;
  U_grad(1) = sum(part1 % theta_obs) + sum(part2 % theta_miss) + (lambda(1) - lam1_0) / sigma2_lam1;
}

void HMC_lam(vec &lambda_new, double &log_r, const vec &lambda, const vec &theta_miss, const vec &theta_obs, const vec &p_rnorm,
             const double &lam0_0, const double &lam1_0, const double &sigma2_lam0, const double &sigma2_lam1,
             const double &epsilon, const int &num_step) {
  //The kinetic energy has the simplest form sum(p^2)/2
  vec U_grad(2);

  lambda_new = lambda;
  vec p_new = p_rnorm;

  // a half step for momentum at the beginning
  U_lam_grad(U_grad, lambda_new, theta_miss, theta_obs, lam0_0, lam1_0, sigma2_lam0, sigma2_lam1);
  p_new -= epsilon *  U_grad / 2.0;
  // full steps for position and momentum
  for (int i = 0; i < (num_step-1); i++) {
    lambda_new += epsilon * p_new;

    U_lam_grad(U_grad, lambda_new, theta_miss, theta_obs, lam0_0, lam1_0, sigma2_lam0, sigma2_lam1);
    p_new -= epsilon * U_grad;
  }
  lambda_new += epsilon * p_new;

  U_lam_grad(U_grad, lambda_new, theta_miss, theta_obs, lam0_0, lam1_0, sigma2_lam0, sigma2_lam1);
  p_new -= epsilon * U_grad / 2.0;

  p_new = -p_new;
  log_r = U_lam(lambda, theta_miss, theta_obs, lam0_0, lam1_0, sigma2_lam0, sigma2_lam1) -
    U_lam(lambda_new, theta_miss, theta_obs, lam0_0, lam1_0, sigma2_lam0, sigma2_lam1) +
    sum(p_rnorm % p_rnorm) / 2.0 - sum(p_new % p_new) / 2.0;
}

void update_lambda(vec &lambda0_t, vec &lambda1_t, const mat &theta_t, const mat &ind_zero, const int &G, const unsigned long &seed,
                   const double &lam0_0=2, const double &lam1_0=-2, const double &sigma2_lam0=0.25, const double &sigma2_lam1=0.25,
                   const double &epsilon = 0.01, const int &num_step = 50) {
  int g;

#pragma omp parallel shared(lambda0_t, lambda1_t) private(g) num_threads(24)
{
  trng::yarn2 rx;
  rx.seed(seed);
  int size = omp_get_num_threads();
  int rank = omp_get_thread_num();
  rx.split(size, rank);
  trng::uniform_dist<> r_unif(0,1);
  trng::normal_dist<> r_norm(0,1);
#pragma omp for schedule(auto)
  for (g = 0; g < G; g++) {
    vec lambda_g(2);
    lambda_g(0) = lambda0_t(g);
    lambda_g(1) = lambda1_t(g);

    vec theta_g = theta_t.col(g);
    uvec ind_0= find(ind_zero.col(g) == true);
    uvec ind_obs = find(ind_zero.col(g) == false);

    vec theta_g_obs = theta_g(ind_obs);

    vec theta_g_miss = theta_g(ind_0);

    vec p_rnorm(2);
    p_rnorm(0) = r_norm(rx);
    p_rnorm(1) = r_norm(rx);

    double tmp_unif = r_unif(rx);

    vec lambda_g_new;
    double log_r;
    HMC_lam(lambda_g_new, log_r, lambda_g, theta_g_miss, theta_g_obs, p_rnorm, lam0_0, lam1_0, sigma2_lam0, sigma2_lam1, epsilon, num_step);

    if (tmp_unif < exp(log_r) && lambda_g_new(1) < 0) {
      lambda0_t(g) = lambda_g_new(0);
      lambda1_t(g) = lambda_g_new(1);
    }
  }
}
}

vec rDirichlet(const vec &alpha_vec) {
  vec tmp(alpha_vec.n_elem, fill::zeros);
  for (unsigned int i = 0; i < alpha_vec.n_elem; i++) {
    tmp(i) = log(randg(distr_param(alpha_vec(i), 1.0)));
  }
  tmp = tmp - max(tmp);
  tmp = exp(tmp);
  tmp = tmp/sum(tmp);
  return tmp;
}

void update_pi(vec &pi_t, const vec &group_t, const vec &gam, const int &K) {
  vec tmp0(K, fill::zeros);
  for (int k = 0; k < K; k++) {
    uvec ind_j = find(group_t == k);
    tmp0(k) = ind_j.n_elem;
  }
  pi_t = rDirichlet(tmp0+gam);
}

void update_group(vec &group_t, const cube &invcov_t, const mat &theta_t, const mat &mu_t, const vec &pi_t, const int &G, const int &N, const int &K, const unsigned long &seed) {
  uvec pi_n0 = find(pi_t > 0);
  int i;
  vec invcov_logdet(K);
  for (int k = 0; k < K; k++) {
    invcov_logdet(k) = log(det(invcov_t.slice(k)))/2.0;
  }
#pragma omp parallel shared(group_t, pi_n0, invcov_logdet) private(i) num_threads(24)
{
  trng::lcg64 rx;
  rx.seed(seed);
  int size = omp_get_num_threads();
  int rank = omp_get_thread_num();
  rx.split(size, rank);

#pragma omp for schedule(auto)
  for (i = 0; i < N; i++) {
    vec theta_i = theta_t.row(i).t();
    vec tmp(K);
    tmp.fill(- datum::inf);
    for (unsigned int k_pi = 0; k_pi < pi_n0.n_elem; k_pi++) {
      vec tmp_k_pi = invcov_logdet(pi_n0(k_pi)) -
        (theta_i - mu_t.col(pi_n0(k_pi))).t() * invcov_t.slice(pi_n0(k_pi)) * (theta_i - mu_t.col(pi_n0(k_pi)))/2.0;
      tmp(pi_n0(k_pi)) = tmp_k_pi(0);
    }
    tmp.replace(datum::inf, pow(10, 308));
    vec tmp_new = tmp - max(tmp);
    tmp_new = exp(tmp_new);
    vec prob = tmp_new % pi_t;
    prob = prob / sum(prob);

    trng::discrete_dist dist_C(prob.begin(), prob.end());
    group_t(i) = dist_C(rx);

  }
}
}

// [[Rcpp::export]]
List MCMC_with_difficulty_based_matching(
  const int num_iter, const int num_save,
  mat theta_t, mat ind_zero,
  mat mu_treat, mat mu_ctrl,
  cube invcov_treat, cube cov_treat,
  cube invcov_ctrl, cube cov_ctrl,
  cube edge_treat, cube edge_ctrl,
  vec group_treat, vec group_ctrl,
  vec lambda0_t, vec lambda1_t,
  vec pi_treat, vec pi_ctrl,
  vec gam, mat cell_t,
  const int G, const int N_int, const int K,
  double ssp_v0, double ssp_v1, double ssp_l, double ssp_xi,
  double epsilon_theta = 0.2, int num_step_theta = 20,
  double eta_mu = 0.0, double tau_sq_mu = 1.0,
  double lam0_0 = 2.0, double lam1_0 = -2.0,
  double sigma2_lam0 = 0.25, double sigma2_lam1 = 0.25,
  double edge_prob_threshold = 0.6, double edge_prob_adj = 0.75,
  double epsilon_lam = 0.01, int num_step_lam = 10) {
  group_treat = group_treat - 1;
  group_ctrl = group_ctrl - 1;
  int save_start = num_iter - num_save;
  const uword N = static_cast<uword>(N_int);
  const uword G_u = static_cast<uword>(G);
  const uword K_u = static_cast<uword>(K);

  uvec treat_idx = find(cell_t.col(0) == 1);
  uvec ctrl_idx = find(cell_t.col(0) == 0);
  const uword n_treat = treat_idx.n_elem;
  const uword n_ctrl = ctrl_idx.n_elem;

  mat theta_treat_t = theta_t.rows(treat_idx);
  mat theta_ctrl_t  = theta_t.rows(ctrl_idx);
  mat ind_zero_treat = ind_zero.rows(treat_idx);
  mat ind_zero_ctrl  = ind_zero.rows(ctrl_idx);

  Rcpp::Rcout << "=== MCMC Parameter Debug ===\n";
  Rcpp::Rcout << "G: " << G << ", N: " << N << ", K: " << K << "\n";
  Rcpp::Rcout << "treat_idx length: " << treat_idx.n_elem << ", ctrl_idx length: " << ctrl_idx.n_elem << "\n";
  Rcpp::Rcout << "group_treat range: " << min(group_treat) << " ~ " << max(group_treat) << "\n";
  Rcpp::Rcout << "group_ctrl range: " << min(group_ctrl) << " ~ " << max(group_ctrl) << "\n";
  Rcpp::Rcout << "Total iterations: " << num_iter << ", Save iterations: " << num_save << "\n";
  Rcpp::Rcout << "============================\n";
  Rcpp::Rcout << "MCMC started...\n";

  cube invcov_treat_save(G_u, G_u, K_u, fill::zeros);
  cube invcov_ctrl_save(G_u, G_u, K_u, fill::zeros);
  cube edge_treat_save(G_u, G_u, K_u, fill::zeros);
  cube edge_ctrl_save(G_u, G_u, K_u, fill::zeros);
  mat mu_treat_save(G_u, K_u, fill::zeros);
  mat mu_ctrl_save(G_u, K_u, fill::zeros);

  mat theta_treat_save(n_treat, G, fill::zeros);
  mat theta_ctrl_save(n_ctrl, G, fill::zeros);
  mat theta_save(N, G, fill::zeros);
  mat group_treat_save(K_u, n_treat, fill::zeros);
  mat group_ctrl_save(K_u, n_ctrl, fill::zeros);

  mat group_iter(N, num_save, fill::zeros);
  cube mu_treat_iter(G_u, K_u, num_save, fill::zeros);
  cube mu_ctrl_iter(G_u, K_u, num_save, fill::zeros);
  arma::field<cube> invcov_treat_iter(num_save);
  arma::field<cube> invcov_ctrl_iter(num_save);
  for (int s = 0; s < num_save; s++) {
    invcov_treat_iter(s) = cube(G_u, G_u, K_u, fill::zeros);
    invcov_ctrl_iter(s) = cube(G_u, G_u, K_u, fill::zeros);
  }
  //save the results in the iterations for the first column of each matrix
  mat edge_treat_iter(G_u, num_save, fill::zeros);
  mat edge_ctrl_iter(G_u, num_save, fill::zeros);
  mat theta_treat_iter(n_treat, num_save, fill::zeros);
  mat theta_ctrl_iter(n_ctrl, num_save, fill::zeros);
  mat theta_iter(N, num_save, fill::zeros);
  mat lam0_iter(G_u, num_save, fill::zeros);
  mat lam1_iter(G_u, num_save, fill::zeros);
  mat pi_treat_iter(K_u, num_save, fill::zeros);
  mat pi_ctrl_iter(K_u, num_save, fill::zeros);

  vec group_t(N, fill::zeros);
  // treatment group：group_iter[i, j] < K
  // control group：group_iter[i, j] >= K
  group_t(treat_idx) = group_treat;
  group_t(ctrl_idx) = group_ctrl + static_cast<double>(K);

  long seed = randi(distr_param(1,1000000));

  for (int t_iter = 0; t_iter < num_iter; t_iter++) {
    double progress = (t_iter + 1.0) / num_iter * 100.0;
    Rcpp::Rcout << "\rMCMC Progress: [" << t_iter + 1 << "/" << num_iter << "] " << fixed << setprecision(1) << progress << "%";
    Rcpp::Rcout.flush();
    // STEP 1: Update theta
    // treatment group
    seed = randi(distr_param(1, 1000000));
    update_theta(theta_treat_t, ind_zero_treat, mu_treat, invcov_treat, cov_treat,
                 lambda0_t, lambda1_t, group_treat, n_treat, seed, epsilon_theta, num_step_theta);
    // control group
    seed = randi(distr_param(1, 1000000));
    update_theta(theta_ctrl_t, ind_zero_ctrl, mu_ctrl, invcov_ctrl, cov_ctrl,
                 lambda0_t, lambda1_t, group_ctrl, n_ctrl, seed, epsilon_theta, num_step_theta);
    // STEP 2: Update lambda
    seed = randi(distr_param(1,1000000));
    update_lambda(lambda0_t, lambda1_t, theta_t, ind_zero, G, seed, lam0_0, lam1_0, sigma2_lam0, sigma2_lam1, epsilon_lam, num_step_lam);
    // STEP 3: Update mu
    // treatment group
    update_mu(mu_treat, theta_treat_t, invcov_treat, group_treat, G, K, eta_mu, tau_sq_mu);
    // control group
    update_mu(mu_ctrl, theta_ctrl_t, invcov_ctrl, group_ctrl, G, K, eta_mu, tau_sq_mu);
    // STEP 4: Update invcov
    // treatment group
    seed = randi(distr_param(1, 1000000));
    update_invcov(invcov_treat, cov_treat, edge_treat, theta_treat_t, mu_treat, group_treat,
                  ssp_v0, ssp_v1, ssp_l, G, K, seed);
    // control group
    seed = randi(distr_param(1, 1000000));
    update_invcov(invcov_ctrl, cov_ctrl, edge_ctrl, theta_ctrl_t, mu_ctrl, group_ctrl,
                  ssp_v0, ssp_v1, ssp_l, G, K, seed);
    seed = randi(distr_param(1,1000000));
    // STEP 5: Update edge
    seed = randi(distr_param(1, 1000000));
    update_edge_with_matching_by_difficulty(
      edge_treat, edge_ctrl, invcov_treat, invcov_ctrl, ssp_v0, ssp_v1, ssp_xi, edge_prob_threshold, edge_prob_adj,
      G, K, N, static_cast<unsigned long>(seed),
      theta_t, cell_t, group_treat, group_ctrl, treat_idx, ctrl_idx
    );
    // STEP 6: Update group
    // treatment group
    seed = randi(distr_param(1, 1000000));
    update_group(group_treat, invcov_treat, theta_treat_t, mu_treat, pi_treat, G, n_treat, K, seed);
    // control group
    seed = randi(distr_param(1, 1000000));
    update_group(group_ctrl, invcov_ctrl, theta_ctrl_t, mu_ctrl, pi_ctrl, G, n_ctrl, K, seed);
    group_t(treat_idx) = group_treat;
    group_t(ctrl_idx) = group_ctrl + static_cast<double>(K);
    // STEP 7: Update pi
    // treatment group
    update_pi(pi_treat, group_treat, gam, K);
    // control group
    update_pi(pi_ctrl, group_ctrl, gam, K);

    if (t_iter >= save_start) {
      int save_i = t_iter - save_start;
      for (uword i = 0; i < n_treat; i++) {
        int k_t = group_treat(i);
        group_treat_save(k_t, i) += 1;
      }
      for (uword i = 0; i < n_ctrl; i++) {
        int k_c = group_ctrl(i);
        group_ctrl_save(k_c, i) += 1;
      }

      mu_treat_save   += mu_treat;
      mu_ctrl_save    += mu_ctrl;
      invcov_treat_save += invcov_treat;
      invcov_ctrl_save += invcov_treat;
      edge_treat_save  += edge_treat;
      edge_ctrl_save   += edge_ctrl;

      theta_treat_save += theta_treat_t;
      theta_ctrl_save  += theta_ctrl_t;
      theta_save       += theta_t;

      group_iter.col(save_i)      = group_t;
      mu_treat_iter.slice(save_i) = mu_treat;
      mu_ctrl_iter.slice(save_i)  = mu_ctrl;
      invcov_treat_iter(save_i)   = invcov_treat;
      invcov_ctrl_iter(save_i)    = invcov_ctrl;
      edge_treat_iter.col(save_i) = edge_treat.slice(0).col(0);
      edge_ctrl_iter.col(save_i)  = edge_ctrl.slice(0).col(0);
      lam0_iter.col(save_i)       = lambda0_t;
      lam1_iter.col(save_i)       = lambda1_t;
      pi_treat_iter.col(save_i)   = pi_treat;
      pi_ctrl_iter.col(save_i)    = pi_ctrl;

      theta_treat_iter.col(save_i) = theta_treat_t.col(0);
      theta_ctrl_iter.col(save_i)  = theta_ctrl_t.col(0);
      theta_iter.col(save_i)       = theta_t.col(0);
    }
  }
  Rcpp::Rcout << "\nMCMC completed!\n";
  mu_treat_save   /= num_save;
  mu_ctrl_save    /= num_save;
  invcov_treat_save /= num_save;
  invcov_ctrl_save /= num_save;
  edge_treat_save  /= num_save;
  edge_ctrl_save   /= num_save;
  theta_treat_save /= num_save;
  theta_ctrl_save  /= num_save;
  theta_save       /= num_save;

  vec group_treat_est(n_treat);
  vec group_ctrl_est(n_ctrl);
  for (uword i = 0; i < n_treat; i++) {
    group_treat_est(i) = group_treat_save.col(i).index_max() + 1;
  }
  for (uword i = 0; i < n_ctrl; i++) {
    group_ctrl_est(i) = group_ctrl_save.col(i).index_max()+ 1;
  }
  vec group_est(N, fill::zeros);
  group_est(treat_idx) = group_treat_est;
  group_est(ctrl_idx)  = group_ctrl_est + static_cast<double>(K);
  List final_matching = network_propensity_matching_by_difficulty(
    theta_t, cell_t, group_treat, group_ctrl, treat_idx, ctrl_idx,
    invcov_treat, invcov_ctrl,
    G, K);
  return List::create(
    Named("group_est")        = group_est, Named("group_treat_est")  = group_treat_est, Named("group_ctrl_est")   = group_ctrl_est,
    Named("mu_treat_est")     = mu_treat_save, Named("mu_ctrl_est")      = mu_ctrl_save,
    Named("invcov_treat_est") = invcov_treat_save, Named("invcov_ctrl_est")  = invcov_ctrl_save,
    Named("edge_treat_est")   = edge_treat_save, Named("edge_ctrl_est")    = edge_ctrl_save,
    Named("theta_treat_est")  = theta_treat_save, Named("theta_ctrl_est")  = theta_ctrl_save, Named("theta_est")        = theta_save,
    Named("group_iter")       = group_iter,
    Named("mu_treat_iter")    = mu_treat_iter, Named("mu_ctrl_iter")     = mu_ctrl_iter,
    Named("invcov_treat_iter")= invcov_treat_iter, Named("invcov_ctrl_iter") = invcov_ctrl_iter,
    Named("edge_treat_iter")  = edge_treat_iter, Named("edge_ctrl_iter")   = edge_ctrl_iter,
    Named("theta_treat_iter") = theta_treat_iter, Named("theta_ctrl_iter")  = theta_ctrl_iter, Named("theta_iter")       = theta_iter,
    Named("lam0_iter")        = lam0_iter, Named("lam1_iter")        = lam1_iter, Named("pi_treat_iter")    = pi_treat_iter, Named("pi_ctrl_iter")     = pi_ctrl_iter,
    Named("group_treat_save") = group_treat_save, Named("group_ctrl_save")  = group_ctrl_save,
    Named("matching_results")   = final_matching
  );
}

// [[Rcpp::export]]
mat update_mu_R(mat theta_t, cube invcov_t, vec group_t, int G, int K, double eta_mu = 0, double tau_sq_mu = 1) {
  mat mu_updated(G, K);
  mat I_tau_sq = eye(G,G);
  I_tau_sq /= tau_sq_mu;
  for (int k = 0; k < K; k++) {
    uvec ind_k = find(group_t == k);
    vec tmp1 = sum(theta_t.rows(ind_k), 0).t();

    mat invcov_k = invcov_t.slice(k);

    mat COV = inv(ind_k.n_elem * invcov_k + I_tau_sq);
    vec MU = COV * (invcov_k * tmp1 + eta_mu/tau_sq_mu);

    vec mu_k = mvnrnd(MU, COV, 1);
    mu_updated.col(k) = mu_k;
  }
  return mu_updated;
}

//' Update pi parameter in scNPM model
//'
//' @param group_t Cluster assignment vector
//' @param gam Prior parameter
//' @param K Number of clusters
//' @return Dirichlet sampling result for pi
//' @export 
// [[Rcpp::export]]
vec update_pi_R(vec group_t, vec gam, int K) {
  vec pi_new(K);
  vec tmp0(K, fill::zeros);
  for (int k = 0; k < K; k++) {
    uvec ind_j = find(group_t == k);
    tmp0(k) = ind_j.n_elem;
  }
  pi_new = rDirichlet(tmp0+gam);
  return pi_new;
}

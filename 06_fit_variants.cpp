// author: haim krupkin+claude
// date: 11.08.2026
// description: compiled kernel of step 6 (06_cis_scan.R): the weighted regression for every cis variant of one gene. this was mostly implemented for speed and looks linda robust.
// Same math as the R version it replaced (results identical to ~1e-15), ~10x faster, so the scan can run genome-wide.
// For variant v (genotype g, centred) and the gene's scaled data:
//   Z     = A[, a] * g for every column a of A (6 cell types, then the death variables)   -> samples x p
//   Zr    = Z - Q Q'Z                       (remove the base design: proportions, death x proportions, covariates)
//   S     = Zr'Zr,  beta = S^-1 Zr'y_r,  e = y_r - Zr beta
//   h     = h_base + diag(Zr S^-1 Zr')      (leverage of each sample in the full model), capped at 0.99
//   V     = S^-1 [Zr' diag(e^2 / (1 - h)^2) Zr] S^-1        (robust HC3 covariance)
// Returned per variant (one row): the 6 cell-type effects, their 6 SEs, the heterogeneity Wald statistic (are the 6
// effects equal, 5 df), the neutrophil-minus-mean-lymphoid contrast and its SE. Variants whose S cannot be inverted
// get NA. p-values are computed in R.
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::export]]
arma::mat fit_variants_cpp(const arma::mat& G, const arma::mat& A, const arma::mat& Q,
                           const arma::vec& y_r, const arma::vec& h_base,
                           const arma::mat& het_C, const arma::vec& c_nl) {
  const uword nv = G.n_cols;
  mat out(nv, 15);
  out.fill(datum::nan);
  for (uword v = 0; v < nv; v++) {
    // genotype columns, then remove the base design
    mat Z = A.each_col() % G.col(v);
    Z -= Q * (Q.t() * Z);
    // normal equations
    mat S = Z.t() * Z;
    mat S_inv;
    if (!inv_sympd(S_inv, S)) continue;
    vec beta = S_inv * (Z.t() * y_r);
    vec e = y_r - Z * beta;
    // leverage and robust (HC3) covariance
    vec h = h_base + sum((Z * S_inv) % Z, 1);
    h = clamp(h, -datum::inf, 0.99);
    vec w = e / (1.0 - h);
    mat Zw = Z.each_col() % w;
    mat V = S_inv * (Zw.t() * Zw) * S_inv;
    // the 6 cell-type effects and their contrasts
    vec b6 = beta.head(6);
    mat Vb = V.submat(0, 0, 5, 5);
    vec se = sqrt(Vb.diag());
    vec d = het_C * b6;
    mat H = het_C * Vb * het_C.t();
    vec Hd;
    if (!solve(Hd, H, d)) continue;
    double het_stat = dot(d, Hd);
    double nl_est = dot(c_nl, b6);
    double nl_se = std::sqrt(as_scalar(c_nl.t() * Vb * c_nl));
    // one output row
    out.submat(v, 0, v, 5) = b6.t();
    out.submat(v, 6, v, 11) = se.t();
    out(v, 12) = het_stat;
    out(v, 13) = nl_est;
    out(v, 14) = nl_se;
  }
  return out;
}

// 11.10: the same regression with k cell-type columns (k = the cell types kept by the phi gate, 1..6; MODEL phi).
// A = [k cell-type columns | death columns]; het_C is (k-1) x k (empty when k = 1); c_nl has length k (all 0: no contrast).
// Returned per variant: k effects, k SEs, het statistic (NaN when k = 1), nl estimate, nl SE  -> 2k + 3 columns.
// [[Rcpp::export]]
arma::mat fit_variants_k_cpp(const arma::mat& G, const arma::mat& A, const arma::mat& Q,
                             const arma::vec& y_r, const arma::vec& h_base,
                             const arma::mat& het_C, const arma::vec& c_nl, const int k) {
  const uword nv = G.n_cols;
  mat out(nv, 2 * k + 3);
  out.fill(datum::nan);
  for (uword v = 0; v < nv; v++) {
    mat Z = A.each_col() % G.col(v);
    Z -= Q * (Q.t() * Z);
    mat S = Z.t() * Z;
    mat S_inv;
    if (!inv_sympd(S_inv, S)) continue;
    vec beta = S_inv * (Z.t() * y_r);
    vec e = y_r - Z * beta;
    vec h = h_base + sum((Z * S_inv) % Z, 1);
    h = clamp(h, -datum::inf, 0.99);
    vec w = e / (1.0 - h);
    mat Zw = Z.each_col() % w;
    mat V = S_inv * (Zw.t() * Zw) * S_inv;
    vec bk = beta.head(k);
    mat Vb = V.submat(0, 0, k - 1, k - 1);
    vec se = sqrt(Vb.diag());
    out.submat(v, 0, v, k - 1) = bk.t();
    out.submat(v, k, v, 2 * k - 1) = se.t();
    if (k > 1) {
      vec d = het_C * bk;
      mat H = het_C * Vb * het_C.t();
      vec Hd;
      if (solve(Hd, H, d)) out(v, 2 * k) = dot(d, Hd);
    }
    if (any(c_nl != 0)) {
      out(v, 2 * k + 1) = dot(c_nl, bk);
      out(v, 2 * k + 2) = std::sqrt(as_scalar(c_nl.t() * Vb * c_nl));
    }
  }
  return out;
}
#math not fully verified by hk 11/14/2026 as its hard to verifiy

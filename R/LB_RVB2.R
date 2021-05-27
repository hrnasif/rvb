#' Evidence lower bound for RVB2
#'
#' @param mu Vector. Mean in variational approximation
#' @param C Lower triangular matrix. Cholesky decomposition of variance
#' @param y List. Responses per cluster
#' @param X List. Covariates per cluster for fixed effects
#' @param Z List. Covariates per cluster for random effects
#' @param Zt List. Output from transformZi
#' @param Zy Numeric or vector. Z^T y
#' @param etahat List. Contains per cluster MLEs under RVB1 regularized model
#' @param vbeta0 Numeric. Prior variance for fixed effects
#' @param Sinv Matrix. Inverse of prior covariance matrix for random effects
#' @param model Character. Either "poisson" or "binomial"
#' @param m Integer. Number of trials if "binomial"
#' @param n Integer. Number of clusters
#' @param p Integer. Number of fixed effects
#' @param r Integer. Number of random effects
#' @param d Integer. Total number of variables
#' @param N Integer. Number of iterations for ELBO average
#'
#' @return Value of the evidence lower bound based on RVB2 optimization
#'
#' @export
LB_RVB2 <- function(mu, C, y, X, Z, Zt, Zy, etahat, vbeta0, Sinv, model, m, n, p, r, d, N = 1000){

  LB = 0
  for (i in 1:N){
    s = stats::rnorm(d)
    ttheta = mu + C %*% s
    LB = LB + (logpy_RVB2(ttheta, y, X, Z, Zt, Zy, etahat, vbeta0, Sinv, model, m, n, p, r) +
      sum(log(diag(C))) + 0.5* crossprod(s))/N
  }
}

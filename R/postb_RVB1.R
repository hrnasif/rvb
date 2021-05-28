#' Sample from posterior under RVB1
#'
#' @param RVB List. Output from Alg_RVB1
#' @param y List. Responses per cluster
#' @param X List. Covariates per cluster for fixed effects
#' @param Z List. Covariates per cluster for random effects#'
#' @param etahat List. Estimate of canonical parameter about which to approximate
#' @param N Integer. Sample size
#' @param L Integer. Number of clusters
#' @param model Character. Either "poisson" or "binomial"
#' @param m Positive integer. Number of trials in binomial. If poisson, keep as m = 1.
#'
#' @return Matrix of posterior samples
#'
#' @export
postb_RVB1_b1 <- function(RVB, y, X, Z, etahat, N, L, model, m = 1){
  n = length(y)
  p = ncol(X[[1]])
  r = ncol(Z[[1]])

  mu = RVB$mu
  C = RVB$C
  d = length(mu)

  if (m == 1){rep(1, n)}
  Zg_ZH <- .ZgH(y, Z, model, etahat, m, n, r)
  Zg = Zg_ZH$Zg
  ZH = Zg_ZH$ZH

  b = matrix(0, nrow = r*L, ncol = N)
  llim = (0:(n-1))*r + 1; rlim = (1:n)*r

  set.seed(572)
  for (j in 1:N){
    s = stats::rnorm(d)
    ttheta = C %*% s + mu
    beta = ttheta[(n*r+1):(n*r+p)]

    if (r == 1){
      W = exp(ttheta[d])
      Omega = W^2

      for (i in 1:L){
        Lambdai = 1/(Omega + crossprod(ZH[[i]], Z[[i]]))
        Lit = sqrt(Lambdai)
        lambda_i = Lambdai*(Zg[[i]] - crossprod(ZH[[i]], X[[i]] %*% beta))
        b[llim[i]:rlim[i], j] = Lit*ttheta[i] + lambda_i
      }
    }
    else{
      omega = ttheta[(n*r+p+1):d]
      W = make_lower_tri_matrix(r, omega)
      diag(W) = exp(diag(W))
      Omega = tcrossprod(W)

      for (i in 1:L){
        Lambdai = (Omega + ZH[[i]] %*% Z[[i]]) %>% chol() %>% chol2inv()
        Lit = chol(Lambdai)
        lambda_i = Lambdai %*% (Zg[i,] - ZH[[i]] %*% X[[i]] %*% beta)
        b[llim[i]:rlim[i], j] = crossprod(Lit, ttheta[((i - 1)*r + 1): (i*r)]) + lambda_i
      }
    }
  }
  return(b)
}

#' Log joint density under RVB2
#'
#' @param ttheta Numeric vector. All global and local variables.
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
#'
#' @return Numeric. Log joint density for RVB2

#' @export
logpy_RVB2 <- function(ttheta, y, X, Z, Zt, Zy, etahat, vbeta0, Sinv, model, m,
                       n, p, r){

  global <- global_var_components(ttheta, n, p, r)
  t1 = global$t1; beta = global$beta; Omega = global$Omega;

  L = t1 - 0.5*sum(Sinv * Omega) - 0.5 * crossprod(beta)/vbeta0

  if (r == 1){
    for (i in 1:n){
      Xibeta = X[[i]] %*% beta
      bhati = fisherscoring(etahat[[i]], Zy[i], Z[[i]], Zt[[i]], Xibeta, Omega, model, m[i])
      etahati = Xibeta + Z[[i]] %*% bhati
      Lit = 1/sqrt(Omega + crossprod(h2(model, etahati, m[i]) * Z[[i]], Z[[i]]))
      bi = Lit * ttheta[i] + bhati
      etai = Xibeta + Z[[i]] %*% bi
      L = L + crossprod(y[[i]], etai) - sum(h0(model, etai, m[i])) + log(Lit) - 0.5*Omega*bi^2
    }
  }

  else{
    for (i in 1:n){
      Xibeta = X[[i]] %*% beta
      bhati = fisherscoring(etahat[[i]], Zy[,i], Z[[i]], Zt[[i]], Xibeta, Omega, model, m[i])
      etahati = Xibeta + Z[[i]] %*% bhati
      Zhi = diag(h2(model, etahati, m[i])) * Z[[i]]
      Lambdai = (Omega + crossprod(Zhi, Z[[i]])) %>% chol() %>% chol2inv()
      Lit = chol(Lambdai)
      bi = crossprod(Lit, ttheta[((i - 1)*r + 1): (i*r)]) + bhati
      etai = Xibeta + Z[[i]] %*% bi
      L = L + crossprod(y[[i]], etai) - sum(h0(model, etai, m[i])) + .logdet(Lit) -
        0.5*crossprod(bi, Omega %*% bi)
    }
  }

  return(L %>% drop())
}

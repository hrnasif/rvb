#' Log joint density and gradient for RVB2
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
#' @param d Integer. Total number of variables
#'
#' @return List. Containing log joint density and the gradient vector
#' for the log joint density under RVB2
#'
#' @export
Lgrad_RVB2 <- function(ttheta, y, X, Z, Zt, Zy, etahat, vbeta0, Sinv, model, m, n, p, r, d){
  global <- global_var_components(ttheta, n, p, r)
  t1 <- global$t1; beta <- global$beta; rseq <- global$rseq;
  W <- global$W; Omega <- global$Omega

  L = t1 - 0.5*sum(Sinv*Omega) - 0.5*crossprod(beta)/vbeta0 # log joint density
  g = numeric(d) # gradient vector
  S1 = -beta/vbeta0

  if (r == 1){
    for (i in 1:n){
      tbi = ttheta[i]
      Xibeta = X[[i]]%*%beta
      bhati = fisherscoring(etahat[[i]], Zy[[i]], Z[[i]], Zt[[i]], Xibeta, Omega, model, m[i])
      etahati = Xibeta + Z[[i]] %*% bhati
      ZHi = h2(model, etahati, m[i]) * Z[[i]]
      Lambdai = 1/(Omega + crossprod(ZHi, Z[[i]]))
      Lit = sqrt(Lambdai)
      bi = Lit*tbi + bhati
      etai = Xibeta + Z[[i]] %*% bi
      ai = -Omega * bi
      L = L + crossprod(y[[i]], etai) - sum(h0(model, etai, m[i])) + log(Lit) + 0.5*bi*ai
      gi = y[[i]] - h1(model, etai, m[i])
      ai = ai + crossprod(Z[[i]], gi)
      gbi = Lit * ai
      t1 = (Lambdai*(1 + gbi*tbi)) %>% drop()
      t2 = t1*(Z[[i]]^2) * h3(model, etahati, m[i])
      ai = ai - 0.5*crossprod(Z[[i]], t2)
      ai = (Lambdai*ai) %>% drop()
      S1 = S1 + crossprod(X[[i]], gi - ai*ZHi - 0.5*t2) # for beta
      Sinv = Sinv + 2*ai*bhati + bi^2 + t1
      g[i] = gbi
    }
    g[(n*r + 1): (n*r+p)] = S1
    g[d] = rseq - Sinv * W^2
  }

  else{
    for (i in 1:n){
      tbi = ttheta[((i - 1)*r + 1):(i*r)]
      Xibeta = X[[i]] %*% beta
      bhati = fisherscoring(etahat[[i]], Zy[,i], Z[[i]], Zt[[i]], Xibeta, Omega, model, m[i])
      etahati = Xibeta + Z[[i]] %*% bhati
      ZHi = diag(h2(model, etahati, m[i])) * Z[[i]]
      Lambdai = (Omega + crossprod(ZHi, Z[[i]])) %>% chol() %>% chol2inv()
      Lit = chol(Lambdai)
      bi = crossprod(Lit, tbi) + bhati
      etai = Xibeta + Z[[i]] %*% bi
      ai = -Omega %*% bi
      L = L + crossprod(y[[i]], etai) - sum(h0(model, etai, m[i])) + .logdet(Lit) + 0.5*crossprod(bi,ai)
      gi = y[[i]] - h1(model, etai, m[i])
      ai = ai + crossprod(Z[[i]], gi)
      gbi = Lit %*% ai
      gbi_tbi <- tcrossprod(gbi, tbi)
      gbi_tbi[upper.tri(gbi_tbi)] = gbi_tbi[lower.tri(gbi_tbi)] # Make symmetric
      t1 = Lambdai + crossprod(Lit, gbi_tbi %*% Lit)
      t2 = 0.5*diag(tcrossprod(Z[[i]] %*% t1, Z[[i]])) * h3(model, etahati, m[i])
      ai = ai - crossprod(Z[[i]], t2)
      ai = Lambdai %*% ai
      S1 = S1 + crossprod(X[[i]], gi - ZHi %*% ai - t2)
      ci = tcrossprod(ai, bhati)
      Sinv = Sinv + ci + t(ci) + tcrossprod(bi) + t1
      g[((i - 1)*r + 1): (i*r)] = gbi
    }
    g[(n*r + 1): (n*r+p)] = S1
    A = -get_lower_tri_vector(Sinv %*% W)
    A[.diag_locs(r)] = A[.diag_locs(r)] * diag(W) + rseq
    g[(n*r + p + 1):d] = A
  }

  return(list(L = (L %>% drop()),
              g = g))
}

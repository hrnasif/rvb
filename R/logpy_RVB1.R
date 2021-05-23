logpy_RVB1 <- function(ttheta, y, X, Z, Zg, ZH, vbeta0, Sinv, model, m, n, p, r){
  global <- global_var_components(ttheta, n, p, r)

  if (r == 1){
    L = global$t1 - 0.5 * Sinv * global$Omega -
      0.5 * crossprod(global$beta)/vbeta0
    for (i in 1:n){
      Xibeta = X[[i]] %*% global$beta
      Lambdai = 1/(global$Omega + crossprod(ZH[[i]], Z[[i]]))
      Lit = sqrt(Lambdai)
      bi = Lambdai * (Zg[i] - crossprod(ZH[[i]], Xibeta))
      bi = bi + Lit*ttheta[i*r]
      etai = Xibeta + Z[[i]]*bi
      L = L + crossprod(y[[i]], etai) - sum(h0(model, etai, m[i])) + log(Lit) -
        0.5*global$Omega*bi^2
    }
  }
  else{
    L = global$t1 - 0.5*sum(Sinv*global$Omega) -
      0.5*crossprod(global$beta)/vbeta0
    for (i in 1:n){
      Xibeta = X[[i]] %*% global$beta
      Lambdai = (global$Omega + ZH[[i]] %*% Z[[i]]) %>% chol() %>% chol2inv()
      Lit = chol(Lambdai)
      bi = Lambdai %*% (Zg[i,] - (ZH[[i]] %*% Xibeta))
      bi = bi + crossprod(Lit, ttheta[((i - 1)*r+1):(i*r)])
      etai = Xibeta + Z[[i]] %*% bi
      L = L + crossprod(y[i], etai) - sum(h0(model, etai, m[i])) + .logdet(Lit) -
        0.5*crossprod(bi, global$Omega %*% bi)
    }
  }
  return(L %>% drop())
}

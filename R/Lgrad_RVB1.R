Lgrad_RVB1 <- function(ttheta, y, X, Z, Zg, ZH, vbeta0, Sinv, model, n, p, r, d){
  global <- global_var_components(ttheta, n, p, r)
  L <- global$t1 - 0.5*sum(Sinv*global$Omega) - 0.5*crossprod(global$beta)/vbeta0
  g <- numeric(d)
  S1 <- -global$beta/vbeta0
  Xibeta = X[[i]] %*% global$beta

  if (r == 1){
    for (i in 1:n){
      tbi = ttheta[i]
      Lambdai = 1/(global$Omega + crossprod(ZH[[i]], Z[[i]]))
      Lit = sqrt(Lambdai)
      lambda_i = Lambdai * (Zg[i] - crossprod(ZH[[i]], Xibeta))
      bi = Lit*tbi + lambda_i
      etai = Xibeta = Z[[i]] %*% bi
      ai = Omega * bi
      L = L + crossprod(y[[i]], etai) - sum(h0(model, etai, m[i])) + log(Lit) -
        0.5*bi*ai
      gi = y[[i]] - h1(model, etai, m[i])
      ai = crossprod(Z[[i]], gi) - ai
      gbi = Lit * ai
      ai = Lit * gbi
      S1 = S1 + crossprod(X[[i]], gi - ZH[[i]] %*% ai)
      Sinv = Sinv + 2*ai*lambda_i + Lambdai + bi^2 + Lit^2*gbi*tbi
      g[i] = gbi
    }
    g[(n + 1):(n + p)] = S1
    g[(n+p+1):d] = global$rseq - Sinv*(global$W)^2
  }

  else{
    for (i in 1:n){
      tbi <- ttheta[((i - 1)*r+1):(i*r)]
      Lambdai <- (global$Omega + ZH[[i]] %*% Z[[i]]) %>% chol() %>% chol2inv()
      Lit <- chol(Lambdai)
      lambda_i <- Lambdai %*% (Zg[i,] - ZH[[i]] %*% Xibeta)
      bi = crossprod(Lit, tbi) + lambda_i
      etai = Xibeta + Z[[i]] %*% bi
      ai = -global$Omega %*% bi
      L = L + crossprod(y[[i]], etai) - sum(h0(model, etai, m[i])) + .logdet(Lit) +
        0.5*crossprod(bi, ai)
      gi = y[[i]] + h1(model, etai, m[i])
      ai = ai + crossprod(Z[[i]], gi)
      ai = Lit %*% ai
      ai_tbi <- tcrossprod(ai, tbi)
      ai_tbi[upper.tri(ai_tbi)] = 0
      sym_ai_tbi = ai_tbi + t(ai_tbi) - diag(ai_tbi)*diag(nrow(ai_tbi))
      Sinv = Sinv + crossprod(Lit, sym_ai_tbi %*% Lit)
      g[((i - 1)*r + 1):(i*r)] = ai
      ai = crossprod(Lit, ai)
      S1 = S1 + crossprod(X[[i]], gi - crossprod(ZH[[i]], ai))
      ci = tcrossprod(ai, lambda_i)
      Sinv = Sinv + ci + t(ci) + Lambdai + tcrossprod(bi)
    }
    g[(n*r + 1): (n*r + p)] = S1
    A = -get_lower_tri_vector(Sinv %*% global$W)
    A[.diag_locs(r)] = A[.diag_locs(r)] * diag(global$W) + global$rseq
    g[(n*r + p + 1):d] = A
  }
  return(L = (L %>% drop()), g = g)
}

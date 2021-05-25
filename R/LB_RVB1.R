LB_RVB1 <- function(mu, C, y, X, Z, Zg, ZH, vbeta0, Sinv, model, m, n, p, r, d,
                    N = 1000){
  LB = 0
  for (i in 1:N){
    s = stats::rnorm(d)
    ttheta = mu + C %*% s
    LB = LB + logpy_RVB1(ttheta, y, X, Z, Zg, ZH, vbeta0, Sinv, model, m, n, p, r) +
      sum(log(diag(C))) + 0.5*crossprod(s)/N
  }
  return(LB)
}

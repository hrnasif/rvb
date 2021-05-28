#' RVB1 Algorithm implementation
#'
#' @param y List. Responses per cluster
#' @param X List. Covariates per cluster for fixed effects
#' @param Z List. Covariates per cluster for random effects
#' @param Wprior List. Wishart prior for random effect covariance
#' @param etahat List. Estimate of canonical parameter about which to approximate
#' @param model Character. Either "poisson" or "binomial"
#' @param m Integer. Number of trials if model is binomial. Keep m = 1 if model is "poisson"
#'
#' @return List containing posterior covariance C, posterior mean mu, the ELBO values per
#' each 1000 iterations, the run duration, and the final ELBO value.
#'
#' @export
Alg_RVB1 <- function(y, X, Z, Wprior, etahat, model, m = 1){
  start_time = Sys.time()

  Interval = 1000
  tol = 6
  Tr = 1.5e5 # max number of iterations
  n = length(y)
  p = ncol(X[[1]])
  r = ncol(Z[[1]])
  G = p + 0.5*r*(r+1)
  d = n*r + G
  if (m == 1){ m = rep(1, n)}

  indices <-  indices_for_posterior_covariance(n,r,p) # indices for nonzero elmnts in C
  I = indices$row_indices
  J = indices$col_indices
  Cdiag = which(I == J) # indices for diag elmnts in Cstar
  nz = length(I) # Number of non-zero elements in C
  vbeta0 = 100 # prior for beta
  Zg_ZH = .ZgH(y, Z, model, etahat, m, n, r)
  Sinv = Wprior$Sinv
  mu = numeric(d)
  C = diag(d)
  gbal = (n*r + 1):d
  diag(C)[gbal] = 0.1
  Cstar = numeric(nz)
  find_id = length(Cdiag) - match(gbal, rev(I[I == J])) + 1
  id = Cdiag[find_id[!is.na(find_id)]]
  Cstar[id] = log(0.1)
  par = numeric(Tr/Interval)

  be1 = 0.9; be2 = 0.999; alpha = 0.001; epsilon = 1e-8;
  mtmu = numeric(d); vtmu = numeric(d); # First and second moment vectors for Adam
  mtCstar = numeric(nz); vtCstar = numeric(nz); # Same for Cstar
  meanLB = 0; it = 0; gradlr = 1; j = 0;

  while ((it < Tr) & (gradlr > 0)){
    it = it + 1
    s = stats::rnorm(d)
    ttheta = mu + C %*% s
    L_and_grad = Lgrad_RVB1(ttheta, y, X, Z, Zg_ZH$Zg, Zg_ZH$ZH, vbeta0, Sinv, model, n, p, r, d)
    meanLB = meanLB + (L_and_grad$L + 0.5*crossprod(s) + sum(log(diag(C))))/Interval
    L_and_grad$g = L_and_grad$g + solve(t(C), s)

    mtmu = be1*mtmu + (1 - be1)*L_and_grad$g
    vtmu = be2*vtmu + (1 - be2)*(L_and_grad$g)^2
    mtmuhat = mtmu/(1 - be1^(it))
    vtmuhat = vtmu/(1 - be2^(it))
    mu = mu + alpha*mtmuhat/(sqrt(vtmuhat) + epsilon)

    Cstargrad = L_and_grad$g[I] * s[J]
    Cstargrad[Cdiag] = diag(C) * Cstargrad[Cdiag]

    mtCstar = be1*mtCstar + (1 - be1)*Cstargrad
    vtCstar = be2*vtCstar + (1 - be2)*(Cstargrad^2)
    mtCstarhat = mtCstar/(1 - be1^(it))
    vtCstarhat = vtCstar/(1 - be2^(it))
    Cstar = Cstar + alpha*mtCstarhat/(sqrt(vtCstarhat) + epsilon)

    C = Matrix::sparseMatrix(i = I, j = J, x = Cstar, dims = c(d,d)) %>% as.matrix()
    diag(C) = exp(diag(C))

    if (it %% Interval == 0){
      j = it/Interval
      par[j] = meanLB
      gradlr = GLB(j, tol, par)
      cat("Iteration:", it, " meanLB=", round(meanLB, 3), " gradlr=", round(gradlr,3))
      meanLB = 0
      cat(" mean=",round(mu[(n*r+1):d], 2),"\n")
    }
  }
  end_time = Sys.time()
  dur = end_time - start_time
  par = par[1:j]
  LB = LB_RVB1(mu, C, y, X, Z, Zg_ZH$Zg, Zg_ZH$ZH, vbeta0, Sinv, model, m, n, p, r, d)
  cat("Final LB:", LB, " Duration:", dur)

  return(list(C = C, mu = mu, par = par, dur = dur, LB = LB))
}

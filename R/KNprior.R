#' Kass and Natarajan prior for random effects precision matrix
#'
#' @param model Character. Either "poisson" or "binomial"
#' @param pred Vector. GLM predictions in y-scale
#' @param Z List. Covariates for random effects
#' @param Integer. Number of trials if model = "binomial". If "poisson", leave as 1.
#'
#' @return List containing Wishart prior parameter nu, S-inverse, and S
#'
#' @export
KNprior <- function(model, pred, Z, m = 1){
  n = length(Z)
  r = ncol(Z[[1]])
  vni = sapply(Z, nrow) %>% unname()
  Sm = matrix(0, r, r)
  startindex = c(1, 1+cumsum(vni)[1:n-1])
  endindex = cumsum(vni)
  if (model == "binomial"){
    if(length(m) == 1){pred = rep(m, times = vni)*pred*(1 - pred)}
    else{pred = m * pred * (1 - pred) }
  }

  for(i in 1:n){
    Sm = Sm + crossprod(Z[[i]], pred[startindex[i]:endindex[i]]*Z[[i]])
  }
  Sm = Sm/n
  R = solve(Sm)
  if (r == 1){
    Wprior = list(nu = r, Sinv = R[1], S = Sm[1])
  }
  else{
    Wprior = list(nu = r+1, Sinv = (r+1)*R, S = Sm/(r + 1))
  }
  return(Wprior)
}

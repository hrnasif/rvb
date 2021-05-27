#' Conditional density for local variables bi
#'
#' @param bi Numeric or vector. Local variables for cluster i
#' @param Zyi Numeric or vector. Zi^T yi
#' @param Zi Vector. Covariates for local variables
#' @param Xibeta Vector. Xi^T beta
#' @param Omega Matrix. Local variable prior precision matrix
#' @param model Character. Either "poisson" or "binomial"
#' @param mi Numeric. Number of trials if binomial
#'
#' @return Conditional density of bi
#'
#' @export
biConditional <- function(bi, Zyi, Zi, Xibeta, Omega, model, mi){
  if (length(bi) == 1){
    cond <- (Zyi - 0.5*Omega*bi)*bi - sum(h0(model, Xibeta + Zi %*% bi, mi))
  }
  else{
    cond <- crossprod(Zyi - 0.5*Omega %*%  bi, bi) - sum(h0(model, Xibeta + Zi %*% bi, mi))
  }
  return(cond)
}

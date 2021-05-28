#' Fisher scoring for RVB2
#'
#' @param etahati Vector. MLE under the regularized model from RVB1
#' @param Zyi Numeric or vector. Zi^T yi
#' @param Zi Vector. Covariates for local variables
#' @param Zti Vector. Output from transformZi
#' @param Xibeta Vector. Xi^T beta
#' @param Omega Matrix. Local variable prior precision matrix
#' @param model Character. Either "poisson" or "binomial"
#' @param mi Numeric. Number of trials if binomial
#'
#' @return A numeric or vector representing the mode of the conditional posterior
#' density of the local variables bi given the global variables and data.
#'
#' @export
fisherscoring <- function(etahati, Zyi, Zi, Zti, Xibeta, Omega, model, mi){
  dif = 1; it = 0;
  if (ncol(Zti) == 1){
    bihat = crossprod(Zti, etahati - Xibeta)
    Ciold = biConditional(bihat, Zyi, Zi, Xibeta, Omega, model, mi)
    while (dif > 1e-4 & it <= 8){
      it = it+1
      etahati = Xibeta + Zi %*% bihat
      Lambdai_inv = Omega + crossprod(h2(model, etahati, mi)*Zi, Zi)
      bihat = bihat + (Zyi - crossprod(Zi, h1(model, etahati,mi)) - Omega*bihat)/Lambdai_inv
      Ci = biConditional(bihat, Zyi, Zi, Xibeta, Omega, model, mi)
      dif = abs((Ci - Ciold)/Ciold)
      Ciold = Ci
    }
  }
  else{
    bihat = Zti %*% (etahati - Xibeta)
    Ciold = biConditional(bihat, Zyi, Zi, Xibeta, Omega, model, mi)
    while (dif > 1e-4 & it <= 8){
      it = it+1
      etahati = Xibeta + Zi %*% bihat
      Zhi = diag(h2(model, etahati, mi))*Zi
      Lambdai_inv = Omega + crossprod(Zhi, Zi)
      Lambdai_inv = chol(Lambdai_inv) %>% chol2inv()
      bihat = bihat + Lambdai_inv %*% (Zyi - crossprod(Zi, h1(model, etahati,mi)) -
                                        Omega %*% bihat)
      Ci = biConditional(bihat, Zyi, Zi, Xibeta, Omega, model, mi)
      dif = abs((Ci - Ciold)/Ciold)
      Ciold = Ci
    }
  }
  return(bihat)
}

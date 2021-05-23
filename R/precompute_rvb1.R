#' Precompute for RVB1.
#'
#' See \lambda_i in equation 7 of Tan (2021).
#'
#' @param y List of numeric vectors. y[[i]] is the response vector for cluster i.
#' @param Z List of vectors (r = 1) or matrices. Z[[i]] is the Z matrix for cluster i.
#' @param model Character. Binomial or Poisson model
#' @param etahat List of numeric vectors. etahat[[i]] is the MLE for cluster i
#' @param m Numeric. Number of trials in Binomial model.
#' @param n Numeric. Number of clusters.
#' @param r Numeric. Number of local variables (dimensionality of random effects).
#'
#' @return List of precomputed vals
#'
#' @keywords internal
.ZgH <- function(y, Z, model, etahat, m, n, r){
  ZH <- vector(mode = "list", length = n)
  if (r == 1){
    Zg = numeric(length = n)
    for (i in 1:n){
      ZH[[i]] <- Z[[i]] * h2(model,etahat[[i]], m[i])
      Zg[i] <- (crossprod(Z[[i]], y[[i]] - h1(model, etahat[[i]], m[i])) +
                  crossprod(ZH[[i]], etahat[[i]])) %>% drop()
    }
  }
  else{
    Zg <- matrix(nrow = n, ncol = r)
    for (i in 1:n){
      ZH[[i]] = (Z[[i]] * h2(model, etahat[[i]], m[i])) %>% t()
      Zg[i,] = crossprod(Z[[i]], y[[i]] - h1(model, etahat[[i]], m[i])) +
                  ZH[[i]] %*% etahat[[i]]
    }
  }
  return(list(Zg = Zg, ZH = ZH))
}

#' Gradient of least squares regression line
#'
#' @param j Iteration in thousands.
#' @param tol Iteration tolerance
#' @param par parameter
#'
#' @return Gradient of least squares regression line fitted to the last means.

#' @export
GLB <- function(j, tol, par){
  gradlr = 1
  if (1 < j & j <= tol){
    linreg <- lm(par[1:j] ~ c(1:j)) %>% coef() %>% unname()
    gradlr = linreg[2]
  }
  else if(j > tol){
    linreg <- lm(par[(j - tol + 1):j] ~ c((j - tol + 1):j)) %>% coef() %>% unname()
    gradlr = linreg[2]
  }
  return(gradlr)
}

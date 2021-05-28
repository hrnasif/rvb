#' Transform Zi for RVB2
#'
#' @param Zi A vector or matrix
#'
#' @return (Zi^T Zi)^{-1} Zi
#'
#'@keywords internal
.transformZi <- function(Zi){
  if(ncol(Zi) == 1){
    Zi = Zi/(crossprod(Zi) %>% drop())
  }
  else{
    Zi = solve(crossprod(Zi), t(Zi))
  }
  return(Zi)
}

#' Vector from lower triangular matrix entries
#'
#' Extract the lower triangular matrix entries into a vector, stacking columns left to right.
#'
#' @param W an n x n matrix
#'
#' @return Vector of size n(n+1)/2 containing lower triangular matrix entries.
#'
#' @export
get_lower_tri_vector <- function(W){
  return(W[lower.tri(W, diag = T)])
}

#' Matrix from lower triangular entries vector
#'
#' Construct a lower triangular matrix from a vector of entries
#'
#' @param n Integer. Size of the output matrix
#' @param omega Numeric vector of size n(n+1)/2. Entries to make matrix.
#'
#' @return An n x n lower triangular matrix with `omega` entries column-wise
#'
#' @export
make_lower_tri_matrix <- function(n,omega){
  W = matrix(0, nrow = n, ncol = n)
  W[lower.tri(W,diag = T)] = omega
  return(W)
}


#' Indices of diagonal elements when in lower triangular vector
#'
#' @param n An integer. Dimenson of matrix
#'
#' @keywords internal
.diag_locs <- function(n){
  c(1, cumsum(n:2)+1)
}




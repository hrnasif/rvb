#' Indices for posterior covariance
#'
#' Finds the indices for the non-zero elements in the Cholesky decomposition
#' of the posterior covariance matrix. The posterior covariance is assumed to
#' to be a block matrix of corresponding to the local and global variables
#'
#' @param n Integer. Number of subjects/clusters.
#' @param r Integer. Dimensionality of random effects
#' @param p Integer. Dimensionality of fixed effects
#'
#' @return List containing row and column indices for non-zero elements.
#'
#' @importFrom magrittr %>%
#'
#' @export
indices_for_posterior_covariance <- function(n, r, p){
  f <- 0.5*r*(r+1)
  G <- p + f
  e <- 0.5 * G * (G+1)
  nz <- e + f*n

  I = integer(nz)
  J = integer(nz)

  # Obtaining indices for local variable blocks in the full matrix
  # First get indices in one block
  local_block_indices <- diag(r) %>% lower.tri(diag = T) %>% which(arr.ind = T)
  local_block_row_indices <- local_block_indices[,1]
  local_block_col_indices <- local_block_indices[,2]

  # Iterate over entire matrix
  for (i in 1:n){
    I[((i-1) * f + 1):(i*f)] = local_block_row_indices + r*(i-1)
    J[((i-1) * f + 1):(i*f)] = local_block_col_indices + r*(i-1)
  }

  # Now for the global variable block

  global_block_indices <- diag(G) %>% lower.tri(diag = T) %>% which(arr.ind = T)

  I[(n*f+1):(n*f+e)] = r*n + global_block_indices[,1]
  J[(n*f+1):(n*f+e)] = r*n + global_block_indices[,2]

  return(list(row_indices = I, col_indices = J))
}

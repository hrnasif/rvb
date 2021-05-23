#' Get global variable density components
#'
#' Obtains log joint density components for global variables, to be used when
#' computing the log joint density. Values are due to the default conjugate
#' Wishart prior from Kass and Natarajan (2006). See Tan (2021) for more details.
#'
#' @param theta A vector of length p + r(r+1)/2 + nr, containing first local then global variables
#' @param n Integer. Number of total observations in data.
#' @param p Integer. Dimension of fixed effects (beta) in a GLMM
#' @param r Integer. Dimension of random effects (b) in a GLMM
#'
#' @return List containing components needed to compute log joint density. In reverse order,
#' includes the precision matrix Omega, its Cholesky decomposition W, the vector of fixed
#' effects beta, the weighted sum of the log diagonal of W, and the sequence used to determine
#' the weights.
#' @export
global_var_components <- function(theta, n, p, r){
    nu <- ifelse(r == 1, r, r + 1) # See Kass and Natarajan prior
    beta <- theta[(n*r+1):(n*r+p)]
    omega <- theta[(n*r+p+1):length(theta)]

    if(r == 1){
      rseq <- n + nu
      log_diag_weighted_sum <- omega * rseq
      W = exp(omega)
      Omega = W^2
    }

    else{
      W <- make_lower_tri_matrix(r, omega)
      rseq <- (n+nu):(n + nu + 1 - r)
      log_diag_weighted_sum <- sum(rseq * diag(W))
      diag(W) <- exp(diag(W))
      Omega <- tcrossprod(W)
    }

    return(list(rseq = rseq,
                t1 = log_diag_weighted_sum,
                beta = beta,
                W = W,
                Omega = Omega))
}

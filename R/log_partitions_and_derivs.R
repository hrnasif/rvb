#' Log partition function
#'
#' @param model Character. Either "binomial" or "poisson"
#' @param x Numeric. Probability for binomial, mean for poisson
#' @param m Integer. Count for binomial
#'
#' @return Log partition function for the exponential family
#'
#' @export
h0 <- function(model, x, m){
  # binomial/bernoulli model
  if(model == "binomial"){
    output <- m*log(1 + exp(x))
  }

  # Poisson model
  else{
    output <- exp(x)
  }
  return(output)
}

#' First derivative of log partition
#'
#' @param model Character. Either "binomial" or "poisson"
#' @param x Numeric. Probability for binomial, mean for poisson
#' @param m Integer. Count for binomial
#'
#' @return First derivative of log partition function.
#'
#' @export
h1 <- function(model, x, m){
  if(model == "binomial"){
    output <- m/(1 + exp(-x))
  }

  else{
    output <- exp(x)
  }

  return(output)
}

#' Second derivative of log partition
#'
#' @param model Character. Either "binomial" or "poisson"
#' @param x Numeric. Probability for binomial, mean for poisson
#' @param m Integer. Count for binomial
#'
#' @return Second derivative of log partition function.
#'
#' @export
h2 <- function(model, x, m){
  if (model == "binomial"){
    output <- m/(2 + exp(x) + exp(-x))
  }

  else{
    output <- exp(x)
  }

  return(output)
}

#' Third derivative of log partition
#'
#' @param model Character. Either "binomial" or "poisson"
#' @param x Numeric. Probability for binomial, mean for poisson
#' @param m Integer. Count for binomial
#'
#' @return Third derivative of log partition function.
#'
#' @export
h3 <- function(model, x, m){
  if (model == "binomial"){
    output <- m*exp(-x)*(exp(-x) - 1)/(1 + exp(-x))^3
  }

  else{
    output <- exp(x)
  }

  return(output)
}

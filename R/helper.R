#' Center data
#'
#' @param x data vector
#'
#' @return data centered at 0
#'
#' @export
center <- function(x){
  return(x - mean(x))
}

#' Standardize data
#'
#' @param x data vector
#'
#' @return data with mean 0 and std 1
#'
#' @export
standardize <- function(x){
  return((x - mean(x))/stats::sd(x))
}

#' Logit transformation
#'
#' @param x data vector
#'
#' @return logit of x
#'
#' @keywords internal
.logit <- function(x){
  return(stats::qlogis(x))
}


#' Log determinant
#'
#' @param x A full rank square matrix
#'
#' @return log determinant of x
#'
#' @keywords internal
.logdet <- function(x){
  return(determinant(x, logarithm = T)[[1]][1])
}

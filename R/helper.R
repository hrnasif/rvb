#' Center data
#'
#' @param x data vector
#'
#' @return data centered at 0
#'
#' @keywords internal
.center <- function(x){
  return(x - mean(x))
}

#' Standardize data
#'
#' @param x data vector
#'
#' @return data with mean 0 and std 1
#'
#' @keywords internal
.standardize <- function(x){
  return((x - mean(x))/sd(x))
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

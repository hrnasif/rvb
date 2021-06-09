#' HERS data
#'
#' Cleaned data set of a portion of the  Heart and Estrogen/Progestin Study
#' (HERS, Hulley et al., 1998) used to demonstrate the RVB algorithm
#'
#' @format A data frame with 9127 rows and 6 variables:
#' \describe{
#'   \item{id}{unique patient id}
#'   \item{response}{Indicator, = 1 if systolic blood pressure > 140}
#'   \item{htn}{Indicator, = 1 if high blood pressure medication was used}
#'   \item{bmi}{Body Mass Index, standardized}
#'   \item{age}{Age at baseline visit, standardized}
#'   \item{visit}{visit number from 0 to 5, zero centered and divided by 2.5}
#' }
#' @source \url{https://regression.ucsf.edu/second-edition/data-examples-and-problems}
"rvb_hers"

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biConditional.R
\name{biConditional}
\alias{biConditional}
\title{Conditional density for local variables bi}
\usage{
biConditional(bi, Zyi, Zi, Xibeta, Omega, model, mi)
}
\arguments{
\item{bi}{Numeric or vector. Local variables for cluster i}

\item{Zyi}{Numeric or vector. Zi^T yi}

\item{Zi}{Vector. Covariates for local variables}

\item{Xibeta}{Vector. Xi^T beta}

\item{Omega}{Matrix. Local variable prior precision matrix}

\item{model}{Character. Either "poisson" or "binomial"}

\item{mi}{Numeric. Number of trials if binomial}
}
\value{
Conditional density of bi
}
\description{
Conditional density for local variables bi
}

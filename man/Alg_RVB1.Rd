% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Alg_RVB1.R
\name{Alg_RVB1}
\alias{Alg_RVB1}
\title{RVB1 Algorithm implementation}
\usage{
Alg_RVB1(y, X, Z, Wprior, etahat, model, m = 1)
}
\arguments{
\item{y}{List. Responses per cluster}

\item{X}{List. Covariates per cluster for fixed effects}

\item{Z}{List. Covariates per cluster for random effects}

\item{Wprior}{List. Wishart prior for random effect covariance}

\item{etahat}{List. Estimate of canonical parameter about which to approximate}

\item{model}{Character. Either "poisson" or "binomial"}

\item{m}{Integer. Number of trials if model is binomial. Keep m = 1 if model is "poisson"}
}
\value{
List containing posterior covariance C, posterior mean mu, the ELBO values per
each 1000 iterations, the run duration, and the final ELBO value.
}
\description{
RVB1 Algorithm implementation
}

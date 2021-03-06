% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc_functions.R
\name{target_tau}
\alias{target_tau}
\title{Conditional Posterior for tau}
\usage{
target_tau(tau, beta, beta_mean, alpha_1, alpha_2, psi, D, g = 10)
}
\arguments{
\item{tau}{current value of precision parameter}

\item{beta}{current value of beta}

\item{beta_mean}{mean of the proposal distribution}

\item{alpha_1}{shape parameter for tau inverse gamma prior.}

\item{alpha_2}{scale parameter for tau inverse gamma prior.}

\item{psi}{current value of psi}

\item{D}{variance covariance matrix of spatial random effects}

\item{g}{value of zellnor's g prior taken as 10}
}
\description{
Conditional Posterior for tau
}

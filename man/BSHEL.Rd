% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BEL_functions.R
\name{BSHEL}
\alias{BSHEL}
\title{BEL spatial function for Moran Basis Prior (POrter et al. 2015)}
\usage{
BSHEL(
  y,
  x,
  n,
  p,
  q,
  var,
  niter,
  beta_init,
  psi_init,
  tau_init,
  M,
  MBM,
  wi,
  sd_psi,
  sd_beta,
  sd_tau
)
}
\arguments{
\item{y}{observed response from data}

\item{x}{design matrix including a column of 1 for intercept and the other from dataset}

\item{n}{number of observations/small areas in the dataset}

\item{p}{number of regression paramters or number of columns in design matrix x}

\item{q}{number of positive eigen values of MBM}

\item{var}{residual variance}

\item{niter}{number of iterations}

\item{beta_init}{initial values for beta}

\item{psi_init}{initial values of spatial random effect psi}

\item{tau_init}{initial value of precision paramter tau}

\item{M}{n × q matrix with the columns being the eigenvectors corresponding to the q largest nonzero eigenvalues of the matrix PBP,P be the latent process space and B be the adjacency matrix}

\item{MBM}{M(B_plus-B)M, where B_plus is a diagonal matrix having number of neighbours for each area in its diagonal}

\item{wi}{current EL weights}

\item{sd_psi}{value of sd for proposal distribution for psi}

\item{sd_beta}{value of sd for proposal distribution for beta}

\item{sd_tau}{value of sd for proposal distribution for tau}
}
\description{
BEL spatial function for Moran Basis Prior (POrter et al. 2015)
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/comp_f_state.R
\name{comp_f_state}
\alias{comp_f_state}
\title{Construct companion matrix for VAR(p) model}
\usage{
comp_f_state(phi, Q, H, con, n, k)
}
\arguments{
\item{phi}{Diagoanl matrix of dimension k x k with vector autoregressive
coefficients.}

\item{Q}{A matrix of .}

\item{H}{A matrix of lambadas from Multi-Move Gibbs sampler in companion
form.}

\item{con}{A scalar, where con = 1 for model with intercept, con = 0 for
model without intercept.}

\item{n}{Number of variables.}

\item{k}{Number of states (number of factors).}
}
\value{
A List with phi, Q, H in companion form.
}
\description{
\code{create_inventory} constructs companion matrix for VAR(p) model.
}
\details{
The function runs inside the Multi-Move Gibbs sampler.
}
\examples{
yt <- as.matrix(t(Xmat))
k <- 2
n <- dim(yt)[1]
lambdasim <- matrix(rep(rnorm(n,0,1)*0.1,k),
nrow = n, ncol = k, byrow = TRUE)
diag(lambdasim) <- 1
lambdasim[upper.tri(lambdasim)] <- 0
lambda <- lambdasim
phi <- diag(rnorm(k,0,1))
Q <- as.matrix(diag(0.1,k))
const <- 0
matcomp <- comp_f_state(phi,Q,lambda,const)
}

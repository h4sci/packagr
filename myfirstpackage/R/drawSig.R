#' Draw variance of idiosyncratic components from inverse gamma distribution
#'
#' function  draws variance of idiosyncratic components from inverse gamma distribution.
#'
#'
#' @param yt A matrix of demeaned and standardized time series data.
#' @param ft A matrix of the factors.
#' @param lambda A vector of dimension n x k of the factor loadings.
#' @param Ttq Number of high-frequency periods minus lag length for state equation.
#' @param nu0 shape (degrees of freedom) of IG distribution.
#' @param s0 Prior for the rate (reciprocal of scale) of IG distribution (1/100=0,01 as benchmark residuals of AR(1) from OLS).
#' @return Matrix R.
#' @import stats
#' @examples
#' q <- 1
#' Tt <- dim(yt)[2]
#' Ttq <- Tt-q
#' V_lam <- diag(k)
#'
#' nu0 <- k + 3
#' s0 <- 100

#' @export
#'
drawSig <- function(yt,lambda,ft,Ttq,nu0,s0){

  nubar <- Ttq/2 + nu0# Posterior shape
  s2bar <- 1/s0 + (yt-lambda%*%ft)%*%t(yt-lambda%*%ft)/2# Posterior scale
  R <- 1/rgamma(1,shape=nubar,rate=s2bar)
  return(R)}


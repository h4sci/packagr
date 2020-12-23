#' Draw factor loadings from their conditional posterior
#'
#' function that allows to draw from posterior of factor loadings.
#'
#' The function conducts a normalization by demeaning and dividing by the
#' standard deviation of the individual time series.
#'
#' @param yt A matrix of demeaned and standardized time series data.
#' @param ft A matrix of the factors.
#' @param R Diagonal matrix of dimension n of idiosyncratic component.
#' @param lam0 Matrix of mean of initial factor loadings.
#' @param V_lam Diagonal matrix of dimension k of Variance-Covariance matrix of
#'   initial factor loadings.
#'
#' @return Element of a vector, coefficient of lambda.
#' @importFrom MASS mvrnorm
#' @examples
#' k <- 2
#' lam0 <- matrix(0,k,1)
#' V_lam <- diag(k)
#' @export
#'
draw_Lam <- function(yt,ft,R,lam0,V_lam){

  # This function is suited for equation-by-equation loop
  D_lam <- solve(ft%*%t(ft)/R + solve(V_lam))
  d_lam <- ft%*%yt/R+ solve(V_lam)%*%lam0

  # Draw lambda from conditional normal distribution
  lambda <- mvrnorm(1, D_lam%*%d_lam, D_lam)

  return(lambda)
}

#' Draws posterior of VAR(p) Model with non-informative prior
#'
#' \code{create_inventory} constructs companion matrix for VAR(p) model.
#'
#' This function draws from the posterior of a VAR(p) model with non-informatibve prior
#'
#' @param Yts n x T matrix of data.
#' @param p the lag lenght of the VAR.
#' @param const A scalar, where const = 1 for model with intercept, const = 0 for model without intercept.
#'
#' @return A List with:
#' \itemize{
#' \item betar_all n x np+1 x R-burn matrices of coefficients.
#' \item Sigr_all n x n x R-burn matrices of variance covariance matrices.
#' }
#' @importFrom MASS mvrnorm
#' @importFrom MCMCpack riwish
#' @examples
#' param <- BVAR_Jeff(ft,q,0)
#' @export
#'
BVAR_Jeff <- function(Yts,p,const){

  # Number of variables
  n <- dim(Yts)[1]

  # Length of time series
  Tt <- dim(Yts)[2]
  Ttp <- Tt-p

  # Data for regression (conditional on first p observations)
  Yt <- t(Yts[,-seq(1,p)])

  # Create X matrix for each t
  if (const==1){
    Xt <-  lapply(seq(p+1,Tt),function(tx){t(c(1,sapply(seq(1,p),function(x){t(Yts[,tx-x])})))})
  }else{
    Xt <-  lapply(seq(p+1,Tt),function(tx){t(c(sapply(seq(1,p),function(x){t(Yts[,tx-x])})))})
  }
  # Stack observations over t
  X <- Reduce(rbind, Xt)


  # Get OLS estimate for beta
  beta_bar <- solve(t(X)%*%X)%*%t(X)%*%Yt

  # Posterior shape and scale
  nubar <- Ttp
  Sbar <- Reduce("+",lapply(seq(1,Ttp),function(tx){t(Yt[tx,]-Xt[[tx]]%*%beta_bar)%*%(Yt[tx,]-Xt[[tx]]%*%beta_bar)})) # is the transpose here wrong as well?

  # Draw Sigma from inverse wishart distribution
  Sigr <- riwish(nubar, Sbar)


  # Draw beta from conditional normal distribution
  betar <- mvrnorm(1, c(beta_bar), Sigr%x%(solve(t(X)%*%X)) )

  # Bring betas back in shape
  betar_m <- t(matrix(betar,n*p+const,n))


  output <- list(FF=betar_m, Q= Sigr)

  return(output)

}

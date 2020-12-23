#' Draws factor conditional on parameters
#'
#' The Multi-Move Gibbs sampler applies the kalman filter by forward and
#' backwards filtering.
#'
#' @inheritParams compFstate
#' @param phi Diagoanl matrix of dimension k x k with vector autoregressive
#'   coefficients.
#' @param lambda A vector of dimension n x k of the factor loadings.
#' @param const A scalar, where const = 1 for model with intercept, const = 0
#'   for model without intercept.
#' @param yt A matrix of demeaned and standardized time series data.
#' @param Tt Number of high-frequency periods.
#' @param q lag length for state equation (adjust starting value of phi
#'   accordingly).
#' @param alpha_0 Vector of dimension m x 1 (Initial conditions for Kalman
#'   filter).
#' @param P_0 Diagonal matrix of dimension m (Initial conditions for Kalman
#'   filter).
#' @param R Diagonal matrix of dimension n of idiosyncratic component.
#' @importFrom MASS mvrnorm
#' @return A Vector of factors conditional parameters.
#' @examples
#' q <- 1
#' yt <- as.matrix(t(Xmat))
#' n <- dim(yt)[1]
#' Tt <- dim(yt)[2]
#' k <- 2
#' m <- k*q
#' alpha_0 <- matrix(0,m,1)
#' P_0 <- diag(m)
#' const <- 0
#' R <- as.matrix(diag(n)*0.01)
#' phi <- diag(rnorm(k,0,1))
#' lambdasim <- matrix(rep(rnorm(n,0,1)*0.1,k), nrow = n, ncol = k, byrow = TRUE)
#' diag(lambdasim) <- 1
#' lambdasim[upper.tri(lambdasim)] <- 0
#' lambda <- lambdasim
#' Q <- as.matrix(diag(0.1,k))
#' matcomp <- compFstate(phi,Q,lambda,const)
#' @export
#'
multimoveGibbs <- function(yt,phi,Q,lambda,const,Tt,q,alpha_0,P_0,R){

  k <- dim(Q)[1]

  # Pre-allocate space within multi-move sampler
  alpha_s <-vector("list",Tt)
  alpha_u <- vector("list", Tt+1)
  P_u <-  vector("list", Tt+1)
  alpha_f <-  vector("list", Tt)
  P_f <-  vector("list", Tt)

  # Initial conditions
  alpha_u[[1]] <- alpha_0
  P_u[[1]] <- P_0

  # Get matrices in companion form
  matcomp <- compFstate(phi,Q,lambda,const)
  FFcom <- matcomp$phicom
  Qcom <- matcomp$Qcom
  Hcom <- matcomp$Hcom
  FF <- phi
  H <- lambda

  for(t in seq(1,Tt)){

    ###### Forward filter iterations
    #Prediction steps
    alpha_f[[t]] <- FFcom%*%alpha_u[[t]] # predicting the state
    P_f[[t]] <- FFcom%*%P_u[[t]]%*%t(FFcom) + Qcom #variance-covariance matrix of prediction
    nut <- yt[,t] - Hcom%*%alpha_f[[t]] # predition error
    nuvart <- Hcom%*%P_f[[t]]%*%t(Hcom) + R # variance-covariance of predicion error # example 5) + diag(R)

    # Updating steps
    Kt <- P_f[[t]]%*%t(Hcom)%*%solve(nuvart)#Kalman gain
    alpha_u[[t+1]] <- alpha_f[[t]] + Kt%*%nut # Update alpha_{t|t-1}
    P_u[[t+1]] <- P_f[[t]] -  Kt%*%Hcom%*%P_f[[t]] #Update P_{t|t-1}
  }
  ###### Backward sampling
  # First step is to sample from N(alpha_{T|T},P{T,T}) using the last step of Kalman filter
  alpha_s[[Tt]] <- mvrnorm(1,alpha_u[[Tt+1]], P_u[[Tt+1]])

  for(t in seq(Tt-1,1)){
    # Because Qcom is singular use only those equations that do not corresponds to the identities of the companion form
    Ktt <- P_u[[t+1]]%*%t(FF)%*%solve(FF%*%P_u[[t+1]]%*%t(FF) + Q)
    alpha_stt <- alpha_u[[t+1]] + Ktt%*%(alpha_s[[t+1]][1:k]- FF%*%alpha_u[[t+1]]) #Update alpha_{t|t}
    P_st <- P_u[[t+1]] - Ktt%*%FF%*%P_u[[t+1]] #Update P_{t|t}

    # Draw alpha for T-1,...,1
    alpha_s[[t]]  <- mvrnorm(1,alpha_stt, P_st)

  }

  # Select those elements from orignial specification
  sele <- cbind(diag(k),matrix(0,k,k*(q-1)))
  alphat <-  sele%*%Reduce(cbind,alpha_s)

  return(alphat)
}

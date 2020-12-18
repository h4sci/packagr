#' Get nowcast from the bayesian DFM
#'
#' function extends the data matrix by the forecasting period, draws extended factor conditional on posterior parameters and makes forecast.
#'
#'
#' @param yt A matrix of demeaned and standardized time series data.
#' @param k Number of states (number of factors)
#' @param q Lag length for state equation.
#' @param m number of states (number of factors) x lag length for state equation (nxq)
#' @param n Number of variables
#' @param Tt Number of high-frequency periods.
#' @param Ttq Number of high-frequency periods minus lag length for state equation
#' @param const A scalar, where const = 1 for model with intercept, const = 0 for model without intercept.
#' @param target Defines the target variable.
#'
#' @return A list with the mean and quantiles of the parameters, the forecast and the number of burns, draws and saves.
#' @importFrom stats rnorm
#' @examples
#'
#' out <- run_model(yt,k,q,m,n,Tt,Ttq,const,target)
#' @export
#'
run_model <- function(yt,k,q,m,n,Tt,Ttq,const,target){

###### PRIORS

# IG distribution
nu0 <- k + 3 # shape (degrees of freedom) of IG distribution
s0 <- 100 # Prior for the rate (reciprocal of scale) of IG distribution (1/100=0,01 as benchmark residuals of AR(1) from OLS)

# Factor loadings
lam0 <- matrix(0,k,1) # mean = 0
V_lam <- diag(k)

# Initial conditions for Kalman filter
alpha_0 <- matrix(0,m,1)
P_0 <- diag(m)

# Initialize Lambda
lambdasim <- matrix(rep(rnorm(n,0,1)*0.1,k), nrow = n, ncol = k, byrow = TRUE) # Must be a nxk Vector! # Could replace this with principal components???
diag(lambdasim) <- 1
lambdasim[upper.tri(lambdasim)] <- 0
lambda <- lambdasim

# VAR coefficient
#phi <- t(c(0.4)) # if q=1, k=1
#phi <-  t(c(0.4,0.2)) # if q=2, k=1 etc.
phi <- diag(rnorm(k,0,1)) # eg. matrix(c(0.2,0,0,0.2),2,2) # vector autoregressive coefficients       # if q=1, k=2
#betas2 <- matrix(c(0.2,0.01,0.01,0.05),2,2)# vector autoregressive coefficient       # phi2 usw.

# Variance-Covariance Matrix Error Terms
Q <- as.matrix(diag(0.1,k)) # How to choose this? # Assume indepenence of factors
R <- as.matrix(diag(n)*0.01) # Assume no serial correlation in the idiosyncratic component


####### RUN ESTIMATION


# Gibbs sampling preliminaries
ndraws <- 500 #5000
burn <- 100 #1000
nsave <- ndraws - burn


# Pre-allocate space
ft_all <- vector("list",nsave)#replicate((R-burn)/thin,matrix(0,n,n*p+1))
phi_all <- vector("list",nsave)
lambda_all <- vector("list", nsave)
R_all <- vector("list", nsave)
Q_all <- vector("list", nsave)


gx <- 1

for(r in seq(1,ndraws)){
  if(r%%100==0){print(r)}

  # Draw factor conditional on parameters
  ft <-  multimoveGibbs(yt,phi,Q,lambda,const,Tt,q,alpha_0,P_0,R)

  # Draw (V)AR parameters of factor equation
  param <- BVAR_Jeff(ft,q,0)
  phi <- param$FF
  Q <- param$Q

  # Draw R from inverse gamma
  for(ix in seq(1,n)){
    # Given that errors independent we can draw them equation-by-equation
    R[ix,ix] <-  drawSig(yt[ix,],lambda[ix,],ft,Ttq,nu0,s0)

    lambda[ix,] <- drawLam(yt[ix,],ft,R[ix,ix],lam0,V_lam)
  }

  # Set the first value of lambda to 1 (this is not really efficient, but within our context okay to do)
  # Also watch out for additional rotational indeterminacy in case you want to increase number of factors to k>1!
  diag(lambda) <- 1
  lambda[upper.tri(lambda)] <- 0

  # # Also try with first quadrat as Identity
  # lambda_head <- lambda[1:dim(lambda)[2],1:dim(lambda)[2]]
  # lambda_head[lower.tri(lambda[1:dim(lambda)[2],1:dim(lambda)[2]])] <- 0
  # lambda[1:dim(lambda)[2],1:dim(lambda)[2]] <- lambda_head

  # Save draws after burn in
  if(r>burn){

    lambda_all[[gx]] <- lambda
    phi_all[[gx]] <- phi
    Q_all[[gx]] <- Q
    R_all[[gx]] <- R
    ft_all[[gx]] <- ft

    gx <- gx + 1}
}


#save results in list
out <- list()
out$ft_all <- ft_all
out$phi_all <- phi_all
out$lambda_all <- lambda_all
out$R_all <- R_all
out$Q_all <- Q_all


###Compute posterior mean for the parameters
out$lambda_mean <- apply(simplify2array(out$lambda_all), 1:2, mean)
out$phi_mean <- apply(simplify2array(out$phi_all), 1:2, mean)
#out$phi_mean <- mean(simplify2array(out$phi_all))
out$Q_mean <- apply(simplify2array(out$Q_all), 1:2, mean)
#out$Q_mean <- mean(simplify2array(out$Q_all))
out$R_mean <- apply(simplify2array(out$R_all), 1:2, mean)

# Compute posterior mean for the factor
out$ft_mean <- apply(simplify2array(out$ft_all), 1:2, mean)


# Compute quantiles for the factor
out$ft_Q <- apply(simplify2array(ft_all), 1:2, quantile, prob = c(0.05,0.5, 0.95))

out$lambda_Q <- apply(simplify2array(lambda_all), 1:2, quantile, prob = c(0.05,0.5, 0.95))

out$lambda_05 <- out$lambda_Q[1, , ]  ## 0.05 quantile
out$lambda_50 <- out$lambda_Q[2, , ]  ## 0.5 quantile (median)
out$lambda_95 <- out$lambda_Q[3, , ]  ## 0.95 quantile


####### NOWCAST

# define parameters
freq_diff <- max(inventory$freq)/min(inventory$freq) # Fraction of high-frequency periods in lowest frequency
s <- 2*(k - 1) # Number of periods for aggregation rule in forumla (3)


# Make prediction
out$ncst <- get_nowcast(Xmat = Xmat, # n, k?
                        s = s,
                        q = q,
                        alpha_0 = alpha_0,
                        P_0 = P_0,
                        inventory = inventory,
                        target = target,
                        flows = flows,
                        lambda_mean = out$lambda_mean,
                        phi_mean = out$phi_mean,
                        R_mean = out$R_mean,
                        Q_mean = out$Q_mean,
                        const = const)

out$nsave <- nsave
out$ndraws <- ndraws
out$burn <- burn

return(out)

}


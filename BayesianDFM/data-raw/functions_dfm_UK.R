

# Create Inventory --------------------------------------------------------

create_inventory <- function(flows, stocks) {

  # construct inventory of time series
  inventory <- rbind(data.frame("key" = as.character(names(flows)),
                                "type" = factor("flow", levels = c("stock","flow")),
                                "freq" = sapply(flows, frequency),
                                "mean" = sapply(flows, mean, na.rm=T),
                                "sd" = sapply(flows, sd, na.rm=T),
                                stringsAsFactors = F,
                                row.names = NULL),
                     data.frame("key" = names(stocks),
                                "type" = factor("stock", levels = c("stock","flow")),
                                "freq" = sapply(stocks, frequency),
                                "mean" = sapply(stocks, mean, na.rm=T),
                                "sd" = sapply(stocks, sd, na.rm=T),
                                stringsAsFactors = F,
                                row.names = NULL))

  # remove NULL entires
  if(length(which(inventory$key == "") > 0)) inventory <- inventory[-which(inventory$key == ""),]

  return(inventory)

}


# Prepare data ------------------------------------------------------------

prepare_data <- function(flows, stocks, inventory, target) {

  data <- c(flows, stocks)

  # standardize data
  data_std <- lapply(inventory$key, function(ix){

    (data[[ix]] - inventory[which(inventory$key == ix),"mean"])/
      inventory[which(inventory$key == ix),"sd"]

  }); names(data_std) <- inventory$key

  # adjust time periods
  freq_max <- max(sapply(data, frequency))

  data_adj <- lapply(inventory$key, function(ix) {

    tim <- time(data_std[[ix]]) + (freq_max/frequency(data_std[[ix]]) - 1)/freq_max
    out <- as.data.frame(cbind(round(tim,5), data_std[[ix]]))
    colnames(out) <- c("time",ix)

    return(out)

  }); names(data_adj) <- inventory$key

  # coerce to sparse matrix
  idx = data.frame("time" = round(seq(1900, 2100, 1/freq_max),5))

  dfs <- lapply(data_adj, function(fx){

    tab <- left_join(idx, fx, by = "time")
    tab[,-1]

  })

  dfs <- do.call(cbind,dfs)
  colnames(dfs) <- inventory$key
  dfts <- ts(dfs, start = 1900, frequency = max(sapply(data, function(x) frequency(x))))
  dfts <- na.trim(dfts, is.na = "all")

  dfts[which(is.na(dfts))] <- 0

  return(dfts)

}


# Information Criterion --------------------------------------------------------

get_ic <- function(Xmat) {
  IC <- list()
  #  See Ahn, Horenstein (2013) - Eigenvalue Ratio Test for the Number of Factors, page 1207
  #Xmat_spline <- na.approx(Xmat, na.rm = TRUE)
  Xmat_spline <- na.spline(Xmat)

  pc <- princomp(Xmat_spline, cor=TRUE)
  summary(pc, loadings=TRUE)
  pdf("text/pc.pdf")
  screeplot(pc, type="lines", main = "Principal Components")
  dev.off()
  kmax <-  10
  kmax = min(kmax,length(pc$sdev)); # Set kmax to ensure a maximum number of possible factors

  # Eigenvalue ratio test
  ER <- pc$sdev[1:kmax-1]/pc$sdev[2:kmax]
  IC$ER_number <- which.max(ER)

  # Growth ratio test
  m <- min(dim(Xmat_spline),length(pc$sdev)) # define m for calculation of V(k)
  V_k <- rev(cumsum(rev(pc$sdev[m:1])))
  GR = log(V_k[1:(kmax-2)]/V_k[2:(kmax-1)])/log(V_k[2:(kmax-1)]/V_k[3:kmax])
  IC$GR_number <- which.max(GR)

  return(IC)
}


# companion form of our VAR(p) model --------------------------------------

# function that computes the companion form of our VAR(p) model for the state space model

comp_f_state <- function(phi,Q,H,con) {
  # This function constructs companion matrix for VAR(p) model
  # input: - n x n*p + 1 matrix of coefficients betam
  #        -  scalar const, where const=1 model with intercept, const=0 model without intercept

    dim1 <- dim(phi)[1]
  p <- dim(phi)[2]/dim1

  if(con==1){
    conm <- phi[,1]
    A <- phi[,-1]}else{A <- phi}
  if(p>1){
    dimdiff <- dim(A)[2]-dim(A)[1]
    Acom <- rbind(A,cbind(diag(dimdiff),matrix(0,dimdiff,dim(A)[1])))
    Hcom <- cbind(H,matrix(0,n,k*(q-1)))
    Qcom <- cbind(rbind(Q,matrix(0,dim1*(p-1),dim1)), matrix(0,dim1*p,dim1*(p-1)))
  }else{
    Acom <- A
    Hcom <- H
    Qcom <- Q}


  output <- list(phicom=Acom,Qcom = Qcom,Hcom = Hcom)
  return(output)
}


# Multi-Move Gibbs sampler ------------------------------------------------

# function that implements the Gibbs sampler outlined in the lecture(forward filtering and backward sampling)
multimove_gibbs <- function(yt,phi,Q,lambda,const,Tt,q,alpha_0,P_0,R) {

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
  matcomp <- comp_g_state(phi,Q,lambda,const)
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


# VAR(p) Model with non-informative prior ---------------------------------

bvar_jeff <- function(Yts,p,const) {
  # This function draws from the posterior of a VAR(p) model with non-informative prior

  # input: - Yts n x T matrix of data
  #         - p the lag lenght of the VAR
  #         - con for choosing intercept (0/1)

  # ouput: - betar_all n x np+1 x R-burn matrices of coefficients
  #        -   Sigr_all n x n x R-burn matrices of variance covariance matrices
  # Start timer

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
  # Stack obeservations over t
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
#why not apply this?:
#For nowcasting purposes and especially in high-frequency time series, it is
#common to assume some persistence in the factors. It is therefore desirable to shrink the
#autoregressive parameters towards a random walk setup. As a result, the model uses a
#normal prior with mean a0 = 1, 0, . . . , 0 and covariance matrix A0 = Ip/t. This shrinks
#the process of the factors towards a random walk.


# Function to draw factor loadings ----------------------------------------

# Write function that allows to draw from posterior of factor loadings
draw_lam <- function(yt,ft,R,lam0,V_lam) {
  # Draw factor loadings from their conditional posterior
  #(see Exercise 2 and see also lecture slide "Gibbs Sampling in a Regression Model")
  # This function is suited for equation-by-equation loop
  D_lam <- solve(ft%*%t(ft)/R + solve(V_lam))
  d_lam <- ft%*%yt/R+ solve(V_lam)%*%lam0

  # Draw lambda from conditional normal distribution
  lambda <- mvrnorm(1, D_lam%*%d_lam, D_lam)

  return(lambda)
}



# Function to draw variances of idiosyncratic components ------------------

draw_sig <- function(yt,lambda,ft,Ttq,nu0,s0) {
  # Draw variance of idiosyncratic components from inverse gamma distribution
  nubar <- Ttq/2 + nu0# Posterior shape
  s2bar <- 1/s0 + (yt-lambda%*%ft)%*%t(yt-lambda%*%ft)/2# Posterior scale
    R <- 1/rgamma(1,shape=nubar,rate=s2bar)
  return(R)}

# It is assumed that the dynamic factor accounts for a majority of the cross-correlation in
#the data. Sigma is therefore assumed to be diagonal and we can draw the error variances sigma_i
#equation-by-equation from an inverse Gamma distribution.

##################### Function for drawing from multivariate normal without worrying too much (but much slower thab mvrnorm) ######################################
# When is this better than mvrnorm?
mvrnC <- function(mu,Sigma) {
  m <- dim(Sigma)[1]
  C <- chol(nearPD(Sigma)$mat)
  #create a matrix of random standard normals
  Z <- matrix(rnorm(m), m)
  #multiply the standard normals by the transpose of the Cholesky and add mean
  X <- as.matrix(mu + t(C) %*% Z)

  return(X)
}


# Function to get nowcast -------------------------------------------------

get_nowcast <- function(Xmat, s, q, alpha_0, P_0, inventory,
                        target, flows, lambda_mean, phi_mean, R_mean, Q_mean, const) {

  # determine nowcast date
  #gdp_ncst <- as.numeric(tail(time(flows[[target]]),1)) + (s+1)/frequency(Xmat)
  gdp_ncst <- as.numeric(tail(time(flows[[target]]),1) + 3) + (s+1)/frequency(Xmat)

  # extend data with zeros and stochastic volatility with random walk if necessary
  if(tail(time(Xmat),1)[1] < gdp_ncst){

    # data
    Xmat_ext <- window(Xmat,
                       end = gdp_ncst,
                       extend = T)
    Xmat_ext[which(is.na(Xmat_ext))] <- 0

    # # stochastic volatility
    # h_ext <- window(ts(h_out,
    #                    start = time(Xmat)[1] - s/frequency(Xmat),
    #                    frequency = frequency(Xmat)),
    #                 end = gdp_ncst,
    #                 extend = T)
    # h_ext[is.na(h_ext)] <- tail(h_out,1)

  } else {

    Xmat_ext <- Xmat
    # h_ext <- ts(h_out,
    #             start = time(Xmat)[1]-s/frequency(Xmat),
    #             frequency = frequency(Xmat))

  }
  yt_ext <- as.matrix(t(Xmat_ext))
  Tt_ext <- dim(yt_ext)[2]

  # Draw extended factor conditional on posterior parameters
  ft_nc <-  multimove_gibbs(yt_ext,phi_mean,Q_mean,lambda_mean,
                            const,Tt_ext,q,alpha_0,P_0,R_mean)

  # Make prediction
  X_fit <- lambda_mean %*% ft_nc
  target_fit <- X_fit[which(inventory$key==target),]
  target_fit <- target_fit[which(target_fit!=0)]

  # rescale
  target_rescaled <- (target_fit * inventory[which(inventory$key == target), "sd"]) +
    inventory[which(inventory$key == target), "mean"]

  out <- ts(target_rescaled,
            start = time(Xmat[,1])[1],
            frequency = 12) #frequency(flows[[target]]))

  return(out)

}


# Run Model ---------------------------------------------------------------

run_model <- function(yt,k,q,m,n,Tt,Ttq,const,target) {


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
lambdasim <- matrix(rep(rnorm(n,0,1)*0.1,k), nrow = n, ncol = k,
                    byrow = TRUE) # Must be a nxk Vector!
# Could replace this with principal components???
diag(lambdasim) <- 1
lambdasim[upper.tri(lambdasim)] <- 0
lambda <- lambdasim

# VAR coefficient
#phi <- t(c(0.4)) # if q=1, k=1
#phi <-  t(c(0.4,0.2)) # if q=2, k=1 etc.
phi <- diag(rnorm(k,0,1)) # eg. matrix(c(0.2,0,0,0.2),2,2) # vector autoregressive coefficients
# if q=1, k=2
#betas2 <- matrix(c(0.2,0.01,0.01,0.05),2,2)# vector autoregressive coefficient
# phi2 usw.

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
  ft <-  multimove_gibbs(yt,phi,Q,lambda,const,Tt,q,alpha_0,P_0,R)

  # Draw (V)AR parameters of factor equation
  param <- BVAR_Jeff(ft,q,0)
  phi <- param$FF
  Q <- param$Q

  # Draw R from inverse gamma
  for(ix in seq(1,n)){
    # Given that errors independent we can draw them equation-by-equation
    R[ix,ix] <-  draw_sig(yt[ix,],lambda[ix,],ft,Ttq,nu0,s0)

    lambda[ix,] <- draw_lam(yt[ix,],ft,R[ix,ix],lam0,V_lam)
  }

  # Set the first value of lambda to 1
  # (this is not really efficient, but within our context okay to do)
  # Also watch out for additional rotational indeterminacy
  # in case you want to increase number of factors to k>1!
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
freq_diff <- max(inventory$freq)/min(inventory$freq)
# Fraction of high-frequency periods in lowest frequency
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


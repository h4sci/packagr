#' Get nowcast from the bayesian DFM
#'
#' function  extends the data matrix by the forecasting period, draw extended factor conditional on posterior parameters and makes forecast.
#'
#'
#' @param Xmat Demeaned and standardized matrix of time series data.
#' @param s Number of periods for aggregation rule.
#' @param q Lag length for state equation.
#' @param alpha_0 Vector of dimension m x 1 (Initial conditions for Kalman filter).
#' @param P_0 Diagonal matrix of dimension m (Initial conditions for Kalman filter).
#' @param inventory Inventory of the corresponding data.
#' @param target Defines the target variable.
#' @param flows A list of time series variable (flow variables).
#' @param lambda_mean A vector of dimension n x k of the posterior mean factor loadings.
#' @param phi_mean Diagonal matrix of dimension k x k with vector of posterior mean autoregressive coefficients.
#' @param R_mean Diagonal matrix of dimension n of posterior mean idiosyncratic component.
#' @param Q_mean matrices of posterior mean variance covariance matrices.
#' @param const A scalar, where const = 1 for model with intercept, const = 0 for model without intercept.
#'
#' @return A time series of the target variable including forecast.
#' @import stats
#' @importFrom utils tail
#' @examples
#' yt <- as.matrix(t(Xmat))
#' k <- 2
#' n <- dim(yt)[1]
#' q <- 1
#' const <- 0
#' m <- k*q
#' s <- 2*(k - 1)
#' alpha_0 <- matrix(0,m,1)
#' P_0 <- diag(m)
#' inventory <- create_inventory(flows = data$flows, stocks = data$stocks)
#' target <- c("UKGDPM.YQ")
#' flows = data$flows
#' out$lambda_mean <- apply(simplify2array(out$lambda_all), 1:2, mean)
#' out$phi_mean <- apply(simplify2array(out$phi_all), 1:2, mean)
#' out$Q_mean <- apply(simplify2array(out$Q_all), 1:2, mean)
#' out$R_mean <- apply(simplify2array(out$R_all), 1:2, mean)
#' out$ncst <- get_nowcast(Xmat = Xmat, s = s, q = q, alpha_0 = alpha_0,
#' P_0 = P_0, inventory = inventory, target = target, flows = flows,
#' lambda_mean = out$lambda_mean, phi_mean = out$phi_mean, R_mean = out$R_mean,
#' Q_mean = out$Q_mean, const = const)
#' @export
#'
get_nowcast <- function(Xmat, s, q, alpha_0, P_0, inventory, target, flows, lambda_mean, phi_mean, R_mean, Q_mean, const){

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
  ft_nc <-  multimoveGibbs(yt_ext,phi_mean,Q_mean,lambda_mean,const,Tt_ext,q,alpha_0,P_0,R_mean)

  # Make prediction
  X_fit <- lambda_mean %*% ft_nc
  target_fit <- X_fit[which(inventory$key==target),]
  target_fit <- target_fit[which(target_fit!=0)]

  # rescale
  target_rescaled <- (target_fit * inventory[which(inventory$key == target), "sd"]) + inventory[which(inventory$key == target), "mean"]

  out <- ts(target_rescaled,
            start = time(Xmat[,1])[1],
            frequency = 12) #frequency(flows[[target]]))

  return(out)

}

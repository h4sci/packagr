
cat("\014")
rm(list = ls())

# NOTES -------------------------------------------------------------------
# Run this file to read in the data and to run the estimation.
# The data is first transformed by read_dta_UK.R
# To get the plots, run plots_UK.R
# All the functions are in functions_dfm_UK.R


#Estimate a dynamic factor model
#Model: y(i,t) = lambda(i) f(t) + epsilon(i,t), with epsilon(i,t) ~ N(0,R)
#f(t) = phi_1 f(t-1) + ... + phi_q f(t-q) + nu(t), with nu(t) ~ N(0,Q)


# PACKAGES AND FUNCTIONS --------------------------------------------------

# packages
library(MASS) # For drawing from a multvariate normal distribution
library(ggplot2) # For plotting nice histogramms
library(openxlsx) # To open xlsx files
library(reshape2) # For reshaping the dataset
library(MCMCpack) # For drawing from inverse Wishart
library(Matrix) # For efficient and convenient matrix manipulations
library(dplyr) # For join in prepare_data function
library(zoo) # For na.trim in prepare_data function
library(tidyr)

source("code/functions_dfm_UK.R")

# IMPORT DATA -------------------------------------------------------------

# load data table created in read_data_UK.R
load("code/out/data_UK.Rda")



# DEFINE TARGET VARIABLE --------------------------------------------------

target <- c("UKGDPM.YQ")


# PREPARE DATA ------------------------------------------------------------

# Create an inventory of the variables and their main summary statistics
inventory <- create_inventory(flows = data$flows, stocks = data$stocks)
save(inventory, file = "code/out/inventory_UK.Rda")

# De-mean & standardize data and bring variables with different frequency in one matrix
Xmat <- prepare_data(flows = data$flows,
                     stocks = data$stocks,
                     inventory,
                     target = target)
save(Xmat, file = "code/out/Xmat_UK.Rda")


yt <- as.matrix(t(Xmat))

# plot the time series as a check for the user
plots <- T
flows = data$flows
stocks = data$flows

if(plots == T){
  pdf("text/var_hist.pdf")
  tsl <- c(stocks,flows)
  par(mfrow = c(length(unique(inventory$freq)), 1))
  for(x in unique(inventory$freq)){
    plot.ts(scale(do.call(cbind,tsl[inventory$key[inventory$freq == x]])),
            xlab = NULL,
            ylab = paste("frequency: ",x),
            ylim = c(-15,15),
            plot.type="single")
  }
  par(mfrow = c(1,1))

}
dev.off()

# MODEL SPECIFICATION -----------------------------------------------------

# Check information criteria
IC <- get_IC(Xmat)

#factors <- list()
#for (k in 2:10){
k <- 2 # number of states (number of factors)

q <- 1 # lag length for state equation (adjust starting value of phi accordingly)
m <- k*q

n <- dim(yt)[1] # Number of variables
Tt <- dim(yt)[2] # Number of high-frequency periods
Ttq <- Tt-q
const <- 0 # choose constant in the state equation (we choose no constant)

# Set seed for reproducibility
set.seed(111)



# RUN ESTIMATION ----------------------------------------------------------

out <- run_model(yt,k,q,m,n,Tt,Ttq,const,target,inventory)


# Evaluation --------------------------------------------------------------

# get nowcast
#out$X_pred <- ((1+tail(out$ncst,1)/100)^4-1)*1e+4
out$X_pred <- tail(out$ncst,1)

# Test
# Run the Estimation including calculation of prediction interval
#output <- get_interval(yt,k,q,m,const,target, Xmat, inventory, flows)

## Is correlation matrix represented well using applied number of factors?
est <- out$lambda_mean %*% t(out$lambda_mean) + out$R_mean
round(est, 3)
dat <- round(cor(t(yt)), 3)
dif <- round(dat - est, 3)
dif

## Check whether factors are approx. mean = 0
for (i in 1:dim(out$ft_mean)[1]){
  print(mean(out$ft_mean[i,]))
}

## Check standard deviation of factors
for (i in 1:dim(out$ft_mean)[1]){
  print(sd(out$ft_mean[i,]))
}

# Save results ------------------------------------------------------------

#save(out, output, k, m, q, Tt, Ttq, target, file = "code/out/results_UK.Rda")
save(out, k, m, q, Tt, Ttq, target, file = "code/out/results_UK.Rda")



<!-- README.md is generated from README.Rmd. Please edit that file -->

# BayesianDFM

<!-- badges: start -->

<!-- badges: end -->

The myfirstpackage package provides relevant functions for running a
Bayesian Dynamic Factor Model.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r

devtools::install_github("h4sci/packagr")
```

## Example

This is a basic example which shows you how to solve a common problem:

Load and prepare the input data

``` r
load("data_UK.Rda")

target <- c("UKGDPM.YQ") # Define target variable

# De-mean & standardize data and bring variables with different frequency in one matrix
Xmat <- prepare_data(flows = data$flows,
                     stocks = data$stocks,
                     inventory,
                     target = target)

yt <- as.matrix(t(Xmat))
```

Make a in-sample evaluation of the optimal number of factors

``` r
IC <- get_IC(Xmat) # Check information criteria
```

Define the parameters for the model

``` r
k <- 2 # number of states (number of factors)

q <- 1 # lag length for state equation (adjust starting value of phi accordingly)
m <- k*q

n <- dim(yt)[1] # Number of variables
Tt <- dim(yt)[2] # Number of high-frequency periods
Ttq <- Tt-q
const <- 0 # choose constant in the state equation (we choose no constant)
```

Having all the parameters, run the model

``` r
out <- run_model(yt,k,q,m,n,Tt,Ttq,const,target)
```

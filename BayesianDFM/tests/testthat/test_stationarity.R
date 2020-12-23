context("Stationarity Test")
library(BayesianDFM)

rowVars <- function(x, na.rm=F) {
  # Vectorised version of variance filter
  rowSums((x - rowMeans(x, na.rm=na.rm))^2, na.rm=na.rm) / (ncol(x) - 1)
}

yt <- as.matrix(t(Xmat))

rowSd <- sqrt(rowVars(yt))

rowM <- rowMeans(yt)

test_that("indicators are already demeaned and normalized such that stationarity is given", {
  expect_equal(unname((rowM)), rep(0,length(rowM)), tolerance=1e-10)
})


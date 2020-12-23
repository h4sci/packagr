context("Stationarity Test")
library(BayesianDFM)

row_vars <- function(x, na.rm=F) {
  # Vectorised version of variance filter
  rowSums((x - rowMeans(x, na.rm=na.rm))^2, na.rm=na.rm) / (ncol(x) - 1)
}

yt <- as.matrix(t(Xmat))

row_sd <- sqrt(row_vars(yt))

row_m <- rowMeans(yt)

test_that("indicators are already demeaned and normalized such that stationarity is given", {
  expect_equal(unname((row_m)), rep(0,length(row_m)), tolerance=1e-10)
})


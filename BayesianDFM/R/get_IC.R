#' Find optimal number of factors in a factor model.
#'
#' \code{create_inventory} constructs the Eigenvalue Ratio Test and Growth Ratio
#' Test according to Ahn, Horenstein (2013), page 1207
#'
#' The function calculates the principal components of the data set. It creates
#' a screeplot of the principal components. Further, it calculates the optimal
#' number of factors of the factor model according to the Eigenvalue Ratio Test
#' and Growth Ratio Test of Ahn, Horenstein (2013), page 1207. The maximam
#' number of factors is set to 10.
#'
#' @param Xmat A data set containing all the time series.
#'
#' @return A List: \itemize{ \item ER_number Optimal number of factors according
#'   to Eigenvalue Ratio Test. \item GR_number Optimal number of factors
#'   according to Growth Ratio Test. }
#' @importFrom zoo na.spline
#' @import stats
#' @importFrom grDevices dev.off
#' @examples
#' IC <- get_IC(Xmat)
#' @export
#'
get_IC <- function(Xmat){

  IC <- list()
  #  See Ahn, Horenstein (2013) - Eigenvalue Ratio Test for the Number of Factors, page 1207
  #Xmat_spline <- na.approx(Xmat, na.rm = TRUE)
  Xmat_spline <- na.spline(Xmat)

  pc <- princomp(Xmat_spline, cor=TRUE)
  summary(pc, loadings=TRUE)
  #screeplot(pc, type="lines", main = "Principal Components")
  #dev.off()
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


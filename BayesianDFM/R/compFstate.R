#' Construct companion matrix for VAR(p) model
#'
#' \code{create_inventory} constructs companion matrix for VAR(p) model.
#'
#' The function runs inside the Multi-Move Gibbs sampler.
#'
#' @param phi Diagoanl matrix of dimension k x k with vector autoregressive coefficients.
#' @param Q A matrix of .
#' @param H A matrix of lambadas from Multi-Move Gibbs sampler in companion form.
#' @param con A scalar, where con = 1 for model with intercept, con = 0 for model without intercept.
#' @param n Number of variables.
#' @param k Number of states (number of factors).
#'
#' @return A List with phi, Q, H in companion form.
#' @examples
#' yt <- as.matrix(t(Xmat))
#' k <- 2
#' n <- dim(yt)[1]
#' lambdasim <- matrix(rep(rnorm(n,0,1)*0.1,k), nrow = n, ncol = k, byrow = TRUE)
#' diag(lambdasim) <- 1
#' lambdasim[upper.tri(lambdasim)] <- 0
#' lambda <- lambdasim
#' phi <- diag(rnorm(k,0,1))
#' Q <- as.matrix(diag(0.1,k))
#' const <- 0
#' matcomp <- compFstate(phi,Q,lambda,const)
#' @export
#'
compFstate <- function(phi,Q,H,con,n,k){

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

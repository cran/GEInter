#' Make predictions from a "bic.PTReg" object
#'
#'
#' This function makes predictions from a BIC PTReg model, using the stored \code{"bic.PTReg"}
#' object.  This function makes it easier to use the results of BIC to make
#' a prediction.
#' @param object Fitted \code{"bic.PTReg"} object.
#' @param newE Matrix of new values for \code{E} at which predictions are to be
#' made.
#' @param newG Matrix of new values for \code{G} at which predictions are to be
#' made.
#' @param \dots Not used. Other arguments to predict.
#'
#' @return The object returned depends on the \dots{} argument which is passed
#' on to the \code{predict} method for \code{PTReg} objects.
#' @seealso  \code{bic.PTReg}, and \code{coef}, and \code{plot} methods, and \code{PTReg}.
#' @references Yaqing Xu, Mengyun Wu, Shuangge Ma, and Syed Ejaz Ahmed.
#' \emph{Robust gene-environment interaction analysis using penalized trimmed regression. Journal
#' of Statistical Computation and Simulation, 88(18):3502-3528, 2018.}
#' @method predict bic.PTReg
#'
#' @export
#' @export predict.bic.PTReg
predict.bic.PTReg <- function(object,newE,newG, ...){
  if((missing(newE))|(missing(newG))){
    stop("You need to supply a value for 'newE' and 'newG'")
  }
  newE=as.matrix(newE)
  newG=as.matrix(newG)
  temp=coef.bic.PTReg(object)
  beta_estimate=temp$beta
  p=dim(beta_estimate)[2]
  q=dim(beta_estimate)[1]-1
  alpha_estimate=temp$alpha
  intercept_estimate=temp$intercept
  b_estimate=matrix(beta_estimate,p*(q+1),1)
  n=dim(newE)[1]
  W=matrix(0,n,p*(q+1))
  W[,seq(from=1,to=p*(q+1),by=(q+1))]=newG
  for (i in 1:n){
    temp3=matrix(newE[i,],q,1)%*%newG[i,]
    ggg=setdiff(seq(from=1,to=p*(q+1),by=1),seq(from=1,to=p*(q+1),by=q+1))
    W[i,ggg]=matrix(temp3,(p*q),1)
  }

  y_predict=newE%*%alpha_estimate+W%*%b_estimate+intercept_estimate
  return(y_predict)
}

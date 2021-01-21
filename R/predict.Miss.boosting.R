#' Make predictions from a "Miss.boosting" object
#'
#'
#' This function makes predictions from a  Miss.boosting model, using the stored
#' \code{"Miss.boosting"} object.
#'
#' @param object Fitted \code{"Miss.boosting"} object.
#' @param newE Matrix of new values for \code{E} at which predictions are to be
#' made.
#' @param newG Matrix of new values for \code{G} at which predictions are to be
#' made.
#' @param \dots Not used. Other arguments to predict.
#'
#' @return The object returned depends on the \dots{} argument which is passed
#' on to the \code{predict} method for \code{Miss.boosting} objects.
#' @seealso \code{Miss.boosting}, and \code{plot} methods.
#' @references Mengyun Wu and Shuangge Ma.
#' \emph{Robust semiparametric gene-environment interaction analysis using sparse boosting.
#' Statistics in Medicine, 38(23):4625-4641, 2019.}
#'
#' @method predict Miss.boosting
#' @export
#' @export predict.Miss.boosting
predict.Miss.boosting<-function(object,newE,newG,...){
  if((missing(newE))|(missing(newG))){
    stop("You need to supply a value for newE and newG")
  }
  E_test=as.matrix(newE)
  G_test=as.matrix(newG)

  p=dim(G_test)[2]
  q=dim(E_test)[2]

  temp=object


  intercept=temp$intercept;
  unique_variable=temp$unique_variable
  unique_coef=temp$unique_coef
  unique_knots=temp$unique_knots
  unique_Boundary.knots=temp$unique_Boundary.knots
  unique_vtype=temp$unique_vtype


  id=length(unique_coef)


  if (id==1){
    unique_variable=matrix(unique_variable,1,2)
  }

  y_predict=0

  for (t in 1:id){
    if (unique_variable[t,1]==0){
      x_temp=G_test[,unique_variable[t,2]]
    } else if (unique_variable[t,2]==0){
      x_temp=E_test[,unique_variable[t,1]]
    } else {
      x_temp=cbind(E_test[,unique_variable[t,1]],G_test[,unique_variable[t,2]])
    }

    N=design.matrix(x=x_temp,knots=unique_knots[[t]],Boundary.knots=unique_Boundary.knots[[t]],degree=object$degree,v_type=unique_vtype[t],NorM=object$NorM)$X
    if(dim(N)[1]==1){
      temp=sum(N[,-1]*unique_coef[[t]])
    }else{
    temp=as.matrix(N[,-1])%*%unique_coef[[t]]
    }
    y_predict=y_predict+temp
  }
  y_predict=y_predict+intercept
  return(y_predict)
}

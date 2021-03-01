#' Extract coefficients from a "RobSBoosting" object
#'
#' This function extracts coefficients from a RobSBoosting model, using the stored
#' \code{"RobSBoosting"} object.
#' @method coef RobSBoosting
#' @param object Fitted \code{"RobSBoosting"} model object.
#' @param \dots Not used. Other arguments to get coefficients.
#'
#' @return
#' \item{intercept}{The intercept estimate.}
#' \item{unique_variable}{A matrix with two columns that represents the variables that are selected
#' for the model after removing the duplicates, since the \code{loop_time} iterations of the method
#' may produce variables that are repeatedly selected into the model. Here, the first and second
#' columns correspond to the indexes of environmental (E) factors and genetic (G) factors. For
#' example, (1, 0) represents that this variable is the first E factor, and (1,2) represents that the
#' variable is the interaction between the first E factor and second G factor.}
#' \item{unique_coef}{Coefficients corresponding to \code{unique_variable}. Here, the coefficients
#' are simple regression coefficients for the linear effect (discrete E factor, G factor, and their
#' interaction), and B spline coefficients for the nonlinear effect (continuous E factor, and
#' corresponding G-E interaction).}
#' \item{unique_knots}{A list of knots corresponding to \code{unique_variable}. Here, when the type
#' of \code{unique_variable} is discrete E factor, G factor, or their interaction, knot will be NULL,
#' and knots will be B spline otherwise.}
#' \item{unique_Boundary.knots}{A list of boundary knots corresponding to \code{unique_variable}.}
#' \item{unique_vtype}{A vector representing the variable type of \code{unique_variable}. Here, "EC"
#' stands for continuous E effect, "ED" for discrete E effect, "G" for G effect, "EC-G" for the
#' interaction between "EC" and "G", and "ED-G" for the interaction between "ED" and "G".}
#' \item{estimation_results}{A list of estimation results for each variable. Here, the first
#' \code{q} elemnets are for the E effects, the (\code{q+1}) element
#' is for the first G effect and the (\code{q+2}) to (\code{2q+1}) elements are for the interactions
#' corresponding to the first G factor, and so on.}
#' @seealso \code{RobSBoosting}, and \code{predict}, and \code{plot} methods.
#' @references Mengyun Wu and Shuangge Ma.
#' \emph{Robust semiparametric gene-environment interaction analysis using sparse boosting.
#' Statistics in Medicine, 38(23):4625-4641, 2019.}
#' @importFrom splines bs
#' @export
#' @export coef.RobSBoosting

coef.RobSBoosting<-function(object,...){
  v_type=object$v_type
  p=sum(v_type=="G")
  q=sum(v_type=="EC"|v_type=="ED")


  variable_pair=unique(object$variable_pair[object$variable[1:object$id],],MARGIN=1)
  if (!is.matrix(variable_pair)){
    variable_pair=matrix(variable_pair,1,2)
  }

  alpha0=matrix(0,q,1)
  G_main=matrix(0,p,1)
  GE_interaction=matrix(0,q,p)
  temp=variable_pair[which(variable_pair[,2]==0),1]
  alpha0[temp]=1
  temp=variable_pair[which(variable_pair[,1]==0),2]
  G_main[temp]=1

  temp=variable_pair[((variable_pair[,1]!=0) & (variable_pair[,2]!=0)),]

  GE_interaction[temp]=1

  beta0=matrix(0,p+p*q,1)

  for (j in 1:p) {
    beta0[(j-1)*(q+1)+1]=G_main[j]
    beta0[((j-1)*(q+1)+2):(j*(q+1))]=GE_interaction[,j]
  }




  unique_temp=unique(object$variable[1:object$id])
  unique_variable=object$variable_pair[unique_temp,]
  unique_coef=vector('list',length(unique_temp))
  unique_knots=vector('list',length(unique_temp))
  unique_Boundary.knots=vector('list',length(unique_temp))
  for (i in 1:length(unique_temp)){
    unique_coef[[i]]=0
  }


  v=object$v
  a=0
  for (i in 1:object$id){
    id_temp=which(unique_temp==object$variable[i])
    unique_coef[[id_temp]]=unique_coef[[id_temp]]+object$spline_result[[i]]$estimates[-1]*v
    a=a+object$spline_result[[i]]$estimates[1]*v
    if ((!is.null(object$spline_result[[i]]$knots))){
      unique_knots[[id_temp]]=object$spline_result[[i]]$knots
      unique_Boundary.knots[[id_temp]]=object$spline_result[[i]]$Boundary.knots
    }
  }



  unique_vtype=object$v_type[unique_temp]
  #
  #
  #
  # unique_set=list(intercept=a,unique_variable=unique_variable,unique_coef=unique_coef,unique_knots=unique_knots,unique_Boundary.knots=unique_Boundary.knots,unique_vtype=unique_vtype)
  id=length(unique_coef)

  if (id==1){
    unique_variable=matrix(unique_variable,1,2)
  }
  #
  #
  # xx_dot=seq(from=0,to=1,by=0.0001)
  #
  # estimation_results=vector('list',p+q+p*q) ## E+G+E*G
  #
  #
  #
  #
  # id_temp1=which(unique_set$unique_vtype=='EC')
  # id_temp2=unique_set$unique_variable[id_temp1,1]+unique_set$unique_variable[id_temp1,2]*(q+1)
  # if(length(id_temp1>0)){
  # for (i in 1:length(id_temp1)){
  #
  #   xx=splines::bs(xx_dot, knots=unique_set$unique_knots[[id_temp1[i]]], intercept=TRUE, degree=object$degree,Boundary.knots = unique_set$unique_Boundary.knots[[id_temp1[i]]])
  #   xx=xx[,-1]
  #   xx=xx-(matrix(1,length(xx_dot),1)%*%object$NorM)
  #   bs_predict=xx%*%unique_set$unique_coef[[id_temp1[i]]]
  #
  #
  #   estimation_results[[id_temp2[i]]]=bs_predict
  # }
  #
  # id_temp1=which(unique_set$unique_vtype=='G-EC')
  # id_temp2=unique_set$unique_variable[id_temp1,1]+unique_set$unique_variable[id_temp1,2]*(q+1)
  # id_temp3=unique_set$unique_variable[id_temp1,2]
  #
  # for (i in 1:length(id_temp1)){
  #
  #   xx=splines::bs(xx_dot, knots=unique_set$unique_knots[[id_temp1[i]]], intercept=TRUE, degree=object$degree,Boundary.knots = unique_set$unique_Boundary.knots[[id_temp1[i]]])
  #   xx=xx[,-1]
  #   xx=xx-(matrix(1,length(xx_dot),1)%*%object$NorM)
  #   bs_predict=xx%*%unique_set$unique_coef[[id_temp1[i]]]
  #
  #
  #
  #   estimation_results[[id_temp2[i]]]=bs_predict
  # }
  # }
  #
  #
  #
  # id_temp1=which(((unique_set$unique_vtype=='ED') | (unique_set$unique_vtype=='G') | (unique_set$unique_vtype=='G-ED'))!=0)
  # id_temp2=unique_set$unique_variable[id_temp1,1]+unique_set$unique_variable[id_temp1,2]*(q+1)
  #
  #
  # for (i in 1:length(id_temp1)){
  #   estimation_results[[id_temp2[i]]]=unique_set$unique_coef[[id_temp1[i]]]
  # }
  #
  # intercept=unique_set$intercept


  id=length(unique_coef)

  if (id==1){
    unique_variable=matrix(unique_variable,1,2)
  }

  estimation_results=object$estimation_results
  return(list(intercept=a,unique_variable=unique_variable,unique_coef=unique_coef,unique_knots=unique_knots,unique_Boundary.knots=unique_Boundary.knots,unique_vtype=unique_vtype,estimation_results=estimation_results))
}

#' Robust semiparametric gene-environment interaction analysis using sparse boosting
#'
#' Robust semiparametric gene-environment interaction analysis using sparse boosting. Here a
#' semiparametric model is assumed to accommodate nonlinear effects, where we model continuous
#' environmental (E) factors in a nonlinear way, and discrete E factors and all genetic (G)
#' factors in a linear way. For estimating the nonlinear functions, the B spline expansion is
#' adopted. The Huber loss function and Qn estimator are adopted to accommodate long-tailed
#' distribution/data contamination. For model estimation and selection of relevant variables, we
#' adopt an effective sparse boosting approach, where the strong hierarchy is respected.

#' @param G Input matrix of \code{p} genetic measurements consisting of \code{n} rows. Each row
#' is an observation vector.
#' @param E Input matrix of \code{q} environmental risk factors, each row is an observation
#' vector.
#' @param Y Response variable. A quantitative vector for \code{family="continuous"}. For
#' \code{family="survival"}, \code{Y} should be a two-column matrix with the first column being
#' the log(survival time) and the second column being the censoring indicator. The indicator is a
#' binary variable, with "1" indicating dead, and "0" indicating right censored.
#' @param loop_time Number of iterations of the sparse boosting.
#' @param num.knots Numbers of knots for the B spline basis.
#' @param Boundary.knots The boundary of knots for the B spline basis.
#' @param degree Degree for the B spline basis.
#' @param v The step size used in the sparse boosting process. Default is 0.1.
#' @param family Response type of \code{Y} (see above).
#' @param knots List of knots for the B spline basis. Default is NULL and knots can be generated
#' with the given \code{num.knots}, \code{degree} and \code{Boundary.knots}.
#' @param E_type A vector indicating the type of each E factor, with "ED" representing discrete E factor, and "EC" representing continuous E factor.
#' @return An object with S3 class \code{"RobSBoosting"} is returned, which is a list with the following components.
#' \item{call}{The call that produced this object.}
#' \item{max_t}{The stopping iteration time of the sparse boosting.}
#' \item{spline_result}{A list of length \code{max_t} that includes the estimation results of
#' each iteration.}
#' \item{BIC}{A vector of length max_t that includes Bayesian Information Criterion based on the
#' Huber's prediction error.}
#' \item{variable}{A vector of length max_t that includes the index of selected variable in each
#' iteration.}
#' \item{id}{The iteration time with the smallest BIC.}
#' \item{variable_pair}{A matrix with two columns that include the set of variables that can
#' potentially enter the regression model at the stopping iteration time. Here, the first and
#' second columns correspond to the indexes of E factors and G factors. For example, (1, 0)
#' represents that this variable is the first E factor, and (1,2) represents that the variable is
#' the interaction between the first E factor and second G factor. }
#' \item{v_type}{A vector whose length is the number of rows of \code{variable_pair}, with each
#' element representing the variable type of the corresponding row of \code{variable_pair}. Here,
#' "EC" stands for continuous E effect, "ED" for discrete E effect, and "G"  for G effect, "EC-G"
#' for the interaction between "EC" and "G", "ED-G" for the interaction between "ED" and "G".}
#' \item{family}{The same as input \code{family}.}
#' \item{degree}{Degree for the B spline basis.}
#' \item{v}{The step size used in the sparse boosting process.}
#' \item{NorM}{The values of B spline basis.}
#'
#' @seealso  \code{coef}, \code{predict}, and \code{plot} methods, and \code{Miss.boosting}
#' method.
#' @references Mengyun Wu and Shuangge Ma.
#' \emph{Robust semiparametric gene-environment interaction analysis using sparse boosting.
#' Statistics in Medicine, 38(23):4625-4641, 2019.}
#' @importFrom stats quantile
#' @importFrom stats mad
#' @importFrom pcaPP qn
#' @importFrom splines bs
#' @examples
#' data(Rob_data)
#' G=Rob_data[,1:20];E=Rob_data[,21:24]
#' Y=Rob_data[,25];Y_s=Rob_data[,26:27]
#' knots = list();Boundary.knots = matrix(0, 24, 2)
#' for(i in 1:4) {
#'   knots[[i]] = c(0, 1)
#'   Boundary.knots[i, ] = c(0, 1)
#'   }
#'
#' #continuous
#' fit1= RobSBoosting(G,E,Y,loop_time = 80,num.knots = 2,Boundary.knots=Boundary.knots,
#' degree = 2,family = "continuous",knots = knots,E_type=c("EC","EC","ED","ED"))
#' coef1 = coef(fit1)
#' predict1=predict(fit1,newE=E[1:2,],newG=G[1:2,])
#' plot(fit1)
#'
#' \donttest{
#' #survival
#' fit2= RobSBoosting(G,E,Y_s,loop_time = 200, num.knots = 2, Boundary.knots=Boundary.knots,
#' family = "survival", knots = knots,E_type=c("EC","EC","ED","ED"))
#' coef2 = coef(fit2)
#' predict2=predict(fit2,newE=E[1:2,],newG=G[1:2,])
#' plot(fit2)
#' }

#' @export
#' @export RobSBoosting


RobSBoosting<-function(G,E,Y,loop_time,num.knots=NULL,Boundary.knots=NULL,degree=1,v=0.1,family=c("continuous","survival"),knots=NULL,E_type){
  if(is.null(knots)){
    if((is.null(Boundary.knots))|(is.null(num.knots))|(is.null(degree)))
    stop("You need to supply 'Boundary.knots', 'degree' and 'num.knots' since you do not input 'knots")
  }

  # get call and family, E_type, Method
  thisCall = match.call()
  family = match.arg(family)
  Method = "Robust"
  y=Y
  n<-dim(E)[1]
  q<-dim(E)[2]
  p<-dim(G)[2]

  # E_type=vector()
  # for (j in 1:q){
  #   if (length(unique(E[,j]))<n/2){
  #     E_type[j]='ED'
  #   }else E_type[j]='EC'
  # }




  n=dim(E)[1]
  if (family=='survival'){
    w=kmw(y[,1],y[,2])
    delta=y[,2]
    y=y[,1]
  } else {
    w=NULL
    delta=matrix(1,n,1)
  }

  kmweight=w[delta==1]
  G=G[delta==1,]
  E=E[delta==1,]
  y=y[delta==1]

  n=dim(E)[1]
  p=dim(G)[2]
  q=dim(E)[2]

  if (!is.null(knots)){
    ss<-seq(from=1/(num.knots+1),to=num.knots/(num.knots+1),by=1/(num.knots+1))
    knots_temp=vector('list',p+q)
    for (ii in 1:q){
      knots_temp[[ii]]=stats::quantile(knots[[ii]],ss)
    }
    knots=knots_temp
    xx_dots=seq(from=0,to=1,by=0.0001)
    NorM = splines::bs(xx_dots, knots=knots[[1]],degree=degree,Boundary.knots = Boundary.knots[1,])
    NorM=colMeans(NorM)
  } else {
    xx_dots=seq(from=0,to=1,by=0.0001)
    NorM = splines::bs(xx_dots, df=num.knots+degree,degree=degree)
    NorM=colMeans(NorM)
  }



  variable<-BIC<-loglike<-c()
  result=list()

  variable_pair=matrix(0,p+q,2)
  variable_pair[,1]=c(1:q,rep(0,p))
  variable_pair[,2]=c(rep(0,q),1:p)

  colnames(variable_pair)<-c('E','G')

  GE=cbind(E,G)

  v_type=c(E_type,rep('G',p))

  designmat=vector('list',dim(GE)[2])


  if (is.null(knots)){
    for (j in 1:dim(GE)[2]){
      designmat[[j]]=design.matrix(GE[,j],num.knots,knots,Boundary.knots,degree,v_type[j],model_type='nonlinear',NorM)
    }
  } else {
    for (j in 1:dim(GE)[2]){
      designmat[[j]]=design.matrix(GE[,j],num.knots,knots[[j]],Boundary.knots[j,],degree,v_type[j],model_type='nonlinear',NorM)
    }
  }


  f=0
  u=y-f


  t=1
  f_temp=f


  while (t<=loop_time) {
    pp=length(designmat)
    res=matrix(1e+20,pp,1)
    hatmat_set=vector('list',pp)


    for (j in 1:pp){

      if (Method=='Robust'){
        temp=robust.spline(designmat[[j]],u,y-f,kmweight)
      } else if (Method=='LS') {
        temp=ls.spline(designmat[[j]],u,kmweight)
      }

      y_predict=temp$y_hat
      f1=y_predict+f_temp
      hatmat_set[[j]]=temp

      if (Method=='Robust'){
        ptemp=dim(temp$X)[2]

        if (ptemp==2){
          sdd1=sqrt(mean((temp$X[,2:ptemp]-mean(temp$X[,2:ptemp]))^2))
          sdd2=pcaPP::qn(u)

          RSS=n*(sdd2^2-(temp$estimate[2]*sdd1)^2)

        } else {
          xtemp=temp$X[,2:ptemp]
          xxtemp=xtemp-matrix(1,n,1)%*%colMeans(xtemp)
          robustcov=(t(xxtemp)%*%xxtemp)/n
          beta_temp=matrix(temp$estimates[2:ptemp],1,ptemp-1)
          RSS=n*(pcaPP::qn(u)^2-beta_temp%*%robustcov%*%t(beta_temp))
        }
      } else if (Method=='LS') {
        RSS=sum((y-f1)^2)
      }
      if (RSS<0){
        RSS=1e+100
      }

      df=length(unique(c(variable[1:(t-1)],j)))

      res[j]=log(RSS)+df*log(n)/n

    }



    id=which.min(res)
    temp=hatmat_set[[id]]
    y_predict=temp$y_hat
    coef=temp$coef

    result[[t]]=temp
    variable[t]=id
    f=f+v*y_predict
    f_temp=f


    u=y-f

    id_unique=unique(variable[1:t])


    if (length(id_unique)>1){ #### have at least two main effects
      if (variable[t]<=p+q){ #### the selected variable is a main effect
        if (!is.element(variable[t],variable[1:(t-1)])){ # the selected variable is a new main effect
          id_E=id_unique[which(id_unique<=q)]
          id_G=id_unique[which((id_unique<=(p+q)) & (id_unique>q))]-q
          if ((length(id_E)>0) & (length(id_G)>0)) {# have at least one E and one G effects
            if (variable[t]<=q) {

              temp_pair=matrix(0,length(id_G),2)
              temp_pair[,1]=rep(variable[t],length(id_G))
              temp_pair[,2]=id_G
              vtype_add=matrix(0,1,length(id_G))
              ppp=pp+length(id_G)
              for (ii in 1:length(id_G)){
                if (E_type[variable[t]]=='ED'){
                  vtype_add[ii]='G-ED'
                } else {
                  vtype_add[ii]='G-EC'
                }
                temp_data=cbind(E[,variable[t]],G[,id_G[ii]])



                if (is.null(knots)){
                  designmat[[pp+ii]]=design.matrix(temp_data,num.knots,knots,Boundary.knots,degree,vtype_add[ii],model_type='nonlinear',NorM)
                } else {
                  designmat[[pp+ii]]=design.matrix(temp_data,num.knots,knots[[variable[t]]],Boundary.knots[variable[t],],degree,vtype_add[ii],model_type='nonlinear',NorM)
                }

              }
            } else {

              temp_pair=matrix(0,length(id_E),2)
              temp_pair[,1]=id_E
              temp_pair[,2]=rep(variable[t]-q,length(id_E))
              vtype_add=matrix(0,1,length(id_E))
              ppp=pp+length(id_E)
              for (ii in 1:length(id_E)){
                if (E_type[id_E[ii]]=='ED'){
                  vtype_add[ii]='G-ED'
                } else {
                  vtype_add[ii]='G-EC'
                }
                temp_data=cbind(E[,id_E[ii]],G[,variable[t]-q])
                if (is.null(knots)){
                  designmat[[pp+ii]]=design.matrix(temp_data,num.knots,knots,Boundary.knots,degree,vtype_add[ii],model_type='nonlinear',NorM)
                } else {

                  designmat[[pp+ii]]=design.matrix(temp_data,num.knots,knots[[id_E[ii]]],Boundary.knots[id_E[ii],],degree,vtype_add[ii],model_type='nonlinear',NorM)


                }
              }
            }

            v_type=c(v_type,vtype_add)

            variable_pair=rbind(variable_pair,temp_pair)
          }
        }
      }
    }



    if (Method=='Robust') {
      absy=abs(y-f)
      MAD_v=stats::mad(y-f)
      c = 1.345*MAD_v
      Huber=absy^2
      Huber[absy>c]=2*c*(absy[absy>c]-c/2)
      RSS=sum(Huber)
    } else if (Method=='LS') {
      RSS=sum((y-f)^2)
    }

    df=length(unique(variable[1:t]))
    BIC[t]=log(RSS)+df*log(n)/(n)


    if (t>1000) {
      if (abs((BIC[t]-BIC[t-1])/BIC[t-1])<1e-8){
        break
      }
    }
    t=t+1
  }



  BIC=BIC[1:t]
  variable=variable[1:t]
  loglike=loglike[1:t]

  id=which.min(BIC)

  output=list(call=thisCall,max_t=t,spline_result=result,BIC=BIC,variable=variable,id=id,variable_pair=variable_pair,v_type=v_type,family=family,degree=degree,v=v,NorM=NorM)
  class(output) ="RobSBoosting"
  return(output)
}


#' Plot coefficients from a "Miss.boosting" object
#'
#' Draw plots for estimated parameters in a fitted
#' \code{"Miss.boosting"} object, including a heatmap for discrete
#' environmental (E) effects, and selected genetic (G) effects and G-E interactions, and plots for each
#' of selected continuous E (EC) effect and interactions between EC and G.
#'
#' @param x Fitted \code{"Miss.boosting"} model.
#'
#' @param \dots Other graphical parameters to plot.
#' @return A heatmap for estimated coefficients.
#' @seealso \code{Miss.boosting}, and \code{predict} methods.
#' @references Mengyun Wu and Shuangge Ma.
#' \emph{Robust semiparametric gene-environment interaction analysis using sparse boosting.
#' Statistics in Medicine, 38(23):4625-4641, 2019.}
#' @method plot Miss.boosting
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom graphics plot
#' @export
#' @export plot.Miss.boosting
plot.Miss.boosting=function(x,...){
  fit=x
  q=length(fit$alpha0)
  p=length(fit$beta0)/(q+1)
  E_type=fit$E_type[1:q]############differ
  v_type=fit$unique_vtype
  loc=matrix(0:(q+p+q*p),(q+1),(p+1))
  variable_pair=fit$unique_variable
  iin=which(is.element(E_type[1:q],c("ED")))
  nolinear_c=setdiff(c(1:q),iin)


  #####change
  cnames=rownames(fit$beta0)
  temp=seq(1,length(fit$beta0),by=(q+1))
  names_G=cnames[temp]

  cnames=rownames(fit$alpha0)
  names_E=cnames
  colnames(loc)=c("G",names_G)
  rownames(loc)=c("E",names_E)

  linear_id_G=as.matrix(loc[1,-1],ncol=1)


  linear_id1=loc[(iin+1),]
  linear_id2=matrix(loc[(iin+1),],(length(iin)*(p+1)),1)


  linear_id=rbind(linear_id_G,linear_id2)

  # estimation_results=coef.Miss.boosting(fit)$estimation_results
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

  estimation_results=vector('list',p+q+p*q) ## E+G+E*G

  xx_dot=seq(from=0,to=1,by=0.0001)


  id_temp1=which(v_type=='EC')
  id_temp2=fit$unique_variable[id_temp1,1]+fit$unique_variable[id_temp1,2]*(q+1)

  if(length(id_temp1)>0){
  for (i in 1:length(id_temp1)){

    xx=splines::bs(xx_dot, knots=fit$unique_knots[[id_temp1[i]]], intercept=TRUE, degree=fit$degree,Boundary.knots =fit$unique_Boundary.knots[[id_temp1[i]]])
    xx=xx[,-1]
    xx=xx-(matrix(1,length(xx_dot),1)%*%fit$NorM)
    bs_predict=xx%*%fit$unique_coef[[id_temp1[i]]]


    estimation_results[[id_temp2[i]]]=bs_predict
  }


  id_temp1=which(v_type=='G-EC')
  id_temp2=fit$unique_variable[id_temp1,1]+fit$unique_variable[id_temp1,2]*(q+1)
  id_temp3=fit$unique_variable[id_temp1,2]

  for (i in 1:length(id_temp1)){

    xx=splines::bs(xx_dot, knots=fit$unique_knots[[id_temp1[i]]], intercept=TRUE, degree=fit$degree,Boundary.knots = fit$unique_Boundary.knots[[id_temp1[i]]])
    xx=xx[,-1]
    xx=xx-(matrix(1,length(xx_dot),1)%*%fit$NorM)
    bs_predict=xx%*%fit$unique_coef[[id_temp1[i]]]



    estimation_results[[id_temp2[i]]]=bs_predict
  }


  id_temp1=which(((fit$unique_vtype=='ED') | (fit$unique_vtype=='G') | (fit$unique_vtype=='G-ED'))!=0)
  id_temp2=fit$unique_variable[id_temp1,1]+fit$unique_variable[id_temp1,2]*(q+1)


  for (i in 1:length(id_temp1)){
    estimation_results[[id_temp2[i]]]=fit$unique_coef[[id_temp1[i]]]
  }
  }


  temp=matrix(0,(p+1),(length(iin)+1))
  if(length(iin)==0) {
  colnames(temp)= "G"
  }else{
    iinE=names_E[iin]
    colnames(temp)=c("G",iinE)
  }
  rownames(temp)=c("E",colnames(loc)[-1])
  ttt=lapply(estimation_results, function(x) ifelse(length(x)==0,0,x))
  temp[1,]=c(0,unlist(ttt[iin]))
  tempp=estimation_results[linear_id_G]
  names(tempp)= names_G
  hi=sapply(tempp,length)
  te=which(hi!=0)



  if(length(te)>1)
    temp[(te+1),1]=unlist(tempp[te])#########################change
  mm=estimation_results[linear_id2]
  pp=as.numeric(sapply(mm,function(x) ifelse(length(x)==0,0,x)))#####change
  t22=matrix(pp,nrow=(p+1),byrow = TRUE)

  temp[1:(p+1),-1]=t22
  x=temp
  xx=as.matrix(x[rowSums(x)!=0,])
  if(rowSums(x)[1]==0){ xx=rbind(x[1,],xx); rownames(xx)[1]="E"}
  x1=xx
  if(dim(x1)[1]>0)
    rownames(x1)[1]="E"
  colnames(x1)=colnames(temp)
  colnames(x1)[1]="G"

  tt=as.matrix(x1[,dim(x1)[2]:1])
  colnames(tt)=rev(colnames(x1))
  rownames(tt)=rownames(x1)
  if(dim(tt)[1]>0){
  data <- as.data.frame(t(tt))
  data$ID <- rownames(data)
  data_m <- reshape2::melt(data, id.vars=c("ID"))
  variable=data_m$variable
  ID=data_m$ID
  value=data_m$value
  da<-factor(data_m$ID)
  tttt=sapply(levels(da),function(x) substr(x,2,nchar(x)))

  cnames=paste("E",1:q,sep="")
  if(sum(rownames(loc)==c("E",cnames))==(q+1)){
    levels_order=levels(da)[order(as.numeric(tttt))]
    temp1=temp=rev(levels_order)
  }else{
    temp1=temp=levels(da)
    temp[1]="G"
    temp[-1]=setdiff(levels(da),temp[1])
    temp1=temp
  }

  if(length(temp)==1){
    temp1=temp
    }else{
  temp1[1:(length(temp)-1)]=temp[-1]
  temp1[(length(temp))]=temp[1]
    }
  figure <- ggplot2::ggplot(data_m, ggplot2::aes(x=variable,y=ID)) +
    ggplot2::xlab("Selected Genetic factors")+ggplot2::ylab("Environment variables")  +  ggplot2::theme_classic() + ggplot2::theme(axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank()) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank()) +
    ggplot2::theme(legend.key=ggplot2::element_blank())  +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=0,hjust=1, vjust=1)) +
    ggplot2::theme(legend.position="top") +
    ggplot2::geom_tile(ggplot2::aes(fill=value)) +
    ggplot2::scale_fill_gradient2("Value",
      low = "grey",
      high = "red",
      mid = "green")+ggplot2::scale_x_discrete(position = "top")  + ggplot2::scale_y_discrete(limits = temp1)


  nolinear_c=setdiff(c(1:q),iin)

  cnames=paste("G",1:p,sep="")


  nonlinear_id1=loc[(nolinear_c+1),]
  nonlinear_id2=matrix(loc[(nolinear_c+1),],(length(nolinear_c)*(p+1)),1)


  # temp=matrix(0,(p),(length(nolinear_c)))
  # colnames(temp)=c(paste("E",nolinear_c,sep=""))
  # rownames(temp)=setdiff(colnames(loc),"G")
  #



  xx_dot=seq(from=0,to=1,by=0.0001)
  for(i in nonlinear_id2){
    if(!(is.null(estimation_results[[i]]))){
      index_G=colnames(loc)[which(loc==i,arr.ind = T)[2]]
      index_E=rownames(loc)[which(loc==i,arr.ind = T)[1]]
      if(index_G=="G"){plot(xx_dot,estimation_results[[i]],col="red",lty=2,lwd=1.8,main=paste0("estimation results for ",index_E),ylab="coefficients")
      }else{
        plot(xx_dot,estimation_results[[i]],col="red",lty=2,lwd=1.8,main=paste0("estimation results for ",index_E," and ", index_G),ylab="coefficients")}
    }
  }
  figure
  }
}


#' Plot coefficients from a "PTReg" object
#'
#' Draw a heatmap for estimated coefficients in a fitted
#' \code{"PTReg"} object.
#'
#' @param x Fitted \code{"PTReg"} model.
#' @param \dots Other graphical parameters to plot.
#' @return A heatmap for estimated coefficients.
#' @seealso \code{PTReg}, and \code{predict}, and \code{coef}
#' methods.
#' @references Yaqing Xu, Mengyun Wu, Shuangge Ma, and Syed Ejaz Ahmed.
#' \emph{Robust gene-environment interaction analysis using penalized trimmed regression. Journal
#' of Statistical Computation and Simulation, 88(18):3502-3528, 2018.}
#' @method plot PTReg
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom graphics plot
#' @export
#' @export plot.PTReg
plot.PTReg=function(x,...){
  object=x
  alpha=c(0,object$alpha)
  q=dim(object$alpha)[1]
  beta=t(object$beta)
  index_G=which(beta[,1]!=0)
  x=rbind(alpha,beta[index_G,])
  #####change
  temp=rownames(object$alpha)
  cnames=temp
  cnames=c("G",cnames)
  colnames(x)=cnames

  temp=colnames(object$beta)[index_G]
  cnames=temp
  cnames=c("E",cnames)
  rownames(x)=cnames

  x1=x
  tt=x1[,dim(x1)[2]:1]
  data <- as.data.frame(t(tt))
  data$ID <- rownames(data)
  data_m <- reshape2::melt(data, id.vars=c("ID"))
  variable=data_m$variable
  ID=data_m$ID
  value=data_m$value
  da<-factor(data_m$ID)
  tttt=sapply(levels(da),function(x) substr(x,2,nchar(x)))
  cnames=paste("E",1:q,sep="")
  if(sum(colnames(x)==c("G",cnames))==(q+1)){
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
  figure
}

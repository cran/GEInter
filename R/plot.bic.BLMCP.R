#' Plot coefficients from a "bic.BLMCP" object
#'
#' Draw a heatmap for estimated coefficients in a fitted
#' \code{"bic.BLMCP"} object.
#'
#' @param x Fitted \code{"bic.BLMCP"} model.
#' @param \dots Other graphical parameters to plot.
#' @return A heatmap for estimated coefficients.
#' @seealso \code{predict}, \code{coef} and \code{BLMCP} methods.
#' @references Mengyun Wu, Yangguang Zang, Sanguo Zhang, Jian Huang, and Shuangge Ma.
#' \emph{Accommodating missingness in environmental measurements in gene-environment interaction
#' analysis. Genetic Epidemiology, 41(6):523-554, 2017.}\cr  Jin Liu, Jian Huang, Yawei Zhang,
#' Qing Lan, Nathaniel Rothman, Tongzhang Zheng, and Shuangge Ma.
#' \emph{Identification of gene-environment interactions in cancer studies using penalization.
#' Genomics, 102(4):189-194, 2013.}
#' @method plot bic.BLMCP
#' @importFrom reshape2 melt
#' @import ggplot2
#' @importFrom graphics plot
#' @export
#' @export plot.bic.BLMCP
plot.bic.BLMCP=function(x,...){
  object=x
  q=length(object$alpha_estimate)
  p=length(object$beta_estimate)/(q+1)
  alpha=c(0,object$alpha_estimate)
  beta=t(matrix(object$beta_estimate,q+1,p))
  index_G=which(beta[,1]!=0)
  x=rbind(alpha,beta[index_G,])

  #####change
  temp=rownames(object$alpha_estimate)
  cnames=temp
  cnames=c("G",cnames)
  colnames(x)=cnames

  temp=colnames(object$beta_estimate)[index_G]
  cnames=temp
  cnames=c("E",cnames)
  rownames(x)=cnames
  # rc <- rainbow(nrow(x), start = 0, end = .3)
  # cc <- rainbow(ncol(x), start = 0, end = .3)
  # hv <- stats::heatmap(x, col = cm.colors(256), scale = "column",
  #   RowSideColors = rc, ColSideColors = cc, margins = c(5,10),
  #   xlab = "Environment variables", ylab =  "Selected Genetic factors",
  #   main = "heatmap for the estimated parameters")
  x1=x
  tt=x1[,dim(x1)[2]:1]
  data <- as.data.frame(t(tt))
  data$ID <- rownames(data)
  data_m <- melt(data, id.vars=c("ID"))
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

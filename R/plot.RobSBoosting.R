#' Plot coefficients from a "RobSBoosting" object
#'
#' Draw a plot for estimated parameters in a fitted
#' \code{"RobSBoosting"} object, where a heat map is displaced for selected discrete
#' environmental (E) effects, genetic (G) effects and their interactions, and a plot for
#' each of selected continuous E (EC) effect and interactions between EC and G.
#'
#' @param x Fitted \code{"RobSBoosting"} model.
#'
#' @param \dots Other graphical parameters to plot.
#' @return Plots for estimated coefficients.
#' @seealso \code{RobSBoosting}, \code{predict} and \code{coef}
#' methods.
#' @references Mengyun Wu and Shuangge Ma.
#' \emph{Robust semiparametric gene-environment interaction analysis using sparse boosting.
#' Statistics in Medicine, 38(23):4625-4641, 2019.}
#' @method plot RobSBoosting
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom graphics plot
#' @export
#' @export plot.RobSBoosting
plot.RobSBoosting=function(x,...){
   fit=x
   v_type=fit$v_type
   p=sum(v_type=="G")
   q=sum(v_type=="EC"|v_type=="ED")
loc=matrix(0:(q+p+q*p),(q+1),(p+1))
variable_pair=fit$variable_pair
iin=which(is.element(fit$v_type[1:q],c("ED")))
nolinear_c=setdiff(c(1:q),iin)

cnames=paste("G",1:p,sep="")
colnames(loc)=c("G",cnames)
cnames=paste("E",1:q,sep="")
rownames(loc)=c("E",cnames)
linear_id_G=as.matrix(loc[1,-1],ncol=1)


linear_id1=loc[(iin+1),]
linear_id2=matrix(loc[(iin+1),],(length(iin)*(p+1)),1)


linear_id=rbind(linear_id_G,linear_id2)

estimation_results=coef.RobSBoosting(fit)$estimation_results




temp=matrix(0,(p+1),(length(iin)+1))
colnames(temp)=c("E",paste("E",iin,sep=""))
rownames(temp)=colnames(loc)
ttt=lapply(estimation_results, function(x) ifelse(length(x)==0,0,x))
temp[1,]=c(0,unlist(ttt[iin]))
tempp=estimation_results[linear_id_G]
names(tempp)= paste("G",1:p,sep="")
hi=sapply(tempp,length)
te=which(hi!=0)

temp[2:(length(te)+1),1]=unlist(tempp[te])
# gg=apply(linear_id1, 1, names)
mm=estimation_results[linear_id2]
pp=sapply(mm,function(x) ifelse(length(x)==0,0,x))
temp[1:(p+1),-1]=pp
x=temp
xx=as.matrix(x[rowSums(x)!=0,])
x1=xx
rownames(x1)[1]="E"
colnames(x1)=colnames(temp)
colnames(x1)[1]="G"
# cnames=paste("E",1:(length(alpha)-1),sep="")
# cnames=c("G",cnames)
# colnames(x)=cnames
# cnames=paste("G",index_G,sep="")
# cnames=c("E",cnames)
# rownames(x)=cnames
# rc <- rainbow(nrow(x1), start = 0, end = .3)
# cc <- rainbow(ncol(x1), start = 0, end = .3)
# hv <- stats::heatmap(x1, col = cm.colors(256), scale = "column",
#   RowSideColors = rc, ColSideColors = cc, margins = c(5,10),
#   xlab = "Environment variables", ylab =  "Selected Genetic factors",
#   main = "heatmap for the estimated parameters")

tt=as.matrix(x1[,dim(x1)[2]:1])
colnames(tt)=rev(colnames(x1))
rownames(tt)=rownames(x1)
data <- as.data.frame(t(tt))
data$ID <- rownames(data)
data_m <- reshape2::melt(data, id.vars=c("ID"))
variable=data_m$variable
ID=data_m$ID
value=data_m$value
da<-factor(data_m$ID)
tttt=sapply(levels(da),function(x) substr(x,2,nchar(x)))
levels_order=levels(da)[order(as.numeric(tttt))]

temp1=temp=rev(levels_order)
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


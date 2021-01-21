###weighted OLS
kmw<-function(y,delta){
  y_s=y
  delta_s=delta
  kmweight<-c()
  nw<-length(y)
  
  comb<-cbind(y,delta)
  oo=order(y)
  ocomb<-comb[oo,]
  y<-ocomb[,1]
  delta<-ocomb[,2]
  kmweight[1]<-delta[1]/nw
  for(ind in 2:nw){
    tmp<-c()
    for(ind2 in 1:(ind-1)){
      tmp[ind2]<-((nw-ind2)/(nw-ind2+1))^delta[ind2]
    }
    kmweight[ind]<-delta[ind]/(nw-ind+1)*prod(tmp)
  }
  kmweight1=matrix(0,nw,0)
  kmweight1[oo]=kmweight
  return(kmweight=kmweight1)
}

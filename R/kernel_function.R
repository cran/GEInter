kernel_function<-function(x,id,h){
  n<-dim(x)[1]
  q<-dim(x)[2]
  k<-matrix(0,n,q)
  k[,id==1]=exp(-x[,id==1]^2/(2*h[id==1]^2))
  k[,id==0]=h[id==0]^(abs(x[,id==0]))
  kernel_r=matrix(0,n,1)
  for (i in 1:n){
    kernel_r[i]<-prod(k[i,])
  }
  return(kernel_r)
}
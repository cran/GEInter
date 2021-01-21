which.column<-function(a,b){

  if (!is.matrix(b)){
    b=matrix(b,1,2)
  }

  id1=is.element(b[,1],a[1])
  id2=is.element(b[,2],a[2])
  id=which(id1*id2==1)

  return(id)
}
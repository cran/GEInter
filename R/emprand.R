emprand<-function(dist,n){


  x = c(dist)
  #Remove missing observations indicated by NaN's.
  id = !is.na(x)
  x = x[id]

  # Compute empirical cumulative distribution function (cdf)
  xlen = length(x)
  x = sort(x)
  p = seq(from=1,to=xlen,by=1)
  p = p/xlen

  # Generate uniform random number between 0 and 1
  ur =  stats::runif(n,0,1)

  # Interpolate ur from empirical cdf and extraplolate for out of range  values.
  xr = Hmisc::approxExtrap(p,x,ur,rule=2)$y
  return(xr)
}

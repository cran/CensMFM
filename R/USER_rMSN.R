rMSN <-
function(n, mu, Sigma, shape){
  #n: the number of generating data;
  #mu: Must be a vector with p entries;
  #Sigma: Must be a matrix with dimension pxp
  #shape: Must be a vector with p entries, with shape is zero we have the normal case.
  
  p <- length(mu)
  delta <- shape/(sqrt(1 + as.vector(t(shape)%*%shape)))
  y <- matrix(0,n,p)
  for (i in 1:n) y[i,] <- mu + sqrtm(Sigma)%*%(delta*abs(rnorm(1)) + sqrtm(diag(p) - delta%*%t(delta))%*%as.vector(rmvnorm(1, mean = rep(0,p), sigma = diag(p))))
  return(y)
}

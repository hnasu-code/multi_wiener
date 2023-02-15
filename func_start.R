## Density of simple Wiener process without time scale transformation-------------------------------------------------------
simpleWiener = function(mu,sigma){
  nmu = mu*deltat
  nsigma = sigma*deltat
  density = dnorm(deltay,mean=nmu,sd=sqrt(nsigma),log=T)
  return(-sum(density))
}

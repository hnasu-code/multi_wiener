## Data generation for model M0-------------------------------------------------
data_generation_M0 = function(seed,n,p,m,interval,gamma,mu_a,mu_b,
                              Sigma_a,Sigma_b,Sigma_ab,sigma_D){
  # generate random effects: y0 and b
  set.seed(seed)
  re_a = mvrnorm(n,mu=mu_a,Sigma=Sigma_a)
  while(sum(re_a<0)>0){
    set.seed(seed*7+12345678)
    re_a = mvrnorm(n,mu=mu_a,Sigma=Sigma_a)
    seed = seed+1
  }
  re_b = matrix(rep(0,n*p),ncol=p)
  Sigma_b_a = round(Sigma_b-t(Sigma_ab)%*%solve(Sigma_a)%*%Sigma_ab,digits=7)
  for (i in 1:n){
    mu_b_a = mu_b + t(Sigma_ab)%*%solve(Sigma_a)%*%(re_a[i,]-mu_a)
    re_b[i,] = mvrnorm(1,mu=mu_b_a,Sigma=Sigma_b_a)
    while(sum(re_b[i,]<0)>0){
      set.seed(seed*7+12345678)
      re_b[i,] = mvrnorm(1,mu=mu_b_a,Sigma=Sigma_b_a)
      seed = seed+1
    }
  }
  re = cbind(re_a,re_b)
  
  # generate degradation paths
  dy = matrix(rep(0,n*m*p),ncol=m*n)
  y = matrix(rep(0,n*(m+1)*p),ncol=(m+1)*n)
  y0 = matrix(rep(0,n*p),ncol=n)
  dt = matrix(rep(interval,n*m*p),ncol=m*n)
  t = matrix(rep(0,n*(m+1)*p),ncol=(m+1)*n)
  
  
  for (i in 1:n){
    for (j in 1:p){
      t[j,((i-1)*m+i+1):(i*m+i)] = cumsum(t(dt[j,((i-1)*m+1):(i*m)]))
      for (k in 1:m){
        dtau = t[j,(i-1)*m+k+i]^gamma[j] - t[j,(i-1)*m+k+i-1]^gamma[j]
        dy[j,k+m*(i-1)] = rnorm(1,mean=dtau*re[i,p+j],sd=sqrt(dtau*sigma_D[j]))
      }
      y[j,((i-1)*m+i):(i*m+i)] = cumsum(t(c(re[i,j],dy[j,((i-1)*m+1):(i*m)])))
      y0[j,i] = y[j,i*(m+1)-m]
    }
  }
  
  return(list(n=n,p=p,m=m,dy=dy,y=y,y0=y0,dt=dt,t=t))
} 
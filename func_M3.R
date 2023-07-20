## EM algorithm=======================================================

## The notations used in the code v.s. notations used in the manuscript
##         Used in code     --     Used in manuscript
##               n          --            n
##               p          --            p
##               m          --            m, identical across all units
##            mu_a.hat      --       \hat{\mu}_a 
##            mu_b.hat      --       \hat{\mu}_b 
##           sigma_a.hat    --       \hat{\sigma}^2_a 
##           sigma_b.hat    --       \hat{\sigma}^2_b 
##           rho_ab.hat     --       \hat{\rho}_ab 
##           sigma_D.hat    --       \hat{\sigma}^2
##           gamma.hat      --       \hat{\gamma} assuming tau=t_{ik}^gamma-t_{i,k-1}^gamma

## Define E_step function---------------------------------------------
E_step_M3 = function(theta.hat,data,mu_a.hat,sigma_a.hat,gammayes){
  # extract estimated parameters
  mu_b.hat = theta.hat$mu_b.hat
  sigma_b.hat = theta.hat$sigma_b.hat
  rho_ab.hat = theta.hat$rho_ab.hat
  sigma_D.hat = theta.hat$sigma_D.hat
  if (gammayes == 0){
    gamma.hat = 1
  }else{
    gamma.hat = theta.hat$gamma.hat
  }
  
  # extract data stored
  n = data$n
  m = data$m
  dy = data$dy # degradation increments
  y = data$y # degradation measurements
  y0 = data$y0 # initial measurements
  dt = data$dt # time interval before transformation
  t = data$t # inspection time points
  
  # calculate dtau
  if (gammayes == 0){
    dtau = dt # transformed time interval  
  }else{
    dtau = dt
    for (i in 1:n){
      for (k in 1:m){
        dtau[i,k] = t[i,k+1]^gamma.hat - t[i,k]^gamma.hat
      }
    }
  }
  
  # calculate up-to-date (transformed) time interval and degradation increments  
  dtau_total = rep(0,n)
  dy_total = rep(0,n)
  for (i in 1:n){
    dtau_total[i] = sum(dtau[i,1:m])
    dy_total[i] = sum(dy[i,1:m])
  }
  
  # calculate check.sigma and check.mu
  check.sigma = rep(0,n)
  check.mu = rep(0,n)
  for (i in 1:n){
    check.sigma[i] = 1/(1/sigma_b.hat+dtau_total[i]/sigma_D.hat+
      rho_ab.hat^2/(sigma_b.hat*(1-rho_ab.hat^2)))
    temp = rho_ab.hat*(y0[i]-mu_a.hat+sqrt(sigma_a.hat/sigma_b.hat)*rho_ab.hat*mu_b.hat)/
      (1-rho_ab.hat^2)/sqrt(sigma_a.hat*sigma_b.hat)
    check.mu[i] = check.sigma[i]*(mu_b.hat/sigma_b.hat+dy_total[i]/sigma_D.hat+temp)
  }
  
  # calculate Q function
  Q = 0
  for (i in 1:n){
    part3 = log(2*pi) + log(sigma_a.hat) + (y0[i]-mu_a.hat)^2/sigma_a.hat
    part2.1 = (1-rho_ab.hat^2)*sigma_b.hat
    part2.2 = (check.mu[i]-(mu_b.hat+sqrt(sigma_b.hat/sigma_a.hat)*rho_ab.hat*(y0[i]-mu_a.hat)))^2+check.sigma[i]
    part2 = log(2*pi) + log(part2.1)+part2.2/part2.1
    Q = Q + part2 + part3
    for (k in 1:m){
        part1 = log(2*pi) + log(dtau[i,k]*sigma_D.hat) + 
          (check.mu[i]^2+check.sigma[i]-2*dy[i,k]*check.mu[i]/dtau[i,k]+
                     (dy[i,k]/dtau[i,k])^2)/(sigma_D.hat/dtau[i,k])
        Q = Q + part1
    }
  }
  Q = -Q/2
  
  if (gammayes == 0){
    return(list(check.mu=check.mu,check.sigma=check.sigma,Q=Q))
  }else{
    return(list(check.mu=check.mu,check.sigma=check.sigma,
                gamma.hat=gamma.hat,Q=Q))
  }
}

## Define M_step function---------------------------------------------
M_step_M3 = function(E.result,data,mu_a.hat,sigma_a.hat,gammayes){
  # extract elements in the result of E_step
  check.mu = E.result$check.mu
  check.sigma = E.result$check.sigma
  
  if (gammayes == 0){
    gamma.hat = 1
  }else{
    gamma.hat = E.result$gamma.hat
  }
  
  # extract data stored
  n = data$n
  m = data$m
  dy = data$dy # degradation increments
  y = data$y # degradation measurements
  y0 = data$y0 # initial measurements
  dt = data$dt # time interval before transformation
  t = data$t # inspection time points
  
  # calculate dtau
  if (gammayes == 0){
    dtau = dt # transformed time interval  
  }else{
    dtau = dt
    for (i in 1:n){
      for (k in 1:m){
        dtau[i,k] = t[i,k+1]^gamma.hat - t[i,k]^gamma.hat
      }
    }
  }
  
  # update mu_b.hat
  mu_b.hat = mean(check.mu) 
  
  # update sigma_D.hat
  sigma_D.sum = 0
  for (i in 1:n){
    for (k in 1:m){
      sigma_D.sum = sigma_D.sum +
          dtau[i,k]*check.mu[i]^2+dtau[i,k]*check.sigma[i]-
          2*dy[i,k]*check.mu[i]+dy[i,k]^2/dtau[i,k]
    }
  }
  sigma_D.hat = sigma_D.sum/(n*m)
  
  # update sigma_ab.hat
  cross.sum = 0
  for (i in 1:n){
    cross.sum = cross.sum + (y0[i]-mu_a.hat)*(check.mu[i]-mu_b.hat)
  }
  sigma_ab.hat = cross.sum/(n-1)
  
  # update sigma_b.hat and rho_ab.hat
  Big.sum = 0
  for (i in 1:n){
    Big.sum = Big.sum + (check.sigma[i]+
                            (check.mu[i]-mu_b.hat-sigma_ab.hat*(y0[i]-mu_a.hat)/sigma_a.hat)^2)
  }
  sigma_b.hat = sigma_ab.hat^2/sigma_a.hat+Big.sum/n
  rho_ab.hat = sigma_ab.hat/sqrt(sigma_b.hat*sigma_a.hat)
  
  # update gamma.hat
  if (gammayes == 0){
    return(list(mu_b.hat=mu_b.hat,sigma_b.hat=sigma_b.hat,
                rho_ab.hat=rho_ab.hat,sigma_D.hat=sigma_D.hat))
  }else{
    gamma.func = function(gamma){
      tau = currentt^gamma - lastt^gamma
      tau_d = (currentt^gamma)*log(currentt) - (lastt^gamma)*log(lastt)
      index.NaN = which(is.nan(tau_d))
      tau_d[index.NaN] = (currentt[index.NaN]^gamma)*log(currentt[index.NaN])
      eqn = tau_d/tau + (tau_d*(check.mu.long^2+check.sigma.long)-tau_d*dy.long^2/tau^2)/sigma_D.hat.long
      
      return(sum(eqn))
    }
    currentt = NULL
    lastt = NULL
    check.mu.long = NULL
    check.sigma.long = NULL
    sigma_D.hat.long = NULL
    dy.long = NULL
    for (i in 1:n){
      currentt = append(currentt,t[i,2:(m+1)])
      lastt = append(lastt,t[i,1:m])
      check.mu.long = append(check.mu.long,rep(check.mu[i],m))
      check.sigma.long = append(check.sigma.long,rep(check.sigma[i],m))
      dy.long = append(dy.long,dy[i,])
    }
    sigma_D.hat.long = rep(sigma_D.hat,n*m)
    root = UnirootAll(gamma.func,interval=c(1e-08,5))
    index = which.min(abs(root-gamma.hat))
    gamma.hat = root[index]
    
    return(list(mu_b.hat=mu_b.hat,sigma_b.hat=sigma_b.hat,
                rho_ab.hat=rho_ab.hat,sigma_D.hat=sigma_D.hat,
                gamma.hat=gamma.hat))
  }
}

## Define the implementation function-------------------------------------------
implementEM_M3 = function(theta.hat,data,gammayes) {
  if (gammayes == 0){
    theta.last = list(mu_b.hat = 0,sigma_b.hat = 0,
                      rho_ab.hat = 0,sigma_D.hat = 0)
  }else{
    theta.last = list(mu_b.hat = 0,sigma_b.hat = 0,
                      rho_ab.hat = 0,sigma_D.hat = 0,
                      gamma.hat=0)
  }
  max.iter = 1e5
  cur.iter = 0
  theta.progress = as.data.frame(t(unlist(theta.last)))
  theta.progress$Q = 0
  theta.progress = theta.progress[-1,]

  while(max(abs(unlist(theta.hat) - unlist(theta.last))) > 1e-5 &&  cur.iter <= max.iter) {
    theta.last = theta.hat
    E.result = E_step_M3(theta.last,data,mu_a.hat,sigma_a.hat,gammayes)
    check.mu=E.result$check.mu
    Q = E.result$Q
    theta.progress = rbind(theta.progress,as.data.frame(cbind(t(unlist(theta.last)),Q=Q)))
    
    theta.hat = M_step_M3(E.result,data,mu_a.hat,sigma_a.hat,gammayes)
    progress(cur.iter,max.value = max.iter)
    cur.iter = cur.iter + 1
  }
  E.result = E_step_M3(theta.hat,data,mu_a.hat,sigma_a.hat,gammayes)
  check.mu = E.result$check.mu
  Q = E.result$Q
  theta.progress = rbind(theta.progress,as.data.frame(cbind(t(unlist(theta.hat)),Q=Q)))
  
  return(list(theta.progress=theta.progress,theta.hat=theta.hat,check.mu=check.mu,Q=Q))
}
















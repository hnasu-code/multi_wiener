## EM algorithm=======================================================

## The notations used in the code v.s. notations used in the manuscript
##         Used in code     --     Used in manuscript
##               n          --            n
##               p          --            p
##               m          --            m, identical across all units
##            mu_a.hat      --       \hat{\mu}_a in bold
##            mu_b.hat      --       \hat{\mu}_b in bold
##           sigma_a.hat    --       \hat{\sigma}_a^2 in bold
##           sigma_b.hat    --       \hat{\sigma}_b^2 in bold
##           rho_ab.hat     --       \hat{\rho}_ab in bold
##           sigma_D.hat    --       \hat{\sigma}^2 in bold
##           gamma.hat      --       \hat{\gamma} in bold, assuming tau=t_{ik}^gamma-t_{i,k-1}^gamma

## Define E_step function---------------------------------------------
E_step_M4 = function(theta.hat,data,gammayes,
                     mu_a.hat,sigma_a.hat){
  
  # extract estimated parameters
  mu_b.hat = theta.hat$mu_b.hat
  sigma_b.hat = theta.hat$sigma_b.hat
  rho_ab.hat = theta.hat$rho_ab.hat
  sigma_D.hat = theta.hat$sigma_D.hat
  if (gammayes == 0){
    gamma.hat = rep(1,p)
  }else{
    gamma.hat = theta.hat$gamma.hat
  }
  
  # extract data stored
  n = data$n
  p = data$p
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
      for (j in 1:p){
        for (k in 1:m){
          dtau[j,m*(i-1)+k] = t[j,m*(i-1)+i+k]^gamma.hat[j] - 
            t[j,m*(i-1)+i+k-1]^gamma.hat[j]
        }
      }
    }
  }
  
  # calculate up-to-date (transformed) time interval and degradation increments  
  dtau_total = matrix(rep(0,n*p),ncol=n)
  dy_total = matrix(rep(0,n*p),ncol=n)
  for (i in 1:n){
    for (j in 1:p){
      dtau_total[j,i] = sum(dtau[j,((i-1)*m+1):(i*m)])
      dy_total[j,i] = sum(dy[j,((i-1)*m+1):(i*m)])
    }
  }
  
  # calculate check.sigma
  sigma_b.hat.inv = 1/sigma_b.hat #part1
  sigma_part2 = rep(0,n) #part2
  for (i in 1:n){
    for (j in 1:p){
      sigma_part2[i] = sigma_part2[i] + dtau_total[j,i]*mu_b.hat[j]^2/sigma_D.hat[j]
    }
  }
  sigma_part3 = p*(rho_ab.hat^2)/(sigma_b.hat-sigma_b.hat*rho_ab.hat^2) #part3
  check.sigma = rep(0,n)
  for(i in 1:n){
    check.sigma[i] = 1/(sigma_b.hat.inv+sigma_part2[i]+sigma_part3)
  }
  
  # calculate check.mu
  mu_part2 = rep(0,n)  #part2
  mu_part3 = rep(0,n)  #part3
  mu_part3.1 = rho_ab.hat/(sqrt(sigma_a.hat*sigma_b.hat)*(1-rho_ab.hat^2))
  mu_part3.2 = rep(0,n)
  check.mu = rep(0,n)
  for (i in 1:n){
    for (j in 1:p){
      mu_part2[i] = mu_part2[i] + mu_b.hat[j]*dy_total[j,i]/sigma_D.hat[j]
      mu_part3.2[i] = mu_part3.2[i] + (y0[j,i]/mu_a.hat[j]-1+rho_ab.hat*sqrt(sigma_a.hat/sigma_b.hat))
    }
    mu_part3[i] = mu_part3.1*mu_part3.2[i]
    check.mu[i] = check.sigma[i]*(sigma_b.hat.inv+mu_part2[i]+mu_part3[i])
  }
  
  # calculate Q function
  Q = 0
  sigma_b_y0.hat = sigma_b.hat*(1-rho_ab.hat^2)
  for (i in 1:n){
    for (j in 1:p){
      part3 = log(2*pi) + log(sigma_a.hat*mu_a.hat[j]^2)+(y0[j,i]-mu_a.hat[j])^2/(sigma_a.hat*mu_a.hat[j]^2)
      part2.1 = (check.mu[i]-(1+rho_ab.hat*sqrt(sigma_b.hat/sigma_a.hat)*(y0[j,i]/mu_a.hat[j]-1)))^2
      part2 = log(2*pi) + log(sigma_b_y0.hat)+(part2.1+check.sigma[i])/sigma_b_y0.hat
      Q = Q + part2 + part3
      for (k in 1:m){
        part1 = log(2*pi) + log(dtau[j,m*(i-1)+k]*sigma_D.hat[j]) + 
          (check.mu[i]^2+check.sigma[i]-2*dy[j,m*(i-1)+k]*check.mu[i]/(dtau[j,m*(i-1)+k]*mu_b.hat[j])+
                     (dy[j,m*(i-1)+k]/dtau[j,m*(i-1)+k]/mu_b.hat[j])^2)/(sigma_D.hat[j]/(dtau[j,m*(i-1)+k]*mu_b.hat[j]^2))
        Q = Q + part1
      }
    }
  }
  Q = -Q/2
  
  if (gammayes == 0){
    return(list(check.mu=check.mu,check.sigma=check.sigma,Q=Q))
  }else{
    return(list(gamma.hat=gamma.hat,check.mu=check.mu,check.sigma=check.sigma,Q=Q))
  }
}

## Define M_step function---------------------------------------------
M_step_M4 = function(E.result,data,gammayes,
                     mu_a.hat,sigma_a.hat){
  
  # extract elements in the result of E_step
  check.mu = E.result$check.mu
  check.sigma = E.result$check.sigma
  if (gammayes == 0){
    gamma.hat = rep(1,p)
  }else{
    gamma.hat = E.result$gamma.hat
  }
  
  # extract data stored
  n = data$n
  p = data$p
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
      for (j in 1:p){
        for (k in 1:m){
          dtau[j,m*(i-1)+k] = t[j,m*(i-1)+i+k]^gamma.hat[j] - 
            t[j,m*(i-1)+i+k-1]^gamma.hat[j]
        }
      }
    }
  }
  
  # calculate up-to-date (transformed) time interval and degradation increments  
  dtau_total = matrix(rep(0,n*p),ncol=n)
  dy_total = matrix(rep(0,n*p),ncol=n)
  for (i in 1:n){
    for (j in 1:p){
      dtau_total[j,i] = sum(dtau[j,((i-1)*m+1):(i*m)])
      dy_total[j,i] = sum(dy[j,((i-1)*m+1):(i*m)])
    }
  }
  
  # update mu_b.hat
  up.mu.b = rep(0,p)
  down.mu.b = rep(0,p)
  mu_b.hat = rep(0,p)
  for (j in 1:p){
    for (i in 1:n){
      up.mu.b[j] = up.mu.b[j] + dy_total[j,i]*check.mu[i]
      down.mu.b[j] = down.mu.b[j] + dtau_total[j,i]*(check.mu[i]^2+check.sigma[i])
    }
    mu_b.hat[j] = up.mu.b[j]/down.mu.b[j]
  }
  
  # update rho_ab.hat
  sum.rho.ab = 0
  for (i in 1:n){
    for (j in 1:p){
      sum.rho.ab = sum.rho.ab+(y0[j,i]/mu_a.hat[j]-1)*(check.mu[i]-1)
    }
  }
  sigma_ab.hat = sum.rho.ab/(n*p-1)
  
  # update sigma_b.hat
  sum.sigma.b = 0
  for (i in 1:n){
    for (j in 1:p){
      sum.sigma.b = sum.sigma.b+check.sigma[i]+
        (check.mu[i]-1-sigma_ab.hat*(y0[j,i]/mu_a.hat[j]-1)/sigma_a.hat)^2
    }
  }
  sigma_b.hat = (sigma_ab.hat^2/sigma_a.hat)+sum.sigma.b/(n*p)
  rho_ab.hat = sigma_ab.hat/sqrt(sigma_a.hat*sigma_b.hat)
  
  # update sigma_D.hat
  sigma_D.sum = rep(0,p)
  for (i in 1:n){
    for (j in 1:p){
      for (k in 1:m){
        sigma_D.sum[j] = sigma_D.sum[j] +
          dtau[j,m*(i-1)+k]*mu_b.hat[j]^2*(check.mu[i]^2+check.sigma[i])-
             2*dy[j,m*(i-1)+k]*check.mu[i]*mu_b.hat[j]+dy[j,m*(i-1)+k]^2/dtau[j,m*(i-1)+k]
      }
    }
  }
  sigma_D.hat = sigma_D.sum/(n*m)
  
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
      eqn = tau_d/tau + (tau_d*mu_b.hat.long^2*(check.mu.long^2+check.sigma.long)-tau_d*dy.long^2/tau^2)/sigma_D.hat.long
      
      return(sum(eqn))
    }
    
    for (j in 1:p){
      currentt = NULL
      lastt = NULL
      mu_b.hat.long = NULL
      check.mu.long = NULL
      check.sigma.long = NULL
      sigma_D.hat.long = NULL
      for (i in 1:n){
        currentt = append(currentt,t[j,(m*(i-1)+i+1):(m*(i-1)+i+m)])
        lastt = append(lastt,t[j,(m*(i-1)+i):(m*(i-1)+i+m-1)])
        check.mu.long = append(check.mu.long,rep(check.mu[i],m))
        check.sigma.long = append(check.sigma.long,rep(check.sigma[i],m))
      }
      mu_b.hat.long = rep(mu_b.hat[j],n*m)
      sigma_D.hat.long = rep(sigma_D.hat[j],n*m)
      dy.long = dy[j,]
      root = UnirootAll(gamma.func,interval=c(1e-08,5))
      index = which.min(abs(root-gamma.hat[j]))
      gamma.hat[j] = root[index]
    }
    return(list(mu_b.hat=mu_b.hat,sigma_b.hat=sigma_b.hat,
                rho_ab.hat=rho_ab.hat,sigma_D.hat=sigma_D.hat,
                gamma.hat=gamma.hat))
  }
}

## Define the implementation function-------------------------------------------
implementEM_M4 = function(theta.hat,data,gammayes) {
  if (gammayes == 0){
    theta.last = list(mu_b.hat = rep(0,p),sigma_b.hat=0,
                      rho_ab.hat=0,sigma_D.hat = rep(0,p))
  }else{
    theta.last = list(mu_b.hat = rep(0,p),sigma_b.hat=0,
                      rho_ab.hat=0,sigma_D.hat = rep(0,p),
                      gamma.hat=rep(0,p))
  }
  max.iter = 1e5
  cur.iter = 0
  theta.progress = as.data.frame(t(unlist(theta.last)))
  theta.progress$Q = 0
  theta.progress = theta.progress[-1,]

  while(max(abs(unlist(theta.hat) - unlist(theta.last))) > 1e-5 &&  cur.iter <= max.iter) {
    theta.last = theta.hat
    E.result = E_step_M4(theta.hat,data,gammayes,
                         mu_a.hat,sigma_a.hat)
    check.mu=E.result$check.mu
    Q = E.result$Q
    theta.progress = rbind(theta.progress,as.data.frame(cbind(t(unlist(theta.last)),Q=Q)))
    
    theta.hat = M_step_M4(E.result,data,gammayes,
                          mu_a.hat,sigma_a.hat)
    progress(cur.iter,max.value = max.iter)
    cur.iter = cur.iter + 1
  }
  E.result = E_step_M4(theta.hat,data,gammayes,
                       mu_a.hat,sigma_a.hat)
  check.mu=E.result$check.mu
  Q = E.result$Q
  theta.progress = rbind(theta.progress,as.data.frame(cbind(t(unlist(theta.hat)),Q=Q)))
  
  return(list(theta.progress=theta.progress,theta.hat=theta.hat,check.mu=check.mu,Q=Q))
}
















## EM algorithm=======================================================

## The notations used in the code v.s. notations used in the manuscript
##         Used in code     --     Used in manuscript
##               n          --            n
##               p          --            p
##               m          --            m, identical across all units
##            mu_a.hat      --       \hat{\mu}_a in bold
##            mu_b.hat      --       \hat{\mu}_b in bold
##           Sigma_a.hat    --       \hat{\Sigma}_a in bold
##           Sigma_b.hat    --       \hat{\Sigma}_b in bold
##           Sigma_ab.hat   --       \hat{\Sigma}_ab in bold
##           sigma_D.hat    --       \hat{\sigma}^2 in bold
##           gamma.hat      --       \hat{\gamma} in bold, assuming tau=t_{ik}^gamma-t_{i,k-1}^gamma

## Define E_step function---------------------------------------------
E_step_M0 = function(theta.hat,data,gammayes,
                     mu_a.hat,Sigma_a.hat,Sigma_a.hat.inv,
                     sigma_a.hat,rho_a.hat){
  
  # extract estimated parameters
  mu_b.hat = theta.hat$mu_b.hat
  Sigma_b.hat = theta.hat$Sigma_b.hat
  Sigma_ab.hat = theta.hat$Sigma_ab.hat
  Sigma_ba.hat = t(Sigma_ab.hat)
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
  
  # calculate check.Sigma
  Sigma_b.hat.inv = solve(Sigma_b.hat) #part1
  D.hat.inv = solve(diag(sigma_D.hat))
  tD.hat.inv = as.matrix(kronecker(t(rep(1,n)),diag(p))) #part2
  for (i in 1:n){
    tD.hat.inv[1:p,((i-1)*p+1):(i*p)] = diag(dtau_total[1:p,i])%*%D.hat.inv
  }
  Sigma_Y0_b.hat.inv = solve(Sigma_a.hat-Sigma_ab.hat%*%Sigma_b.hat.inv%*%Sigma_ba.hat)
  Sigma_part3 = Sigma_b.hat.inv%*%Sigma_ba.hat%*%Sigma_Y0_b.hat.inv%*%Sigma_ab.hat%*%Sigma_b.hat.inv #part3
  
  check.Sigma = as.matrix(kronecker(t(rep(1,n)),diag(p)))
  for(i in 1:n){
    check.Sigma[1:p,((i-1)*p+1):(i*p)] = solve(Sigma_b.hat.inv+tD.hat.inv[1:p,((i-1)*p+1):(i*p)]+Sigma_part3)
  }
  
  # calculate check.mu
  mu_part1 = Sigma_b.hat.inv%*%mu_b.hat #part1
  mu_part2 = matrix(rep(0,n*p),ncol=n)  #part2
  mu_part3 = matrix(rep(0,n*p),ncol=n)  #part3
  check.mu = matrix(rep(0,n*p),ncol=n)
  for (i in 1:n){
    mu_part2[1:p,i] = D.hat.inv%*%dy_total[1:p,i]
    mu_part3[1:p,i] = Sigma_b.hat.inv%*%Sigma_ba.hat%*%Sigma_Y0_b.hat.inv%*%
      (y0[1:p,i]-mu_a.hat+Sigma_ab.hat%*%Sigma_b.hat.inv%*%mu_b.hat)
    check.mu[1:p,i] = check.Sigma[1:p,((i-1)*p+1):(i*p)]%*%(mu_part1+mu_part2[1:p,i]+mu_part3[1:p,i])
  }
  
  # calculate Q function
  Q = 0
  Sigma_b_y0.hat = Sigma_b.hat-Sigma_ba.hat%*%Sigma_a.hat.inv%*%Sigma_ab.hat
  Sigma_b_y0.hat.inv = solve(Sigma_b_y0.hat)
  for (i in 1:n){
    mu_b_y0.hat = mu_b.hat+Sigma_ba.hat%*%Sigma_a.hat.inv%*%(y0[1:p,i]-mu_a.hat)
    part3 = p*log(2*pi) + log(det(Sigma_a.hat)) + as.numeric(t(y0[1:p,i]-mu_a.hat)%*%Sigma_a.hat.inv%*%(y0[1:p,i]-mu_a.hat))
    part2 = p*log(2*pi) + log(det(Sigma_b_y0.hat)) + matrix.trace(Sigma_b_y0.hat.inv%*%check.Sigma[1:p,((i-1)*p+1):(i*p)])+
      as.numeric(t(check.mu[1:p,i]-mu_b_y0.hat)%*%
                           Sigma_b_y0.hat.inv%*%(check.mu[1:p,i]-mu_b_y0.hat))
    Q = Q + part2 + part3
    for (j in 1:p){
      for (k in 1:m){
        part1 = log(2*pi) + log(dtau[j,m*(i-1)+k]*sigma_D.hat[j]) + 
          (check.mu[j,i]^2+check.Sigma[j,(i-1)*p+j]-2*dy[j,m*(i-1)+k]*check.mu[j,i]/dtau[j,m*(i-1)+k]+
                     (dy[j,m*(i-1)+k]/dtau[j,m*(i-1)+k])^2)/(sigma_D.hat[j]/dtau[j,m*(i-1)+k])
        Q = Q + part1
      }
    }
  }
  Q = -Q/2
  
  if (gammayes == 0){
    return(list(check.mu=check.mu,check.Sigma=check.Sigma,Q=Q))
  }else{
    return(list(gamma.hat=gamma.hat,check.mu=check.mu,check.Sigma=check.Sigma,Q=Q))
  }
}

## Define M_step function---------------------------------------------
M_step_M0 = function(E.result,data,gammayes,
                     mu_a.hat,Sigma_a.hat,Sigma_a.hat.inv,
                     sigma_a.hat,rho_a.hat){
  
  # extract elements in the result of E_step
  check.mu = E.result$check.mu
  check.Sigma = E.result$check.Sigma
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
  
  # update mu_b.hat
  mu_b.hat = rowMeans(check.mu) 
  
  # update Sigma_ab.hat
  Sigma_ab.hat = cov(t(y0),t(check.mu))
  Sigma_ba.hat = t(Sigma_ab.hat)
  
  # update Sigma_b.hat
  Sigma_b.sum = matrix(rep(0,p*p),ncol=p)
  for (i in 1:n){
    Sigma_b.sum = Sigma_b.sum + check.Sigma[1:p,((i-1)*p+1):(i*p)] + 
      (check.mu[1:p,i]-mu_b.hat-Sigma_ba.hat%*%Sigma_a.hat.inv%*%(y0[1:p,i]-mu_a.hat))%*%
      t(check.mu[1:p,i]-mu_b.hat-Sigma_ba.hat%*%Sigma_a.hat.inv%*%(y0[1:p,i]-mu_a.hat))
  }
  Sigma_b.hat = Sigma_ba.hat%*%Sigma_a.hat.inv%*%Sigma_ab.hat + Sigma_b.sum/n  
  sigma_b.hat = diag(Sigma_b.hat)
  rho_b.hat = NULL
  for (j in 1:(p-1)){
    rho_b.hat=append(rho_b.hat,cov2cor(Sigma_b.hat)[(j+1):p,j])
  }
  rho_ab.hat = NULL
  for (j in 1:p){
    for (i in 1:p){
      rho_ab.hat = append(rho_ab.hat,Sigma_ab.hat[i,j]/sqrt(sigma_a.hat[i]*sigma_b.hat[j]))  
    }
  }
  
  # update sigma_D.hat
  sigma_D.sum = rep(0,p)
  for (i in 1:n){
    for (j in 1:p){
      for (k in 1:m){
        sigma_D.sum[j] = sigma_D.sum[j] +
          dtau[j,m*(i-1)+k]*check.mu[j,i]^2+dtau[j,m*(i-1)+k]*check.Sigma[j,(i-1)*p+j]-
             2*dy[j,m*(i-1)+k]*check.mu[j,i]+dy[j,m*(i-1)+k]^2/dtau[j,m*(i-1)+k]
      }
    }
  }
  sigma_D.hat = sigma_D.sum/(n*m)
  
  # update gamma.hat
  if (gammayes == 0){
    return(list(mu_b.hat=mu_b.hat,sigma_b.hat=sigma_b.hat,rho_b.hat=rho_b.hat,
                rho_ab.hat=rho_ab.hat,sigma_D.hat=sigma_D.hat,
                Sigma_b.hat = Sigma_b.hat,Sigma_ab.hat=Sigma_ab.hat))
  }else{
    gamma.func = function(gamma){
      tau = currentt^gamma - lastt^gamma
      tau_d = (currentt^gamma)*log(currentt) - (lastt^gamma)*log(lastt)
      index.NaN = which(is.nan(tau_d))
      tau_d[index.NaN] = (currentt[index.NaN]^gamma)*log(currentt[index.NaN])
      eqn = tau_d/tau + (tau_d*(check.mu.long^2+check.Sigma.long)-tau_d*dy.long^2/tau^2)/sigma_D.hat.long
      
      return(sum(eqn))
    }
    
    for (j in 1:p){
      currentt = NULL
      lastt = NULL
      check.mu.long = NULL
      check.Sigma.long = NULL
      sigma_D.hat.long = NULL
      for (i in 1:n){
        currentt = append(currentt,t[j,(m*(i-1)+i+1):(m*(i-1)+i+m)])
        lastt = append(lastt,t[j,(m*(i-1)+i):(m*(i-1)+i+m-1)])
        check.mu.long = append(check.mu.long,rep(check.mu[j,i],m))
        check.Sigma.long = append(check.Sigma.long,rep(check.Sigma[j,(i-1)*p+j],m))
      }
      sigma_D.hat.long = rep(sigma_D.hat[j],n*m)
      dy.long = dy[j,]
      root = UnirootAll(gamma.func,interval=c(1e-08,5))
      index = which.min(abs(root-gamma.hat[j]))
      gamma.hat[j] = root[index]
    }
    return(list(mu_b.hat=mu_b.hat,sigma_b.hat=sigma_b.hat,rho_b.hat=rho_b.hat,
                rho_ab.hat=rho_ab.hat,sigma_D.hat=sigma_D.hat,gamma.hat=gamma.hat,
                Sigma_b.hat = Sigma_b.hat,Sigma_ab.hat=Sigma_ab.hat))
  }
}

## Define the implementation function-------------------------------------------
implementEM_M0 = function(theta.hat,data,gammayes) {
  if (gammayes == 0){
    theta.last = list(mu_b.hat = rep(0,p),sigma_b.hat=rep(0,p),
                      rho_b.hat=rep(0,p*(p-1)/2),rho_ab.hat=rep(0,p^2),
                      sigma_D.hat = rep(0,p),
                      Sigma_b.hat = matrix(rep(0,p*p),ncol=p),
                      Sigma_ab.hat = matrix(rep(0,p*p),ncol=p))
  }else{
    theta.last = list(mu_b.hat = rep(0,p),sigma_b.hat=rep(0,p),
                      rho_b.hat=rep(0,p*(p-1)/2),rho_ab.hat=rep(0,p^2),
                      sigma_D.hat = rep(0,p),gamma.hat=rep(0,p),
                      Sigma_b.hat = matrix(rep(0,p*p),ncol=p),
                      Sigma_ab.hat = matrix(rep(0,p*p),ncol=p))
  }
  len = length(theta.last)
  max.iter = 1e5
  cur.iter = 0
  theta.progress = as.data.frame(t(unlist(theta.last)))
  theta.progress$Q = 0
  theta.progress = theta.progress[-1,]

  while(max(abs(unlist(theta.hat[-c(len-1,len)]) - unlist(theta.last[-c(len-1,len)]))) > 1e-5 &&  cur.iter <= max.iter) {
    theta.last = theta.hat
    E.result = E_step_M0(theta.last,data,gammayes,
                         mu_a.hat,Sigma_a.hat,Sigma_a.hat.inv,
                         sigma_a.hat,rho_a.hat)
    check.mu = E.result$check.mu
    check.Sigma = E.result$check.Sigma
    Q = E.result$Q
    theta.progress = rbind(theta.progress,as.data.frame(cbind(t(unlist(theta.last)),Q=Q)))
    
    theta.hat = M_step_M0(E.result,data,gammayes,
                          mu_a.hat,Sigma_a.hat,Sigma_a.hat.inv,
                          sigma_a.hat,rho_a.hat)
    progress(cur.iter,max.value = max.iter)
    cur.iter = cur.iter + 1
  }
  E.result = E_step_M0(theta.hat,data,gammayes,
                       mu_a.hat,Sigma_a.hat,Sigma_a.hat.inv,
                       sigma_a.hat,rho_a.hat)
  check.mu = E.result$check.mu
  check.Sigma = E.result$check.Sigma
  Q = E.result$Q
  theta.progress = rbind(theta.progress,as.data.frame(cbind(t(unlist(theta.hat)),Q=Q)))
  
  return(list(theta.progress=theta.progress,theta.hat=theta.hat,
              check.mu=check.mu,check.Sigma=check.Sigma,Q=Q))
}



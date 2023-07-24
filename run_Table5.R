cat("\014")
rm(list=ls())
# ensure all the packages below are installed
library(statmod)
library(stats4)
library(reshape2)
library(matrixcalc)
library(DescTools)
library(svMisc)
library(MASS)
library(DataCombine)
library(dplyr)
library(Matrix)
library(matrixStats)
library(gridExtra)
library(cubature)
library(plyr)
library(tidyr)
library(mvtnorm)
library(rstudioapi)

# import data=========================================================
# Set working directory to source file location
setwd(dirname(getActiveDocumentContext()$path))   
LED = read.csv(file="IRLED.csv",header=T)
LED$Unit = as.factor(LED$Unit)
LED$PC = as.factor(LED$PC) 
LED$BurnIn = as.factor(LED$BurnIn)
LED_burnin = LED[LED$BurnIn=='No',]
rownames(LED_burnin) = NULL
n = as.numeric(length(unique(LED_burnin$Unit)))
p = as.numeric(length(unique(LED_burnin$PC)))
m = as.numeric(length(unique(LED_burnin$Current_T)))-1
LED_burnin[,2] = LED_burnin[,2] - 6
for (i in 1:(n*p)){
  LED_burnin[((i-1)*(m+1)+2):(i*(m+1)),3] = LED_burnin[((i-1)*(m+1)+1):(i*(m+1)-1),2]
  LED_burnin[(i-1)*(m+1)+1,c(3,4,6,7)] = c(NA,NA,NA,NA)
}
dy = NULL
y = NULL
y0 = NULL
dt = NULL
t = NULL
for (i in 1:n){
  dy_temp = rbind(t(na.omit(LED_burnin$Delta_Y[LED_burnin$PC=="PC1" & LED_burnin$Unit == i])),
                  t(na.omit(LED_burnin$Delta_Y[LED_burnin$PC=="PC2" & LED_burnin$Unit == i])))
  y_temp = rbind(t(LED_burnin$Current_Y[LED_burnin$PC=="PC1" & LED_burnin$Unit == i]),
                 t(LED_burnin$Current_Y[LED_burnin$PC=="PC2" & LED_burnin$Unit == i]))
  y0_temp = rbind(t(LED_burnin$Current_Y[LED_burnin$PC=="PC1" & LED_burnin$Current_T==0 & LED_burnin$Unit == i]),
                  t(LED_burnin$Current_Y[LED_burnin$PC=="PC2" & LED_burnin$Current_T==0 & LED_burnin$Unit == i]))
  dt_temp = rbind(as.numeric(t(na.omit(LED_burnin$Delta_T[LED_burnin$PC=="PC1" & LED_burnin$Unit == i]))/10),
                  as.numeric(t(na.omit(LED_burnin$Delta_T[LED_burnin$PC=="PC2" & LED_burnin$Unit == i]))/10))
  t_temp = rbind(as.numeric(t(LED_burnin$Current_T[LED_burnin$PC=="PC1" & LED_burnin$Unit == i])/10),
                 as.numeric(t(LED_burnin$Current_T[LED_burnin$PC=="PC2" & LED_burnin$Unit == i])/10))
  dy = cbind(dy,dy_temp)
  y = cbind(y,y_temp)
  y0 = cbind(y0,y0_temp)
  dt = cbind(dt,dt_temp)
  t = cbind(t,t_temp)
}
data = list(n=n,p=p,m=m,dy=dy,y=y,y0=y0,dt=dt,t=t)

# parameter estimation with model M0===========================================
mu_a.hat = rowMeans(data$y0)
Sigma_a.hat = cov(t(data$y0))
Sigma_a.hat.inv = solve(Sigma_a.hat)
sigma_a.hat = diag(Sigma_a.hat)
rho_a.hat = NULL
for (j in 1:(p-1)){
  rho_a.hat=append(rho_a.hat,cov2cor(Sigma_a.hat)[(j+1):p,j])
}

## fit simple Wiener processes and get starting values
source("func_start.R")
est = matrix(rep(0,n*p*3),ncol=p*3)
colnames(est) = c("mu_b_PC1","sigma_D_PC1","gam_PC1",
                  "mu_b_PC2","sigma_D_PC2","gam_PC2")
for (i in 1:n){
  for (j in 1:p){
    ct = data$t[j,(m*(i-1)+i+1):(m*(i-1)+i+m)]
    lt = data$t[j,(m*(i-1)+i):(m*(i-1)+i+m-1)]
    deltay = data$dy[j,((i-1)*m+1):(i*m)]
    fit = suppressWarnings(mle(simpleWiener_tran,start=list(mu=2,sigma=3,gam=1),
                               method="Nelder-Mead"))
    est[i,((j-1)*3+1):(j*3)] = as.numeric(coef(fit))
  }
}
start.mu_b = as.numeric(colMeans(est[,c(1,4)]))
start.sigma_D = as.numeric(colMeans(est[,c(2,5)]))
start.gam = as.numeric(colMeans(est[,c(3,6)]))
start.Sigma_b = unname(as.matrix(nearPD(cov(est[,c(1,4)]))$mat))
start.Sigma_ab = unname(as.matrix(cov(t(data$y0),est[,c(1,4)])))
start.sigma_b = diag(start.Sigma_b)
start.rho_b = NULL
for (j in 1:(p-1)){
  start.rho_b=append(start.rho_b,cov2cor(start.Sigma_b)[(j+1):p,j])
}
start.rho_ab = NULL
for (j in 1:p){
  for (i in 1:p){
    start.rho_ab = append(start.rho_ab,start.Sigma_ab[i,j]/sqrt(sigma_a.hat[i]*start.sigma_b[j]))  
  }
}

## fit the proposed model
source("func_M0.R")
theta.hat = list(mu_b.hat=start.mu_b,sigma_b.hat=start.sigma_b,rho_b.hat=start.rho_b,
                 rho_ab.hat=start.rho_ab,sigma_D.hat=start.sigma_D,gamma.hat=start.gam,
                 Sigma_b.hat = start.Sigma_b,Sigma_ab.hat=start.Sigma_ab)
len = length(theta.hat)
theta = implementEM_M0(theta.hat,data,1)
log_likelihood = theta$Q
AIC = 2*(13+5)-2*log_likelihood

cat("Model M0:","log-likehood:",log_likelihood,"AIC:", AIC)

# parameter estimation with model M1===========================================
mu_a.hat = rowMeans(y0)
Sigma_a.hat = cov(t(y0))
Sigma_a.hat.inv = solve(Sigma_a.hat)
sigma_a.hat = diag(Sigma_a.hat)
rho_a.hat = NULL
for (j in 1:(p-1)){
  rho_a.hat=append(rho_a.hat,cov2cor(Sigma_a.hat)[(j+1):p,j])
}
data = list(n=n,p=p,m=m,dy=dy,dt=dt,t=t)

## fit simple Wiener processes and get starting values
source("func_start.R")
est = matrix(rep(0,n*p*3),ncol=p*3)
colnames(est) = c("mu_PC1","sigma_PC1","gam_PC1",
                  "mu_PC2","sigma_PC2","gam_PC2")
for (i in 1:n){
  for (j in 1:p){
    ct = data$t[j,(m*(i-1)+i+1):(m*(i-1)+i+m)]
    lt = data$t[j,(m*(i-1)+i):(m*(i-1)+i+m-1)]
    deltay = data$dy[j,((i-1)*m+1):(i*m)]
    fit = suppressWarnings(mle(simpleWiener_tran,start=list(mu=2,sigma=3,gam=1),method="Nelder-Mead"))
    est[i,((j-1)*3+1):(j*3)] = as.numeric(coef(fit))
  }
}
start.mu_b = as.numeric(colMeans(est[,c(1,4)]))
start.sigma_D = as.numeric(colMeans(est[,c(2,5)]))
start.gam = as.numeric(colMeans(est[,c(3,6)]))
start.Sigma_b = unname(as.matrix(nearPD(var(est[,c(1,4)]))$mat))
start.sigma_b = diag(start.Sigma_b)
start.rho_b = NULL
for (j in 1:(p-1)){
  start.rho_b=append(start.rho_b,cov2cor(start.Sigma_b)[(j+1):p,j])
}

## fit the proposed model
source("func_M2.R")
theta.hat = list(mu_b.hat=start.mu_b,sigma_b.hat=start.sigma_b,
                 rho_b.hat=start.rho_b,sigma_D.hat=start.sigma_D,
                 gamma.hat=start.gam,Sigma_b.hat=start.Sigma_b)
theta = implementEM_M2(theta.hat,data,gammayes=1)
Q=0
for (i in 1:n){
  Q = Q + dmvnorm(y0[,i],mean=mu_a.hat,sigma=Sigma_a.hat,log=T)
}
log_likelihood = theta$Q + Q
AIC = 2*(13+5-4)-2*log_likelihood

cat("Model M1:","log-likehood:",log_likelihood,"AIC:", AIC)

# parameter estimation with model M3===========================================
data1 = list(n=n,m=m,dy=t(matrix(dy[1,],nrow=m)),y=t(matrix(y[1,],nrow=m+1)),
             y0=y0[1,],dt=t(matrix(dt[1,],nrow=m)),t=t(matrix(t[1,],nrow=m+1)))
data2 = list(n=n,m=m,dy=t(matrix(dy[2,],nrow=m)),y=t(matrix(y[2,],nrow=m+1)),
             y0=y0[2,],dt=t(matrix(dt[2,],nrow=m)),t=t(matrix(t[2,],nrow=m+1)))

## fit simple Wiener processes and get starting values
source("func_start.R")
est1 = matrix(rep(0,n*3),ncol=3)
est2 = matrix(rep(0,n*3),ncol=3)
colnames(est1) = c("mu_b","sigma_D","gam")
colnames(est2) = c("mu_b","sigma_D","gam")
for (i in 1:n){
  ct = data1$t[i,2:(m+1)]
  lt = data1$t[i,1:m]
  deltay = data1$dy[i,1:m]
  est1[i,1:3] = suppressWarnings(as.numeric(coef(mle(simpleWiener_tran,start=list(mu=5,sigma=1,gam=1),
                                    method = "Nelder-Mead"))))
  
  ct = data2$t[i,2:(m+1)]
  lt = data2$t[i,1:m]
  deltay = data2$dy[i,1:m]
  est2[i,1:3] = suppressWarnings(as.numeric(coef(mle(simpleWiener_tran,start=list(mu=5,sigma=1,gam=1),
                                    method = "Nelder-Mead"))))
}
start.mu_b1 = as.numeric(mean(est1[,1]))
start.mu_b2 = as.numeric(mean(est2[,1]))
start.sigma_D1 = as.numeric(mean(est1[,2]))
start.sigma_D2 = as.numeric(mean(est2[,2]))
start.gam1 = as.numeric(mean(est1[,3]))
start.gam2 = as.numeric(mean(est2[,3]))
start.sigma_b1 = var(as.numeric(est1[,1]))
start.sigma_b2 = var(as.numeric(est2[,1]))
start.rho_ab1 = cov2cor(cov(cbind(data1$y0,as.numeric(est1[,1]))))[1,2]
start.rho_ab2 = cov2cor(cov(cbind(data2$y0,as.numeric(est2[,1]))))[1,2]

## fit the proposed model
source("func_M3.R")
theta.hat1 = list(mu_b.hat = start.mu_b1,sigma_b.hat = start.sigma_b1,
                  rho_ab.hat = start.rho_ab1,sigma_D.hat = start.sigma_D1,
                  gamma.hat=start.gam1)
mu_a.hat = mean(data1$y0)
sigma_a.hat = var(data1$y0)
theta1 = implementEM_M3(theta.hat1,data1,gammayes=1)

theta.hat2 = list(mu_b.hat = start.mu_b2,sigma_b.hat = start.sigma_b2,
                  rho_ab.hat = start.rho_ab2,sigma_D.hat = start.sigma_D2,
                  gamma.hat=start.gam2)
mu_a.hat = mean(data2$y0)
sigma_a.hat = var(data2$y0)
theta2 = implementEM_M3(theta.hat2,data2,gammayes=1)
log_likelihood = theta1$Q + theta2$Q
AIC = 2*(length(theta.hat1)+length(theta.hat2)+4)-2*log_likelihood

cat("Model M3:","log-likehood:",log_likelihood,"AIC:", AIC)

# parameter estimation with model M4===========================================
data = list(n=n,p=p,m=m,dy=dy,y=y,y0=y0,dt=dt,t=t)
mu_a.hat = rowMeans(data$y0)
sum.sigma.a = 0
for (i in 1:n){
  for (j in 1:p){
    sum.sigma.a = sum.sigma.a+(data$y0[j,i]/mu_a.hat[j]-1)^2
  }
}
sigma_a.hat = sum.sigma.a/(n*p-1)

## fit simple Wiener processes and get starting values
source("func_start.R")
est = matrix(rep(0,n*p*3),ncol=p*3)
colnames(est) = c("mu_b_PC1","sigma_D_PC1","gam_PC1",
                  "mu_b_PC2","sigma_D_PC2","gam_PC2")
for (i in 1:n){
  for (j in 1:p){
    ct = data$t[j,(m*(i-1)+i+1):(m*(i-1)+i+m)]
    lt = data$t[j,(m*(i-1)+i):(m*(i-1)+i+m-1)]
    deltay = data$dy[j,((i-1)*m+1):(i*m)]
    est[i,((j-1)*3+1):(j*3)] = suppressWarnings(as.numeric(coef(mle(simpleWiener_tran,
                                                   start=list(mu=2,sigma=3,gam=1),method = "Nelder-Mead"))))
  }
}
start.mu_b = as.numeric(colMeans(est[,c(1,4)]))
start.sigma_D = as.numeric(colMeans(est[,c(2,5)]))
start.gam = as.numeric(colMeans(est[,c(3,6)]))
start.b = t(est[,c(1,4)])
sum.sigma.b = 0
for (i in 1:n){
  for (j in 1:p){
    sum.sigma.b = sum.sigma.b+(start.b[j,i]/start.mu_b[j]-1)^2
  }
}
start.sigma_b = as.numeric(sum.sigma.b)/(n*p-1)

sum.rho.ab = 0
for (i in 1:n){
  for (j in 1:p){
    sum.rho.ab = sum.rho.ab+(data$y0[j,i]/mu_a.hat[j]-1)*(start.b[j,i]/start.mu_b[j]-1)
  }
}
start.rho_ab = as.numeric(sum.rho.ab)/(n*p-1)/sqrt(sigma_a.hat*start.sigma_b)

## fit the proposed model
source("func_M4.R")
theta.hat = list(mu_b.hat=start.mu_b,sigma_b.hat=start.sigma_b,
                 rho_ab.hat=start.rho_ab,sigma_D.hat=start.sigma_D,
                 gamma.hat = start.gam)
theta = implementEM_M4(theta.hat,data,gammayes=1)
log_likelihood = theta$Q
AIC = 2*11-2*log_likelihood

cat("Model M4:","log-likehood:",log_likelihood,"AIC:", AIC)

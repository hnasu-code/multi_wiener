cat("\014")
rm(list=ls())
library(statmod)
library(stats4)
library(reshape2)
library(matrixcalc)
library(svMisc)
library(DescTools)
library(MASS)
library(mvtnorm)
library(MBESS)
library(MTS)
library(Matrix)
library(matrixStats)
library(ggplot2)
library(tidyr)
library(latex2exp)
library(ggalt)
library(gridExtra)
## model set-up=================================================================
n = 10
p = 3 
m = 15
interval = 1
gamma = c(1,1,1)
mu_a = c(3,4,5)
mu_b = c(5,4,3)
sigma_D = c(4,5,6)

Vech_Sigma_a = c(1, 0.5, 0.5, 1, 0.5, 1)
Sigma_a = VechM(Vech_Sigma_a)
sigma_a = diag(Sigma_a)
rho_a = NULL
for (j in 1:(p-1)){
  rho_a=append(rho_a,cov2cor(Sigma_a)[(j+1):p,j])
}

Vech_Sigma_b = c(1, 0.8, 0.8, 1, 0.8, 1)
Sigma_b = VechM(Vech_Sigma_b)
sigma_b = diag(Sigma_b)
rho_b = NULL
for (j in 1:(p-1)){
  rho_b=append(rho_b,cov2cor(Sigma_b)[(j+1):p,j])
}

Vec_Sigma_ab = c(0.3, 0.1, 0.1, 0.1, 0.3, 0.1, 0.1, 0.1, 0.3)
Sigma_ab = matrix(Vec_Sigma_ab,nrow=3)
rho_ab = NULL
for (j in 1:p){
  for (i in 1:p){
    rho_ab = append(rho_ab,Sigma_ab[i,j]/sqrt(sigma_a[i]*sigma_b[j]))  
  }
}

source("func_simu.R")

## plot degradation paths=======================================================
data = data_generation_M0(12,n,p,m,interval,gamma,mu_a,mu_b,Sigma_a,Sigma_b,Sigma_ab,sigma_D)
data_plot = data.frame(matrix(ncol = 4, nrow = 0))
colnames(data_plot) = c("Unit","PC","Current_T","Current_Y")
for (i in 1:n){
  for (j in 1:p){
    Unit = rep(i,m+1)
    PC = rep(paste(c("PC", j), collapse = ""),m+1)
    Current_T = data$t[j,((i-1)*m+i):(i*m+i)]
    Current_Y = data$y[j,((i-1)*m+i):(i*m+i)]
    
    add = data.frame(Unit,PC,Current_T,Current_Y)
    data_plot = rbind(data_plot,add)
  }
}
data_plot$Unit = as.factor(data_plot$Unit)
data_plot$PC = as.factor(data_plot$PC)

theme_set(theme_bw())
a = ggplot(data_plot,aes(x=Current_T,y=Current_Y,group=Unit:PC,colour=Unit,cex.lab=2))+
  geom_line()+geom_point(size=0.5,shape=2)+
  labs(y="Responses",x="Time")+
  facet_grid(.~PC)+theme(panel.spacing.x = unit(1.1, "lines"),
                         legend.position="bottom",legend.box = "horizontal")+
  guides(colour = guide_legend(nrow = 1))+
  theme(text = element_text(size = 18))
a

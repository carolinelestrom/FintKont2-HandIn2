setwd("/Users/carolinelestrom/Documents/KU/FinKont2/Handins/HandIn2")

library(png) ### Plots
library(knitr) ### Data manipulations
library(tidyverse) ### Data manipulations
library(dtplyr) ### Data manipulations - merge datasets
library(ggplot2) ### Plots
library(gridExtra) ### Plots






# Q2 - Options on Coupon-Bearing Bonds

## Q2.c

### INITIALIZE
### =========

r0 <- 0.02
thetaQ <- 0.05
kappa <- 0.1
sigma <- 0.015
timetomat <- 1
N <- 5
alpha <- rep(1, N)  ### alpha_i = 1 for all i
Ti <- seq(1, N) + timetomat  ### Ti = i + 1
K <- 4.5


### FUNCTIONS
### =========

### Vasicek B-function in "drift = kappa(thetaQ - x)"-notation
Bvasicek <- function(t, T) {
  return((1 - exp(-kappa * (T - t))) / kappa)
}

### Vasicek B-function in "drift = kappa(thetaQ - x)"-notation
Avasicek <- function(t, T) {
  return((thetaQ - 0.5 * (sigma / kappa) ^ 2) * (Bvasicek(t, T) - (T - t)) - (sigma ^ 2) / (4 * kappa) * Bvasicek(t, T) ^ 2)
}

### Vasicek ZCB prices for ZCB's all starting at time T=1 at maturing at time Ti>T
ZCBpricevasicek <- function(T, Ti, r_star) {
  B <- Bvasicek(T, Ti)
  A <- Avasicek(T, Ti)
  return(sum(exp(A - B * r_star)))
}

### Equation to solve for r*
r_star_solution <- function(r_star, K) {
  return(ZCBpricevasicek(timetomat, Ti, r_star) - K)
}

### Vasicek ZCB prices
ZCBvasicek<-function(r,tau,Pparam,riskpremium=0)
{ thetaQ<-Pparam[1]+riskpremium
kap<-Pparam[2]
sig<-Pparam[3]

Btau<-(1-exp(-kap*tau))/kap
Atau<-((Btau-tau)*(kap^2*thetaQ-0.5*sig^2)/kap^2 - sig^2*Btau^2/(4*kap))

ZCBvasicek<-exp(-r*Btau+Atau) 

}

### Vasicek time-0 put price
PUTvasicek<-function(r,tauPUT,tauZCB,strike,Pparam,riskpremium=0){ 
  kap<-Pparam[2]
  sig<-Pparam[3]
  
  PtauPUT<-ZCBvasicek(r,tauPUT,Pparam,riskpremium)
  PtauZCB<-ZCBvasicek(r,tauZCB,Pparam,riskpremium)
  
  sigZCB<-sig*(1-exp(-kap*(tauZCB-tauPUT)))*sqrt((1-exp(-2*kap*(tauPUT)))/(2*kap))/kap
  
  h<-log(PtauZCB/(PtauPUT*strike))/sigZCB+0.5*sigZCB
  
  PUTvasicek<-strike*PtauPUT*pnorm(-h+sigZCB)-PtauZCB*pnorm(-h)
  
}

### Vasicek time-0 call price
CALLvasicek<-function(r,tauCALL,tauZCB,strike,Pparam,riskpremium=0){
  PtauCALL<-ZCBvasicek(r,tauCALL,Pparam,riskpremium)
  PtauZCB<-ZCBvasicek(r,tauZCB,Pparam,riskpremium)
  CALLvasicek<-PUTvasicek(r,tauCALL,tauZCB,strike,Pparam,riskpremium)+PtauZCB-strike*PtauCALL
}




### CALCULATIONS
### ============


### Solve for r*
r_star <- uniroot(r_star_solution, interval = c(-0.1, 0.2), K = K)$root

### Calculate adjusted strikes using r*
AdjustedStrikes <- rep(NA, N)
for (i in 1:N){
  T[i] <- i+1
  AdjustedStrikes[i] <- ZCBpricevasicek(timetomat, T[i], r_star)
}

### Alternative calculation of adjusted strikes using r*
Ki <- exp(Avasicek(timetomat, Ti) - Bvasicek(timetomat, Ti) * r_star)



### Calculate time-0 Vasicek call prices
call <- rep(NA, length(Ti))
for (i in 1:length(Ti)){
  call[i] <- CALLvasicek(r0, tauCALL = timetomat, tauZCB = Ti[i], 
                         AdjustedStrikes[i], c(thetaQ, kappa, sigma), riskpremium=0)
}

### Calculate price of call on coupon-bearing bond
CallPrice <- sum(alpha * call)








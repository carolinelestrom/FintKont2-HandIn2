setwd("/Users/carolinelestrom/Documents/KU/FinKont2/Handins/HandIn2")

library(png) ### Plots
library(knitr) ### Data manipulations
library(tidyverse) ### Data manipulations
library(dtplyr) ### Data manipulations - merge datasets
library(ggplot2) ### Plots
library(gridExtra) ### Plots






# Q3 - Linking Local and Implied Volatility with Dupire

## Q3.e


### FUNCTIONS
### =========

### Bachelier Pricing
BachelierFormula <- function(spot, strike, timetomat, t, sigma){
  d1 = (spot-strike)/(sigma*sqrt(timetomat-t))
  pi_bach <- (spot-strike)*pnorm(d1) + sigma*sqrt(timetomat-t) * dnorm(d1)
  return(pi_bach) 
}

### Black-Scholes Pricing
BlackScholesFormula <- function(spot, strike, t, r, div, timetomat, sigma){
  d1 = (1/(sigma*sqrt(timetomat-t))) * (log(spot/strike) + (r + 0.5*sigma^2)*(timetomat-t))
  d2 = d1 - sigma * sqrt(timetomat-t)
  pi_BS = spot * exp(-div*timetomat)*pnorm(d1) - exp(-r*(timetomat-t))*strike*pnorm(d2)
  return(pi_BS)
}

### Implied Volatility
IV <- function(obsprice, K){ ### Takes input: Observed Call Price and Strike
  difference <- function(sigmaIMP){ ### sigmaIMP = implied volatility
    abs(BlackScholesFormula(S0, K, 0, 0, 0, T, sigmaIMP) - obsprice)
  }
  optimize(difference, c(10^(-6), 10), maximum = FALSE)$minimum
}

#Now we define function for the benchmark and the approximations
### Function for Benchmark approximation
BenchmarkApprox <- function(strike){
  sigma/strike  
}

### Function for zero-sigma approximation
ZeroSigmaApprox <- function(strike){
  sqrt(((sigma/S0)^2 + (sigma/strike)^2)/2)
}

LampertiTransformation <- function(u){
  1/(u*(sigma/u))
}

LampertiApprox <- function(strike){
  result <- integrate(LampertiTransformation, strike, S0)
  if (result$message != "OK") {
    return(NA) ### return NA if the integral fails to converge
  }
  return (log(S0/strike)/result$value)
}



### INITIALIZE
### ==========

S0 = 100
T = 1
sigma = 15




### CALCULATIONS
### ============



strikes <- seq(50, 150, length.out = 100)
IVBach <- rep(0, length(strikes))
BM <- rep(0, length(strikes))
ZS <- rep(0, length(strikes))
L <- rep(0, length(strikes))


for (i in 1:length(strikes)){
  K <- i+min(strikes)-1
  price <- BachelierFormula(S0, K, T, 0, sigma)
  IVBach[i] <- IV(price, K)
  BM[i] <- BenchmarkApprox(K)
  ZS[i] <- ZeroSigmaApprox(K)
  L[i] <- LampertiApprox(K)
}




### PLOT
### ====


data <- as.data.frame(cbind(strikes, IVBach, BM, ZS, L))
colors <- c("True Implied Volatility" = "#901a1E", "Benchmark" = "steelblue", "Zero-Sigma Approximation" = "#39641c", "Lamperti Approximation" = "hotpink")



p1 <- data %>%
  ggplot(aes(x = strikes, y = IVBach)) +
  geom_line(aes(color = "True Implied Volatility"), linewidth = 1.7) +
  theme_bw() +
  xlab("Strike") + ylab("Implied volatilities") +
  ggtitle("Implied volatility across strikes in the Bachelier model") +
  theme(plot.title = element_text(size=20, hjust=0)) +
  theme(axis.title = element_text(size=17)) +
  theme(axis.text.x = element_text(size=13, angle = 0, vjust = 0.7),
        axis.text.y = element_text(size=13)) +
  scale_x_continuous(breaks=c(50, 75, 100, 125, 150), limits = c(50, 150), 
                     labels = c(50, 75, "100 (ATM)", 125, 150)) +
  scale_y_continuous(breaks=c(0.12, 0.14, 0.16, 0.18, 0.20), limits = c(0.117, 0.21)) +
  geom_line(aes(y = BM, color = "Benchmark"), linewidth = 1.7, linetype = "dotted", alpha = 1) +
  geom_line(aes(y = ZS, color = "Zero-Sigma Approximation"), linewidth = 1.7, linetype = "dotted", alpha = 1) +
  geom_line(aes(y = L, color = "Lamperti Approximation"), linewidth = 1.7, linetype = "dotted", alpha = 1) +
  labs(color = "Approximations") +
  scale_color_manual(values = colors) +
  theme(legend.position = c(0.71, 0.83), 
        legend.text = element_text(colour="black", size=17),
        legend.title = element_text(colour="black", size=17,face="bold"),
        legend.background = element_rect(fill="white",size=0.5, linetype="solid", colour ="black"))



p2 <- data %>%
  ggplot(aes(x = strikes)) +
  theme_bw() +
  xlab("Strike") + ylab("Error") +
  ggtitle("Approximation error") +
  theme(plot.title = element_text(size=20, hjust=0)) +
  theme(axis.title = element_text(size=17)) +
  theme(axis.text.x = element_text(size=13, angle = 0, vjust = 0.7),
        axis.text.y = element_text(size=13)) +
  scale_x_continuous(breaks=c(50, 75, 100, 125, 150), limits = c(50, 150), 
                     labels = c(50, 75, "100 (ATM)", 125, 150)) +
  scale_y_continuous(breaks=c(-0.10, -0.05, 0.00, 0.05, 0.10), limits = c(-0.10, 0.10)) +
  geom_line(aes(y = (BM - IVBach)/IVBach, color = "Benchmark"), linewidth = 1.7, linetype = "dotted", alpha = 1) +
  geom_line(aes(y = (ZS - IVBach)/IVBach, color = "Zero-Sigma Approximation"), linewidth = 1.7, linetype = "dotted", alpha = 1) +
  geom_line(aes(y = (L - IVBach)/IVBach, color = "Lamperti Approximation"), linewidth = 1.7, linetype = "dotted", alpha = 1) +
  labs(color = "Approximations") +
  scale_color_manual(values = colors) +
  theme(legend.position = c(0.73, 0.83), 
        legend.text = element_text(colour="black", size=17),
        legend.title = element_text(colour="black", size=17,face="bold"),
        legend.background = element_rect(fill="white",size=0.5, linetype="solid", colour ="black"))





#png(file="2b.png")


grid.arrange(p1, p2, ncol = 2)


#dev.off()




### BASIC R PLOT
### ============


### Vectors for colors and legend
colourvec <- c("black", "deepskyblue", "hotpink", "chartreuse")
legendvec <- c("True IV", "Benchmark Approximation", "Zero-Sigma Approximation", "Lamperti Approximation")




### Plots
par(mfrow = c(1,2))
plot(strikes, IVBach, type = "l", lwd = 3, col = "black", xlab = "Strike", ylab = "Implied volatility ", ylim = c(0.1, 0.3))

lines(strikes, BM, type = "l", lwd = 3, lty = 2, col = "deepskyblue")
lines(strikes, ZS, type = "l", lwd = 3, lty = 2, col = "hotpink")
lines(strikes, L, type = "l", lwd = 3, lty = 2, col = "chartreuse")

grid()
legend(x= 100, y= 0.3, legend = legendvec, col = colourvec, lty = 1, lwd = 2, cex = 0.8)
title(main = "Bachelier market prices with B/S IV")

plot(strikes, (BM - IVBach)/IVBach, type = "l", lwd = 3, col = "deepskyblue", xlab = "Strike", ylab = "Error", ylim = c(-0.1, 0.1))
lines(strikes, (ZS - IVBach)/IVBach, type = "l", lwd = 3, lty = 2, col = "hotpink")
lines(strikes, (L - IVBach)/IVBach, type = "l", lwd = 3, lty = 2, col = "chartreuse")

grid()
legend(x= 100, y= 0.1, legend = legendvec[2:4], col = colourvec[2:4], lty = 1, lwd = 2, cex = 0.8)
title(main = "Approximation error")
















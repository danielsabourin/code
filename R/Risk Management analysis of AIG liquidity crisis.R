#############################
# ACTSC 445 Assignment 1    #
# Daniel Matheson: 20270871 #
#      May 27th 2017        #
#############################

### Question 2 Code

# Import data
setwd("C:/Users/Daniel/Dropbox/Spring 2017/ACTSC 445/Assignment 1")

# Take losses values only
losses_data <- read.csv("MC_Output_A1.csv",header=T)
losses <- losses_data$Simulated.Losses.using.Monte.Carlo

# Sorted losses for estimating VaR/CTE
sorted_losses <- sort(losses)

N <- length(losses) # number of observations

### Question 2(a):
alpha2a <- 0.97
q2a <- 0.97 * N
VaR2a <- sorted_losses[q2a]

### Question 2(b):
alpha2b <- 0.97
q <- 0.95      # 95% confidence interval
qc <- (1+q)/2  # centered quantile; as in equal tails
k <- qnorm(qc, 0,1) * sqrt(N * alpha2b * (1 - alpha2b))
k <- ceiling(k)
CI2b <- c(sorted_losses[N * alpha2b - k], sorted_losses[N * alpha2b + k])

### Question 2(c)
alpha2c <- 0.95 
cte_values <- sorted_losses[(N * alpha2c + 1):N]
cte2c <- (1/(N * (1-alpha2c))) * sum(cte_values)

### Question 2(d)
alpha2d <- 0.95
q2d <- alpha2d * N
VaR2d <- sorted_losses[q2d]
S_n <- (1/(N*(1 - alpha2d) -1)) * sum((cte_values - cte2c)^2)
varcte <- (1/(N*(1 - alpha2d))) *(S_n + alpha2d *(cte2c - VaR2d)^2)
ctese <- sqrt(varcte)

### Question 2(f)
n <- 1000
VaR2f <- rep(0,10)
VaRCI2f <- rep(0,20)
dim(VaRCI2f) <- c(10,2)
CTE2f <- rep(0,10)
alphaVaR <- 0.97
qVaR <- n * alphaVaR
alphaCTE <- 0.95
qCTE <- n * alphaCTE
losses_temp <- rep(0,n)

for(ii in 1:10){ 
  losses_temp <- losses[((ii-1)*n+1):(ii*n)]
  losses_temp <- sort(losses_temp)
  VaR2f[ii] <- losses_temp[qVaR]
  cte_values2f <- losses_temp[(qCTE + 1):n]
  CTE2f[ii] <- (1/(n * (1- alphaCTE))) * sum(cte_values2f)
}

VaR2fmean <- mean(VaR2f)
CTE2fmean <- mean(CTE2f)
VaR2fsd <- sd(VaR2f)
CTE2fsd <- sd(CTE2f)

### Question 4(a)
S_0 <- 100
r <- 0.04
sigma <- 0.3
T <- 10/256
K <- 87
d_1 <- (log(S_0/K) + (r + 1/2 * sigma^2) * T)/(sigma * sqrt(T))
d_2 <- d_1 - sigma * sqrt(T)
phi1 <- pnorm(d_1,0,1)
phi2 <- pnorm(d_2,0,1)
C <- S_0 * phi1 - K * exp(-r * T) * phi2

phi1put <- pnorm(-d_1,0,1)
phi2put <- pnorm(-d_2,0,1)
P <- K * exp(-r * T) * phi2put - S_0 * phi1put

### Question 4(b)
mu <- 0.05
S_T = c(90,89,88,87,86,85, 87.293)
obs <- ((log(S_T/S_0)) - mu * T)/(sigma * sqrt(T))
p_obs <- pnorm(obs,0,1)

C1_val <- pmax(S_T - 87,0)
C2_val <- - pmax(S_T - 91,0)
P1_val <- - pmax(87 - S_T,0)
M <- 239798539

portfolio_value <- 1000000 * (C1_val + C2_val) + M * P1_val
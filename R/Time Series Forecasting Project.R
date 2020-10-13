setwd("C:/Users/Daniel/Dropbox/Spring 2017/STAT 443/Final Project")
#################################   
#    STAT 443: Final Project    #
#   Daniel Matheson: 20270871   #
#         Spring 2017           #
#################################
library(tseries)
library(ggplot2)
###############################
####      Question 1       ####
###############################
input.Data <- read.csv("GDP_CONS_CANADA.csv")
w1t <- as.numeric(as.character(input.Data$GDP))
w2t <- as.numeric(as.character(input.Data$CONS))

plot(w1t, xlab = "Time", ylab = "Raw Values", 
     ylim = c(min(min(w1t), min(w2t)) - 500, max(max(w1t),max(w2t)) + 500),
     type = "l", col = "blue", lwd = 1)
lines(w2t, col = "red", lwd = 1)
legend(20, 1200000, legend = c("GDP", "CONS"), lty = c(1,1),lwd=c(2.5,2.5),
       col=c("blue","red"), cex = 0.7)
# CONS = Canadian Personal Expenditures on Consumer Goods and Service

###############################
####      Question 2       ####
###############################
xt <- log(w1t)
## Trend Stationary Method
tSeq <- 1:length(xt)
TSmodel <- lm(xt ~ tSeq + 1)
TSyt <- TSmodel$residuals

## Difference Stationary Method
xt.diff <- diff(xt)
DSmodel <- lm(xt.diff ~ 1)
DSyt <- DSmodel$residuals
plot(TSyt, type = "l", xlab = "Time", 
     ylab = "Cycle Yt", main = "Difference Stationary Model")
plot(DSyt, type = "l", xlab = "Time", 
     ylab = "Cycle Yt", main = "Trend Stationary Model")

###############################
####      Question 3       ####
###############################
BIC <- function(res, k, N) {
  bic = log(sum(res^2) / N)
  bic = bic + log(N) * k / N
  bic}

# DS Method
Q3DSBIC <- rep(0,10)
for (k in 1:10){
  Q3DSmodel <- arima(DSyt, order = c(k,0,0), include.mean = F)
  Q3DS.resid <- Q3DSmodel$residuals
  Q3DSBIC[k] <- BIC(Q3DS.resid, k, length(DSyt))}

# TS Method
Q3TSBIC <- rep(0,10)
for (k in 1:10){
  Q3TSmodel <- arima(TSyt, order = c(k,0,0), include.mean = F)
  Q3TS.resid <- Q3TSmodel$residuals
  Q3TSBIC[k] <- BIC(Q3TS.resid, k, length(TSyt))}

indices <- 1:10
BICmins <- c(indices[Q3DSBIC == min(Q3DSBIC)], indices[Q3TSBIC == min(Q3TSBIC)])
# DS Method: p = 1 (set to p = 2)
# TS Method: p = 2

# Final Q3 model, TS with p = 2:
TSyt.N<- length(TSyt)
Q3model <-  lm(TSyt[-(1:2)] ~TSyt[-c(1,TSyt.N)] + TSyt[-((TSyt.N-1):TSyt.N)] -1)
Q3phi <- Q3model$coef

# rho(k) and psi_k for k = 0,1,2,3,4,5,6:
Q3rho <- ARMAacf(ar = Q3phi, pacf = F, lag.max = 6)
Q3psi <- rep(NA,7)
Q3psi[1] <- 1
Q3psi[2] <- Q3phi[1] * Q3psi[1] 
Q3psi[3] <- Q3phi[1] * Q3psi[2] + Q3phi[2] * Q3psi[1]
for( j in 3:7) {
  Q3psi[j] <- Q3psi[c((j-1),(j-2))] %*% Q3phi[1:2]}
# gamma(0)^{1\2}
Q3.resid <- Q3model$residuals 
Q3sigma.sq <- sum((Q3.resid)^2)/length(TSyt)
Q3gamma.zero <- Q3sigma.sq/(1 - Q3phi[1] * Q3rho[2] - Q3phi[2] * Q3rho[3])
Q3gamma.zero.sqrt <- sqrt(Q3gamma.zero)

###############################
####      Question 4       ####
###############################
# DS Method: under estimation that Y_t ~ AR(2)
DSyt.N = length(DSyt)
Q4DS.model <- lm(DSyt[-(1:2)] ~DSyt[-c(1,DSyt.N)] + DSyt[-((DSyt.N-1):DSyt.N)] -1)
Q4DS.mu <- DSmodel$coefficients["(Intercept)"]
Q4DS.phi <- Q4DS.model$coef
Q4DS.N <- length(xt.diff)
Q4DS.E <- rep(0,9) # E_t[X_t+k], k from 0 to 8
Q4DS.Y <- rep(0,9)  # E_t[Y_t+k], k from 0 to 8
Q4DS.Y[1] <- DSyt[Q4DS.N]
Q4DS.Y[2] <- Q4DS.phi %*% DSyt[c(Q4DS.N,(Q4DS.N-1))] # Y_T+1

for (j in 3:9){
  Q4DS.Y[j] <- Q4DS.phi %*% Q4DS.Y[c(j-1,j-2)]}

Q4DS.E[1] <- xt.diff[Q4DS.N]
Q4DS.E[2:9] <- Q4DS.mu + Q4DS.Y[2:9]

Q4DS.psi <- rep(0,9)
Q4DS.psi[1] <- 1
Q4DS.psi[2] <- Q4DS.phi[1] * Q4DS.psi[1] 
for( j in 3:9) {
  Q4DS.psi[j] <- Q4DS.phi %*% Q4DS.psi[c((j-1),(j-2))]}

Q4DS.sigma.sq <- sum((Q4DS.model$residuals)^2)/(Q4DS.N - 2)
Q4DS.var <- rep(0,9) # Var(delta X_t+k) for k = 0 ... 8
Q4DS.var[1] <- 0
for (k in 2:9){
Q4DS.var[k] <- Q4DS.sigma.sq * sum((Q4DS.psi[c(1:(k-1))])^2)}

Q4DS.CI.U <- rep(0,9) # 95% confidence interval bounds
Q4DS.CI.L <- rep(0,9)
Q4DS.CI.U <- Q4DS.E + 1.96 * sqrt(Q4DS.var)
Q4DS.CI.L <- Q4DS.E - 1.96 * sqrt(Q4DS.var)
Q4DS.CI <- cbind(round(Q4DS.CI.L,5), round(Q4DS.CI.U,5))
colnames(Q4DS.CI) <- c("Lower Bound", "Upper Bound")

Q4DS.E.vec <- rep(0,20)
Q4DS.E.vec[13:20] <- Q4DS.E[-1]
Q4DS.E.vec[1:12] <- xt.diff[(length(xt.diff)-11):length(xt.diff)]
Q4DS.CI.L.vec <- rep(0,20)
Q4DS.CI.L.vec[1:12] <- rep(NA,12)
Q4DS.CI.L.vec[13:20] <- Q4DS.CI.L[-1]
Q4DS.CI.U.vec <- rep(0,20)
Q4DS.CI.U.vec[1:12] <- rep(NA,12)
Q4DS.CI.U.vec[13:20] <- Q4DS.CI.U[-1]

plot(Q4DS.E.vec, type = "l", ylim = c(min(Q4DS.CI.L)*1.1, max(Q4DS.CI.U)*1.1),
     xlab = "Time (Final Obs. = T = 12)", ylab = "Delta Xt",
     main = "Difference Stationary Delta Xt Forecast")
lines(rep(Q4DS.mu,20), col = "red", lty =3)
lines(Q4DS.CI.L.vec, col = "blue")
lines(Q4DS.CI.U.vec, col = "blue")
abline(v = 13, col = "green", lty = 6)
legend(x = "topleft", legend = c("Delta Xt", "95% CI Bounds", "Delta Xt Mu", "k = 1"),
       lty = c(1,1,3,6),lwd= c(0.5,0.5,0.5,0.5),col= c("black","blue","red", "green"), cex = 0.8,
       seg.len = c(1,1,1,1), text.width = 3.2)
  
# TS Method
Q4TS.model <- lm(TSyt[-(1:2)] ~TSyt[-c(1,TSyt.N)] + TSyt[-((TSyt.N-1):TSyt.N)] -1)
Q4TS.mu <- TSmodel$coefficients["tSeq"]
Q4TS.phi <- Q4TS.model$coef
Q4TS.N <- length(xt.diff)
Q4TS.E <- rep(0,9)   # E_t[X_t+k], k from 0 to 8
Q4TS.Y <- rep(0,10)  # E_t[Y_t+k], k from -1 to 8

Q4TS.Y[1] <- TSyt[(Q4TS.N)]   #Yt[184]
#Yt[185] : last Yt since length(delta X_t) < length(Yt) for TS:
Q4TS.Y[2] <- TSyt[(Q4TS.N+1)] 

for (j in 3:10){
  Q4TS.Y[j] <- Q4TS.phi %*% Q4TS.Y[c((j-1),(j-2))]}

Q4TS.E[1] <- xt.diff[Q4TS.N] # k = 0
Q4TS.E[2:9] <- Q4TS.mu + Q4TS.Y[-c(1,2)] - Q4TS.Y[2:9]
Q4TS.psi <- rep(0,8)
Q4TS.psi[1] <- 1
Q4TS.psi[2] <- Q4TS.phi[1] * Q4TS.psi[1] 
Q4TS.psi[3] <- Q4TS.phi[1] * Q4TS.psi[2] + Q4TS.phi[2] * Q4TS.psi[1]
for( j in 3:8) {
  Q4TS.psi[j] <- Q4TS.psi[c((j-1),(j-2))] %*% Q4TS.phi[1:2]}

Q4TS.sigma.sq <- sum((Q4TS.model$residuals^2))/(Q4TS.N-2)
Q4TS.var <- rep(0,9) # Var(delta X_t+k) for k = 1 ... 8
Q4TS.var[1] <- 0
Q4TS.var[2] <- Q4TS.sigma.sq
for (k in 3:9){
  Q4TS.var[k] <- Q4TS.sigma.sq * sum( (Q4TS.psi[2:(k-1)] - Q4TS.psi[1:(k-2)])^2)}

Q4TS.CI.U <- rep(0,9) # 95% confidence interval bounTS
Q4TS.CI.L <- rep(0,9)
Q4TS.CI.U <- Q4TS.E + 1.96 * sqrt(Q4TS.var)
Q4TS.CI.L <- Q4TS.E - 1.96 * sqrt(Q4TS.var)
Q4TS.CI <- cbind(round(Q4TS.CI.L,5), round(Q4TS.CI.U,5))
colnames(Q4TS.CI) <- c("Lower Bound", "Upper Bound")
Q4TS.E.vec <- rep(0,20)
Q4TS.E.vec[13:20] <- Q4TS.E[-1]
Q4TS.E.vec[1:12] <- xt.diff[(length(xt.diff)-11):length(xt.diff)]
Q4TS.CI.L.vec <- rep(0,20)
Q4TS.CI.L.vec[1:12] <- rep(NA,12)
Q4TS.CI.L.vec[13:20] <- Q4TS.CI.L[-1]
Q4TS.CI.U.vec <- rep(0,20)
Q4TS.CI.U.vec[1:12] <- rep(NA,12)
Q4TS.CI.U.vec[13:20] <- Q4TS.CI.U[-1]
plot(Q4TS.E.vec, type = "l", ylim = c(min(Q4TS.CI.L)*1.1, max(Q4TS.CI.U)*1.1),
     xlab = "Time (Final Obs. = T = 12)", ylab = "Delta Xt",
     main = "Trend Stationary Delta Xt Forecast")
lines(rep(Q4TS.mu,20), col = "red", lty =3)
lines(Q4TS.CI.L.vec, col = "blue")
lines(Q4TS.CI.U.vec, col = "blue")
abline(v = 13, col = "green", lty = 6)
legend(x = "topleft", legend = c("Delta Xt", "95% CI Bounds", "Delta Xt Mu", "k = 1"),
       lty = c(1,1,3,6),lwd= c(0.5,0.5,0.5,0.5),col= c("black","blue","red", "green"),
       cex = 0.8, seg.len = c(1,1,1,1), text.width = 3.2)

###############################
####      Question 5       ####
###############################
# Here the H0 is that the series is difference stationary.
Q5ADF <- adf.test(xt,k=5)
# p = 0.2865
# We do not reject H_0: x_t is difference stationary.

###############################
####      Question 6       ####
###############################
# DS Method
DSkmax <- 20
Q6DSrho <- acf(DSyt, lag.max = DSkmax, plot = F)$acf
Q6DSpsi <- rep(0,DSkmax+1)
Q6DSpsi[1] <- 1
for (k in 1:DSkmax){
  Q6DSmodels <- arima(DSyt, order = c(k,0,0), include.mean = F)
  Q6DSpsi[(k+1)] <- Q6DSmodels$coef[length(Q6DSmodels$coef)]}

# TS Method
TSkmax <- 20
Q6TSrho <- acf(TSyt, lag.max = TSkmax, plot = F)$acf
Q6TSpsi <- rep(0,TSkmax+1)
Q6TSpsi[1] <- 1
for (k in 1:TSkmax){
  Q6TSmodels <- arima(TSyt, order = c(k,0,0), include.mean = F)
  Q6TSpsi[(k+1)] <-  Q6TSmodels$coef[length(Q6TSmodels$coef)]}

################################################
# NOTE ON IDENTIFICATION OF CUT-OFF:
# k is ahead by 1 since psi_0/rho(0) are included
# then we take the index previous to the one which
# passes H_0, so this is index (k-2)
# with max(k-2,0) in case the cut-off happens immediately 
################################################
#### Identifying p,q in ARMA(p,q) for DS:
Q6DSp <- 0 # will return 0 if H_0: psi_k = 0 is rejected for every k
for (k in 1:(DSkmax+1)){
  if(abs(Q6DSpsi[k]) < 2/sqrt(length(DSyt))){
    Q6DSp <- max(k-2,0) 
    break}}
Q6DSq <- 0 # will return 0 if H_0: rho_k = 0 is rejected for every k
for (k in 1:(DSkmax+1)){
  if(abs(Q6DSrho[k]) < 2/sqrt(length(DSyt))){
    Q6DSq <- max(k-2,0)  
    break}}

#### Identifying p,q in ARMA(p,q) for TS:
Q6TSp <- 0 # will return 0 if H_0: psi_k = 0 is rejected for every k
for (k in 1:TSkmax){
  if(abs(Q6TSpsi[k]) < 2/sqrt(length(TSyt))){
    Q6TSp <- max(k-2,0)  
    break}}
Q6TSq <- 0 # will return 0 if H_0: rho_k = 0 is rejected for every k
for (k in 1:length(Q6TSrho)){
  if(abs(Q6TSrho[k]) < 2/sqrt(length(TSyt))){
    Q6TSq <- max(k-2,0) 
    break}}

Q6matrix <- matrix(c(Q6DSp, Q6TSp, Q6DSq, Q6TSq), nrow = 2, ncol = 2)
rownames(Q6matrix) <- c("DS", "TS")
colnames(Q6matrix) <- c("p","q")
##
#### Recall: rho(k) determines q, psi_k determines p
#### Plots to show cut-off properties:
##
DSplot.x <- rep(0,21)
DSplot.x[2:21] <- 1:20
par(mfrow = c(1,2), mar = c(2,2,1,1))
plot(x = DSplot.x, y =Q6DSpsi, type = "l", main = "Difference Stationary Phi_kk")
abline(h = 0, col = "red", lty = 2)
abline(v = Q6DSp, col = "blue", lty = 2)
abline(h = -1.96*sqrt(1/length(DSyt)), lty = 6, col = "orange")
abline(h = 1.96*sqrt(1/length(DSyt)), lty = 6, col = "orange")
legend(x = "topright", legend = c("Phi_kk", "Cut-off k = 1", "+/- 1.96SE(Phi_kk)"),
       lty = c(1,2,2), col = c("black", "blue", "orange"), cex = 0.7)
plot(x = DSplot.x, y =Q6DSrho, type = "l", main = "Difference Stationary Rho(k)",
     ylim = c(-1.96*sqrt(1/length(DSyt)),max(Q6DSrho)*1.1))
abline(h = 0, col = "red", lty = 2)
abline(v = Q6DSq, col = "blue", lty =2)
abline(h = -1.96*sqrt(1/length(DSyt)), lty = 6,col = "orange")
abline(h = 1.96*sqrt(1/length(DSyt)), lty = 6,col = "orange")
legend(x = "topright", legend = c("Rho(k)", "Cut-off k = 3", "+/- 1.96SE(Rho(k))"),
       lty = c(1,2,2), col = c("black", "blue", "orange"), cex = 0.7)

TSplot.x <- rep(0,21)
TSplot.x[2:21] <- 1:20
par(mfrow = c(1,2), mar = c(2,2,1,1))
plot(x = TSplot.x , y = Q6TSpsi, type = "l", main = "Trend Stationary Phi_kk")
abline(h = 0, col = "red", lty = 2)
abline(v = Q6TSp, col = "blue", lty = 2)
abline(h = -1.96*sqrt(1/length(TSyt)), lty = 6,col = "orange")
abline(h = 1.96*sqrt(1/length(TSyt)), lty = 6,col = "orange")
legend(x = "topright", legend = c("Phi_kk", "Cut-off k = 2", "+/- 1.96SE(Phi_kk)"),
       lty = c(1,2,2), col = c("black", "blue", "orange"), cex = 0.7)
plot(x = TSplot.x, y= Q6TSrho, type = "l", main = "Trend Stationary Rho(k)", xlab = "k")

#### ARMA models
Q6DS.model1 <- arima(DSyt, order = c(1,0,0), include.mean = F,
                     optim.control = list(maxit = 10000))
Q6DS.model2 <- arima(DSyt, order = c(0,0,3), include.mean = F,
                     optim.control = list(maxit = 10000))
Q6TS.model <- arima(TSyt, order = c(2,0,0), include.mean = F,
                    optim.control = list(maxit = 10000))
# AR(1) DS
Q6DS.resid1 <- Q6DS.model1$residuals
Q6DS.sigmasq1 <- sum((Q6DS.resid1)^2)/length(DSyt)
Q6DS.std.resid1 <- Q6DS.resid1/sqrt(Q6DS.sigmasq1)
# MA(3) DS
Q6DS.resid2 <- Q6DS.model2$residuals
Q6DS.sigmasq2 <- sum((Q6DS.resid2)^2)/length(DSyt)
Q6DS.std.resid2 <- Q6DS.resid2/sqrt(Q6DS.sigmasq2)
Q6TS.resid <- Q6TS.model$residuals
Q6TS.sigmasq <- sum((Q6TS.resid)^2)/length(TSyt)
Q6TS.std.resid <- Q6TS.resid/sqrt(Q6TS.sigmasq)

####### Plots of std. resid
# Means of residuals are extremely close to 0, so ignore when plotting
par(mfrow = c(1,2), mar = c(4,4,2,2))
hist(Q6DS.std.resid1, breaks = 30, freq = F, 
     xlab = "Standardized Residuals", main = "Difference Stationary AR(1)",
     cex.main = 0.8)
curve(dnorm(x),col = "red", add = T)
legend(x = "topright", legend = "N(0,1) PDF",
       lty = 1,lwd=0.5,col="red", cex = 0.6,
       seg.len = 0.5, text.width = 1.35, x.intersp = 0.3)
hist(Q6DS.std.resid2, breaks = 30, freq = F,
     xlab = "Standardized Residuals", main = "Difference Stationary MA(3)",
     cex.main = 0.8)
curve(dnorm(x),col = "red", add = T)
legend(x = "topright", legend = "N(0,1) PDF",
       lty = 1,lwd=0.5,col="red", cex = 0.6,
       seg.len = 0.5, text.width =1.35, x.intersp = 0.3)

par(mfrow= c(1,1), mar = c(4,4,2,2))
hist(Q6TS.std.resid, breaks = 30, freq = F,
     xlab = "Standardized Residuals", main = "Trend Stationary AR(2)",
     cex.main = 0.8)
curve(dnorm(x),col = "red", add = T)
legend(x = "topright", legend = "N(0,1) PDF",
       lty = 1,lwd=0.5,col="red", cex = 0.6,
       seg.len = 0.5, text.width =1.35, x.intersp = 0.3)

####### Box Pierce Tests
# DS AR(1)
# sqrt(length(Q6DS.resid1)) = 13.56 -> round to 14
DS1.BoxPierce <- Box.test(x = Q6DS.resid1, type = "Box-Pierce",lag=14) 
# p = 0.005105; reject H_0  -> residuals correlated
# DS MA(3)
DS2.BoxPierce <- Box.test(x = Q6DS.resid2, type = "Box-Pierce",lag=14)
# p = 0.09984; do not reject H_0 -> residuals uncorrelated (but close)
# TS
TS.BoxPierce <- Box.test(x = Q6TS.resid, type = "Box-Pierce",lag=14) 
# p = 0.003365; reject H_0 -> residuals correlated


####### Overfitting with r = 4
Q6DS.model1of <- arima(DSyt, order = c(5,0,0), include.mean = F,
                       optim.control = list(maxit = 10000))
Q6DS.model2of <- arima(DSyt, order = c(0,0,7), include.mean = F,
                       optim.control = list(maxit = 10000))
Q6TS.modelof <- arima(TSyt, order = c(6,0,0), include.mean = F,
                      optim.control = list(maxit = 10000))
Q6DS.sigmasq1.of <- Q6DS.model1of$sigma2
Q6DS.sigmasq2.of <- Q6DS.model2of$sigma2
Q6TS.sigmasq.of <- Q6TS.modelof$sigma2
Q6DS.Lambda1 <- (length(DSyt) - 5) * log(Q6DS.sigmasq1/Q6DS.sigmasq1.of)
Q6DS.Lambda2 <- (length(DSyt) - 7) * log(Q6DS.sigmasq2/Q6DS.sigmasq2.of)
Q6TS.Lambda <- (length(TSyt) - 5) * log(Q6TS.sigmasq/Q6TS.sigmasq.of)
of.crit <- qchisq(0.95, 4)

# if true then reject H_0:
Q6DS.of.test1 <- (Q6DS.Lambda1 > of.crit)
Q6DS.of.test2 <- (Q6DS.Lambda2 > of.crit)
Q6TS.of.test <- (Q6TS.Lambda > of.crit)
# all false -> do not reject any of H_0. overfitting fails.

####### Jarque Bera Test
# can fail this test in 2 ways:
# 1. normality is rejected completely
# 2. the model is normal but didn't pass JB test because of
# the existence of significant outliers or structure breaks

############################################################
# if an ARMA model passes all diagnostics except for JB due
# to outliers, it's still a good model
############################################################

# DS AR(1)
DS.resid1 <- Q6DS.resid1
DS.sigma.sq1 <- Q6DS.sigmasq1
DS.Z1 <- (DS.resid1 - mean(DS.resid1))/sqrt(DS.sigma.sq1)
DS.n1 <- length(DS.Z1)
DS.k31 <- 1/DS.n1 * sum((DS.Z1)^3)
DS.k41 <- 1/DS.n1 * sum((DS.Z1)^4)
DS.J1 <- DS.n1 * ( (DS.k31^2)/6 + (DS.k41-3)^2/24)
chisq.crit <- qchisq(0.95,2)
DS.J.test1 <- (DS.J1 > chisq.crit) # if true, reject H_0

# DS MA(3)
DS.resid2 <- Q6DS.resid2
DS.sigma.sq2 <- Q6DS.sigmasq2
DS.Z2 <- (DS.resid2 - mean(DS.resid2))/sqrt(DS.sigma.sq2)
DS.n2 <- length(DS.Z2)
DS.k32 <- 1/DS.n2 * sum((DS.Z2)^3)
DS.k42 <- 1/DS.n2 * sum((DS.Z2)^4)
DS.J2 <- DS.n2 * ( (DS.k32^2)/6 + (DS.k42-3)^2/24)
DS.J.test2 <- (DS.J2 > chisq.crit) # if true, reject H_0

# TS 
TS.resid <- Q6TS.resid
TS.sigma.sq <- Q6TS.sigmasq
TS.Z <- (TS.resid - mean(TS.resid))/sqrt(TS.sigma.sq)
TS.n <- length(TS.Z)
TS.k3 <- 1/TS.n * sum((TS.Z)^3)
TS.k4 <- 1/TS.n * sum((TS.Z)^4)
TS.J <- TS.n * ( (TS.k3^2/6) + (TS.k4-3)^2/24)
TS.J.test <- (TS.J > chisq.crit) # if true, reject H_0

####### ARCH(6) test for independence of resid
# DS AR(1)
ARCH.q <- 6
N1 <- length(Q6DS.resid1)
Q6DS.resid1 <- Q6DS.resid1^2
ARCH.DS.model1 <- lm(Q6DS.resid1[-(1:6)]~Q6DS.resid1[-c((1:5), N1)] 
                     +  Q6DS.resid1[-c(1:4,(N1-1),N1)] + Q6DS.resid1[-c(1:3, (N1-2):N1)] 
                     + Q6DS.resid1[-c(1:2, (N1-3):N1)] + Q6DS.resid1[-c(1, (N1-4):N1)] 
                     + Q6DS.resid1[-((N1-5):N1)] - 1)
DS.Rsq1 <- summary(ARCH.DS.model1)$r.squared
DS.ARCH.test1 <- N1 * DS.Rsq1
DS.ARCH.crit <- qchisq(0.95,6)
DS.ARCH.H01 <- (DS.ARCH.test1 > DS.ARCH.crit) # if true reject H_0

# DS MA(3)
N2 <- length(Q6DS.resid2)
Q6DS.resid2 <- Q6DS.resid2^2
ARCH.DS.model2 <- lm(Q6DS.resid2[-(1:6)]~Q6DS.resid2[-c((1:5), N2)] 
                     + Q6DS.resid2[-c(1:4,(N2-1),N2)] + Q6DS.resid2[-c(1:3, (N2-2):N2)] 
                     + Q6DS.resid2[-c(1:2, (N2-3):N2)] + Q6DS.resid2[-c(1, (N2-4):N2)] 
                     + Q6DS.resid2[-((N2-5):N2)] - 1)
DS.Rsq2 <- summary(ARCH.DS.model2)$r.squared
DS.ARCH.test2 <- N2 * DS.Rsq2
DS.ARCH.crit <- qchisq(0.95,6)
DS.ARCH.H02 <- (DS.ARCH.test2 > DS.ARCH.crit) # if true reject H_0

# TS
N3 <- length(Q6TS.resid)
Q6TS.resid <- Q6TS.resid^2
ARCH.TS.model <- lm(Q6TS.resid[-(1:6)]~Q6TS.resid[-c((1:5), N3)] 
                    + Q6TS.resid[-c(1:4,(N3-1),N3)] + Q6TS.resid[-c(1:3, (N3-2):N3)] 
                    + Q6TS.resid[-c(1:2, (N3-3):N3)] + Q6TS.resid[-c(1, (N3-4):N3)] 
                    + Q6TS.resid[-((N3-5):N3)] - 1)
TS.Rsq <- summary(ARCH.TS.model)$r.squared
TS.ARCH.test <- N3 * TS.Rsq
TS.ARCH.crit <- qchisq(0.95,6)
TS.ARCH.H0 <- (TS.ARCH.test > TS.ARCH.crit) # if true reject H_0

###############################
####      Question 8       ####
###############################
Q8data <- read.csv("SP_data_for_Q8.csv")
P <- Q8data$P
lnP <- log(P)
N <- length(lnP)
Q8.model <- lm(lnP[-1] ~ 1 + lnP[-N])
at <- Q8.model$residuals
at.N <- length(at)

at.rho <- acf(at, lag.max = at.N, plot = F)$acf
at.psi <- rep(0,10)
for (k in 1:10){
  at.models <- arima(at, order = c(k,0,0), include.mean = F)
  at.psi[k] <-  at.models$coef[k]}
# sqrt(at.N) ~ 25
at.BoxPierce <- Box.test(x = at, type = "Box-Pierce",lag=25) 
# p-value = 0.6529; do not reject H_0; uncorrelated

## Checking if a_t is normally distributed:
at.sigmasq <- sum((Q8.model$residuals)^2)/(at.N-1)
at.std <- (at - mean(at))/sqrt(at.sigmasq)

# Jarque-Bera Test for Normality
at.k3 <- 1/at.N * sum((at.std)^3)
at.k4 <- 1/at.N * sum((at.std)^4)
at.J <- at.N * ( (at.k3^2/6) + (at.k4-3)^2/24)
chisq.crit <- qchisq(0.95,2)


at.J.test <- (at.J > chisq.crit) # if true, reject H_0
# true -> reject H_0 -> residuals not normal, or there are outliers:
# at.J = 547

# outliers are indices 16, 585:
ato.std <- at.std[-c(16,585)]
ato.N <- at.N - 2

# plots
par(mfrow= c(1,2), mar=c(4,4,2,2))
hist(at.std, breaks = 50, freq = F, 
     xlab = "Standardized Residuals", main = expression("Standardized Residuals a"[t]))
curve(dnorm(x),col = "red", add = T)
legend(x = "topleft", legend = "N(0,1) PDF",
       lty = 1,lwd=0.5,col="red", cex = 0.6,
       seg.len = 0.5, text.width = 3)
hist(ato.std, breaks = 50, freq = F, 
     xlab = "Standardized Residuals", main = expression("2 Outliers Removed"))
curve(dnorm(x),col = "red", add = T)
# looks normal

# JB test without outliers:
ato.k3 <- 1/ato.N * sum((ato.std)^3)
ato.k4 <- 1/ato.N * sum((ato.std)^4)
ato.J <- ato.N * ( (ato.k3^2/6) + (ato.k4-3)^2/24)

# if true, reject H_0
ato.J.test <- (ato.J > chisq.crit) # value 1.398
# false: do not reject H_0, residuals are normally distributed
# when only 2 of 647 observations removed
# -> will assume normal then

# autocorrelation function for a_t^2:
at.sq <- at^2
par(mfrow = c(1,1), mar = c(4,4,4,2))
at.sq.rho <- acf(at.sq, lag.max = 20, plot = T,
                 main = expression("Autocorrelation Function of squared a"[t]))$acf
# none of them are significant other than rho(0)
# remove rho(0) to fit into a smaller space and for visibility

### Fitting to GARCH(1,1):
require(fGarch)
# Here Yt is the appropriate time series you are trying to fit.
at.GARCH <- garchFit(formula = ~garch(1, 1), data = at) 

################################################################################
##                         STAT 331 Final Project                             ##
##                            Daniel Matheson                                 ##
##                                20270871                                    ##
################################################################################

options(scipen=3) # sets the threshold for scientific notation
strikedata <- read.csv("strikes_clean.csv") # import data
## strikedata:
# Country: 18 countries
# Year: every year from 1951 to 1985
# Strike: number of days in the year lost per 1,000 workers due to strikes
# Unemp: Unemployment rate (%)
# Infl: Inflation rate (%)
# Demo: A democracy index, defined as proportion of left-party parliamentary 
# representation
# Centr: Measure of union centralization, which refers to "the authority that 
# union condeferations have over their members". The higher this value,
# the more powerful the union. The measure is aggregated over all years in a 
# given country.
# Dens: Trade union density, which is the fraction of wage earners in the 
# country who belong to a trade union.

################################################################################
########################### VARIABLE MANIPULATION ##############################
# Creating the Region co-variate:

countries_temp <- strikedata$Country 
wEurope <- c("Austria","Belgium","France","Germany","Ireland","Italy",
             "UnitedKingdom","Switzerland") 
strikedata$Region[ (countries_temp %in% wEurope)] <- "Europe"

north_america <- c("Canada","UnitedStates")
strikedata$Region[(countries_temp %in% north_america)] <- "NorthAmerica"

ausnz <- c("Australia", "NewZealand")
strikedata$Region[(countries_temp %in% ausnz)] <- "AusNZ"

strikedata$Region[(countries_temp == "Japan")] <- "Japan"

scandanavia <- c("Denmark","Finland","Netherlands","Norway","Sweden")
strikedata$Region[(countries_temp %in% scandanavia)] <- "Scandanavia"

# Taking log(Strike + 1)
strikedata$Strike <- log(strikedata$Strike + 1)

################################################################################
#######################  AUTOMATED MODEL SELECTION ############################# 
M0 <- lm(Strike ~ 1, data = strikedata)           # Minimal Model
Mfull <- lm(Strike ~ (.)^2, data = strikedata)    # Maximal Model

## Foward Model
Mfwd <- step(object = M0, scope = list(lower = M0, upper = Mfull),
             direction = "forward", trace = F)
## Backward Model
Mback <- step(object = Mfull, scope = list(lower = M0, upper = Mfull),
              direction = "backward", trace = F)
## Stepwise Model
Mstart <- lm(Strike ~ ., data = strikedata)       # Stepwise starting model
Mstep <- step(object = Mstart, scope = list(lower = M0, upper = Mfull),
              direction = "both", trace = F)

models <- c(Mfwd$call, Mback$call, Mstep$call)    # Summarize Results
names(models) <- c("FWD", "BACK","STEP")
 
## Model which has the co-variates common to Mfwd, Mback, Mstep
## Used to perform three ANOVA tests (We will omit these as they are long)
Mtest1 <- lm(Strike ~ Country + Year + Infl + Dens + Country:Year,
             data = strikedata)

################################################################################
########################## MODEL DEFINITIONS ###################################
Mquant <- lm(Strike ~ Country + Year + Infl + Unemp + Dens + Country:Year
             + Country:Unemp + Year:Dens + Unemp:Dens, data = strikedata)
Mqual <- lm(Strike ~ Country + Infl + Dens:Centr + Year - 1, data = strikedata)
Mqualregion<- lm(Strike ~ Region + Year + Infl + Dens:Centr + Year -1,
                 data = strikedata) # Not used

################################################################################
######################## 2. MODEL DIAGNOSTICS ##################################

### PRESS residuals
Mquant.h <- hatvalues(Mquant)                    # Hat values
Mqual.h <- hatvalues(Mqual)

Mquant.press <- resid(Mquant)/(1-Mquant.h)       # PRESS residuals
Mqual.press <- resid(Mqual)/(1-Mqual.h)

Mquant.press.sq <- sum(Mquant.press^2)           # Sum-of-squared Press 
Mqual.press.sq <- sum(Mqual.press^2)               # residuals

press.list <- c(Mquant.press.sq, Mqual.press.sq)     # Combined for nice output

### DFFITS residuals
Mquant.DFFITS <- dffits(Mquant)                  # DFFITS residuals
Mqual.DFFITS <- dffits(Mqual)

Mquant.DFFITS.sq <- sum(Mquant.DFFITS^2)         # Sum-of-squared DFFITS
Mqual.DFFITS.sq <- sum(Mqual.DFFITS^2)            # residuals

DFFITS.list <- c(Mquant.DFFITS.sq, Mqual.DFFITS.sq)  # Combined for nice output

#### Akaike Information Criterion (AIC)
Mquant.AIC <- AIC(Mquant)                        # AIC
Mqual.AIC <- AIC(Mqual)
AIC.list <- c(Mquant.AIC, Mqual.AIC)             # Combined for nice output

#### R^2 
Mquant.R2 <- summary(Mquant)$r.squared           # R squared values
Mqual.R2 <- summary(Mqual)$r.squared
R2.list <- c(Mquant.R2, Mqual.R2)                # Combined for nice output

#### Final Output as a Matrix (Bottom of Page 6) ###
diagnost.matrix <- rbind(press.list,DFFITS.list,AIC.list, R2.list)
row.names(diagnost.matrix) <- c("Sum-of-squared Press", "Sum-of-squared DFFITS",
                                "Akaike Information Criterion", "R^2")
colnames(diagnost.matrix) <- c("Quantitative Model", "Qualitative Model")

################################################################################
##################### MODEL SELECTION: Cross Validation ########################
M1 <- Mquant     # models                      
M2 <- Mqual
nreps <- 8000                           # number of replications
ntot <- nrow(strikedata)                # total number of observations
ntrain <- 600                           # size of training set
ntest <- ntot-ntrain                    # size of test set
sse1 <- rep(NA, nreps)                  # sum-of-square errors for
sse2 <- rep(NA, nreps)                     #  each CV replication
Lambda <- rep(NA, nreps)        # likelihod ratio statistic for each replication
system.time({     # measures time taken to perform analysis
  for(ii in 1:nreps) {
    if(ii%%400 == 0) message("ii = ", ii)  # progress indicator
    # randomly select training observations
    train.ind <- sample(ntot, ntrain) # training observations
    # predict from training observations
    M1.cv <- update(M1, subset = train.ind)
    M2.cv <- update(M2, subset = train.ind)
    # test models with trained estimates
    M1.res <- strikedata$Strike[-train.ind] - 
      predict(M1.cv, newdata = strikedata[-train.ind,])
    M2.res <- strikedata$Strike[-train.ind] - 
      predict(M2.cv, newdata = strikedata[-train.ind,])
    # total sum of square errors (applying exp to return to original scale)
    sse1[ii] <- sum((exp(M1.res))^2)
    sse2[ii] <- sum(exp((M2.res))^2)
    # testing likelihood ratio
    M1.sigma <- sqrt(sum(resid(M1.cv)^2)/ntrain) # MLE of sigma
    M2.sigma <- sqrt(sum(resid(M2.cv)^2)/ntrain)
    Lambda[ii] <- sum(dnorm(M1.res, mean = 0, sd = M1.sigma, log = TRUE))
    Lambda[ii] <- Lambda[ii] - sum(dnorm(M2.res, mean = 0, sd = M2.sigma,
                                         log = TRUE))
  }})

# SSE.list <- c(SSE1 = mean(sse1), SSE2 = mean(sse2)) # not used right now

######## FIGURE 8: Cross validation with Histogram #######

par(mfrow = c(1,2), mar = c(4.5, 4.5, .1, .1))     # set up graph frame
# boxplots:
boxplot(x = list(sse1, sse2), names = c("Quantitative","Qualitative"), cex = .7,
        ylab = expression(SS[err]^{test}), col = c("yellow", "orange"))
# histogram:
hist(Lambda, breaks = 50, freq = FALSE, xlab = expression(Lambda^{test}),
     main = "", cex = .7)
abline(v = mean(Lambda), col = "red") # overlays average Lambda value
legend("topleft", legend = c(paste("mean(Lambda) = ", 
                                   as.character(round(mean(Lambda),2)))),
       cex = 0.8, text.col = "red")   # legend indicating Lambda

################################################################################
####### 2.1: Paramater Estimates and Confidence Intervals Of Final Model #######

Mqual.coefdata <- summary(Mqual)$coefficients # Pulls estimates, sigmas, F-stats 
                                              # and p-values from summary
final.estimates <- Mqual.coefdata[,1]      # estimates
final.sigmas <- Mqual.coefdata[,2]         # sigmas
final.CIs <- paste("[", round(final.estimates - 1.96 * abs(final.sigmas), 3),
                   ",", round(final.estimates + 1.96 * abs(final.sigmas), 3), 
                   "]")                    # creates 95% confidence intervals 
                                            # for all betas
final.estimates <- round(final.estimates, 3)  # rounds estimates for output
final.matrix <- data.frame(as.numeric(final.estimates),final.CIs, 
                           stringsAsFactors = FALSE)  # arrange for output
colnames(final.matrix) <- c("Estimate", "95% Confidence Interval") 
rownames(final.matrix) <- rownames(Mqualtest.coefdata)

################################################################################
################################ FIGURES #######################################


   #### FIGURE 1 & 2: Pairs Plots
# Figure 1 before log(Strike), Figure 2 after
pairs(~ Strike + Year + Infl + Unemp + Demo + Centr + Dens, data = strikedata)

   #### FIGURE 3 & 7: Residual Fit and Histogram
M <- Mquant # Mquant for Figure 3, Mqual for Figure 7.

M.stand.res <- resid(M)/summary(M)$sigma        # Standardized Residuals
M.stud.res <- M.stand.res/sqrt(1-hatvalues(M))  # Studentized Residuals
M.press <- resid(M)/(1-hatvalues(M))            # Press Residuals
M.dffits <- dffits(M)                           # DFFITS residuals

cex = 0.9 # dot size
par(mfrow = c(1,2), mar = c(4,4,2,2))           # sets up graph frame
plot(x = predict(M), y = rep(0,length(predict(M))), type = "n", # empty plot
     ylim = range(M.stand.res,M.stud.res,M.press,M.dffits),      # for all the
     xlab = "Predicted Values", ylab = "Residuals")              # points
segments(x0 = predict(M),
         y0 = pmin(M.stand.res, M.stud.res, M.press, M.dffits), # lines between
         y1 = pmax(M.stand.res, M.stud.res, M.press, M.dffits),  # residuals
         lty = 2)
points(predict(M), M.stand.res, pch = 21, bg = "black", cex = cex) # draws
points(predict(M), M.stud.res, pch = 22, bg = "blue", cex = cex)    # the
points(predict(M), M.press, pch = 23, bg = "red", cex = cex)        # points
points(predict(M), M.dffits, pch = 24, bg = "orange", cex = cex)
legend("topright", legend = c("Standardized", "Studentized", "PRESS", "DFFITS"),
       pch = c(21,22,23,24), pt.bg = c("black", "blue", "red", "orange"),
       title = "Residual Type:", cex = 0.7, pt.cex = 1)            # legend
hist.resid <- M.stand.res # Change residual to be used for Histogram
hist(hist.resid, breaks = 50, freq = F, cex.axis = 0.8,
     xlab = "Standardized Residuals", main = "")  # Histogram
curve(dnorm(x),col = "red", add = T)              # Overlay Normal curve
legend("topright", legend = c("Normal PDF"), lty = 1, # Legend for Normal curve
       col = "red", cex = 0.7)


   #### FIGURE 4: Linear Relationships with log(Strike) 
M.strike.year <- lm(Strike ~ Year, data = strikedata)   # Linear regressions to
M.strike.infl <- lm(Strike ~ Infl, data = strikedata)   # plot the linear 
M.strike.unemp <- lm(Strike ~ Unemp, data = strikedata) # relationships on top
M.strike.centr <- lm(Strike ~ Centr, data = strikedata) # of the scatter plots
M.strike.dens <- lm(Strike ~ Dens, data = strikedata)

cex4 = 1.4
par(mfrow = c(3,2), mar = c(4,4,2,2))                 # set up the graph frame
plot(strikedata$Year, strikedata$Strike, xlab = "Year",    # plots vvvv
     ylab = "log Strike days", cex.axis = cex4, cex.lab = cex4)
abline(M.strike.year , col = "red") # overlays the linear regression line
plot(strikedata$Infl, strikedata$Strike, xlab = "Inflation",
     ylab = "log Strike days", cex.axis = cex4, cex.lab = cex4)
abline(M.strike.infl, col = "red")  # overlays the linear regression line
plot(strikedata$Unemp, strikedata$Strike, xlab = "Unemployment",
     ylab = "log Strike days", cex.axis = cex4, cex.lab = cex4)
abline(M.strike.unemp, col = "red") # overlays the linear regression line
plot(strikedata$Centr, strikedata$Strike, xlab = "Union Centralization",
     ylab = "log Strike days", cex.axis = cex4, cex.lab = cex4)
abline(M.strike.centr, col = "red") # overlays the linear regression line
plot(strikedata$Dens, strikedata$Strike, xlab = "Union Density",
     ylab = "log Strike days", cex.axis = cex4, cex.lab = cex4)
abline(M.strike.dens, col = "red")  # overlays the linear regression line

    #### FIGURE 5: Linear Relationships between co-variates
M.year.infl <- lm(Infl ~ Year, data = strikedata)   # similar as above, linear
M.year.unemp <- lm(Unemp ~ Year, data = strikedata) # regressions to plot
M.year.dens <- lm(Dens ~ Year, data = strikedata)   # over the scatter plots
M.centr.dens <- lm(Centr ~ Dens, data = strikedata)

cex5 = 1.3
par(mfrow = c(2,2), mar = c(4,4,2,2))          # sets up the graph frame
plot(strikedata$Year, strikedata$Infl, xlab = "Year", ylab = "Inflation(%)",
     cex.axis = cex5, cex.lab = cex5)
abline(M.year.infl , col = "red") # overlays the linear regression line
plot(strikedata$Year, strikedata$Unemp, xlab = "Year", ylab = "Unemployment(%)",
     cex.axis = cex5, cex.lab = cex5)
abline(M.year.unemp, col = "red") # overlays the linear regression line
plot(strikedata$Year, strikedata$Dens, xlab = "Year", ylab = "Union Density(%)",
     cex.axis = cex5, cex.lab = cex5)
abline(M.year.dens, col = "red")  # overlays the linear regression line
plot(strikedata$Dens, strikedata$Centr, xlab = "Union Density(%)",
     ylab = "Union Centralization", cex.axis = cex5, cex.lab = cex5)
abline(M.centr.dens, col = "red") # overlays the linear regression line


    #### FIGURE 6: Leverage vs. Influence graphs 
Mquant.h <- hatvalues(Mquant)                    # Hat values (leverages)
Mqual.h <- hatvalues(Mqual)

Mquant.h.bar <- mean(Mquant.h)                   # Average leverages
Mqual.h.bar <- mean(Mqual.h)

Mquant.highlev <- (Mquant.h >= 2*Mquant.h.bar)   # indices that are more than
Mqual.highlev <- (Mqual.h >= 2*Mqual.h.bar)       # twice the average leverage

Mquant.cookD <- cooks.distance(Mquant)           # cook's distance or cook's 
Mqual.cookD <- cooks.distance(Mqual)              # influence measure

Mquant.high.cookD <- (Mquant.cookD >=        # indices of top 15 influence obs
                        quantile(Mquant.cookD, probs = (610/625)))
Mqual.high.cookD <- (Mqual.cookD >=          # indices of top 15 influence obs
                       quantile(Mqual.cookD, probs = (610/625)))

### Figure starts here:
clrs <- rep("black", len = nrow(strikedata))    # empty colors vector
clrs[Mquant.highlev] <- "blue"                  # high leverage in blue
clrs[Mquant.high.cookD] <- "red"                # high influence in red
par(mfrow = c(2,1), mar = c(4,4,2,4), xpd = T)  # set up the graph frame
                                        # xpd = T allows legend outside of plot
plot(Mquant.h, Mquant.cookD, xlab = "Leverage", ylab = "Cook's Distance",
     main = "Quantitative Model", pch = 21, bg = clrs) # plot
abline(v = 2*Mquant.h.bar, col = "grey", lty = 2)   # add 2*mean(leverage) line
legend(x = 0.25, y = 0.061, legend = c("Top 15 Influence Observations",
                                       "Leverage > 2*mean(Leverage)"),
       cex = 0.6, pt.cex = 1.2,                 # legend
       pch = 19, col =c("red", "blue"))
clrs2 <- rep("black", len = nrow(strikedata))   # empty colors vector
clrs2[Mqual.highlev] <- "blue"                  # high leverage in blue
clrs2[Mqual.high.cookD] <- "red"                # high influence in red
plot(Mqual.h, Mqual.cookD, xlab = "Leverage", ylab = "Cook's Distance",
     main = "Qualitative Model", pch = 21, bg = clrs2) # plot
abline(v = 2*Mqual.h.bar, col = "grey", lty = 2) # add 2*mean(leverage) line
legend(x = 0.14, y = 0.06, legend = c("Top 15 Influence Observations",
                                      "Leverage > 2*mean(Leverage)"),
       cex = 0.6, pt.cex = 1.2,                 # legend
       pch = 19, col =c("red", "blue"))

################################################################################
############################### APPENDIX A ##################################### 

#####  IDENTIFYING OUTLIERS 

strike_outliers <- strikedata[(strikedata$Strike > 5000),] # finds indices of
outliers_indices <- as.numeric(row.names(strike_outliers))   # outliers

##### FIGURE 9
M <- Mfwd   

H <- hatvalues(M)                        # Hat values/leverages
M.res <- resid(M)                        # Residuals
M.y.hat <- predict(M)                    # Predicted Values
M.sigma.hat <- summary(M)$sigma          # Sigma hat
M.stand.res <- M.res/M.sigma.hat         # Standardized Residuals
M.stud.res <- M.stand.res/sqrt(1-H)      # Studentized Residuals

color.index <- (strikedata$Strike <= 302.3)       # Setting colors for obs.
pt.col <- rep("black", length(strikedata$Strike)) # with strikes <= mean
pt.col[color.index] <- "red"
pt.col[outliers_indices] <- "blue"                # Blue for 3 large outliers

par(mfrow = c(1,2), mar = c(4,4,2,2))             # Set up graph frame
plot(x = M.y.hat, y = M.stand.res, col = pt.col, xlab = "Predicted Values",
     ylab = "Standardized Residuals")             # plot
legend(x = -150, y = 8, legend = c("Strikes <= 302.3", "Strikes > 302.3",
                                   "Strikes >= 5000"),
       cex = 0.6, y.intersp = 2, xjust = 0.1, pt.cex = 1.2, pch = 19,
       col = c("red", "black","blue"))            # legend
hist(M.stand.res, breaks = 50, freq = FALSE, cex.axis = 0.8,
     xlab = "Standardized Residuals", main = "")  # histogram of residuals
curve(dnorm(x),col = "red", add = TRUE)           # normal curve overlay
legend("topright", legend = c("Normal PDF"), lty = 1, # legend
       col = "red", cex = 0.7)

##### FIGURE 10
cex = 0.9       # point size
par(mfrow = c(1,2), mar = c(4,4,2,2))   # set up graph frame
plot(x = M.y.hat, y= rep(0,length(M.y.hat)), type = "n", # empty plot to put the
     ylim = range(M.stand.res, M.stud.res),               # points onto
     xlab = "Predicted Values", ylab = "Residuals")
segments(x0 = M.y.hat,                                   # lines between the
         y0 = pmin(M.stand.res, M.stud.res),              # residuals
         y1 = pmax(M.stand.res, M.stud.res),
         lty = 2)
points(M.y.hat, M.stand.res, pch = 23, bg = "blue", cex = cex)   # points
points(M.y.hat, M.stud.res, pch = 21, bg = "orange", cex = cex)
legend("topright", legend = c("Standardized", "Studentized"),    # legend
        pch = c(21,23), pt.bg = c("blue", "orange"),
        title = "Residual Type:", cex = 0.7, pt.cex = 1)
hist(M.stand.res, breaks = 50, freq = FALSE, cex.axis = 0.8,
     xlab = "Standardized Residuals", main = "")  # histogram of residuals
curve(dnorm(x),col = "red", add = TRUE)           # normal curve overlay
legend("topright", legend = c("Normal PDF"), lty = 1, # legend
       col = "red", cex = 0.7)

################################################################################
################################ APPENDIX D ####################################

### Models which all add in one co-variate. In this order: Unemp, Dens,
  # Year:Infl, Year:Unemp
Mqualtest1 <- lm(Strike ~ Country + Infl + Dens:Centr + Year - 1 + Unemp,
                 data = strikedata)
Mqualtest2 <- lm(Strike ~ Country + Infl + Dens:Centr + Year - 1 + Dens,
                 data = strikedata)
Mqualtest3 <- lm(Strike ~ Country + Infl + Dens:Centr + Year - 1 + Year:Infl,
                 data = strikedata)
Mqualtest4 <- lm(Strike ~ Country + Infl + Dens:Centr + Year - 1 + Year:Unemp,
                 data = strikedata)

# Get sum-of-squared residuals for all 4 models
test.resid.sq <- c(sum(resid(Mqualtest1)^2),sum(resid(Mqualtest2)^2),
                   sum(resid(Mqualtest3)^2),sum(resid(Mqualtest4)^2))
names(test.resid.sq) <- c("Unemployment", "Union Density", "Year:Inflation",
                          "Year:Unemployment")
# Compare to Mqual's sum-of-squared residuals
test.resid.sq.change <- test.resid.sq - sum(resid(Mqual)^2)

# Get R^2 values for all 4 models
test.rsq <- c(summary(Mqualtest1)$r.squared,summary(Mqualtest2)$r.squared,
             summary(Mqualtest3)$r.squared,summary(Mqualtest4)$r.squared)
# Compare to Mqual's R^2 value
test.rsq.change <- test.rsq - summary(Mqual)$r.squared

# Get sigma hat for all 4 models
test.sigma <- c(summary(Mqualtest1)$sigma,summary(Mqualtest2)$sigma,
                summary(Mqualtest3)$sigma,summary(Mqualtest4)$sigma)
# Compare to Mqual's sigma hat
test.sigma.change <- test.sigma - summary(Mqual)$sigma

# Get p-values for the significance of each additional co-variate in their
# respective models.
test.pval <- c(summary(Mqualtest1)$coefficients[21,4],
               summary(Mqualtest2)$coefficients[21,4],
               summary(Mqualtest3)$coefficients[22,4],
               summary(Mqualtest4)$coefficients[22,4])

### Final Output ###
test.matrix <- rbind(round(test.resid.sq,3), round(test.resid.sq.change,3),
                     round(test.rsq,3), round(test.rsq.change, 3),
                     round(test.sigma,3), round(test.sigma.change,3), 
                     round(test.pval,3))
row.names(test.matrix) <- c("SSQ Resid", "Change", "R^2", "Change","Sigma hat",
                            "Change", "p-value")

################################################################################
############################### APPENDIX E #####################################

Mquant.h <- hatvalues(Mquant)                    # Hat values (leverages)
Mquant.h.bar <- mean(Mquant.h)                   # Average leverage
Mquant.highlev <- (Mquant.h >= 2*Mquant.h.bar)   # High leverage indices
Mquant.cookD <- cooks.distance(Mquant)           # Cook's distance
Mquant.high.cookD <- (Mquant.cookD >=           # top 15 influence observations
                        quantile(Mquant.cookD, probs = (610/625)))
Mqual.high.cookD <- (Mqual.cookD >=             # top 15 influence observations
                       quantile(Mqual.cookD, probs = (610/625)))
# Finding indices which are both high leverage and high influence:
Mquant.remove <- (Mquant.highlev & Mquant.high.cookD)
Mquant.remove.index <- as.numeric(
  row.names(strikedata[(Mquant.highlev & Mquant.high.cookD),]))

# Observations from strikedata which correspond to the indices above
 # i.e high leverage and influence
E.obs <- strikedata[Mquant.remove.index,]    
# The magnitude of the DFFITS residuals of the observations found in E.obs
E.DFFITS <- abs(Mquant.DFFITS[Mquant.remove.index])
# The magnitude of the press residuals of the observations found in E.obs
E.press <- abs(Mquant.press[Mquant.remove.index])
# The unemployment figures for the observations in E.obs
E.unemp <- E.obs[,4]
# The union centralization figures for the observations in E.obs
E.centr <- E.obs[,7]

### Percentiles of DFFITS:
E.DFFITS.CDF <- ecdf(abs(Mquant.DFFITS)) # Empirical CDF function
E.DFFITS.p <- c(E.DFFITS.CDF(E.DFFITS)) 

### Percentiles of PRESS:
E.press.CDF <- ecdf(abs(Mquant.press))   # Empirical CDF function
E.press.p <- c(E.press.CDF(E.press))

### Percentiles of Unemployment 
E.unemp.CDF <- ecdf(strikedata$Unemp)    # Empirical CDF function
E.unemp.p <- c(E.unemp.CDF(E.unemp))

#### FINAL OUTPUT ####
E.matrix <- rbind(E.DFFITS.p, E.press.p, E.centr, E.unemp.p)
row.names(E.matrix) <- c("DFFITS Percentile","Press Percentile",
                         "Union Centralization", "Unemployment Percentile")
colnames(E.matrix) <- c(1,2,3,4,5)

###### FIGURE 12: Pairs Plot with Quantitative Model influencial observations
##### labelled to show that they don't appear to be "special" in any way.
clrs.quant.outliers <- rep("black", nrow(strikedata))
clrs.quant.outliers[Mquant.remove.index] <- "red"
pairs(~ Strike + Year + Infl + Unemp + Dens + Year:Dens + Unemp:Dens,
      col = clrs.quant.outliers, data = strikedata, cex = 1.3)

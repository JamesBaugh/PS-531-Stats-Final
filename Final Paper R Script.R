library(car)
library(margins)
library(ProbMarg)
library(VGAM)
library(haven)
library(glmnet)
library(arsenal)
library(selectiveInference)
library(DeclareDesign)
library(pwr)

ubdata <- read_dta("C:/Users/knigh/Desktop/Stata Work/Underbanking/ubanking.immigration.v2.dta")

# recoding a variable for international remittances
attach(ubdata)
intremit2 <- ifelse(hes130==2 | intremit==0, 0,
                    ifelse(hes130==1 & intremit==1, 1, 99))
# creating the binary variable for underbanking with the new int. remit. variable
ubdata$ubdummy2 <- ifelse(revanyaccount==1 | paydayloan==1 | checkcashing==1 | moneyorder==1 | 
                            pawnshop==1 | autoloan==1 | taxrefund==1 | renttoown==1 | 
                            intremit2==1, 1, ifelse(revanyaccount==0 & paydayloan==0 &
                                                      checkcashing==0 & moneyorder==0 &
                                                      pawnshop==0 & autoloan==0 &
                                                      taxrefund==0 & renttoown==0 &
                                                      intremit2==0, 0, 99))

# removing observations with NA's in the underbanking indicator
NArows <- c(which(is.na(ubdata$ubdummy2), arr.ind=TRUE))
newdata <- ubdata[-c(NArows),]
detach(ubdata)

X <- data.matrix(newdata[, c("fgen_developing", "fgen_developed",
                             "sgen_developed", "sgen_developing",
                             "black", "hispanic1", "asian", "other1",
                             "lessthanhs", "hs", "somecollege",
                             "famincmean",
                             "unemployed", "nilf",
                             "age",
                             "male",
                             "married", "wsd_married",
                             "numchildren2")])

lambdavec <- c(rep(0, 100))

set.seed(123456789)
for(i in 1:100){
  cv.lambda <- cv.glmnet(X, newdata$ubdummy2, alpha=1)
  lambdavec[i] <- cv.lambda$lambda.min
}

hist(lambdavec, main="Histogram of Calculated Penalties")

meanlam <- mean(lambdavec)

lassoreg <- glmnet(X, newdata$ubdummy2, alpha = 1, lambda = meanlam)

lassonames1 <- c("First-generation, developing country",
                 "First-generation, developed country",
                 "Second-generation, developing country",
                 "Second-generation, developed country",
                 "Black", "Hispanic", "Asian", "Other",
                 "Less Than High School", "High School",
                 "Some College",
                 "Annual Income",
                 "Unemployed", "Not in Labor Force",
                 "Age", "Sex",
                 "Married", "W/S/D",
                 "Num. of Children")
cbind(Covariates=lassonames1, Estimates=round(coef(lassoreg)[2:20],5))


X2 <- data.matrix(newdata[, c("fgen_developing",
                              "sgen_developed", "sgen_developing",
                              "black", "hispanic1", "asian", "other1",
                              "lessthanhs", "hs", "somecollege",
                              "famincmean",
                              "unemployed", "nilf",
                              "age",
                              "married", "wsd_married",
                              "numchildren2")])

lambdavec2 <- c(rep(0, 100))

set.seed(123456789)
for(i in 1:100){
  cv.lambda2 <- cv.glmnet(X2, newdata$ubdummy2, alpha=1)
  lambdavec2[i] <- cv.lambda2$lambda.min
}

meanlam2 <- mean(lambdavec2)
cbind("Penalty 1"=round(meanlam,5), "Penalty 2"=round(meanlam2,5))

lassoreg2 <- glmnet(X2, newdata$ubdummy2, alpha = 1, lambda = meanlam2)

lassonames2 <- c("First-generation, developing country",
                 "Second-generation, developing country",
                 "Second-generation, developed country",
                 "Black", "Hispanic", "Asian", "Other",
                 "Less Than High School", "High School",
                 "Some College",
                 "Annual Income",
                 "Unemployed", "Not in Labor Force",
                 "Age", "Married", "W/S/D",
                 "Num. of Children")
cbind(Covariates=lassonames2, Estimates=round(coef(lassoreg2)[2:18],5))

# trying to calculate p-values from the LASSO model directly
n <- length(newdata$ubdummy2)
beta <- coef(lassoreg2,s=meanlam2/n)[-1]

fixedLassoInf(X2, newdata$ubdummy2, beta, meanlam2, family=c("gaussian"), type=c("partial"))

# the linear model
linreg <- lm_robust(ubdummy2~fgen_developing+sgen_developing+sgen_developed+
                      black+hispanic1+asian+other1+
                      lessthanhs+hs+somecollege+
                      famincmean+
                      unemployed+nilf+
                      age+
                      married+ wsd_married+
                      numchildren2,
                    data=newdata)


# 18 predictors including the intercept - 1 = 17
# 34243 observations - 18 predictors = 34225
# f2 = R^2 / (1-R^2), R^2 from linreg is 0.1593

pwr.f2.test(u = 17, # 18 predictors minus the intercept
            v = 34225, # N - 18 predictors
            f2 = 0.1593/0.8407, # R-squared / 1 - R-squared
            power = NULL) # function returns the null option


fdr <- p.adjust(p=linreg$p.value, method="BH")
ptable <- cbind(FPR=round(linreg$p.value, 2), FDR=round(fdr,2))
ptable

# Apologies for the raw variable names
cbind(Estimates=round(coef(linreg), 5))

# rmse
rmse <- function(o, p) {
  sqrt(mean((o-p)^2))
}

# mse
mse <- function(o, p) {
  (mean((o-p)^2))
}

# bias
bias <- function(reg, y){
  x<-mean(fitted(reg))
  y<-mean(y)
  result <- x-y
  rbind("Mean Fitted"=x, "Mean Observed"=y, "Bias"=result)
}

cbind(MSE=mse(newdata$ubdummy2, fitted(linreg)), 
      RMSE=rmse(newdata$ubdummy2, fitted(linreg)))

bias(linreg, newdata$ubdummy2)

predreg <- lm_robust(ubdummy2~famincmean+
                       fgen_developing+sgen_developing+sgen_developed+
                       black+hispanic1+asian+other1+
                       lessthanhs+hs+somecollege+
                       unemployed+nilf+
                       age+
                       married+ wsd_married+
                       numchildren2,
                     data=newdata)
incvec <- c(3,6,9,11,14,18,23,28,33,38,45,55,68,88,125,200)
predgrid <- data.frame(famincmean=rep(incvec, 2),
                       fgen_developing=c(rep(0, 16), rep(1, 16)),
                       sgen_developing=mean(newdata$sgen_developing),
                       sgen_developed=mean(newdata$sgen_developed),
                       black=mean(newdata$black),
                       hispanic1=mean(newdata$hispanic1),
                       asian=mean(newdata$asian),
                       other1=mean(newdata$other1),
                       lessthanhs=mean(newdata$lessthanhs),
                       hs=mean(newdata$hs),
                       somecollege=mean(newdata$somecollege),
                       unemployed=mean(newdata$unemployed),
                       nilf=mean(newdata$nilf),
                       age=mean(newdata$age),
                       married=mean(newdata$married),
                       wsd_married=mean(newdata$wsd_married),
                       numchildren2=mean(newdata$numchildren2))


predline1 <- predict(predreg, predgrid)

plot(predline1[1:16], type="b", lty=2, pch=7,
     xlab="Annual Income (Thousands)",
     ylab="Predicted Values",
     ylim=c(0.1,0.5),
     main="Figure: Predicted Underbanking Values by Income",
     xaxt="n")
lines(predline1[17:32], type="b", lty=1, pch=17)
axis(1,at=seq(1,16,1),labels=incvec,las=0)
legend("bottomleft",legend=c("First Gen., Developing", "Citizen/Other"),
       lty=c(1,2), pch=c(17,7))


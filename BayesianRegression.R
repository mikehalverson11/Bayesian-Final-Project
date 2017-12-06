install.packages("invgamma")
library(invgamma)
library(MASS)
library(tidyverse)
# set up
numberOfSims=1
numberOfVariables=2
specifiedVariables=c("Intercept","Height")
response="REB"
a0=.1
b0=.001
file="C:/Users/Michael/Downloads/nbadata.csv"

# initialization
sds=rinvgamma(numberOfSims,a0,b0)
meanguess=runif(numberOfVariables)
# make sure invertible covariance
#covariance0=matrix(c(2,.5,1,.5,5,1.5,1,1.5,5),numberOfVariables,numberOfVariables)
covariance0=matrix(c(1,0,0,1),numberOfVariables)
#covariance0=cov(data[,specifiedVariables])
betaorig=mvrnorm(numberOfSims,meanguess,sds*covariance0)
data=read_csv(file)
data$Intercept=1
data=data[!(is.na(data[,"Height"])),]

# posteriors
theta1=solve(solve(covariance0)+t(as.matrix(data[specifiedVariables]))%*%as.matrix(data[,specifiedVariables]))%*%
  (solve(covariance0)%*%betaorig+
     t(as.matrix(data[,specifiedVariables]))%*%as.matrix(data[,response]))
covariance1=solve((solve(covariance0)+t(as.matrix(data[specifiedVariables]))
                   %*%as.matrix(data[,specifiedVariables])))
a1=a0+.5*nrow(data)
b1=b0+.5*(t(betaorig)%*%solve(covariance0)%*%betaorig+as.matrix(t(data[,response]))%*%as.matrix(data[,response])-
            t(theta1)%*%solve(covariance1)%*%theta1)

# metrics 
predictions=t(theta1)%*%t(as.matrix(data[,specifiedVariables]))
MSE=(1/nrow(data))*sum((predictions-data[,response])^2)
R2=1-sum((predictions-data[,response])^2)/(sum((data[,response]-mean(as.matrix(data[,response])))^2))

# comparison to actual linear regression
fit<-lm(as.formula(paste(response,"~",paste(specifiedVariables[specifiedVariables!="Intercept"],collapse = "+"),sep="")),data=data)
summary(fit)$r.squared;R2
mean(fit$residuals^2);MSE
fit$coefficients;theta1

# confidence intervals
sims1=10000
sds1=rinvgamma(sims1,a1,b1)
mean(sds1);sort(sds1)[.025*sims1];sort(sds1)[.975*sims1]
Intercept = rep(0, 10000)
MIN = rep(0, 10000)
PTS = rep(0, 10000)
betafinal = data_frame(Intercept, MIN, PTS)
for(j in 1:length(sds1)){
  out = mvrnorm(1,theta1,sds1[j]*covariance1)
  betafinal$Intercept[j] = out[1]
  betafinal$MIN[j] = out[2]
  betafinal$PTS[j] = out[3]
}

for (name in 1:length(colnames(betafinal))){
  message(colnames(betafinal)[name])
  message("Mean: ",round(mean(betafinal[[name]]), 4))
  message("SD: ", round(sd(betafinal[[name]]), 4))
  message("95% Credible Interval: ", "(",round(sort(betafinal[[name]])[.025*sims1], 4),
          ", ",round(sort(betafinal[[name]])[.975*sims1], 4), ")" )
}

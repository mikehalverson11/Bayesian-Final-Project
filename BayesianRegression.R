install.packages("invgamma")
library(invgamma)
library(MASS)
# set up
numberOfSims=1
numberOfVariables=3
specifiedVariables=c("Intercept","MIN","PTS")
response="REB"
a0=.1
b0=.001
file="C:/Users/Michael/Downloads/nbadata.csv"

# initialization
sds=rinvgamma(numberOfSims,a0,b0)
meanguess=runif(numberOfVariables)
# make sure invertible covariance
covariance0=matrix(c(2,0,0,0,5,0,0,0,5),numberOfVariables,numberOfVariables)
#covariance0=cov(data[,specifiedVariables])
betaorig=mvrnorm(numberOfSims,meanguess,sds*covariance0)
data=read_csv(file)
data$Intercept=1
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
summary(fit);R2
mean(fit$residuals^2);MSE
fit$coefficients;theta1

# confidence intervals
sims1=10000
sds1=rinvgamma(sims1,a1,b1)
mean(sds1);sort(sds1)[.025*sims1];sort(sds1)[.975*sims1]
betafinal=mvrnorm(sims1,theta1,sds*covariance1)
for (i in 1:length(specifiedVariables)){
  current=specifiedVariables[i]
  message(current)
  message(mean(betafinal[,current]))
  message(sd(betafinal[,current]))
  message(sort(betafinal[,current])[.025*sims1])
  message(sort(betafinal[,current])[.975*sims1])
}
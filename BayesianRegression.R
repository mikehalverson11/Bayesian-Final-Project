install.packages("invgamma")
library(invgamma)
library(MASS)
numberOfSims=1
numberOfVariables=3
specifiedVariables=c("Intercept","MIN","PTS")
response="REB"
a0=.01
b0=.01
sds=rinvgamma(numberOfSims,.01,.01)
meanguess=rep(0,numberOfVariables)
covariance0=matrix(c(1,0,0,0,1,0,0,0,1),numberOfVariables,numberOfVariables)
betaorig=mvrnorm(numberOfSims,meanguess,sds*covariance0)
data=read_csv("C:/Users/Michael/Downloads/nbadata.csv")
data$Intercept=1
theta1=solve(solve(covariance0)+t(as.matrix(data[specifiedVariables]))%*%as.matrix(data[,specifiedVariables]))%*%
  (solve(covariance0)%*%betaorig+
     t(as.matrix(data[,specifiedVariables]))%*%as.matrix(data[,response]))
covariance1=solve((solve(covariance0)+t(as.matrix(data[specifiedVariables]))
                   %*%as.matrix(data[,specifiedVariables])))
a1=a0+.5*nrow(data)
b1=b0+.5*(t(betaorig)%*%solve(covariance0)%*%betaorig+as.matrix(t(data[,response]))%*%as.matrix(data[,response])-
            t(theta1)%*%solve(covariance1)%*%theta1)
predictions=t(theta1)%*%t(as.matrix(data[,specifiedVariables]))
MSE=(1/nrow(data))*sum((predictions-data[,response])^2)
R2=1-sum((predictions-data[,response])^2)/(sum((data[,response]-mean(as.matrix(data[,response])))^2))

fit<-lm(as.formula(paste(response,"~",paste(specifiedVariables[specifiedVariables!="Intercept"],collapse = "+"),sep="")),data=data)
summary(fit);R2
mean(fit$residuals^2);MSE
fit$coefficients;theta1

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
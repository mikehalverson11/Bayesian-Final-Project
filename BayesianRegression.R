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
betaorig=mvrnorm(numberOfSims,meanguess,covariance0)
data=read_csv("C:/Users/Michael/Downloads/nbadata.csv")
data$Intercept=1
theta1=solve(solve(covariance0)+t(as.matrix(data[specifiedVariables]))%*%as.matrix(data[,specifiedVariables]))%*%
  (solve(covariance0)%*%betaorig+
     t(as.matrix(data[,specifiedVariables]))%*%as.matrix(data[,response]))
covariance1=solve((solve(covariance0)+t(as.matrix(data[specifiedVariables]))
                   %*%as.matrix(data[,specifiedVariables])))
a1=a0+.5*nrow(data)
b1=b0+.5*(t(betaorig)%*%solve(covariance0)%*%betaorig+as.matrix(t(data[,response]))%*%as.matrix(data[,response])-
            t(theta1)%*%solve(covariance0)%*%theta1)
predictions=t(theta1)%*%t(as.matrix(data[,specifiedVariables]))
MSE=(1/nrow(data))*sum((predictions-data[,response])^2)

fit<-lm(as.formula(paste(response,"~",paste(specifiedVariables[specifiedVariables!="Intercept"],collapse = "+"),sep="")),data=data)
summary(fit)
mean(fit$residuals^2);MSE
fit$coefficients;theta1





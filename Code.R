#importing the libraries
library(TSA)
library(tseries)
library(fUnitRoots)
library(lmtest)
library(forecast)
library(FitAR)
source('./Documents/RMIT_Sem_3/Time Series Analysis/Assignment-2/sort.score.R')

#Task 1
#Read the data
eggs <- read.csv("./Documents/RMIT_Sem_3/Time Series Analysis/Assignment-2/eggs.csv")
head(eggs)
eggs$year <- NULL

#convert the data to timeseries format
eggs <- ts(as.vector(eggs), start=1981, end=1996)
#Checking whether the dataset is in timeseries format
is.ts(eggs)

#Plotting the time series plot
plot(eggs,type='o', xlab="Time in Years",ylab='eggs depositions in millions' ,main='Egg depositions from 1981-1996')

#To check the thickness of the ozone layer thickness over the year shown using correlation value
x = eggs
y = zlag(eggs)
index = 2:length(y)
cor(x[index], y[index])

#Scatter plot
plot(y=eggs,x=zlag(eggs),col=c("red"),xlab = "previous years egg depositions",main = "Scatter Plot of Eggs Depositions")



#Plot the acf and pacf plot
par(mfrow = c(1,2))
acf(eggs)
pacf(eggs)
par(mfrow = c(1,1))

#adf test
adf.test(eggs)

#### From the acf and pacf we can conclude that the TimeSeries is Non stationary 
# and there is a trend in the series 

qqnorm(eggs)
qqline(eggs, col = 2)



#Due to lot of variation in the data we perform Box-Cox Transformations
eggs_box = BoxCox.ar(eggs, method = 'yule-walker')

#Confidence Interval
eggs_box$ci

#drawing the plot with box-cox 
lambda = 0.45 
eggs.bc = ((eggs^lambda)-1)/lambda
plot(eggs.bc, type = "o", ylab = "Egg Depositions (in millions)", main = "Transformed - Egg depositions between 1981 and 1996")

qqnorm(eggs.bc)
qqline(eggs.bc, col = 2)
shapiro.test(eggs.bc)

###### There is a little change in the variance


#For making the non stationary timeseries stationary, we are using 1st differencing.
eggs.diff1 <- diff(eggs, differences=1)
plot(eggs.diff1, ylab='Change in eggs depositions', type='o')

# -------- ADF Test --------
#We are using ADF test to check whether the timeseries is stationary or non-stationary.
#Hypothesis (Null value) -- Time Series is Non-stationary.

#if the p-value of the test is less than the significance level (0.05) 
#then you reject the null hypothesis and infer that the time series is indeed stationary.

order  = ar(diff(eggs.diff1))$order # To pass the order to adfTest function
adfTest(eggs.diff1, lags = order,  title = NULL,description = NULL)


#differencing 2
eggs.diff2 <- diff(eggs, differences=2)
plot(eggs.diff2, ylab='Change in eggs depositions', type='o')

order  = ar(diff(eggs.diff2))$order # To pass the order to adfTest function
adfTest(eggs.diff2, lags = order,  title = NULL,description = NULL)


#differencing 3
eggs.diff3 <- diff(eggs, differences=3)
plot(eggs.diff3, ylab='Change in eggs depositions', type='o')

order  = ar(diff(eggs.diff3))$order # To pass the order to adfTest function
adfTest(eggs.diff3, lags = order,  title = NULL,description = NULL)

#differencing 4
eggs.diff4 <- diff(eggs, differences=4)
plot(eggs.diff4, ylab='Change in eggs depositions', type='o')

order  = ar(diff(eggs.diff4))$order # To pass the order to adfTest function
adfTest(eggs.diff4, lags = order,  title = NULL,description = NULL)



##### Hence we confirm that the time series is stationary as the value of p = 0.0443


#For making the non stationary timeseries stationary, we are using 1st differencing.
#eggs.bc.diff2 <- diff(eggs.bc, differences=2)
#plot(eggs.bc.diff2, ylab='Change in eggs depositions', type='l')
#adf.test(eggs.bc.diff2)

#Plot the acf and pacf plot
par(mfrow = c(1,2))
acf(eggs.diff4)
pacf(eggs.diff4)
par(mfrow = c(1,1))

#ARIMA(1,4,1)

#eacf
eacf(eggs.diff4,ar.max = 2, ma.max = 2) 

#Answers - (0,4,1), (0,4,2), (1,4,1), (1,4,2), (2,4,1), (2,4,2)

#BIC Table
par(mfrow=c(1,1))
res3 = armasubsets(y=eggs.diff4,nar=2,nma=3,y.name='test',ar.method='ols')
plot(res3)

#Answers - (1,4,1)
-----------------------------------------------------------
#ARIMA(0,4,1)
model_041_css = arima(eggs,order=c(0,4,1),method='CSS')
coeftest(model_041_css)

model_041_ml = arima(eggs,order=c(0,4,1),method='ML')
coeftest(model_041_ml)

#ARIMA(0,4,2)
model_042_css = arima(eggs,order=c(0,4,2),method='CSS')
coeftest(model_042_css)

model_042_ml = arima(eggs,order=c(0,4,2),method='ML')
coeftest(model_042_ml)

#ARIMA(1,4,1)
model_141_css = arima(eggs,order=c(1,4,1),method='CSS')
coeftest(model_141_css)

model_141_ml = arima(eggs,order=c(1,4,1),method='ML')
coeftest(model_141_ml)

#ARIMA(1,4, 2)
model_142_css = arima(eggs,order=c(1,4,2),method='CSS')
coeftest(model_142_css)

model_142_ml = arima(eggs,order=c(1,4,2),method='ML')
coeftest(model_142_ml)

#ARIMA(2,4,1)
model_241_css = arima(eggs.bc,order=c(2,4,1),method='CSS')
coeftest(model_241_css)

model_241_ml = arima(eggs,order=c(2,4,1),method='ML')
coeftest(model_241_ml)




# AIC and BIC values

sort.score(AIC(model_041_ml,model_042_ml,model_141_ml,model_142_ml,model_241_ml), score = "aic")
sort.score(BIC(model_041_ml,model_042_ml,model_141_ml,model_142_ml,model_241_ml), score = "bic")


#Model Diagnostic

residual.analysis <- function(model, std = TRUE,start = 2, class = c("ARIMA","GARCH","ARMA-GA
RCH")[1]){
  # If you have an output from arima() function use class = "ARIMA"
  # If you have an output from garch() function use class = "GARCH"
  # If you have an output from ugarchfit() function use class = "ARMA-GARCH"
  library(TSA)
  library(FitAR)
  if (class == "ARIMA"){
    if (std == TRUE){
      res.model = rstandard(model)
    }else{
      res.model = residuals(model)
    }
  }else if (class == "GARCH"){
    res.model = model$residuals[start:model$n.used]
  }else if (class == "ARMA-GARCH"){
    res.model = model@fit$residuals
  }else {
    stop("The argument 'class' must be either 'ARIMA' or 'GARCH' ")
  }
  par(mfrow=c(3,2))
  plot(res.model,type='o',ylab='Standardised residuals', main="Time series plot of standardis
       ed residuals")
  abline(h=0)
  hist(res.model,main="Histogram of standardised residuals")
  acf(res.model,main="ACF of standardised residuals")
  pacf(res.model,main="PACF of standardised residuals")
  qqnorm(res.model,main="QQ plot of standardised residuals")
  qqline(res.model, col = 2)
  print(shapiro.test(res.model))
  k=0
  LBQPlot(res.model, lag.max = 10, StartLag = k + 1, k = 0, SquaredQ = FALSE)
}
residual.analysis(model = model_042_css)

#Forecasting

forecasting = Arima(eggs,c(0,4,2))
plot(forecast(forecasting, h=5), ylim =c(-2,3))
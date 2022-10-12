pacman::p_load(tidyverse, lubridate, tseries, forecast, haven, fma, expsmooth, lmtest, zoo, seasonal, ggplot2, seasonalview, 
               aTSA, imputeTS, prophet)
rm(list=ls())

# Bring in data
fulldata1 <- read.csv("hrl_load_metered.csv")
test1 <- read.csv("hrl_load_metered - test1.csv")
test2 <- read.csv("hrl_load_metered - test2.csv")

fulldata<- rbind(fulldata1, test1)

fulldata <- fulldata %>% dplyr::select(datetime_beginning_ept, mw)

fulldata$datetime <- as.POSIXct(fulldata$datetime_beginning_ept, format="%m/%d/%Y %H:%M")
as.Date(electricity_train$datetime)#this grabs just the date portion

# Check for missing values
fulldata %>%
  filter(is.na(mw)) %>% View(.)
# No missing values

# How often is the load 0?
fulldata %>%
  filter(mw == 0) %>% View(.)
# Twice: 3/8/20 at 1:00 and 3/14/21 at 1:00

# plot missing values
ggplot_na_distribution(fulldata$mw) + labs(y='mw')

# imputing missing values. ie if 0s, make them NAs
fulldata$mw[fulldata$mw == 0] <- NA
fulldata$mw <- na_interpolation(fulldata$mw, option='spline')

# fulldata$mw[fulldata$mw==0] <-  median(fulldata$mw)


ggplot(data=fulldata, aes(x=as.Date(datetime), y=mw))+geom_line() + ggtitle("Time Plot of Hourly MW") +
  labs(y="Hourly usage MW", x = "Datetime") #+ scale_x_continuous(breaks = seq(min(fulldata$datetime), max(fulldata$datetime), by = 50))

# Use all the data as train as per prof
electricity_train = fulldata[1:27744, ]
# electricity_val = fulldata[26857:27576, ]

#ASSUME WEEKLY SEASONALITY
# Time series object
electricity_ts <- ts(electricity_train$mw, frequency = 24) #24*7 is too high for holt winter models


# STL Decomposition
electricity_stl <- stl(electricity_ts, s.window = 13)
autoplot(electricity_stl) #figure 3

# ESM models

# SINGLE MODEL: MAPE 0.11
electricity_single <- ses(electricity_ts, initial = "simple", h = 7*24)
autoplot(electricity_single)+
  autolayer(fitted(electricity_single),series="Fitted")+ylab("MW") + geom_vline(xintercept = 8/1/2021,color="orange",linetype="dashed")

single_error = test2$mw - electricity_single$mean
single_MAE <- mean(abs(single_error))
single_MAPE = mean(abs(single_error)/abs(test2$mw))

# HOLT: MAPE 5.98
electricity_holt <- holt(electricity_ts, initial = "optimal", h = 7*24)
#summary(electricity_holt)
plot(electricity_holt)
autoplot(electricity_holt)+
  autolayer(fitted(electricity_holt),series="Fitted")+ylab("MW") + geom_vline(xintercept = 8/1/2021,color="orange",linetype="dashed")

holt_error = test2$mw - electricity_holt$mean
holt_MAE <- mean(abs(holt_error))
holt_MAPE = mean(abs(holt_error)/abs(test2$mw))

# DAMPED HOLT: MAPE 0.31
electricity_dampedholt <- holt(electricity_ts,initial = "optimal", h = 7*24, damped = TRUE)
autoplot(electricity_dampedholt)+
  autolayer(fitted(electricity_dampedholt),series="Fitted")+ylab("MW") + geom_vline(xintercept = 8/1/2021,color="orange",linetype="dashed")

dampedholt_error = test2$mw - electricity_dampedholt$mean
dampedholt_MAE <- mean(abs(dampedholt_error))
dampedholt_MAPE = mean(abs(dampedholt_error)/abs(test2$mw))

# HOLT WINTERS ADDITIVE: MAPE 0.21 #Error frequency is too high if go 24*7
electricity_additive <- hw(electricity_ts, seasonal = "additive", h = 7*24)

autoplot(electricity_additive)+
  autolayer(fitted(electricity_additive),series="Fitted")+ylab("MW") + geom_vline(xintercept = 8/1/2021,color="orange",linetype="dashed")

additive_error = test2$mw - electricity_additive$mean
additive_MAE <- mean(abs(additive_error))
additive_MAPE = mean(abs(additive_error)/abs(test2$mw))

# HOLT WINTERS MULTIPLICATIVE: MAPE 0.12 #Error frequency is too high if go 24*7
electricity_multiplicative <- hw(electricity_ts,seasonal = "multiplicative", h = 7*24)

#Figure 4 in report
autoplot(electricity_multiplicative)+
  autolayer(fitted(electricity_multiplicative),series="Fitted")+ylab("Energy Use (MegaWatts)") +
  geom_vline(color = 'orange', linetype = 'dashed', xintercept = 1157) + xlim(900, 1170) + ggtitle('Forecast of Energy Use with Holt-Winters Model')

mult_error = test2$mw - electricity_multiplicative$mean
mult_MAE <- mean(abs(mult_error))
mult_MAPE = mean(abs(mult_error)/abs(test2$mw))

#ETS method
ets <- ets(electricity_ts)
ets.forecast <-forecast::forecast(ets,h=7*24)
summary(ets)

ets_error = test2$mw - ets.forecast$mean
ets_MAE <- mean(abs(ets_error))
ets_MAPE = mean(abs(ets_error)/abs(test2$mw))

##############################################SEASONAL ARIMA MODELS##################################################

# Canova Hansen Test
# 𝐻0: Deteministic Seasonality Differencing not going to help
# 𝐻a: Stochastic Seasonality Differencing needed

electricity_ts %>%  nsdiffs() #nsdiff results are valid because freq=24 now.

#result from this test not reliable because freq season lenght is 168 (breaks down at 12 for stochastic and 24 for deterministic)
#SO THEN RUN BOTH APPROACHES STOCHASTIC AND DETERMINISTIC AND PICK BEST MAPE
#answer is 1 ie diff needed, use stochansic approach. This means we need differencing to solve. If 0 then no diff needed, use deterministic way

#graphically compare differenced data (weekly) to original
cbind("Hourly MW Usage" = electricity_ts,
      "Annual change in MW Usage" = diff(electricity_ts, 7*24)) %>%
  autoplot(facets=TRUE) +
  xlab("Time") + ylab("") +
  ggtitle("Comparison of Differenced Data to Original")

#test for regular unit root in data for trend
electricity_ts %>% diff(lag = 7*24) %>% ndiffs()
#0 no diff needed. Reason for lags 24*7 is because of weekly seasonality. 

# Seasonal Dummy Variables
Hour <- rep(0, length(electricity_ts))
Hour <- Hour + 1:168

H <- factor(Week)
H <- relevel(H, ref="168")

M.Matrix <- model.matrix(~H)

Trend <- 1:length(electricity_ts)

SD.ARIMA <- auto.arima(electricity_ts, xreg = M.Matrix[,2:168], method="ML", seasonal = FALSE)
summary(SD.ARIMA)
#Code above will take too long, move on to Fourier

#Fourier Transforms (Harmonic Regression) #MAPE is 0.0923
plots <- list()
for (i in seq(6)) {
  fit <- auto.arima(electricity_ts, xreg = fourier(electricity_ts, K = i),
                    seasonal = FALSE, lambda = NULL)
  plots[[i]] <- autoplot(forecast::forecast(fit,
                                            xreg = fourier(electricity_ts, K=i, h=7*24))) +
    xlab(paste("K=",i,"   BIC=",round(fit$bic,2))) +
    ylab("") + ylim(30000,80000)
}
gridExtra::grid.arrange(
  plots[[1]],plots[[2]],plots[[3]],
  plots[[4]],plots[[5]],plots[[6]], nrow=3)

#Pick k=6 based on lowest BIC
F.ARIMA <- auto.arima(electricity_ts, xreg = fourier(electricity_ts, K = 6), seasonal = FALSE)
summary(F.ARIMA) 
checkresiduals(F.ARIMA)

# Series: electricity_ts 
# Regression with ARIMA(5,1,1) errors 


# Box-Ljung test
Box.test(F.ARIMA$residuals, lag = 24, type = "Ljung")

# data:  F.ARIMA$residuals
# X-squared = 3477.5, df = 24, p-value < 2.2e-16

autoplot(forecast::forecast(F.ARIMA, xreg = fourier(electricity_ts, K = 6))) +
  autolayer(fitted(F.ARIMA), series="Fitted") + ylab("Airlines Passengers") +
  geom_vline(xintercept = 2007.25,color="orange",linetype="dashed")

z <- forecast::forecast(F.ARIMA, xreg = fourier(electricity_ts, K = 6))$mean

F.ARIMA.error <- test2$mw - z
F.ARIMA.MAE <- mean(abs(F.ARIMA.error))
F.ARIMA.MAPE <- mean(abs(F.ARIMA.error)/abs(test2$mw))

#SEASONAL ARIMA DIFFERENCING APPROACH - model the difference of weekly lag

electricity_ts %>% diff(lag = 24) %>% ggtsdisplay() #this gives you seasonal ACF and PACF

# Manual Way 
# S.ARIMA <- Arima(electricity_ts, order=c(1,0,0), seasonal=c(1,1,1))


#MULTIPLICATIVE SEASONAL ARIMA by default #MAPE 0.074 
S.ARIMA <- auto.arima(electricity_ts, method="ML", seasonal = TRUE)
summary(S.ARIMA) 
checkresiduals(S.ARIMA) #Figure 1 in report

# Series: electricity_ts 
# ARIMA(0,0,2)(0,1,2)[24] with drift 

#Ljung Box Test
Box.test(S.ARIMA$residuals, lag = 24, type = "Ljung")

# data:  S.ARIMA$residuals
# X-squared = 92739, df = 24, p-value < 2.2e-16

#Figure 2 in Report
autoplot(forecast::forecast(S.ARIMA, h = 7*24)) +
  autolayer(fitted(S.ARIMA), series="Fitted") + ylab("Energy Use (MegaWatts)") +
  geom_vline(color = 'orange', linetype = 'dashed', xintercept = 1157) + xlim(900, 1170) + ggtitle('Forecast of Energy Use with Seasonal ARIMA (0,0,2) (0,1,2) 24')

# ACF and PACF plot 
electricity_ts %>%
  Arima(order=c(0,0,2), seasonal=c(0,1,2)) %>%
  residuals() %>% ggtsdisplay()

g <- forecast::forecast(S.ARIMA, h = 7*24)$mean
S.ARIMA.error <- test2$mw - g
S.ARIMA.MAE <- mean(abs(S.ARIMA.error))
S.ARIMA.MAPE <- mean(abs(S.ARIMA.error)/abs(test2$mw))
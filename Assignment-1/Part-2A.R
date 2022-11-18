if(!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("yukai-yang/FE")
install.packages("kableExtra")
library(FE)
ls( grep("FE", search()) )
library(tidyverse)
library(psych)
library(kableExtra)



################################################################
##### Part 2: Asset Return Predictability and Market Efficiency
################################################################
## A. Testing for asset return predictability

############## Task 1. #################

logreturn_DJ_d = DJ_d$r_Dow_Jones[!is.na(DJ_d$r_Dow_Jones)]
logreturn_DJ_w = DJ_w$r_close[!is.na(DJ_w$r_close)]
describe(DJ_d$r_Dow_Jones)
describe(DJ_w$r_close)


# calculate ACF and PACF

# ACF and PACF values of daily data
acf(logreturn_DJ_d, main = "ACF of daily Dow Jones returns")
pacf(logreturn_DJ_d, main="PACF of daily Dow Jones returns")


# ACF and PACF values of weekly data
acf(logreturn_DJ_w, main = "ACF of Weekly Dow Jones returns")
pacf(logreturn_DJ_w, main="PACF of weekly Dow Jones returns")

### Ljung-Box test function
LB <- function(vx,lag,ip){
  tmp = acf(vx,lag.max=lag,plot=F)$acf
  tmp = tmp[2:(lag+1)]**2
  test = sum(tmp/(length(vx)-1:lag))*length(vx)*(length(vx)+2)
  return(list(test=test, pval=1-pchisq(test,df=lag-ip)))
}

### Ljung-Box test for daily data
LB(vx=logreturn_DJ_d, lag=10, ip=0)
LB(vx=logreturn_DJ_d, lag=50, ip=0)
LB(vx=logreturn_DJ_d, lag=100, ip=0)

### Ljung-Box test for weekly data
LB(vx=logreturn_DJ_w, lag=10, ip=0)
LB(vx=logreturn_DJ_w, lag=50, ip=0)
LB(vx=logreturn_DJ_w, lag=100, ip=0)

## Variance ratio test function
VDR = function(vr, iq){
  iTT = length(vr)
  im = floor(iTT/iq)
  iT = im*iq

  rr = vr[1:iT]
  mu = mean(rr)
  sa2 = var(rr)

  arr = NULL
  for(iter in 1:(iT-iq+1))
    arr = c(arr,sum(rr[iter:(iter+iq-1)]))

  sc2 = sum((arr-mu*iq)**2)/iq/(iT-iq+1)/(1-(1/im))

  VD = sc2 - sa2
  VR = sc2/sa2
  tmp = sqrt(2*(2*iq-1)*(iq-1)/3/iq)

  VD = VD*sqrt(iT)/sa2/tmp
  VR = (VR-1)*sqrt(iT)/tmp

  return(list(VD=VD, VR=VR))
}

# Variance ratio test for daily data
VDR(vr = logreturn_DJ_d, iq = 5)

# Variance ratio test for weekly data
VDR(vr = logreturn_DJ_w, iq = 5)



################# task 2 ################
###### Two-day and two-week log returns

##### generating two-day and two-week returns
DJ_d_even = DJ_d$r_Dow_Jones[seq(2, nrow(DJ_d), 2)]
DJ_d_odd = DJ_d$r_Dow_Jones[seq(1, nrow(DJ_d), 2)]
DJ_w_even = DJ_w$r_close[seq(2, nrow(DJ_w), 2)]
DJ_w_odd = DJ_w$r_close[seq(1, nrow(DJ_w), 2)]

# calculate ACF and PACF

# ACF and PACF values of daily data
acf(DJ_d_even, main = "ACF of two-day Dow Jones returns")
acf(DJ_d_odd, main = "ACF of two-day Dow Jones returns")
pacf(DJ_d_even, main="PACF of two-day Dow Jones returns")
pacf(DJ_d_odd, main="PACF of two-day Dow Jones returns")


# ACF and PACF values of weekly data
acf(DJ_w_even, main = "ACF of two-week Dow Jones returns")
acf(DJ_w_odd, main = "ACF of two-week Dow Jones returns")
pacf(DJ_w_even, main="PACF of two-week Dow Jones returns")
pacf(DJ_w_odd, main="PACF of two-week Dow Jones returns")



### Ljung-Box test for two-day data
LB(vx=DJ_d_even, lag=10, ip=0)
LB(vx=DJ_d_even, lag=20, ip=0)
LB(vx=DJ_d_odd, lag=10, ip=0)
LB(vx=DJ_d_odd, lag=20, ip=0)

### Ljung-Box test for two-week data
LB(vx=DJ_w_even, lag=10, ip=0)
LB(vx=DJ_w_even, lag=50, ip=0)
LB(vx=DJ_w_odd, lag=10, ip=0)
LB(vx=DJ_w_odd, lag=50, ip=0)

# Variance ratio test for daily data
VDR(vr = DJ_d_even, iq = 5)
VDR(vr = DJ_d_odd, iq = 5)

# Variance ratio test for weekly data
VDR(vr = DJ_w_even, iq = 5)
VDR(vr = DJ_w_odd, iq = 5)



######## Higher aggregated log returns ##########
##### Generating Higher aggregated log returns data
DJ_d5 = DJ_d$r_Dow_Jones[seq(1, nrow(DJ_d), 5)]
DJ_d10 = DJ_d$r_Dow_Jones[seq(1, nrow(DJ_d), 10)]
DJ_w4 = DJ_w$r_close[seq(1, nrow(DJ_w), 4)]
DJ_w12 = DJ_w$r_close[seq(1, nrow(DJ_w), 12)]


# calculate ACF and PACF

# ACF and PACF values of 5-day and 10-day data
acf(DJ_d5, main = "ACF of 5-day Dow Jones returns")
acf(DJ_d10, main = "ACF of 10-day Dow Jones returns")
pacf(DJ_d5, main="PACF of 5-day Dow Jones returns")
pacf(DJ_d10, main="PACF of 10-day Dow Jones returns")


# ACF and PACF values of 4 week and 12 week data
acf(DJ_w4, main = "ACF of 4-week Dow Jones returns")
acf(DJ_w12, main = "ACF of 12-week Dow Jones returns")
pacf(DJ_w4, main="PACF of 4-week Dow Jones returns")
pacf(DJ_w12, main="PACF of 12-week Dow Jones returns")


### Ljung-Box test for 5-day and 10-day data
LB(vx=DJ_d5, lag=10, ip=0)
LB(vx=DJ_d5, lag=20, ip=0)
LB(vx=DJ_d10, lag=10, ip=0)
LB(vx=DJ_d10, lag=20, ip=0)

### Ljung-Box test for two-week data
LB(vx=DJ_w4, lag=10, ip=0)
LB(vx=DJ_w4, lag=50, ip=0)
LB(vx=DJ_w12, lag=10, ip=0)
LB(vx=DJ_w12, lag=50, ip=0)

# Variance ratio test for daily data
VDR(vr = DJ_d5, iq = 5)
VDR(vr = DJ_d10, iq = 5)

# Variance ratio test for weekly data
VDR(vr = DJ_w4, iq = 5)
VDR(vr = DJ_w12, iq = 5)

?pnorm

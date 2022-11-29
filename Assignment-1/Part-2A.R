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
LB = function(data, lag, p){

  # time periods
  T = length(data)

  # get ACF coeffs
  rho = acf(data, lag.max=lag, plot=F)$acf

  # remove 0th lag coeffs then square
  rho_sq = rho[2:(lag+1)]**2

  # test statistic
  Q = T * (T+2) * sum(rho_sq / (T-1:lag))

  # return test and p-value
  return(list(test=Q, pval=1-pchisq(Q, df=lag-p)))
}

### Ljung-Box test for daily data
LB(logreturn_DJ_d, lag=10, p=0)
LB(logreturn_DJ_d, lag=50, p=0)
LB(logreturn_DJ_d, lag=100, p=0)

### Ljung-Box test for weekly data
LB(vx=logreturn_DJ_w, lag=10, p=0)
LB(vx=logreturn_DJ_w, lag=50, p=0)
LB(vx=logreturn_DJ_w, lag=100, p=0)

## Variance ratio test function

# key slide lec 3 pg 21
VDR = function(vr, iq){
  iTT = length(vr)
  im = floor(iTT/iq)
  iT = im*iq

  rr = vr[1:iT]
  mu = mean(rr)
  sa2 = var(rr)

  #
  arr = NULL
  for(iter in 1:(iT-iq+1))
    arr = c(arr,sum(rr[iter:(iter+iq-1)]))

  sc2 = sum((arr-mu*iq)**2)/iq/(iT-iq+1)/(1-(1/im))

  VD = sc2 - sa2
  VR = sc2/sa2
  tmp = sqrt(2*(2*iq-1)*(iq-1)/3/iq)

  VD = VD*sqrt(iT)/sa2/tmp
  VR = (VR-1)*sqrt(iT)/tmp

  # returns test stats ~ N(0, 1)
  return(list(VD=VD, VR=VR))
}

# Variance ratio test for daily data
VDR(vr = logreturn_DJ_d, iq = 5)

# Variance ratio test for weekly data
VDR(vr = logreturn_DJ_w, iq = 5)



################# task 2 ################
###### Two-day and two-week log returns

##### generating two-day and two-week returns

### higher aggregate data generation with overlap
data_aggr_overlap = function(data, aggr) {
  iTT = length(data)
  im = floor(iTT/aggr)
  iT = im*aggr
  rr = data[1:iT]

  arr = NULL
  for(iter in 1:(iT-aggr+1))
    arr = c(arr,sum(rr[iter:(iter+aggr-1)]))

  return(arr)
}

### higher aggregate data generation with out overlap
data_aggr = function(data, aggr) {
  iTT = length(data)
  im = floor(iTT/aggr)
  iT = im*aggr
  rr = data[1:iT]

  arr = NULL
  for(iter in seq(1, iT, aggr))
    arr = c(arr, sum(rr[iter:(iter+aggr-1)]))

  return(arr)
}

DJ_d_two = data_aggr(DJ_d$r_Dow_Jones, 2)
DJ_w_two = data_aggr(DJ_w$r_close, 2)

# calculate ACF and PACF

# ACF and PACF values of daily data
acf(DJ_d_two, main = "ACF of two-day Dow Jones returns")
pacf(DJ_d_two, main="PACF of two-day Dow Jones returns")


# ACF and PACF values of weekly data
acf(DJ_w_two, main = "ACF of two-week Dow Jones returns")
pacf(DJ_w_two, main="PACF of two-week Dow Jones returns")



### Ljung-Box test for two-day data
LB(DJ_d_two, lag=10, p=0)
LB(DJ_d_two, lag=20, p=0)

### Ljung-Box test for two-week data
LB(DJ_w_two, lag=10, p=0)
LB(DJ_w_two, lag=50, p=0)

# Variance ratio test for daily data
VDR(vr = DJ_d_two, iq = 5)

# Variance ratio test for weekly data
VDR(vr = DJ_w_two, iq = 5)



######## Higher aggregated log returns ##########
##### Generating Higher aggregated log returns data
DJ_d5 = data_aggr(DJ_d$r_Dow_Jones, 5)
DJ_d10 = data_aggr(DJ_d$r_Dow_Jones, 10)
DJ_w4 = data_aggr(DJ_w$r_close, 4)
DJ_w12 = data_aggr(DJ_w$r_close, 12)


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
LB(DJ_d5, lag=10, p=0)
LB(DJ_d5, lag=20, p=0)
LB(DJ_d10, lag=10, p=0)
LB(DJ_d10, lag=20, p=0)

### Ljung-Box test for two-week data
LB(DJ_w4, lag=10, p=0)
LB(DJ_w4, lag=50, p=0)
LB(DJ_w12, lag=10, p=0)
LB(DJ_w12, lag=50, p=0)

# Variance ratio test for daily data
VDR(vr = DJ_d5, iq = 5)
VDR(vr = DJ_d10, iq = 5)

# Variance ratio test for weekly data
VDR(vr = DJ_w4, iq = 5)
VDR(vr = DJ_w12, iq = 5)



############## Task 3 ################
DJ_d$Year = as.character(DJ_d$Date)
DJ_d$Year = as.integer( str_sub(DJ_d$Year,-2,-1))
DJ_d1 = (DJ_d[DJ_d$Year>= 15 & DJ_d$Year <= 61,])
DJ_d2 = (DJ_d[DJ_d$Year>= 62 & DJ_d$Year <= 90,])

### Line Plot
ggplot()+
  geom_line(mapping = aes( x=1:length(DJ_d1$r_Dow_Jones),
                           y = DJ_d1$r_Dow_Jones), size=.5) +
  labs(y='log returns', x='time horizon') +
  ggtitle("from 1915 to 1955")

par(mfrow)


### acf plot
acf(DJ_d1$r_Dow_Jones, main = "ACF of Dow Jones returns from 1915 to 1955")

LB(DJ_d1$r_Dow_Jones, lag=10, p=0)
VDR(vr = DJ_d1$r_Dow_Jones, iq = 5)

## PS 1 Part 1:A

library(tidyverse)
library(FE)
library(psych)
library(knitr)

View(DJ_d)
View(DJ_w)

daily <- DJ_d$r_Dow_Jones
weekly <- DJ_w$r_close

## No. 1
plot1A1.1 <- ggplot(data=DJ_d) + geom_line(aes(x=1:nrow(DJ_d), y=r_Dow_Jones)) +
  labs(x='Time horizon', y='Log return', title='Dow Jones daily log returns')
plot1A1 + theme_classic()

plot1A1.2 <- ggplot(data=DJ_w) + geom_line(aes(x=1:nrow(DJ_w), y=r_close)) +
  labs(x='Time horizon', y='Log return', title='Dow Jones weekly log returns')
plot1A1.2 + theme_classic()

# library(xtable)
# distr_daily <- describe(DJ_d$r_Dow_Jones)
# distr_weekly <- describe(DJ_w$r_close)
# xtable(distr_daily)

    #``{r, results='asis', message=FALSE, warning=FALSE}
    #library(xtable)

    #print(xtable(dat, caption="My table using xtable package", label="xtabletab"), comment=FALSE,
    #      caption.placement="top", type="html")
    #```


distr_daily <- data.frame(describe(DJ_d$r_Dow_Jones), row.names = "Daily")
distr_weekly <- data.frame(describe(DJ_w$r_close), row.names = "Weekly")

dat.df <- t(rbind.data.frame(distr_daily, distr_weekly))
kable(dat.df)


## No. 2
# a)
std_daily = (daily - mean(daily))/sd(daily) # standardize the daily return
qqnorm(std_daily)
qqline(std_daily, col = "red")

std_weekly = (weekly - mean(weekly))/sd(weekly)
qqnorm(std_weekly)
qqline(std_weekly, col = "red")

# b)
df = 4 # the lower the df the bigger the kurtosis
std_daily_t = std_daily*sqrt(df/(df-2))
qqplot(rt(length(std_daily_t), df = df), std_daily_t)
qqline(std_daily_t, col = "red")

df = 4
std_weekly_t = std_weekly*sqrt(df/(df-2))
qqplot(rt(length(std_weekly_t), df = df), std_weekly_t)
qqline(std_weekly_t, col = "red")

## No. 3
# a) Normal distribution
#daily
b = 30
bin = 1:b/b
r_nor = pnorm(std_daily)
hist(r_nor) # I have smaller bins than Kai; this would have to be uniformly distributed if 0 was true
# (that distribution of daily log returns is normal distribution)

vn = NULL
for(val in bin) vn = c(vn,sum(r_nor <= val)) #try to compute the frequency for the bins
vn = c(vn[1],diff(vn)) # put 1st value because differencing (why?) loses the first observation; vector of probabilities?
chi_test = sum((vn-length(std_daily)/b)**2/(length(std_daily)/b))
cat("test =",chi_test," df =",b-3," p-value =",1-pchisq(chi_test,df=b-3)) # p-value = 0, reject norm. distr.
#diff. output than his (df, test etc)

#weekly
r_nor_w = pnorm(std_weekly)
hist(r_nor_w)

vn_w = NULL
for(val in bin) vn_w = c(vn_w,sum(r_nor_w <= val)) #try to compute the frequency for the bins
vn_w = c(vn_w[1],diff(vn_w)) # put 1st value because differencing (why?) loses the first observation; vector of probabilities?
chi_test_w = sum((vn_w-length(std_daily)/b)**2/(length(std_weekly)/b))
cat("test =",chi_test_w," df =",b-3," p-value =",1-pchisq(chi_test_w,df=b-3))

# b) t distribution
df = 4

#daily
r_t = pt(std_daily_t,df=df)
hist(r_t)
vt = NULL; for(val in bin) vt = c(vt,sum(vq <= val))
vt = c(vt[1],diff(vt))
test = sum((vt-length(std_daily_t)/b)**2/(length(std_daily_t)/b)) # he used the std. normal data, not for t-dist.
cat("test =",test," df =",b-3," p-value =",1-pchisq(test,df=b-3)) # p-value = 0, reject t distr.
#diff. output than his (df, test etc)

#weekly
r_t_w = pt(std_weekly_t,df=df)
hist(r_t_w)
vt_w = NULL; for(val in bin) vt = c(vt_w,sum(vq <= val))
vt_w = c(vt[1],diff(vt_w))
test = sum((vt-length(std_daily_t)/b)**2/(length(std_daily_t)/b))
cat("test =",test," df =",b-3," p-value =",1-pchisq(test,df=b-3))

# c) Mixture of normal distributions

#arbiotrary choice
alpha = 0.5
sigma = 3.7 # this is the var of the second normal distr. and the first one is std_daily?

#daily
mix_daily = std_daily*sqrt((1-alpha) + alpha*sigma**2) #as data is std. and needs to be multiplied with var of mixed distr.
r_mix = (1-alpha)*pnorm(mix_daily) + alpha*pnorm(mix_daily,sd=sigma) #corresp. CDF of the mixture normal distr.
hist(r_mix)
vn_mix = NULL; for(val in bin) vn_mix = c(vn_mix,sum(r_mix <= val))
vn_mix = c(vn_mix[1],diff(vn_mix))
test = sum((vn_mix-length(std_daily)/b)**2/(length(std_daily)/b))
cat("test =",test," df =",b-3," p-value =",1-pchisq(test,df=b-3))


#optimization daily
func <- function(vx){
  alpha = vx[1]
  sigma = vx[2]
  mix_daily_opt = std_daily*sqrt((1-alpha) + alpha*sigma**2)

  r_mix_opt = (1-alpha)*pnorm(mix_daily_opt) + alpha*pnorm(mix_daily_opt,sd=sigma)
  vn_mix_opt = NULL; for(val in bin) vn_mix_opt = c(vn_mix_opt,sum(r_mix_opt <= val))
  vn_mix_opt = c(vn_mix_opt[1],diff(vn_mix_opt))
  return(sum((vn_mix_opt-length(std_daily)/b)**2/(length(std_daily)/b)))
}
#vx just a pair of parameters

##### what's happening here??
func(c(0.15,4))

optim(par=c(0.1,4),fn=func,method="BFGS") # this is not cherry-picking, this is optimizing (similar to OLS)

#alpha = 0.196 irgendwas, sigma = 2.66 irgendwas --> p-value becomes bigger than 1%, best we can get
####

#optimization weekly
func <- function(vx){
  alpha = vx[1]
  sigma = vx[2]
  mix_weekly_opt = std_weekly*sqrt((1-alpha) + alpha*sigma**2)

  r_mix_opt_w = (1-alpha)*pnorm(mix_weekly_opt) + alpha*pnorm(mix_weekly_opt,sd=sigma)
  vn_mix_opt_w = NULL; for(val in bin) vn_mix_opt_w = c(vn_mix_opt_w,sum(r_mix_opt_w <= val))
  vn_mix_opt_w = c(vn_mix_opt_w[1],diff(vn_mix_opt_w))
  return(sum((vn_mix_opt_w-length(std_weekly)/b)**2/(length(std_weekly)/b)))
}

#### adjust
func(c(0.15,4))

optim(par=c(0.1,4),fn=func,method="BFGS")
####


# Part-1B. ####

library(tidyverse)
library(FE)
library(xtable)

# 1. ####

# generate log returns
# log(P_t) - log(P_{t-1}) = log(P_t / P_{t-1}) = log(R_t)
lret = apply(log(index_d), 2, diff)

# summarize data
descr <- describe(lret)
descr <- descr[, c("n", "mean", "sd", "median", "min", "max", "skew", "kurtosis")]

xtable(descr)

# calculate ACF and PACF
# interpretation of ACF/PACF
# https://medium.com/@ooemma83/how-to-interpret-acf-and-pacf-plots-for-identifying-ar-ma-arma-or-arima-models-498717e815b6

acf_tabs <- function(data, lag.max=10, plot=T){

  # create tables for ACF and PACF values
  acf_tab <- matrix(nrow=ncol(data), ncol=lag.max,
                    dimnames=list(index=colnames(data),
                                  lag=1:lag.max
                    )
  )
  pacf_tab <- matrix(nrow=ncol(data), ncol=lag.max,
                     dimnames=list(index=colnames(data),
                                   lag=1:lag.max)
  )

  for (c in colnames(data)){

    # data for each column removing na
    x <- data[!is.na(data[,c]), c]


    if (plot){
      # two columns for ACF and PACF
      par(mfcol=c(1, 2))
    }

    # plot ACf and PACF
    acf_val <- acf(x, main=c, plot=plot)
    pacf_val <- pacf(x, main=c, plot=plot)

    # save coefficients
    acf_tab[c, ] <- acf_val$acf[2:(lag.max+1)]
    pacf_tab[c, ] <- pacf_val$acf[1:lag.max]

  }

  return(list(acf=acf_tab, pacf=pacf_tab))

}

tabs = acf_tabs(lret)

xtable(t(tabs$acf))

xtable(t(tabs$pacf))

# 2. ####

# Ljung-Box test function
# See lec 2 slide 38
LB <- function(data, lag, p){

  # time periods
  T <- length(data)

  # get ACF coeffs
  rho = acf(data, lag.max=lag, plot=F)$acf

  # remove 0th lag coeffs then square
  rho_sq = rho[2:(lag+1)]**2

  # test statistic
  Q = T * (T+2) * sum(rho_sq / (T-1:lag))

  # return test and p-value
  return(list(test=Q, pval=1-pchisq(Q, df=lag-p)))
}


# lag periods to test
lags <- c(10, 50, 100)

# degrees of freedom df = lag - p
# p is number of parameters estimated
p <- 0

# table of LB test statistics
LB_tests <- matrix(nrow=length(lags), ncol=ncol(lret),
                  dimnames=list(lag=lags, index=colnames(lret))
                  )

# table of Chi-square test p-values
LB_pvals <- matrix(nrow=length(lags), ncol=ncol(lret),
                   dimnames=list(lag=lags, index=colnames(lret))
                   )

# run Ljung-Box tests
for (c in colnames(lret)){

  for (i in 1:length(lags)){

    # remove na from data
    x <- lret[!is.na(lret[,c]), c]
    LB_res <- LB(x, lag=lags[i], p)

    # significance level stars
    pval <- signif(LB_res$pval, 3)
    for (qt in c(0.1, 0.05, 0.01)){
      if (LB_res$pval < qt){
        pval <- paste(pval,'*', sep='')
      }
    }

    # save test stat and pvals
    LB_tests[i, c] <- LB_res$test
    LB_pvals[i, c] <- pval

  }
}

comment <- paste('This table presents p-values of',
                 'Ljung-Box tests of no autocorrelation.',
                 'Stars indicate significance levels.',
                 '*: p-value < 0.1,',
                 '**: p-value < 0.5,',
                 '***: p-value <0.01.')


print(xtable(LB_pvals, caption=comment), type='html')

# 3. ####

corr_mat <- function(data, lag=1) {

  # column names
  cols <- colnames(data)

  # correlation matrix
  corr_mat <- matrix(nrow = ncol(data), ncol = ncol(data),
                     dimnames=list(lags=cols, cols))

  for (c1 in cols){
    for (c0 in cols){
      # periods with either are NA
      isna <- (is.na(data[, c1]) | is.na(data[, c0]))

      # non-na returns of pair
      pair_data <- data[!isna, c(c1, c0)]
      T <- nrow(pair_data)

      # add period lag autocorrelation to matrix
      corr_mat[c1, c0] <- cor(pair_data[1:(T-lag),c1], pair_data[(1+lag):T,c0])

      # acf(pair_data)
      # pacf(pair_data)
    }
  }

  # re-label past/current period returns
  rownames(corr_mat) <- paste(cols, '(t-1)', sep='')
  colnames(corr_mat) <- paste(cols, '(t)', sep='')

  return(corr_mat)

}

corr_mat(lret)

# 4.
sqlret = lret^2


acf_tab(sqlret)

#############################################################
## Financial Econometrics 2ST119
## Problem Set 2, solutions
## author: Yukai Yang
## Department of Statistics, Uppsala University
#############################################################

library(FE)

#############################################################
## Part 1: Capital Asset Pricing Model
#############################################################

## A. Size and book-to-market effects

portfolio_m
data("portfolio_m")
?factors_m


# use function Fara_MacBeth
Fama_MacBeth()
EstCAPM()

#############################################################
## Some useful functions
#############################################################

# Fama and Macbeth two-step regression for excess asset returns
# input:
#  mZ, T by N excess returns
#  vZm, vector of excess market returns
# output:
#  alpha
#  beta, estimated betas for the asset returns
#  J1, Wald, mLR, three joint tests and their p-values for alpha=0
#  mum, sample mean of the market excess return
#  gamma0
#  mrprem, market risk premia
#  g0sd, prsd, standard errors of gamma0 and market risk premia
FamaMacbeth0 <- function(mZ, vZm)
{
  iT = nrow(mZ)
  iN = ncol(mZ)

  # First Pass
  mX = cbind(1,vZm)
  pars = chol2inv(chol(crossprod(mX)))%*%crossprod(mX,mZ)
  alpha = c(pars[1,])
  # estimated betas
  beta = c(pars[2,])
  # sample mean of each excess return
  mu = apply(mZ,2,mean)
  # sample mean of the market excess return
  mum = mean(vZm)
  # variance of the market excess return
  sigm2 = c(crossprod(vZm-mum))/iT
  # residules
  me = mZ-mX%*%pars
  # covariance matrix of the residules
  Sigma = crossprod(me)/iT

  # do the test on pp.24 & 25 in slides for lecture 4
  tmp = c(t(alpha)%*%chol2inv(chol(Sigma))%*%alpha)/(1+mum^2/sigm2)

  # Gibbons/Ross/Shanken (1989)
  J1 = tmp*(iT-iN-1)/iN
  jpval = 1-pf(J1,df1=iN,df2=iT-iN-1)

  # Wald
  W = tmp*iT
  wpval = 1-pchisq(W,df=iN)

  # modified LR, Jobson/Korkie (1982)
  mX = matrix(vZm,iT,1)
  pars = chol2inv(chol(crossprod(mX)))%*%crossprod(mX,mZ)
  me = mZ-mX%*%pars
  # covariance matrix of the residules
  SigmaR = crossprod(me)/iT
  LR = (sum(log(eigen(SigmaR)$values))-sum(log(eigen(Sigma)$values)))*(iT-iN/2-2)
  lrpval = 1-pchisq(LR,df=iN)


  # Second Pass
  mX = cbind(1,beta)
  xxinv = chol2inv(chol(crossprod(mX)))
  xxinvx = tcrossprod(xxinv,mX)
  gamma = NULL; gamsd = NULL
  for(iter in 1:iT){
    vy = mZ[iter,]
    tmp = xxinvx%*%vy
    # gamma
    gamma = rbind(gamma, c(tmp))
    s2 = c(crossprod(vy - mX%*%tmp))/iN
    # gammas' standard error
    gamsd = rbind(gamsd, sqrt(diag(s2 * xxinv)))
  }

  # gamma0 and market risk premia
  gamma0 = c(gamma[,1]); g0sd = c(gamsd[,1])
  mrprem = c(gamma[,2]); prsd = c(gamsd[,2])

  # w tests on pp.30 in slides for lecture 4
  # should be compared with student t with T-1 degrees of freedom
  # or with standard normal if T is large
  wgamma0 = mean(gamma0)/sqrt(sum((gamma0-mean(gamma0))**2)/iT/(iT-1))
  wgamma1 = mean(mrprem)/sqrt(sum((mrprem-mean(mrprem))**2)/iT/(iT-1))

  return(list(alpha=alpha,beta=beta,J1=c(J1,jpval),W=c(W,wpval),mLR=c(LR,lrpval),
              gamma0=gamma0,g0sd=g0sd,mrprem=mrprem,prsd=prsd,wgamma0=wgamma0,wgamma1=wgamma1))
}

subsample = 501:700
mR = as.matrix(portfolio_m[subsample,25:124])
Rf = as.matrix(portfolio_m[subsample,'Tbill'])
Rm = as.matrix(portfolio_m[subsample,'Market'])
mZ = sweep(mR,1,Rf)
mR - c(Rf)

vZm = Rm - Rf

tmp = EstCAPM(mZ, vZm)
tmp$beta

Fama_MacBeth(mZ, vZm)
ret = FamaMacbeth0(mZ=mZ, vZm=vZm)
ret

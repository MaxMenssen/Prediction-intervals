#-------------------------------------------------------------------------------
#------- Computation of the six different prediction intervals -----------------
#-------------------------------------------------------------------------------


# Load the packages that will be used for calculations
library(tidyverse)
library(plyr)
library(gamlss)


#-------------------------------------------------------------------------------

# Source the R-files for each interval from github

source("https://raw.githubusercontent.com/MaxMenssen/Prediction-intervals/master/nelson.R")
source("https://raw.githubusercontent.com/MaxMenssen/Prediction-intervals/master/nelsonphi.R")
source("https://raw.githubusercontent.com/MaxMenssen/Prediction-intervals/master/nelsonphi1.R")
source("https://raw.githubusercontent.com/MaxMenssen/Prediction-intervals/master/nelsonphi1_bisec.R")
source("https://raw.githubusercontent.com/MaxMenssen/Prediction-intervals/master/qBB.R")
source("https://raw.githubusercontent.com/MaxMenssen/Prediction-intervals/master/qBB_bisec.R")


#-------------------------------------------------------------------------------


# Example data from the paper
Dead <- c(15, 10, 12, 12, 13, 11, 19, 11, 14, 21)
Surviving <- c(35, 40, 38, 38, 37, 39, 31, 39, 36, 29)

dat <- data.frame(Dead, Surviving)


#-------------------------------------------------------------------------------

### Arguments for the interval functions
# succ: Vector of historical binomial observations (number of success, xk in the paper)
# fail: Number of failures (nk-xk)
# m: Future sample size
# alpha: Type one error rate
# nboot: Number of bootstraps

# Nelson interval eq(9)
nelson(succ=dat$Dead, fail=dat$Surviving, m=50, alpha=0.05)

# Nelsonphi interval eq(12)
nelsonphi(succ=dat$Dead, fail=dat$Surviving, m=50, alpha=0.05)

# Nelsonphi1 interval eq(13)
nelsonphi1(succ=dat$Dead, fail=dat$Surviving, m=50, alpha=0.05)

# Nelsonphi1_bisec interval eq(27)
nelson_bisec(succ=dat$Dead, fail=dat$Surviving, m=50, alpha=0.05, nboot=10000)

# qBB läuft
qBB_noboot(succ=dat$Dead, fail=dat$Surviving,m=50, alpha=0.05)

# qBB_bisec lauft
qBB_bisec(succ=dat$Dead, fail=dat$Surviving, m=50, alpha=0.05, nboot=10000)






#-------------------------------------------------------------------------------
#------------------------ Simple qBB Interval ----------------------------------
#-------------------------------------------------------------------------------


phiBBLui <- function(data){
        
        
        # Adjusting the data if all nk=xk or all xk=0
        if(all(data$succ==0))
        {
                data$succ[1] <- 0.5
                data$fail[1] <- data$fail[1]-0.5
        }
        #
        
        if(all(data$fail==0))
        {
                data$fail[1] <- 0.5
                data$succ[1] <- data$succ[1]-0.5
        }
        # 
        
        
        xk <- data$succ
        nk <- rowSums(data)
        k <- nrow(data)
        pi <- sum(xk)/sum(nk)
        
        
        ### Intra class correlation nach Lui et al 2000
        
        # Between mean squared error
        part1B <- sum(xk^2/nk)
        part2B <- (sum(xk)^2/sum(nk))
        part3B <- k-1
        
        BMS <- (part1B-part2B)/part3B
        
        
        # Within mean squared error
        part1W <- sum(xk)
        part2W <- sum(xk^2/nk)
        part3W <- sum(nk-1)
        
        WMS <- (part1W-part2W)/part3W
        
        
        # mstar
        part1m <- sum(nk)^2
        part2m <- sum(nk^2)
        part3m <- (k-1)*sum(nk)
        
        mstar <- (part1m-part2m)/part3m
        
        
        # Estimated intrasclass correlation
        rhohatLui <- (BMS-WMS)/(BMS+(mstar-1)*WMS)
        
        
        # phiBB
        phiBBLui <- 1+(nk-1)*rhohatLui
        phiBBLui <- max(1.001, phiBBLui[1])
        
        
        # estimated sigma
        sigmacalc <- (phiBBLui-1)/(nk[1]-phiBBLui)
        
        
        return(c(pi=pi, disp=phiBBLui, sigma=sigmacalc))
}
#

#-------------------------------------------------------------------------------



qBB_noboot <- function(succ, fail, m, alpha=0.05){
        
        dat <- data.frame(succ, fail)
        
        # Phi and Pi estimates
        pd <- phiBBLui(data=dat)
        
        # Mu and Sigma
        mu <- as.numeric(pd[1])
        sigma <- as.numeric(pd[3])
        
        # Interval
        l <- qBB(alpha/2, mu=mu, sigma=sigma, bd=m)
        u <- qBB(1-alpha/2, mu=mu, sigma=sigma, bd=m)
        
        
        # Interval as vector
        intcalib <- c(lower=l, upper=u)
        
        return(intcalib)
}






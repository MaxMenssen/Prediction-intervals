#-------------------------------------------------------------------------------
#-------------------------- Alpha Calib qBB ------------------------------------
#-------------------------------------------------------------------------------


# Data generating function 
betbin <- function(anzhist, nhist, exprop, phi){
        
        # eq(5)
        asum <- (phi-nhist)/(1-phi) 
        
        # eq(6)
        a <- exprop*asum 
        b <- asum-a
        
        # Sampling of overdispersed bin. data according to eq(1) 
        pis <- rbeta(n=anzhist, shape1=a, shape2=b) 
        y <- rbinom(n=anzhist, size=nhist, prob=pis) 
        x <- nhist-y 
        
        dat <- data.frame(succ=y, fail=x)
        return(dat)
}
#

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
        
        
        # mu (f?r rBB, qBB)
        sigmacalc <- (phiBBLui-1)/(nk[1]-phiBBLui)
        
        
        return(c(pi=pi, disp=phiBBLui, sigma=sigmacalc))
}
#

#-------------------------------------------------------------------------------

qBB_int <- function(dat, m=rowSums(dat)[1], alpha=0.05)
{
        
        # Phi and Pi estimates
        pd <- phiBBLui(data=dat)
        
        # Mu and Sigma
        mu <- as.numeric(pd[1])
        sigma <- as.numeric(pd[3])
        
        # Interval calculation
        l <- qBB(alpha/2, mu=mu, sigma=sigma, bd=m)
        u <- qBB(1-alpha/2, mu=mu, sigma=sigma, bd=m)
        
        
        # Interval as vector
        intcalib <- data.frame(lower=l, upper=u, alpha=alpha)
        
        
        
        return(intcalib)
}


#-----------------------------------------------------------------------

### Bootstrap function

histdat <- function(data, m=rowSums(data)[1], k=nrow(data), pi1, disp1, nboot=1000){
        
        # BS historical data
        datastar <- replicate(n=nboot, simplify=FALSE, 
                              expr=betbin(anzhist=k, nhist=m, exprop=pi1,
                                          phi=disp1))
        
        # BS future observation
        futvalue <- unlist(replicate(n=nboot, simplify=FALSE, 
                                     expr=betbin(anzhist=1, nhist=m, exprop=pi1,
                                                 phi=disp1)[1]))
        
        hd <- list(datastar, futvalue)
        
        return(hd)
        
}



#-------------------------------------------------------------------------------

### Bisection function

bisection <- function(f, a1=0.00001, b1=0.3, steps = 15, tolerance = 0.003, 
                      alpha=0.05, traceplot=TRUE) {
        
        # Nominal coverage probability
        nomcov <- 1-alpha
        
        
        # Startpoints
        minval <- f(lambda=a1, alpha=alpha)
        maxval <- f(lambda=b1, alpha=alpha)
        
        # Both coverages are higher than 0.95: Return b1 as lambda
        if ((minval > 0) && (maxval > 0)) {
                
                val <- c(minval, maxval)
                ab <- c(a1, b1)
                
                # Traceplot
                if(traceplot==TRUE){
                        plot(x=ab, y=val, type="p", pch=20, 
                             xlab="lambda", ylab=paste("coverage-", nomcov), 
                             main="Both higher than 0")
                        lines(x=ab, y=val, type="s", col="red")
                        abline(a=0, b=0, lty="dashed")
                }
                
                
                return(b1)
                
                
                # Both coverages are smaller than 0.95: Return a1 as lambda       
        } else if ((minval < 0) && (maxval < 0)) {
                
                val <- c(minval, maxval)
                ab <- c(a1, b1)
                
                # Traceplot
                if(traceplot==TRUE){
                        plot(x=ab, y=val, type="p", pch=20, 
                             xlab="lambda", ylab=paste("coverage-", nomcov),
                             main="Both smaller than 0")
                        lines(x=ab, y=val, type="s", col="red")
                        abline(a=0, b=0, lty="dashed")
                }
                
                
                return(a1)
        }
        
        
        # Open vectors for the estimated lambda and coverage probabilities
        lam <- vector()
        cover <- vector()
        
        for (i in 1:steps) {
                
                # Calculate midpoint
                c1 <- (a1 + b1) / 2 
                
                # Calculation of coverdiff based on c1
                runval <- f(lambda=c1, alpha=alpha)
                
                # Assigning c1 and runval into the vectors
                lam[i] <- c1
                cover[i] <- runval
                
                # If runval is 0 or
                # If runval is smaller than the tolerance and positve
                # return the value of c1
                if ((runval == 0)  || near(runval, tolerance) || (runval < tolerance) && sign(runval)==1) {
                        
                        # Traceplot
                        if(traceplot==TRUE){
                                plot(x=lam, y=cover, type="p", pch=20, 
                                     xlab="lambda", ylab=paste("coverage-", nomcov), 
                                     main=paste("Trace with", i, "iterations"))
                                lines(x=lam, y=cover, type="s", col="red")
                                abline(a=0, b=0, lty="dashed")
                        }
                        
                        
                        return(c1)
                }
                
                
                
                # If another iteration is required,
                # check the signs of the function at the points c and a and reassign
                # a or b accordingly as the midpoint to be used in the next iteration.
                if(sign(runval)==1){
                        a1 <- c1}
                
                else if(sign(runval)==-1){
                        b1 <- c1}
        }
        
        # If the full n iterations are done:
        
        # Open a data.frame
        lamdf <- data.frame("lambda"=lam, "coverdiff"=cover)
        
        # Traceplot
        if(traceplot==TRUE){
                plot(x=lamdf$lambda, y=lamdf$coverdiff, type="p", pch=20, 
                     xlab="lambda", ylab=paste("coverage-", nomcov), 
                     main=paste("Trace with", i, "iterations"))
                lines(x=lam, y=cover, type="s", col="red")
                abline(a=0, b=0, lty="dashed")
        }
        
        
        # Take the lambda for wich the coveragedifference lies between 0 and -tolerance
        if(all(lamdf$coverdiff < 0) & (any(lamdf$coverdiff >= -tolerance) | 
                                       any(near(lamdf$coverdiff, -tolerance)))){
                
                # Order lamdf by increasing lambda
                ldfo <- lamdf[order(lamdf$lambda, decreasing = FALSE),]
                
                # Take the lamda which is closest to 0
                c1 <- ldfo[which.max(ldfo$coverdiff),][1,1]
                return(c1)
        }
        
        # All coverdiff are smaller than -tolerance, return a1
        else if(all(lamdf$coverdiff < -tolerance)){
                return(a1)
        }
        
        # Take the lambda with the smallest positive coveragedifference
        else{
                ldf <- filter(lamdf, coverdiff > 0) %>% 
                        arrange(desc(lambda))
                c1 <- ldf$lambda[which.min(ldf$coverdiff)]
                return(c1)
        }
}


#-------------------------------------------------------------------------------

### Calibrated interval

qBB_bisec <- function (succ, fail, alpha, m=rowSums(data)[1], k=nrow(data), nboot=10000,
                       a1=0.00001, b1=0.3, steps = 15, tolerance = 0.003, traceplot=FALSE){
        
        data <- data.frame(succ, fail)
        
        # Parameter Estimates
        
        pd <- phiBBLui(data)
        pi1 <- pd[1]
        disp1 <- pd[2]
        
        # Bootstrapped historical data and future observation
        hd <- histdat(data=data, m=m, k=k, pi1=pi1, disp1=disp1, nboot=nboot)
        
        # Coverage Probability (function is needed in the next step)
        covfun <- function(input=hd, alpha, lambda){

                dint <- map(.x=input[[1]], .f=qBB_int, alpha=lambda)  %>%
                        ldply()

                lambdaint <- split(dint, f=factor(dint$alpha))

                #-----------------------------------------------------------------------

                # Is futvalue an element of each interval?
                futinint <- as.data.frame(map(lambdaint, function(x){x$lower<=hd[[2]]&
                                x$upper>=hd[[2]]}))

                # Coverage probabilities for each lambda
                covprob <- unname(apply(futinint, MARGIN=2, mean))

                # Difference between observed coverage and 0.95
                covdiff <- covprob-(1-alpha)

                return(covdiff)
        }
        
        # Generagtion of the calibrated alpha
        alphacalib <- bisection(f=covfun, a1=a1, b1=b1, steps=steps, tolerance=tolerance, 
                                alpha=alpha, traceplot=traceplot)
        
        # Calculation of the interval with alpha calib
        int <- qBB_int(dat=data, m=m, alpha=alphacalib)
        
        return(int)
}





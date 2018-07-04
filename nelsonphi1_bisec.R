#-------------------------------------------------------------------------------
#-------------------------- Alpha Calib nelson ---------------------------------
#-------------------------------------------------------------------------------


betbin <- function(anzhist, nhist, exprop, phi) 
{
        asum <- (phi-nhist)/(1-phi) 
        a <- exprop*asum 
        b <- asum-a
        pis <- rbeta(n=anzhist, shape1=a, shape2=b) 
        y <- rbinom(n=anzhist, size=nhist, prob=pis) 
        x <- nhist-y 
        dat <- data.frame(succ=y, fail=x)
        return(dat)
}


#-------------------------------------------------------------------------------

pidisp<- function(dat){
        
        if(all(dat$succ==0))
        {
                dat$succ[1] <- 0.5
                dat$fail[1] <- dat$fail[1]-0.5
        }
        # 
        
        if(all(dat$fail==0))
        {
                dat$fail[1] <- 0.5
                dat$succ[1] <- dat$succ[1]-0.5
        }
        # 
        
        
        fit <- glm(cbind(succ, fail)~1, dat,
                   family=quasibinomial())
        
        disp <- max(1.001, summary(fit)$dispersion)
        eta <- as.numeric(coef(fit))
        pi <- exp(eta)/(1+exp(eta))
        return(c(disp=disp, pi=pi))
}


#-------------------------------------------------------------------------------

### Parametric Bootstrap function

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

nelsonphi1b <- function(dat, alpha, m=rowSums(dat)[1]){
        
        # Adjusting the data if all nk=xk or all xk=0
        if(all(dat$succ==0))
        {
                dat$succ[1] <- 0.5
                dat$fail[1] <- dat$fail[1]-0.5
        }
        # 
        
        if(all(dat$fail==0))
        {
                dat$fail[1] <- 0.5
                dat$succ[1] <- dat$succ[1]-0.5
        }
        # 
        
        # estimation of phi and pi
        fit <- glm(cbind(succ, fail)~1, dat, 
                   family=quasibinomial(link="identity"))
        
        # estimated phi
        estphi <- max(1, summary(fit)$dispersion)
        
        # estimated pi
        estprop <- coef(fit)[1] 
        names(estprop) <- NULL
        estq <- 1-estprop
        
        # Extimated variance of pi
        varestprop <- as.vector(vcov(fit))
        names(varestprop) <- NULL
        
        
        # standard normal qunatile
        z <- qnorm(alpha/2, 0, 1)
        
        # future observation
        esty <- m*estprop
        
        # Lower and upper boundaries of the interval
        lower <- esty + z*sqrt(estphi*m*estprop*estq + 
                                       (m^2)*varestprop)
        upper <- esty - z*sqrt(estphi*m*estprop*estq +
                                       (m^2)*varestprop)
        
        nelsonphi1b <- c(lower=lower, upper=upper, alpha=alpha)
        
        return(nelsonphi1b)
}
#


#-------------------------------------------------------------------------------

### Bisection

bisection <- function(f, a1=0.00001, b1=0.3, steps = 15, tolerance = 0.003, 
                      alpha=0.05, traceplot=TRUE) {
        
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
        
        # print(cover)
        # 
        # print(lamdf)
        
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

nelson_bisec <- function (succ, fail, alpha=0.05, m, k=nrow(data), nboot=10000,
                          a1=0.00001, b1=0.3, steps = 15, tolerance = 0.003, traceplot=FALSE){
        
        data <- data.frame(succ, fail)
        
        # Parameter Estimates
        
        pd <- pidisp(data)
        pi1 <- pd[2]
        disp1 <- pd[1]
        
        # Bootstrapped historical data and future observation
        hd <- histdat(data=data, m=m, k=k, pi1=pi1, disp1=disp1, nboot=nboot)
        
        # Coverage Probability (function is needed in the next step)
        covfun <- function(input=hd, alpha, lambda){
                
                dint <- map(.x=input[[1]], .f=nelsonphi1b, alpha=lambda)  %>% 
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
        
        # print(hd)
        
        # Generagtion of the calibrated alpha
        alphacalib <- bisection(f=covfun, a1=a1, b1=b1, steps=steps, tolerance=tolerance, 
                                alpha=alpha, traceplot=traceplot)
        
        
        
        # Calculation of the interval with alpha calib
        int <- nelsonphi1b(dat=data, m=m, alpha=alphacalib)
        
        return(int)
}









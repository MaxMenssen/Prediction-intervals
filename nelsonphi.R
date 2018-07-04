#-------------------------------------------------------------------------------
#--------------------------- Nelsonphi -----------------------------------------
#-------------------------------------------------------------------------------



nelsonphi <- function(succ, fail, m, alpha=0.05){
        
        data <- data.frame(succ, fail)
        
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
        
        # Estimation of phi and pi
        fit <- glm(cbind(succ, fail)~1, data, 
                   family=quasibinomial(link="identity"))
        
        
        # Estimated dispersion parameter
        estphi <- summary(fit)$dispersion
        
        # Estimated proportion
        estprop <- coef(fit)[1] 
        names(estprop) <- NULL
        
        # Expected proportion of no success in the future sample
        estq <- 1-estprop
        
        # Variance of the estimated proportion
        varestprop <- estphi*estprop*estq*(1/(sum(rowSums(data))))
        
        # Expected number of success in the future study
        esty <- m*estprop
        
        # Standard normal quantile
        z <- qnorm(alpha/2, 0, 1)
        
        # Lower boundary of the interval
        lower <- esty + z*sqrt(estphi*m*estprop*estq + 
                                       (m^2)*varestprop)
        
        # Upper boundary of the interval
        upper <- esty - z*sqrt(estphi*m*estprop*estq +
                                       (m^2)*varestprop)
        
        
        nelsonphi <- data.frame(lower=lower, upper=upper)
        
        return(nelsonphi)
}



























#-------------------------------------------------------------------------------
#-------- Data and code about figure 1 of Menssen & Schaarschmidt 2018 ---------
#-------------------------------------------------------------------------------


# Please download the rda file "coverages.rda" from github using the link
# https://github.com/MaxMenssen/Prediction-intervals/blob/master/coverages.rda
# into your working directory.

load("coverages.rda")


varcov <- (0.95*0.05)/5000
secov <- sqrt(varcov)

ggplot(coverages, aes(y=value, x=anzhistnhist)) + 
        # geom_boxplot(outlier.shape = NA) +
        # geom_hline is for horizontal lines
        geom_hline(yintercept=0.95)+
        geom_hline(yintercept=(0.95+secov*2), linetype="dashed", color="gray")+
        geom_hline(yintercept=(0.95+secov), linetype="dashed")+
        geom_hline(yintercept=(0.95-secov), linetype="dashed")+
        geom_hline(yintercept=(0.95-secov*2), linetype="dashed",color="gray")+
        geom_point(aes(shape=pi, color=pi),
                   position=position_jitter(width=0.05, height=0)) +
        facet_grid(phif~variable, labeller=label_parsed)+
        theme_bw()+
        # theme(axis.text.x....) changes the angle of the text
        theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=8))+
        xlab("Number of historical studies (k) : Sample size (n)")+
        ylab(expression("Simulated coverage probability (" ~ hat(Psi)[r] ~ ")"))+
        scale_y_continuous(breaks=c(0.84, 0.87, 0.9, 0.93, 0.95, 0.97, 1.0))+
        scale_color_grey(start=0.8, end=0.2, name=expression(pi))+
        scale_shape_discrete(name=expression(pi))+
        geom_hline(yintercept=1.0, alpha=0.0001)

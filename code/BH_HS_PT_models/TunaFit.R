# Gets data from RAM Legacy Stock assessment database
# RLSADBv3.0 (assessment data only).xlsx
#

#### Get RAM data data ##

library("ggplot2")
library("tidyr")
library("dplyr")
source('Functions_Escapement.R')
dat  = read.csv("RLSADB.v3.0_timeSeries.csv")

#extracting Southern bluefin tuna variables
  ind.tuna = which(dat$stockid == 'SBT')
  n = length(ind.tuna)
  years = dat$year[ind.tuna]
  biomass = dat$TB[ind.tuna]
  Fish = dat$F[ind.tuna]
  TC = dat$TC[ind.tuna]

#escapement
  yield = TC[1:(n-1)]
  esc = biomass[1:(n-1)] - yield
  

#plot production data (biomass vs escapement)
  plot(esc, biomass[2:n], 
       xlim = c(0,max(esc)), ylim = c(0,max(biomass))
       )
  abline(0,1, col = 'red')
  
#fit data to functions
  lowerBound = c( 0, median(biomass) )
  upperBound = c( 1, 1.5*max(biomass) )
  fn = SSE.log
  
  optHS = optim( c( 0.1, max(biomass) ), fn = fn,
                 x = esc, y = biomass[2:n],
                 func = HockeyStick, std = 0.1, 
                 control = list( parscale =  c( 1, max(biomass))),lower = lowerBound, upper = upperBound, method = "L-BFGS-B" ) 

  optBH = optim( c( 0.1, max(biomass) ), fn = fn,
                 x = esc, y =  biomass[2:n],
                 func = BevertonHolt, std = 0.1, 
                 control = list( parscale = c( 1, max(biomass))),lower = lowerBound, upper = upperBound,method = "L-BFGS-B" )

  optPT = optim( c( 0.1, max(biomass), 1 ), fn = fn,
                 x = esc, y =  biomass[2:n],
                 func = PellaTomlinson, std = 0.1,
                 control = list( parscale = c( 1, max(biomass),1)), lower = c(lowerBound,.1), upper = c(upperBound, 1000), method = "L-BFGS-B" )

  
#plot production data (biomass vs escapement)
  pdf('FittedTuna.pdf')
    par(mfrow=c(1,1), mar = c(5,5,5,1))
    unit = 1000; mult = 1.1;
    lw = 2; cexl = 2; cexa=1.2;
    lims = c( min(esc, biomass), mult*max(biomass, optBH$par[2], optPT$par[2], optHS$par[2]))/unit
    plot( esc/unit, biomass[2:n]/unit, 
          xlim = lims, ylim = lims,
          xlab = 'Escapement (1,000 MT)',
          ylab = 'Biomass (1,000 MT)',
          lwd = lw, cex.axis = cexa, cex.lab = cexl)
    x = seq(0, mult*max(biomass, optBH$par[2], optPT$par[2], optHS$par[2]), by = 20)
    lines(x/unit, HockeyStick(x,optHS$par)/unit, col = 'blue', lwd = lw)
    lines(x/unit, BevertonHolt(x,optBH$par)/unit, col = 'black', lwd = lw)
    lines(x/unit, PellaTomlinson(x,optPT$par)/unit, col = 'green', lwd = lw)
    abline(0,1, lty = 3)  
  dev.off()
  
#Optimal escapment calculations
  oe.HS = OptEscHS(array(optHS$par, dim=c(1,1,3)))
  oe.BH = OptEscBH(array(optBH$par, dim=c(1,1,3)))
  oe.PT = OptEscPT(array(optPT$par, dim=c(1,1,3)))

#plot timeseries data
  pdf('TunaBiomass.pdf')
    par(mfrow=c(1,1), mar = c(5,5,5,1))
    plot(years, biomass/unit, 
         xlab = 'years', ylab = 'Biomass (1,000 MT)',
         lwd = lw, cex.axis = cexa, cex.lab = cexl, cex=1.5
    )
  dev.off()
  
  pdf('TunaCatch.pdf')
    par(mfrow=c(1,1), mar = c(5,5,5,1))
    plot(years, TC/unit, 
         xlab = 'years', ylab = 'Total catch (1,000 MT)',
         lwd = lw, cex.axis = cexa, cex.lab = cexl, cex=1.5
    )
  dev.off()
    
#Idea: Generate many scenarios where the fits are the true model 
#and test the optimal escapement strategies produced here 
  
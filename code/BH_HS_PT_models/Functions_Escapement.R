### Pella-Tomlinson Surplus Production functions
# (phi controls strength of weak allee effect phi=1 recovers logistic/schaeffer)
PellaTomlinson = function(s, pars){
  r = pars[1]; k = pars[2]; phi = pars[3];
  y = s + r* ((phi + 1)/phi) *s*(1 - (s/k)^phi)
  y[y < 0] = 1e-12
  return(y)
}

  OptEscPT = function(pars){
    r = pars[,,1]; k = pars[,,2]; phi = pars[,,3];
    opt = ifelse(r<0, 0,  (1/(phi+1))^(1/phi) * k)
    return(opt)
  }

### Hockey-stick functions
HockeyStick = function(s, pars){
  r = pars[1]; k = pars[2]; 
  y = (1+r)*s
  y[y > k] = k
  return(y)
}

  OptEscHS = function(pars){
    r = pars[,,1]; k = pars[,,2]; phi = pars[,,3];
    opt = ifelse(r<0, 0,  k/(1+r))
    return(opt)
  }

### Beverton-Holt functions
BevertonHolt = function(s, pars){
  r = pars[1]; k = pars[2];
  y = (1+r)*s / (1 + r*s/k)
  return(y)
}

  OptEscBH = function(pars){
    r = pars[,,1]; k = pars[,,2]; phi = pars[,,3];
    opt = ifelse(r<0, 0,  k/r*(sqrt((1+r)) -1) ) 
    return(opt)
  }


GenerateDataR = function(func, std, parms, n, n.sims){
  
  #generates y data and x data for a model
  #parms is the parameter values for the escapement function
  #std is the standard deviation of the lognormal random variable
  #n = number of data points
  #n.sims = number of times you do the simulation
  k = parms[2]
  
  #random escapement data (sorted for easier plotting)
    xDat = t( apply(matrix( runif(n.sim*n.dat,0,(1+1*std)*k), n.sim, n), 1, sort) )
  #  SimulateStrategy(func, n, optEsc, std, pars)
  #mean of normal rv such that lognomormal has mean 1
    mu = log( 1 / sqrt(1 + std^2) )
    sig = sqrt( log(1 + std^2) )
  #environmental stochasticisity
    z = matrix( rlnorm(n.sim*n.dat, meanlog = mu, sdlog = sig), n.sim, n) 
  #biomass data
    yDat = z*func(xDat, parms) 
  
    dat = array( NA, dim = c(n.sim, n, 2) )
    dat[,,1] = xDat
    dat[,,2] = yDat
  
  return(dat)
}

GenerateDataS = function(func, std, parms, n, n.sims){
  #function that generates escapement data by simulating
  #a population for n years, n.sims number of times
  
  #generates y data and x data for a model
  #parms is the parameter values for the escapement function
  #std is the standard deviation of the lognormal random variable
  #n = number of data points
  #n.sims = number of times you do the simulation
  k = parms[2]
  xDat = array( NA, dim = c(n.sim, n+1) );
  
  #random escapement data (sorted for easier plotting)
  for(i in 1:n.sims){
    xDat[i,] = SimulateStrategy(func, n, 1000*k, std, pars, output.type='x')[[3]]
  }
  
  dat = array( NA, dim = c(n.sim, n, 2) )
  dat[,,1] = xDat[,1:n]
  dat[,,2] = xDat[,2:(n+1)]
  
  return(dat)
}

SSE = function(p, x, y, func, std){
  
  #Sum squared error between data and the predicted model value
  #func = the escapement model (e.g. Beverton Holt)
  #p = parameters in the model 
  
  yMod = func(x,p)
  
  return( sum( (y - yMod)^2 ) )
  
}

SSE.log = function(p, x, y, func, std){
  
  #Sum squared error between data and the predicted model value
  #func = the escapement model (e.g. Beverton Holt)
  #p = parameters in the model 
  
  yMod = func(x,p)
  yMod[yMod<0] = 1e-300
  mu = log(1/sqrt(1+std^2))
  
  return( sum( ( log(y) - ( log(yMod) - mu ) )^2 ) )
  
}

FitData = function(funcs, dat, p0, type.fit, std){
  #returns the parameters that minimizes sum squared error 
  #between the fitted model and the data
  #funcs = list of escapement models
  #dat = array, dim1 = generating model, dim2 = simulation #,
  #             dim3 = data point #, dim4 = 1 if x, 2 if y
  #output = array, dim1 = generating model, dim2 = fitted model,
  #             dim3 = simulation number, dim4 = parameters
  
  
  d = dim(dat)
  fittedParms = array( NA, dim = c(d[1], d[1], d[2], length(p0)) )
  value = array( NA, dim = c(d[1], d[1], d[2], 1) )
  converge = array( NA, dim = c(d[1], d[1], d[2], 1) )
  
  #fitting minimizing residuals of unlogged or logged data
    if(type.fit=='log') SSEfun = SSE.log else SSEfun = SSE
  
  #perform fitting
    for(mt in 1:d[1]){
      for(mf in 1:d[1]){
      
      #for models with 2 parameters, remove the third from
      #the initial guess, and define the upper and lower bound 
      #for possible parameter guesses
        p0temp = p0
        lowerBound = c(0.0001, 10, 0.01);
        upperBound = c(3, 1000, 10);
        if(mf<3){ 
          p0temp = p0[-3]; 
          lowerBound = lowerBound[-3];
          upperBound = upperBound[-3];
        } 
        lp = length(p0temp)
      
      #do the curve fitting for each simulation
        for(i in 1:d[2]){
  
          optTemp = optim( p0temp, fn = SSEfun,
                           x = dat[mt,i,,1], y = dat[mt,i,,2],
                           func = funcs[[mf]], std = std,
                           lower = lowerBound, upper = upperBound,
                           method = "L-BFGS-B" )
  
          
          fittedParms[mt,mf,i,1:lp] = optTemp$par
          value[mt,mf,i,] = optTemp$value
          converge[mt,mf,i,] = optTemp$convergence
        }
    }
  }
  
  return(fittedParms)
  
}

GenerateZ = function(std, n){
  #generates lognormal random variables for an n term game, with specified std
  mu = log( 1 / sqrt(1 + std^2) )
  sig = sqrt( log(1 + std^2) )
  z = rlnorm(n+1, meanlog = mu, sdlog = sig)
  return(z)
}

SimulateStrategy = function(func, n, optEsc, std, pars, output.type='ty', z=NULL){
  #n = number of turns for harvest
  #y = yield
  #ty = total yield, 
  #pars = true parameters, [1] r [2] k [3] phi
  #p = fitted parameters
  #std is the mean of 
  #z is a vector of random variables, computed if not provided
    k = pars[2]
    
  #initialize variables, yield (y), biomass (x), noise (z)
    y = rep(0,n)
    x = rep(NA,n+1); x[1] = .3*k
    if(is.null(z)) z = GenerateZ(std,n)

  #simulate model
    for(j in 1:n){
      y[j]     = max(0, x[j] - optEsc) #min( x[j], )
      x[j + 1] = z[j + 1] * 
                 func( x[j] - y[j], pars )
    }
    
    ty = sum(y) + x[j + 1]
    
    #return either yield, total yield, or biomass time series
      # if (output.type=='ty') return(ty)
      # else if (output.type=='x') return(x)
      # else if (output.type=='y') return(y)
    return(list(ty, y, x))
} 

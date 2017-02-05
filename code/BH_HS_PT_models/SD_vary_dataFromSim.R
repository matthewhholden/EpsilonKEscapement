###### Parameters


rm(list = ls())


source('Functions_Escapement.R')
source('Function_Analysis.R')
#GenerateData
#outputs a list of data for each model

set.seed(1)

###########################################################################################
###########################################################################################
###########################################################################################

###
#parameters
  r = 1
  k = 100
  phi = 1
  std = .2
  pars = c(r, k, phi)
  escFunctions = c(BevertonHolt, HockeyStick, PellaTomlinson)
  optFunctions = c(OptEscBH, OptEscHS, OptEscPT)
  n.dat = 50
  n.sim = 500
  n.mod = length(escFunctions)
  n.turn = 10000
  sV = seq(0.05, .4, by = .05)
  n.r = length(sV)

#initialize variables to report results, variables that end in T are for the "true model"
#yield (total harvest), 
#fittedParms (parameters, for Bev-Holt and Hockey stick 3rd parameter is NA
# because it is a 2D model)
#y (harvest at each time step), 
#x (biomass at each time step)
#optEsc (optimal escapement rule given the parameters)  
  yield = array(NA, dim = c(n.r, n.mod, n.mod, n.sim) )
  yieldT = array(NA, dim = c(n.r, n.mod,n.sim) )
  fittedPars = array(NA, dim = c(n.r, n.mod, n.mod, n.sim, 3) )
  optEsc = array(NA, dim = c(n.r, n.mod, n.mod, n.sim) )
  optEscT = array(NA, dim = c(n.r, n.mod) )
  yT = array(NA, dim = c(n.r, n.mod, n.sim, n.turn));
  y = array(NA, dim = c(n.r, n.mod, n.mod, n.sim, n.turn));
  xT = array(NA, dim = c(n.r, n.mod, n.sim, n.turn+1));
  x = array(NA, dim = c(n.r, n.mod, n.mod, n.sim, n.turn+1));

  #Run the analysis for each value of std of env. stochasticity and save the results
    tic = Sys.time()
    for(i in 1:n.r){
      print(i)
      tempOpt = RunAnalyis( r=r, k=k, phi=phi, std=sV[i], escFunctions, optFunctions,
                            n.dat=n.dat, n.sim=n.sim, n.mod=n.mod, n.turn=n.turn,
                            type.fit='log', type.dat = 'sim' 
      )
      
      yield[i,,,] = tempOpt[[1]]
      optEsc[i,,,] = tempOpt[[2]]
      fittedPars[i,,,,] = tempOpt[[3]]
      yieldT[i,,] = tempOpt[[4]]
      optEscT[i,] = tempOpt[[5]]
      yT[i,,,] = tempOpt[[6]]
      y[i,,,,] = tempOpt[[7]]
      xT[i,,,] = tempOpt[[8]]
      x[i,,,,] = tempOpt[[9]]
      
    }
    print(Sys.time()-tic)

  #scale the yield from the fitted models by the yield achieved under 
  #the optimal escapement rule with parameters and model known    
    yield.scaled = array(NA, dim = c(n.r, n.mod, n.mod, n.sim) )
    for(i in 1:n.mod){ yield.scaled[,,i,] = yield[,,i,]/yieldT} 
    yieldm = apply(yield.scaled, 1:3, mean)
  
  #plot in plotting window before to PDF for debugging sake
    for(i in 1:n.mod){
      matplot(sV, yieldm[,i,] , ylim = c(0,1.3))
    }

### Plot a graph for each true model with the relative performance of the 
### incorrect model to the fitted true model and save them as seperate PDFs
  
  cx = 2; lw = 3; cl=2;ca=1.5; 
  xlimit = c(0, max(sV)); ylimit = c(0,1.0);
  cols = c('black', 'blue', 'green');
  names = c("Beverton-Holt", "Hockey-Stick", "Pella-Tomlinson")
  for(i in 1:n.mod){
    filename = paste('True_', i, '_relativeYield_sim_high_r.pdf', sep = '')
    addit=FALSE
    pdf(filename)
      par(mfrow=c(1,1), mar = c(5,5,5,1))
      for(j in 1:n.mod){
        if (j!=1) addit = TRUE;
        par(new=addit)
        plot(sV, yieldm[ , i , j ],
             ylim = ylimit, xlim = xlimit,
             xlab = 'Standard deviation of environmental noise',
             ylab = 'Fitted model yield / True model yield',
             main = paste('True model: ', names[i]),
             type = 'b',
             col = cols[j], pch =j, lty = j,
             cex = cx, lwd = lw, cex.lab = cl, cex.axis=ca, cex.main=cl)
      }
      legend('bottom', names, 
             pch=1:n.mod, lty = 1:n.mod, col = cols, lwd = lw, cex=cx, bty='n')
    dev.off()
  }


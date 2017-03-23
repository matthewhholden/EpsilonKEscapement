###### Parameters

# need to run TunaFit.R first

#GenerateData
#outputs a list of data for each model

set.seed(1)
source('Functions_Escapement.R')
source('Function_Analysis.R')


###########################################################################################
###########################################################################################
###########################################################################################

###
#parameters
  r = c( optBH$par[1], optHS$par[1], optPT$par[1] )
  k = c( optBH$par[2], optHS$par[2], optPT$par[2] )
  phi = c( NA, NA, optPT$par[3] )
  std = .1
  pars = cbind(r, k, phi)
  escFunctions = c(BevertonHolt, HockeyStick, PellaTomlinson)
  optFunctions = c(OptEscBH, OptEscHS, OptEscPT)
  n.sim = 100
  n.mod = length(escFunctions)
  n.turn = 10000
  
  #calculate optimal escapements from parameters
    oe.BH = OptEscBH(array(pars[1,], dim=c(1,1,3)))
    oe.HS = OptEscHS(array(pars[2,], dim=c(1,1,3)))
    oe.PT = OptEscPT(array(pars[3,], dim=c(1,1,3)))
    optEscT = c(oe.BH,oe.HS,oe.PT) #true optimal escapements from tuna model


  #Run the analysis for each value of std of env. stochasticity and save the results
  #####
  ## calculate Yield
  #####
  #Use optimal escapement rule in a n turn game (with all fish harvested at time n+1)  
  #mt is the true model index, mf is the fitted model index, i is the simulation #

  #initialize variables
  yield = array(NA, dim = c(n.mod, n.mod, n.sim))
  yieldT = array(NA, dim = c(n.mod, n.sim));
  yT = array(NA, dim = c(n.mod, n.sim, n.turn));
  y = array(NA, dim = c(n.mod, n.mod, n.sim, n.turn));
  xT = array(NA, dim = c(n.mod, n.sim, n.turn+1));
  x = array(NA, dim = c(n.mod, n.mod, n.sim, n.turn+1));
  
  start=Sys.time()
  for(i in 1:n.sim){
    
    z = GenerateZ(std, n.turn)
    
    for(mt in 1:n.mod){
      #yield generated from opt escapement of true model with parameters known
      outPutT = SimulateStrategy(func = escFunctions[[mt]], 
                                 n = n.turn, optEscT[mt], 
                                 std = std, pars = pars[mt,], z=z)
      yieldT[mt, i] = outPutT[[1]]
      yT[mt, i, ] = outPutT[[2]]
      xT[mt, i, ] = outPutT[[3]]
      
      for(mf in 1:n.mod){
        #yield generated from opt escapement of each fitted model 
        outPut = SimulateStrategy(func = escFunctions[[mt]], 
                                  n = n.turn, optEscT[mf], 
                                  std = std, pars = pars[mt,], z=z)
        yield[mt,mf,i] = outPut[[1]]
        y[mt,mf,i, ] = outPut[[2]]
        x[mt,mf,i, ] = outPut[[3]]
      }
    }
  }
  stop=Sys.time()-start; print(stop)

  #scale the yield from the fitted models by the yield achieved under 
  #the optimal escapement rule with parameters and model known    
    yield.scaled = array(NA, dim = c( n.mod, n.mod, n.sim) )
    for(i in 1:n.mod){ yield.scaled[,i,] = yield[,i,]/yieldT} 
    yieldm = apply(yield.scaled, 1:2, mean)
  
print(yieldm)

par(mfrow=c(3,3))
for(i in 1: n.mod){
  for(j in 1:n.mod){
    bins = seq(0,1, by=.025)
    hist(yield.scaled[i,j,], bins)
  }
}


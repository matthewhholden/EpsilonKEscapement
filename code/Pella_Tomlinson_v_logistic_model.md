Recruitment functions
---------------------

``` r
### Pella-Tomlinson Surplus Production function 
  # (phi controls strength of weak allee effect phi=1 recovers logistic/schaeffer)
  PT = function(s, pars){
    r = pars[1]; k = pars[2]; phi = pars[3];
    y = s + r* (phi/(phi + 1)) *s*(1 - (s/k)^phi)
    return(y)
  }
```

Parameters and Data
-------------------

``` r
  #parameters
    rho = 1       #discrete discount factor
    r = 0.25      #growth rate
    k = 100       #carying capacity
    sig = .4      #std of normal rv (for lognormal rvs)
    m = -.5*sig^2 #mean of normal rv such that lognomormal has mean 1
    phi = 10         #power in logistic function
    N=200
    parms = c(r,k,phi)
    
    #data (x = escapement, y = recruitment, z = random variable)
    x = runif(N,0,1.02*k)
    z = rlnorm(N,meanlog = m, sdlog = sig)
    #z = runif(N, .4,1.6) #use this for uniform data
    y = z*PT(x,parms)
    x1 = seq(0,k,length.out=N)
    y1 = PT(x1,parms)
  
  #plot data
    plot(x,y)
    lines(x1, y1)
```

![](Pella_Tomlinson_v_logistic_model_files/figure-markdown_github/parameters-1.png)

Optimal Escapement
------------------

``` r
  OptEscPT = k/r*(sqrt(rho*(1+r)) -1)
  OptEscL = k/(1+r)
  OptEscL = .5*k
  OptEscPTp = (1/(phi+1))^(1/phi) * k
```

Logistic escapement is 50 and Paella Tomilson Escapement is \(\left[\frac{1}{\phi+1}\right]^{1/\phi}k\), with \(\phi=\) 10 escapement is 78.6793442. Increasing \(\phi\) increases optimal escapment

![](Pella_Tomlinson_v_logistic_model_files/figure-markdown_github/PaellaTomilson-1.png)

Data fitting
------------

``` r
  ssePT = function(p,x,y){
    yMod = PT(x,p)
    return( sum( (y-yMod)^2 ) )
  }

  sseL = function(p,x,y){
    yMod = PT(x,c(p,1))
    return( sum( (y-yMod)^2 ) )
  }
  
  parInit = c(r,k)
  optOutPT = optim(c(parInit,phi), fn = ssePT, x=x, y=y)
  optOutL = optim(parInit, fn = sseL, x=x, y=y)
  
  #plot fit
    plot(x,y)
    xSort = sort(x)
    lines( xSort, PT(xSort, optOutPT$par ) )
    lines( xSort, PT(xSort, c(optOutL$par,1) ), col ='red')
```

![](Pella_Tomlinson_v_logistic_model_files/figure-markdown_github/SumSquaredErrorFunctions-1.png)

``` r
### Calculate Optimal Escapements

    
  rPT = optOutPT$par[1];
  kPT = optOutPT$par[2]; 
  phiPT = optOutPT$par[3];
  
  rL = optOutL$par[1];
  kL = optOutL$par[2]; 
  
  OptEscPT = kPT*(1/(1 + phiPT))^(1/phiPT)
  OptEscL = .5*kL
```

Optimal escapement in the fitted PT model is 53.7940424 and fitted Logistic is 44.0143076

Test Escapement Strategies
--------------------------

``` r
## calculate the sum of harvest for fitted PT model    
  xPT=rep(NA,N+1)
  xPT[1]=k    
  objPT=0    
  for(i in 1:N){
    xPT[i+1] = PT(xPT[i], c(r, k, phi))
    h = max(0, xPT[i+1] - OptEscPT)
    xPT[i+1] = xPT[i+1] - h
    objPT = objPT + h
  }    
  ObjPT = objPT + xPT[N+1] #take all fish at eNd

## calculate the sum of harvest for fitted L model   
  xL=rep(NA,N+1)
  xL[1]=k    
  objL=0    
  for(i in 1:N){
    xL[i+1] = PT(xL[i], c(r, k, 1))
    h = max(0, xL[i+1] - OptEscL)
    xL[i+1] = xL[i+1] - h
    objL = objL + h
  }  
  objL = objL + xL[N+1] #take all fish at eNd
  
## calculate the sum of harvest for perfect PT model   
  xPTp=rep(NA,N+1)
  xPTp[1]=k    
  objPTp=0    
  for(i in 1:N){
    xPTp[i+1] = PT(xPTp[i], c(r, k, phi))
    h = max(0, xPTp[i+1] - OptEscPTp)
    xPTp[i+1] = xPTp[i+1] - h
    objPTp = objPTp + h
  }  
  objPTp = objPTp + xPTp[N+1] #take all fish at eNd
```

The objective under perfect knowledge is 3334.9565082 and using the PT function but with parameters estimated from the data 2474.2266318 and using the logistic escapement rule with parameters estimated from the data, 712.9626573. This is despite the SSEs being 1.053478110^{5} and 1.055263310^{5} respectively. Note the above objective is calculated using deterministic dynamics. I will soon include stochastic dynamics below.

CI_hill<-function(data){
  
  #Define the model
  CI_hill<-function(time, state, pars){
    
    with(as.list(c(state, pars)),{
      
      #Fluorescent product/substrate
      dPf<-Vmax*S^n/(Km*(1+Pt/Kic) + S^n)#fluorescent product
      dS<--Vmax*S^n/(Km*(1+Pt/Kic) + S^n)
      dPt<-Vmax*S^n/(Km*(1+Pt/Kic) + S^n)#P-PO4
      
      return(list(c(dPf, dS, dPt)))
                 
    })
  }
  
  #Define the cost function
  cost<-function(x){
    ##the ODE models is run for each initial Substrate and Product concentration
    ##the model runs are stored here
    yhat<-data.frame(Pred=numeric(), Product=numeric(), Substrate=numeric())
    for(i in unique(data$Substrate)){
      for(n in unique(data$InhibitorSRP)){
        out<-as.data.frame(ode(y=c(Pf=0, S=i, Pt=n), 
                               parms = c(Vmax=x[1], Km=x[2], Kic=x[3], n=x[4]), 
                               CI_hill, times=sort(unique((data[(data$Substrate==i & data$InhibitorSRP==n), "time"]))))) 
      
      out<-out[, c("time", "Pf")]
      colnames(out)<-c("time", "Pred")
      outm<-merge(out, data[(data$Substrate==i & data$InhibitorSRP==n), c("time", "Product", "Substrate")], by = c("time"))[,-1]
      yhat<-rbind(yhat, outm)
      }
    }
    
    RMSE<-with(yhat, sum((((Pred-Product)/Substrate)^2), na.rm = T))
    return(RMSE)
  }
  
  #Use MCMC to define ranges of possible model parameters
  par_mcmc<-modMCMC(f=cost, p=c(1e-2, 20, 20, 1), 
                    lower=c(1e-3, 1e-2, 1e-2, 0.01),
                    upper=c(100, 500, 500, 10), niter=10000)
  #lower and upper limits for parameters are extracted
  pl<-as.numeric(summary(par_mcmc)["min",])
  pu<-as.numeric(summary(par_mcmc)["max",])
  
  #these limits are used to find global optimum by DEoptim
  # opt_par<-DEoptim(fn=cost, lower=pl, upper=pu, 
  #                control = c(itermax = 10000, steptol = 50, reltol = 1e-8, 
  #                            trace=FALSE, strategy=3, NP=250))
  
  
  #these limits are used to find global optimum by rgenoud
  #opt_par<-genoud(fn=cost, print.level = 0, pop.size=1e6, max=FALSE, nvars=6, Domains = cbind(pl, pu),
  #boundary.enforcement = 2)
  
  #these limits are used to find global optimum by ABCotpim
  opt_par<-abc_optim(fn=cost, par=as.numeric(summary(par_mcmc)["mean",]), lb=pl, ub=pu, maxCycle = 1e6)
  
  #Calculate goodness of correspondence
  goodness<-function(x){
    yhat<-data.frame(time = numeric(), Pred=numeric(), Product=numeric(), Substrate=numeric(), InhibitorSRP=numeric(),
                     Catchment=character(), Horizon=character(), SRP=numeric())
    
    for(i in unique(data$Substrate)){
      for(n in unique(data$InhibitorSRP)){
        out<-as.data.frame(ode(y=c(Pf=0, S=i, Pt=n), 
                               parms = c(Vmax=x[1], Km=x[2], Kic=x[3], n=x[4]), 
                               CI_hill, times=sort(unique((data[(data$Substrate==i & data$InhibitorSRP==n), "time"]))))) 
        out<-out[, c("time", "Pf", "Pt")]
        colnames(out)<-c("time", "Pred", "SRP")
        outm<-merge(out, data[(data$Substrate==i & data$InhibitorSRP==n), c("time", "Product")], by = c("time"))
        outm$Substrate<-rep(i, times=nrow(outm))
        outm$InhibitorSRP<-rep(n, times=nrow(outm))
        yhat<-rbind(yhat, outm)
      }
    }
    
    yhat$Catchment<-rep(data$Catchment[1], times=nrow(yhat))
    yhat$Horizon<-rep(data$Horizon[1], times=nrow(yhat))
    
    SSres=with(yhat, sum((((Product-Pred)/Substrate)^2), na.rm = T))
    SStot=with(yhat, sum((((Product-mean(Product, na.rm = T))/Substrate)^2), na.rm = T))
    ll=with(yhat, -sum((((Product-Pred)/Substrate)^2), na.rm = T)/2/(sd(Product/Substrate, na.rm = T)^2))
    R2<-1-SSres/SStot
    N<-length(x)
    AIC<-2*N-2*ll
    Gfit<-c(R2=R2, N=N, AIC=AIC, ll=ll, SSres=SSres, SStot=SStot)
    goodness_out<-list(Yhat=yhat, Gfit=Gfit)
    return(goodness_out)
  }
  
  #Parameters<-opt_par$optim$bestmem#DEoptim algorithm
  Parameters<-opt_par$par#genoud/ABC algorithm
  names(Parameters)<-c("Vmax", "Kmf", "Kic", "n")
  
  out_all<-list(Parameters = Parameters,
                #Goodness = goodness(as.numeric(opt_par$optim$bestmem)),#DEoptim algorithm
                Goodness = goodness(as.numeric(opt_par$par)),#genoud/ABC algorithm
                MCMC = par_mcmc)
  
  return(out_all)
}

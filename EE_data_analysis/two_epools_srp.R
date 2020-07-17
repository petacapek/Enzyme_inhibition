two_epools_srp<-function(data){
  
  #Define the model
  two_esrp_model<-function(time, state, pars){
    
    with(as.list(c(state, pars)),{
      
      dPs<-(Vmax1*S/(Km1 + S))+(Vmax2*S/(Km2*(1+Pt/Kic) + S))#fluorescent product
      dPt<-(Vmax1*S/(Km1 + S))+(Vmax2*S/(Km2*(1+Pt/Kic) + S)) + alfa#SRP
      dS<--(Vmax1*S/(Km1 + S))-(Vmax2*S/(Km2*(1+Pt/Kic) + S))
                 
      return(list(c(dPs, dPt, dS)))
                 
    })
  }
  
  #Define the cost function
  cost<-function(x){
    ##the ODE models is run for each initial Substrate and Product concentration
    ##the model runs are stored here
    yhat<-data.frame(Pred=numeric(), Product=numeric(), Substrate=numeric())
    for(i in unique(data$Substrate)){
      for(n in unique(data$InhibitorSRP)){
        out<-as.data.frame(ode(y=c(Ps=0, Pt=n, S=i), parms = c(Vmax1=x[1], Km1=x[2], Vmax2=x[3], Km2=x[4], Kic=x[5],
                                                               alfa=mean(data[(data$Substrate==i & data$InhibitorSRP==n), "slope"])), 
                           two_esrp_model, times=sort(unique((data[(data$Substrate==i & data$InhibitorSRP==n), "time"]))))) 
      
      out<-out[, c("time", "Ps")]
      colnames(out)<-c("time", "Pred")
      outm<-merge(out, data[(data$Substrate==i & data$InhibitorSRP==n), c("time", "Product", "Substrate")], by = c("time"))[,-1]
      yhat<-rbind(yhat, outm)
      }
    }
    
    RMSE<-with(yhat, sum((((Pred-Product)/Substrate)^2), na.rm = T))
    return(RMSE)
  }
  
  #Use MCMC to define ranges of possible model parameters
  par_mcmc<-modMCMC(f=cost, p=c(1e-2, 20, 1e-2, 20, 20), 
                    lower=c(1e-3, 1e-2, 1e-3, 1e-2, 1e-2),
                    upper=c(100, 500, 100, 500, 500), niter=10000)
  #lower and upper limits for parameters are extracted
  pl<-summary(par_mcmc)["min",]
  pu<-summary(par_mcmc)["max",]
  
  #these limits are used to find global optimum by DEoptim
  opt_par<-DEoptim(fn=cost, lower=pl, upper=pu, 
                   control = c(itermax = 10000, steptol = 50, reltol = 1e-8, 
                               trace=FALSE, strategy=3, NP=250))
  
  #Calculate goodness of correspondence
  goodness<-function(x){
    yhat<-data.frame(time = numeric(), Pred=numeric(), Product=numeric(), Substrate=numeric(), InhibitorSRP=numeric(),
                     Catchment=character(), Horizon=character(), Pt=numeric())
    
    for(i in unique(data$Substrate)){
      for(n in unique(data$InhibitorSRP)){
        out<-as.data.frame(ode(y=c(Ps=0, Pt=n, S=i), parms = c(Vmax1=x[1], Km1=x[2], Vmax2=x[3], Km2=x[4], Kic=x[5],
                                                               alfa=mean(data[(data$Substrate==i & data$InhibitorSRP==n), "slope"])), 
                           two_esrp_model, times=sort(unique((data[(data$Substrate==i & data$InhibitorSRP==n), "time"]))))) 
        out<-out[, c("time", "Ps", "Pt")]
        colnames(out)<-c("time", "Pred", "Pt")
        outm<-merge(out, data[(data$Substrate==i & data$InhibitorSRP==n), c("time", "Product")], by = c("time"))
        outm$Substrate<-rep(i, times=nrow(outm))
        outm$InhibitorSRP<-rep(i, times=nrow(outm))
        yhat<-rbind(yhat, outm)
      }
    }
    
    yhat$Catchment<-rep(data$Catchment[1], times=nrow(yhat))
    yhat$Horizon<-rep(data$Horizon[1], times=nrow(yhat))
    
    SSres=with(yhat, sum(((Product-Pred)^2), na.rm = T))
    SStot=with(yhat, sum(((Product-mean(Product, na.rm = T))^2), na.rm = T))
    ll=with(yhat, -sum(((Product-Pred)^2), na.rm = T)/2/(sd(Product, na.rm = T)^2))
    R2<-1-SSres/SStot
    N<-length(x)
    AIC<-2*N-2*ll
    Gfit<-c(R2=R2, N=N, AIC=AIC, ll=ll, SSres=SSres, SStot=SStot)
    goodness_out<-list(Yhat=yhat, Gfit=Gfit)
    return(goodness_out)
  }
  
  Parameters<-opt_par$optim$bestmem
  names(Parameters)<-c("Vmax1", "Km1", "Vmax2", "Km2", "Kic")
  
  out_all<-list(Parameters = Parameters,
                Goodness = goodness(as.numeric(opt_par$optim$bestmem)),
                MCMC = par_mcmc)
  
  return(out_all)
}

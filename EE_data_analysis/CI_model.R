CI_model<-function(data, parameters){
  
  #Define the model
  CI_model<-function(time, state, pars){
    
    with(as.list(c(state, pars)),{
      
      #Fluorescent product/substrate
      dPf<-Vmax*S/(Km*(1+SRP/Kic) + S)#fluorescent product
      dS<--Vmax*S/(Km*(1+SRP/Kic) + S)
      
      #SRP/Porg
      dSRP<-Vmax*S/(Km*(1+SRP/Kic) + S)#SRP
      
      return(list(c(dPf, dS, dSRP)))
      
    })
  }
  
  #Calculate goodness of correspondence
  goodness<-function(x){
    yhat<-data.frame(time = numeric(), Pred=numeric(), Product=numeric(), Substrate=numeric(), InhibitorSRP=numeric(),
                     Catchment=character(), Horizon=character(), SRP=numeric())
    
    for(i in unique(data$Substrate)){
      for(n in unique(data$InhibitorSRP)){
        out<-as.data.frame(ode(y=c(Pf=0, S=i, SRP=n), parms = c(Vmax=x[1], Km=x[2], Kic=x[3]), 
                               CI_model, times=sort(unique((data[(data$Substrate==i & data$InhibitorSRP==n), "time"])))))  
        out<-out[, c("time", "Pf", "SRP")]
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
  
  
  out_all<-goodness(as.numeric(parameters))
  
  return(out_all)
}

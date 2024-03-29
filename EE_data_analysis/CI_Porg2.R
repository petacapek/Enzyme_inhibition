CI_Porg2<-function(data, parameters){
  
  #Define the model
  CI_Porg_model<-function(time, state, pars){
    
    with(as.list(c(state, pars)),{
      
      #Fluorescent product/substrate
      dPf<-Vmax*S/(Kmf*(1+SRP/Kic)*(1 + Porg/Kmorg) + S)#fluorescent product
      dS<--Vmax*S/(Kmf*(1+SRP/Kic)*(1 + Porg/Kmorg) + S)
      
      #SRP/Porg
      dSRP<-Vmax*Porg/(Kmorg*(1+SRP/Kic)*(1 + Pf/Kmf) + Porg) + 
        Vmax*S/(Kmf*(1+SRP/Kic)*(1 + Porg/Kmorg) + S)#SRP
      dPorg<--Vmax*Porg/(Kmorg*(1+SRP/Kic)*(1 + Pf/Kmf) + Porg)
                 
      return(list(c(dPf, dS, dSRP, dPorg)))
                 
    })
  }
  
  #Calculate goodness of correspondence
  goodness<-function(x){
    yhat<-data.frame(time = numeric(), Pred=numeric(), Product=numeric(), Substrate=numeric(), InhibitorSRP=numeric(),
                     Catchment=character(), Horizon=character(), SRP=numeric(), Porg=numeric())
    
    for(i in unique(data$Substrate)){
      for(n in unique(data$InhibitorSRP)){
        out<-as.data.frame(ode(y=c(Pf=0, S=i, SRP=n,
                                   Porg=mean(data[(data$Substrate==i & data$InhibitorSRP==n), "Porg"])), parms = c(Vmax=x[1], Kmf=x[2], Kic=x[3], Kmorg=x[4]), 
                               CI_Porg_model, times=sort(unique((data[(data$Substrate==i & data$InhibitorSRP==n), "time"])))))  
        out<-out[, c("time", "Pf", "SRP", "Porg")]
        colnames(out)<-c("time", "Pred", "SRP", "Porg")
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

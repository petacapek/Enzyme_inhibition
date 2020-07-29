two_epools2<-function(data, parameters){
  
  #Define the model
  two_e_model<-function(time, state, pars){
    
    with(as.list(c(state, pars)),{
      
      dPs<-(Vmax1*S/(Km1 + S))+(Vmax2*S/(Km2*(1+Pt/Kic) + S))#fluorescent product
      dPt<-(Vmax1*S/(Km1 + S))+(Vmax2*S/(Km2*(1+Pt/Kic) + S))#SRP
      dS<--(Vmax1*S/(Km1 + S))-(Vmax2*S/(Km2*(1+Pt/Kic) + S))
                 
      return(list(c(dPs, dPt, dS), Slow=(Vmax1*S/(Km1 + S)), Fast=(Vmax2*S/(Km2*(1+Pt/Kic) + S))))
                 
    })
  }
  
  #Calculate goodness of correspondence
  goodness<-function(x){
    yhat<-data.frame(time = numeric(), Pred=numeric(), Product=numeric(), Substrate=numeric(), InhibitorSRP=numeric(),
                     Fast=numeric(), Slow=numeric())
    
    for(i in unique(data$Substrate)){
      for(n in unique(data$InhibitorSRP)){
        out<-as.data.frame(ode(y=c(Ps=0, Pt=n, S=i), parms = c(Vmax1=x[1], Km1=x[2], Vmax2=x[3], Km2=x[4], Kic=x[5]), 
                           two_e_model, times=sort(unique((data[(data$Substrate==i & data$InhibitorSRP==n), "time"]))))) 
        out<-out[, c("time", "Ps", "Fast", "Slow")]
        colnames(out)<-c("time", "Pred", "Fast", "Slow")
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

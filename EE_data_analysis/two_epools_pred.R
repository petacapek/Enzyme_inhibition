two_epools_pred<-function(data, parameters){
  
  #Define the model
  two_e_model<-function(time, state, pars){
    
    with(as.list(c(state, pars)),{
      
      dPt<-(Vmax1*Porg/(Km1 + Porg))+(Vmax2*Porg/(Km2*(1+Pt/Kic) + Porg))
      #dPorg<--(Vmax1*Porg/(Km1 + Porg))-(Vmax2*Porg/(Km2*(1+Pt/Kic) + Porg))
             
      return(list(c(dPt)))
                 
    })
  }
  
  #Calculate goodness of correspondence
  goodness<-function(x){
    yhat<-data.frame(Time = numeric(), Pred=numeric(), SRP_o=numeric())
    
    for(n in unique(data$Inhibitor)){
        out<-as.data.frame(ode(y=c(Pt=n), 
                               parms = c(Vmax1=x[1], Km1=x[2], Vmax2=x[3], Km2=x[4], Kic=x[5],
                                         Porg=mean(data[(data$Inhibitor==n), "Porg2"], na.rm=T)), 
                               two_e_model, times=sort(unique((data[(data$Inhibitor==n), "Time"]))))) 
        out<-out[, c("time", "Pt")]
        colnames(out)<-c("Time", "Pred")
        outm<-merge(out, data[(data$Inhibitor==n), c("Time", "SRP_o")], by = c("Time"))
        outm$Inhibitor<-rep(n, times=nrow(outm))
        yhat<-rbind(yhat, outm)
      }
    
    
    yhat$Catchment<-rep(data$Catchment[1], times=nrow(yhat))
    yhat$Horizon<-rep(data$Horizon[1], times=nrow(yhat))
    
    SSres=with(yhat, sum((((SRP_o-Pred))^2), na.rm = T))
    SStot=with(yhat, sum((((SRP_o-mean(SRP_o, na.rm = T)))^2), na.rm = T))
    ll=with(yhat, -sum((((SRP_o-Pred))^2), na.rm = T)/2/(sd(SRP_o, na.rm = T)^2))
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

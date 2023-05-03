def CO (y, t, pars):
    #define initial pools
    S=y[0];    P=y[1];	  FP=y[2];  
    #define parameters
    Vmax=pars[0]; 
    Km=pars[1]; 
    n=pars[2];    
    #Fluxes
    decay=Vmax*S**n/(Km + S**n)
        
    #Define derivatives
    dSdt = -decay
    dPdt = decay
    dFPdt = decay
         
    return dSdt, dPdt, dFPdt;

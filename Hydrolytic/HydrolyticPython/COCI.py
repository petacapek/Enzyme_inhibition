def COCI (y, t, pars):
    #define initial pools
    S=y[0];    P=y[1];	  FP=y[2];  
    #define parameters
    Vmax=pars[0]; 
    Km=pars[1];
    Ki=pars[2];
    n=pars[3];     
    #Fluxes
    decay=Vmax*S**n/(Km*(1 + P/Ki) + S**n)
        
    #Define derivatives
    dSdt = -decay
    dPdt = decay
    dFPdt = decay
         
    return dSdt, dPdt, dFPdt;

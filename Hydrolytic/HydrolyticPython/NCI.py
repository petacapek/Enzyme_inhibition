def NCI (y, t, pars):
    #define initial pools
    S=y[0];    P=y[1];	  FP=y[2];  
    #define parameters
    Vmax=pars[0]; 
    Km=pars[1];
    Kic=pars[2];
    Kiu=pars[3];     
    #Fluxes
    decay=Vmax*S/(Km*(1 + P/Kic) + S*(1 + P/Kiu))
        
    #Define derivatives
    dSdt = -decay
    dPdt = decay
    dFPdt = decay
         
    return dSdt, dPdt, dFPdt;

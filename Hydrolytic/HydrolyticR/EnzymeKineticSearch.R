EnzymeKineticSearch <- function(data){
  #================================ODE for different enzyme kinetic 
  source_python("../HydrolyticPython/WI.py")
  source_python("../HydrolyticPython/CI.py")
  source_python("../HydrolyticPython/UCI.py")
  source_python("../HydrolyticPython/NCI.py")
  source_python("../HydrolyticPython/SI.py")
  #===========================================================================================
  #================================Parameter estimation
  source("WIoptim.R")
  source("CIoptim.R")
  source("UCIoptim.R")
  source("NCIoptim.R")
  source("SIoptim.R")
  
  #Results are stored here
  Parms <- data.frame(Model = c("WI", "CI", "UCI", "NCI", "SI"),
                      Vmax = numeric(length = 5),
                      Km = numeric(length = 5),
                      Ki = numeric(length = 5),
                      Kiu = numeric(length = 5))
  Parms[Parms$Model == "WI", c("Vmax", "Km")] <- WIoptim(data)
  Parms[Parms$Model == "CI", c("Vmax", "Km", "Ki")] <- CIoptim(data)
  Parms[Parms$Model == "UCI", c("Vmax", "Km", "Ki")] <- UCIoptim(data)
  Parms[Parms$Model == "NCI", c("Vmax", "Km", "Ki", "Kiu")] <- NCIoptim(data)
  Parms[Parms$Model == "SI", c("Vmax", "Km", "Ki")] <- SIoptim(data)
  #===========================================================================================
  #================================Goodness of fit
  source("WIfit.R")
  source("CIfit.R")
  source("UCIfit.R")
  source("NCIfit.R")
  source("SIfit.R")
  
  #Results are stored here
  Goodness <- data.frame(Model = c("WI", "CI", "UCI", "NCI", "SI"),
                      R2 = numeric(length = 5),
                      R2adj = numeric(length = 5),
                      ll = numeric(length = 5),
                      AIC = numeric(length = 5),
                      Fnorm = numeric(length = 5),
                      n = numeric(length = 5),
                      p = numeric(length = 5))
  Goodness[Goodness$Model == "WI", c("R2", "R2adj", "ll", "AIC", "Fnorm", "n", "p")] <- WIfit(data)
  Goodness[Goodness$Model == "CI", c("R2", "R2adj", "ll", "AIC", "Fnorm", "n", "p")] <- CIfit(data)
  Goodness[Goodness$Model == "UCI", c("R2", "R2adj", "ll", "AIC", "Fnorm", "n", "p")] <- UCIfit(data)
  Goodness[Goodness$Model == "NCI", c("R2", "R2adj", "ll", "AIC", "Fnorm", "n", "p")] <- NCIfit(data)
  Goodness[Goodness$Model == "SI", c("R2", "R2adj", "ll", "AIC", "Fnorm", "n", "p")] <- SIfit(data)
  
  return(list(Parms, Goodness))
}
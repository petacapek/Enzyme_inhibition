EnzymeKineticSearch <- function(data){
  #================================Parameter estimation
  source("WIoptim.R")
  source("CIoptim.R")
  source("UCIoptim.R")
  source("NCIoptim.R")
  source("SIoptim.R")
  source("COoptim.R")
  source("COCIoptim.R")
  
  #Results are stored here
  Parms <- data.frame(Model = c("WI", "CI", "UCI", "NCI", "SI", "CO", "COCI"),
                      Vmax = numeric(length = 7),
                      Km = numeric(length = 7),
                      Ki = numeric(length = 7),
                      Kiu = numeric(length = 7),
                      n = numeric(length = 7))
  Parms[Parms$Model == "WI", c("Vmax", "Km")] <- WIoptim(data)
  Parms[Parms$Model == "CI", c("Vmax", "Km", "Ki")] <- CIoptim(data)
  Parms[Parms$Model == "UCI", c("Vmax", "Km", "Ki")] <- UCIoptim(data)
  Parms[Parms$Model == "NCI", c("Vmax", "Km", "Ki", "Kiu")] <- NCIoptim(data)
  Parms[Parms$Model == "SI", c("Vmax", "Km", "Ki")] <- SIoptim(data)
  Parms[Parms$Model == "CO", c("Vmax", "Km", "n")] <- COoptim(data)
  Parms[Parms$Model == "COCI", c("Vmax", "Km", "Ki", "n")] <- COoptim(data)
  #===========================================================================================
  #================================Goodness of fit
  source("WIfit.R")
  source("CIfit.R")
  source("UCIfit.R")
  source("NCIfit.R")
  source("SIfit.R")
  source("COfit.R")
  source("COCIfit.R")
  
  #Results are stored here
  Goodness <- data.frame(Model = c("WI", "CI", "UCI", "NCI", "SI", "CO", "COCI"),
                      R2 = numeric(length = 7),
                      R2adj = numeric(length = 7),
                      ll = numeric(length = 7),
                      AIC = numeric(length = 7),
                      Fstat = numeric(length = 7),
                      n = numeric(length = 7),
                      p = numeric(length = 7))
  Goodness[Goodness$Model == "WI", c("R2", "R2adj", "ll", "AIC", "Fnorm", "n", "p")] <- WIfit(data, as.numeric(Parms[1, c(2, 3)]))
  Goodness[Goodness$Model == "CI", c("R2", "R2adj", "ll", "AIC", "Fnorm", "n", "p")] <- CIfit(data, as.numeric(Parms[2, c(2, 3, 4)]))
  Goodness[Goodness$Model == "UCI", c("R2", "R2adj", "ll", "AIC", "Fnorm", "n", "p")] <- UCIfit(data, as.numeric(Parms[3, c(2, 3, 4)]))
  Goodness[Goodness$Model == "NCI", c("R2", "R2adj", "ll", "AIC", "Fnorm", "n", "p")] <- NCIfit(data, as.numeric(Parms[4, c(2, 3, 4, 5)]))
  Goodness[Goodness$Model == "SI", c("R2", "R2adj", "ll", "AIC", "Fnorm", "n", "p")] <- SIfit(data, as.numeric(Parms[5, c(2, 3, 4)]))
  Goodness[Goodness$Model == "CO", c("R2", "R2adj", "ll", "AIC", "Fnorm", "n", "p")] <- COfit(data, as.numeric(Parms[6, c(2, 3, 6)]))
  Goodness[Goodness$Model == "COCI", c("R2", "R2adj", "ll", "AIC", "Fnorm", "n", "p")] <- COfit(data, as.numeric(Parms[7, c(2, 3, 4, 6)]))
  
  return(list(Parms = Parms, Goodness = Goodness))
}
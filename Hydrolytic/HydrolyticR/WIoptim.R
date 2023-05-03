WIoptim <- function(data){
  #Define a matrix (Uqs) with unique combinations of initial substrate and product concentrations
  Uqs <- as.matrix((unique(data[order(data$SubstrateInitial, 
                                      data$InhibitorInitial, 
                                      data$ReactionTime), 
                                c("SubstrateInitial", "InhibitorInitial")])))
  #Create a matrix of observations (Obs)
  ##number of columns represent number of unique combinations of initial substrate and inhibitor concentrations defined in Uqs
  ##rows are sorted by time
  Obs <- matrix(as.numeric(data[order(data$SubstrateInitial, 
                                      data$InhibitorInitial, 
                                      data$ReactionTime), "Product"]), ncol = dim(Uqs)[1])
  #Create a matrix of weights (W) with the same dimensions as Obs 
  ##Weights are defined as substrate/inhibitor - specific standard deviations
  W <- matrix(rep(apply(Obs, 2, sd, na.rm = T), each = dim(Obs)[1]), nrow = dim(Obs)[1], ncol = dim(Obs)[2])
  #Create a matrix of reaction  times (RT) with the same dimensions as (Obs)
  RT <- matrix(as.numeric(data[order(data$SubstrateInitial, 
                                     data$InhibitorInitial, 
                                     data$ReactionTime), "ReactionTime"]), ncol = dim(Uqs)[1])/60
  #Define cost function
  Cost <- function(x){
    Yhat <- matrix(unlist(lapply(1:dim(Uqs)[1], function(i){
      scpy$odeint(WI, y0 = c(Uqs[i,1], Uqs[i,2], 0), t = as.numeric(RT[,i]), args = tuple(x))[,3]}))
      , ncol = dim(Uqs)[1])
    return(sum((Obs - Yhat)^2/2/W^2, na.rm = T))
  }
  
  ##First guess of model parameters by MCMC 
  Guess <- ifelse(data$Enzyme == "Beta-Glucosidase",
                 modMCMC(Cost, p = c(1, 10), lower = c(1e-3, 1e-1), upper = c(10, 100), niter = 10000),
                 modMCMC(Cost, p = c(0.01, 10), lower = c(1e-5, 1e-1), upper = c(1, 100), niter = 10000))
  
  ##Estimate
  Optimized <- abc_optim(fn = Cost,
                      par = as.numeric(apply(Guess[[1]], 2, mean)), 
                      lb = as.numeric(apply(Guess[[1]], 2, min)), 
                      ub = as.numeric(apply(Guess[[1]], 2, max)))
  
  return(Optimized$par)
}
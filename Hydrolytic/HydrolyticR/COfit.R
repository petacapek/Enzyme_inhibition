COfit <- function(data, x){
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
  #Run simulation with optimized parameters 
  Yhat <- matrix(unlist(lapply(1:dim(Uqs)[1], function(i){
    scpy$odeint(CO, y0 = c(Uqs[i,1], Uqs[i,2], 0), t = as.numeric(RT[,i]), args = tuple(x))[,3]})), 
    ncol = dim(Uqs)[1])
  
  #Calculate goodness of fit indices
  ##Weights
  What <- matrix(rep(apply(Yhat, 2, sd, na.rm = T), each = dim(Yhat)[1]), nrow = dim(Yhat)[1], ncol = dim(Yhat)[2])
  ##Means
  M <- matrix(rep(apply(Obs, 2, mean, na.rm = T), each = dim(Obs)[1]), nrow = dim(Obs)[1], ncol = dim(Obs)[2])
  Mhat <- matrix(rep(apply(Yhat, 2, mean, na.rm = T), each = dim(Yhat)[1]), nrow = dim(Yhat)[1], ncol = dim(Yhat)[2])
  ##goodness of fit
  ###R2 adjusted
  R2 = 1 - (sum((Obs - Yhat)^2, na.rm = T)/sum((Obs - M)^2, na.rm = T))
  R2adj = 1 - ((1 - R2)*((length(Obs) - 1)/(length(Obs) - 1 - length(x))))
  ###Log-Likelihood
  ll = - sum((Obs - Yhat)^2/2/W^2, na.rm = T)
  ###AIC
  AIC = -2*ll + 2*length(x)
  ###normalized F
  Fnorm = sum((Obs - Yhat)^2, na.rm = T)
                                     
  return(c(R2 = R2, R2adj = R2adj, ll = ll, AIC = AIC, Fnorm = Fnorm, n = length(Obs), p = length(x)))
}
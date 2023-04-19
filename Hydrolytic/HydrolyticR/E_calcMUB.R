E_calcMUB<-function(dataset, MUBconc, APconc, Iconc, Nmeasure, DW, empty){#
  
  #Loop function
  lfunction<-function(sheet){
    
    #Reading data
    ##Careful, sheets are readed from the opposite direction (last measurements represent first sheets) 
    e_0<-read.xlsx(xlsxFile = dataset, sheet=sheet,
                   rows=46:53, cols = 2:13, colNames = F)
    # e_0<-read.xlsx(xlsxFile = "/mnt/580CBE2464C5F83D/pracovni/data_statistika/MalakBC/MalakEnz/Soil2MUF-G25112022.xlsx", sheet=16,
    #               rows=46:53, cols = 2:13, colNames = F)
    #Rearanging the data
    e_0r<-data.frame(Slurry = c(rep("NO", 6), rep("YES", 96-5-6-4), rep("NO", 5)), 
                     Calibration = c(rep("TRUE", times = 12), #define the data set for calibration and measurement
                                     rep("FALSE", times = 80)),
                     C_MUB = c(rep(MUBconc, times = 2), #concentration of MUB in umols per litter after 5 times dilution by soil slurry or water 
                               rep(NA, times = 80)),
                     C_AP = c(rep(NA, times = 12),
                              rep(c(rep(APconc, each = 3)), times = 5),#concentrations of added substrate in umols per gram dry soil
                              APconc),
                     C_I = c(rep(NA, times = 12), #concentrations of added inhibitor in umols per gram dry soil
                             rep(Iconc, each = 15),
                             rep(NA, times=5)),
                     E = c(as.numeric(e_0[1, 1:6]), #Calibration for controls
                           as.numeric(e_0[1, 7:12]), #Calibration for Slurry
                           as.numeric(e_0[2, 1:12]), as.numeric(e_0[3, 1:3]), #I = 0
                           as.numeric(e_0[3, 4:12]), as.numeric(e_0[4, 1:6]), #I1
                           as.numeric(e_0[4, 7:12]), as.numeric(e_0[5, 1:9]), #I2
                           as.numeric(e_0[5, 10:12]), as.numeric(e_0[6, 1:12]), #I3
                           as.numeric(e_0[7, 1:12]), as.numeric(e_0[8, 1:3]), #I4
                           as.numeric(e_0[8, 4:8]) #AP controls
                     ),
                     TimeOff = (c(1:12, 24:13, 25:36, 48:37, 49:60, 72:61, 73:84, 92:85)*33.34/91-(33.34/91))/60 #Time offset - it results from the fact that wells are measured separately one-by-one. The time difference (between the first and last measured well) is 33.34 seconds. Offset id defined in minutes
                     )
    ##Add time stamp
    xdate<-read.xlsx(xlsxFile = dataset, 
                     sheet=sheet,
                     rows=42, cols = 5, colNames = F)
    
    # xdate<-read.xlsx(xlsxFile = "/mnt/580CBE2464C5F83D/pracovni/data_statistika/MalakBC/MalakEnz/Soil1MUF-G21112022.xlsx", 
    #                  sheet=20,
    #                  rows=42, cols = 5, colNames = F)
    
    e_0r$Date<-as.POSIXct(xdate[1,1], format = "%Y-%m-%d %H:%M:%OS")
    
    return(e_0r)
  }
  
  #Go through the loop
  ID=seq(from=Nmeasure-empty, to=1)
  
  #first measurement
  all<-lfunction(sheet = ID[1])
  
  #all other measurements
  for(i in ID[-1]){
    all<-rbind(all, lfunction(sheet = i))
    #print(i)
  }
  
  #Calculating time differences in minutes - from first (time 0) to last time of measurement
  all$TimDiff <- as.numeric(difftime(all$Date, all$Date[1], units = "min")) + all$TimeOff
  
  #divide data set into calibrations, controls and real data
  calsCntr<-subset(all, Calibration == "TRUE" & Slurry=="NO")
  calsSl<-subset(all, Calibration == "TRUE" & Slurry!="NO")
  contr<-subset(all, Calibration != "TRUE" & Slurry=="NO")
  d<-subset(all, Calibration != "TRUE" & Slurry!="NO")
  
  #The slurry is added to substrate by rows so the first time of measurement is very close to start of reaction
  ##This creates another time offset that needs to be accounted for in following loop
  d$TimDiffFirst <- numeric(length = nrow(d))
  for(i in 1:nrow(d)){
    if(d$C_I[i] == Iconc[1] & d$C_AP[i] < APconc[5]){ #Substrate added to second row (B) 
      d$TimDiffFirst[i] <- all$TimDiff[24]
    }else{
      if((d$C_I[i] == Iconc[1] & d$C_AP[i] == APconc[5]) | ((d$C_I[i] == Iconc[2] & d$C_AP[i] < APconc[4]))){#Substrate added to third row (C)
        d$TimDiffFirst[i] <- all$TimDiff[25+92]
      }else{
        if((d$C_I[i] == Iconc[2] & d$C_AP[i] >= APconc[4]) | ((d$C_I[i] == Iconc[3] & d$C_AP[i] < APconc[3]))){#Substrate added to fourth row (D)
          d$TimDiffFirst[i] <- all$TimDiff[48+92*2]
        }else{
          if((d$C_I[i] == Iconc[3] & d$C_AP[i] >= APconc[3]) | ((d$C_I[i] == Iconc[4] & d$C_AP[i] == APconc[1]))){#Substrate added to fifth row (E)
            d$TimDiffFirst[i] <- all$TimDiff[49+92*3]
          }else{
            if((d$C_I[i] == Iconc[4] & d$C_AP[i] >= APconc[2])){#Substrate added to sixth row (F)
              d$TimDiffFirst[i] <- all$TimDiff[72+92*4]
            }else{
              if((d$C_I[i] == Iconc[5])){#Substrate added to last two rows (G and H)
                d$TimDiffFirst[i] <- all$TimDiff[73+92*5]
              }
            }
          }
        }
      }
    }
  }
  
  d$ReactionTime0 <- d$TimDiff-d$TimDiffFirst
  d$ReactionTime <- ifelse(d$ReactionTime0 < 0, 0, d$ReactionTime0)
  
  #Use linear regression to make a calibration line for slurry and controls
  #it is calculated separately for each measurement and coefficients of linear regression are added to respective datasets
  ##controls
  contr$Slope<-NA
  contr$Intercept<-NA
  contr$CalRSQ<-NA
  
  for(i in unique(contr$Date)){
    contr[contr$Date == i, "Intercept"] <- as.numeric(coef(lm(E ~ C_MUB, data = calsCntr[calsCntr$Date == i, ]))[1])
    contr[contr$Date == i, "Slope"] <- as.numeric(coef(lm(E ~ C_MUB, data = calsCntr[calsCntr$Date == i, ]))[2])
    contr[contr$Date == i, "CalRSQ"] <- as.numeric(summary(lm(E ~ C_MUB, data = calsCntr[calsCntr$Date == i, ]))[9])
    #summary(lm(E ~ C_MUB, data = calsCntr[calsCntr$Date == unique(contr$Date)[1], ]))[9]
  }
  
  ##Slurries
  d$Slope<-NA
  d$Intercept<-NA
  d$CalRSQ<-NA
  
  for(i in unique(d$Date)){
    d[d$Date == i, "Intercept"] <- as.numeric(coef(lm(E ~ C_MUB, data = calsSl[calsSl$Date == i, ]))[1])
    d[d$Date == i, "Slope"] <- as.numeric(coef(lm(E ~ C_MUB, data = calsSl[calsSl$Date == i, ]))[2])
    d[d$Date == i, "CalRSQ"] <- as.numeric(summary(lm(E ~ C_MUB, data = calsSl[calsSl$Date == i, ]))[9])
  }
  
  #Calculate product concentration in umols of product per liter
  ##Controls
  contr$Pcontr <- with(contr, (E - Intercept)/Slope)
  contr$Pcontr <- ifelse(contr$Pcontr<0, 0, contr$Pcontr)
  contrMean <- as.data.frame(contr %>% group_by(C_AP) %>% summarize(PcontrMean = mean(Pcontr, na.rm = T)))
  #Slurries
  d$P <- with(d, (E - Intercept)/Slope)
  
  #Add product concentrations to respective measurements in dataset d
  d$PcontrMean <- numeric(length = nrow(d)) 
  for(l in APconc){
    d[d$C_AP == l, "PcontrMean"] <- as.numeric(contrMean[contrMean$C_AP == l, "PcontrMean"])
  }
 
  #Recalculate product concentration to umols of product per gram of dry soil and set negative values to 0
  d$Product <- with(d, (P-PcontrMean)/(10/1.25)/DW)
  #d$Product <- with(d, P/(10/1.25)/DW)
  d$Product <- ifelse(d$Product < 0, 0, d$Product)
  
  #Correct initial concentration of substrate and inhibitor for spontaneous decay of substrate in water
  d$SubstrateInitial <- with(d, C_AP - PcontrMean/(10/1.25)/DW)
  d$InhibitorInitial <- with(d, C_I + PcontrMean/(10/1.25)/DW)
  d$SubstrateInitial <- round(d$SubstrateInitial, 2)
  d$InhibitorInitial <- round(ifelse(d$InhibitorInitial<0, 0, d$InhibitorInitial), 2)
  
  all_out<-list(data=d, 
                contr=contr,
                contrMean=contrMean,
                cal_dataCntr=calsCntr,
                cal_dataSl=calsSl
                )
  return(all_out)
}
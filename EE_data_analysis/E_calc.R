E_calc<-function(dataset, MUBconc, APconc, Iconc, Nmeasure, Times, empty){#
  
  #Loop function
  lfunction<-function(sheet){
    
    #Reading data
    ##Careful, sheets are readed from the opposite direction (last measurements represent first sheets) 
    e_0<-read.xlsx(xlsxFile = dataset, sheet=sheet,
                   rows=46:53, cols = 2:13, colNames = F)
    
    #Rearanging the data
    e_0r<-data.frame(Slurry = c(rep("NO", 6), rep("YES", 96-5-6-4), rep("NO", 5)), 
                     Calibration = c(rep("TRUE", times = 12), #define the data set for calibration and measurment
                                     rep("FALSE", times = 80)),
                     C_MUB = c(rep(MUBconc/5, times = 2), #concentration of MUB
                               rep(NA, times = 80)),
                     C_AP = c(rep(NA, times = 12),
                              rep(c(rep(APconc/5, each = 3)), times = 5),#concentrations of added AP substrate
                              APconc/5),
                     C_I = c(rep(NA, times = 12),
                             rep(Iconc/5, each = 15),
                             rep(NA, times=5)),
                     E = c(as.numeric(e_0[1, 1:6]), #Calibration for controls
                           as.numeric(e_0[1, 7:12]), #Calibration for Slurry
                           as.numeric(e_0[2, 1:12]), as.numeric(e_0[3, 1:3]), #I = 0
                           as.numeric(e_0[3, 4:12]), as.numeric(e_0[4, 1:6]), #I1
                           as.numeric(e_0[4, 7:12]), as.numeric(e_0[5, 1:9]), #I2
                           as.numeric(e_0[5, 10:12]), as.numeric(e_0[6, 1:12]), #I3
                           as.numeric(e_0[7, 1:12]), as.numeric(e_0[8, 1:3]), #I4
                           as.numeric(e_0[8, 4:8]) #AP controls
                     ))
    
    return(e_0r)
  }
  
  #Go through the loop
  ID=seq(from=Nmeasure-empty, to=1)
  
  #first measurement
  all<-lfunction(sheet = ID[1])
  
  #all other measurements
  for(i in ID[-1]){
    all<-rbind(all, lfunction(sheet = i))
  }
  
  #divide data set into calibrations, controls and real data
  calsCntr<-subset(all, Calibration == "TRUE" & Slurry=="NO")
  calsSl<-subset(all, Calibration == "TRUE" & Slurry!="NO")
  contr<-subset(all, Calibration != "TRUE" & Slurry=="NO")
  d<-subset(all, Calibration != "TRUE" & Slurry!="NO")
  
  #~~~~~~~~~~Unused argument from previous version~~~~~~~~~~~~~~~~~~~~#
  # #Subtract blanks
  # cals$Ecorr<-numeric(length = nrow(cals))
  # blanks$Ecorr<-numeric(length = nrow(blanks))
  # d$Ecorr<-numeric(length = nrow(d))
  # ##Calibrations
  # for(i in unique(cals$Slurry)){
  #   cals[cals$Slurry==i, "Ecorr"]<-cals[cals$Slurry==i, "E"]-
  #     mean(cals[(cals$Slurry==i & cals$C_MUB==0), "E"], na.rm=T)
  # }
  # ##blanks
  # blanks$Ecorr<-blanks$E-mean(cals[(cals$Slurry==0 & cals$C_MUB==0), "E"], na.rm=T)
  # ##data
  # d$Time<-Times
  # 
  # if(intens=="TRUE"){
  #   for(i in unique(d$Slurry)){
  #     d[d$Slurry==i, "Ecorr"]<-d[d$Slurry==i, "E"]-
  #       mean(c(as.numeric(cals[(cals$Slurry==i & cals$C_MUB==0), "E"]),
  #              as.numeric(d[(d$Slurry==i & d$Time==0), "E"])), na.rm=T)
  #   }
  # }else{
  #   for(i in unique(d$Slurry)){
  #     d[d$Slurry==i, "Ecorr"]<-d[d$Slurry==i, "E"]-
  #       mean(as.numeric(cals[(cals$Slurry==i & cals$C_MUB==0), "E"]), na.rm=T)
  #   }
  # }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  
  #Use linear regression to make a calibration line for slurry and controls
  #it is calculated separately for each measurement
  ##controls
  contr$P<-NA
  ID2<-seq(from=1, to=(Nmeasure-empty)*6, by=6)
  ID3<-seq(from=1, to=(Nmeasure-empty)*5, by=5)
  
  for(i in ID2){
    for(n in ID3){
      #polynomial regression
      cal_contr<-lm(C_MUB~poly(E, 3, raw=T), data = calsCntr[c(i:(i+5)), ])
      
      contr[c(n:(n+4)), "P"]<-as.numeric(predict(cal_contr, newdata=contr[c(n:(n+4)),]))
    }
    
  }
  
  
  ##Slurries
  d$P<-NA
  ID4<-seq(from=1, to=(Nmeasure-empty)*5*3*5, by=5*3*5)
  
  for(i in ID2){
    for(n in ID4){
      #linear regression
      cal_Sl<-lm(C_MUB~poly(E, 3, raw=T), data = calsSl[c(i:(i+5)), ])
      
      d[c(n:(n+74)), "P"]<-as.numeric(predict(cal_Sl, newdata=d[c(n:(n+74)),]))
    }
    
  }
  
  #add controls to the data frame
  d$Pcontr<-NA
  
  for(i in ID3){
    for(n in ID4){
      d[c(n:(n+74)), "Pcontr"]<-rep(rep(as.numeric(contr[c(i:(i+4)), "P"]), each=3), times=5)
    }
    
  }
  
  #add reaction time
  d$time<-Times
  
  all_out<-list(data=d, 
                contr=contr,
                cal_dataCntr=calsCntr,
                cal_dataSl=calsSl
                )
  return(all_out)
}
###############################################################################################################
###############################################################################################################
############################Biochemical inhibition of hydrolytic enzymes in soil###############################
###############################################################################################################
###############################################################################################################

############################################Loading libraries##################################################
library(readODS)
library(openxlsx)
library(ggplot2)
library(minpack.lm)
library(dplyr)
library(lamW)
library(bbmle)
library(FME)
library(DEoptim)
library(ABCoptim)
library(gridExtra)
library(reshape)
library(segmented)
library(reticulate)
library(set)
#############################################GGPLOT THEME######################################################
theme_min <- readRDS("/mnt/580CBE2464C5F83D/pracovni/helpfull_R_Matlab_script/ggtheme.rds")
#################################Raw data recalculation functions##############################################
source("./E_calcMUB.R")
source("./E_calcAMCL.R")
#========================Plesne catchment - litter layer
#=============Beta-Glucosidase
#1. Without the benzoic acid
#Weighted 49.488 mg of MUB-G to 100 ml and diluted
#DW is 0.34234
PLO_MUBG_0 <- E_calcMUB(dataset = "../HydrolyticRaw/Soil1MUF-G21112022.xlsx",
                    MUBconc = c(1, 5, 10, 25, 50, 0),
                    APconc = c(106.83/100, 106.83/20, 106.83/4, 106.83/2, 106.83),#106.83/100, 106.83/20, 106.83/4, 106.83/2
                    Iconc = c(0, 1, 5, 10, 20),
                    Nmeasure = 20, DW = 0.34234, empty = 0)
#Exporting re-calculated data
Edata <- PLO_MUBG_0$data
Edata$Catchment <- c("Plesne")
Edata$Horizon <- c("Litter layer")
Edata$Enzyme <- c("Beta-Glucosidase")
Edata$BenzoicAcid <- c("False")
#Rounding the initial Substrate and Inhibitor concentrations for nicer plotting
Edata$SubstrateInitialR <- round(Edata$SubstrateInitial, 0)
Edata$InhibitorInitialR <- round(Edata$InhibitorInitial, 0)

#Plotting the progress curves
Edata %>% mutate(Tbins = cut(ReactionTime, 
                             breaks = seq(0, max(Edata$ReactionTime), by = 0.1))) %>% 
  group_by(Tbins, Catchment, Horizon, Enzyme, BenzoicAcid, SubstrateInitialR, InhibitorInitialR) %>% 
  filter(Catchment == "Plesne" & Horizon == "Litter layer" & Enzyme == "Beta-Glucosidase" & BenzoicAcid == "False") %>% 
  summarise(x = mean(ReactionTime/60), y = mean(Product), ySD = sd(Product)) %>% 
  ggplot(aes(x, y)) + geom_point(pch = 21, cex = 6, aes(fill = factor(SubstrateInitialR))) + 
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) +
  facet_wrap(~InhibitorInitialR, scales = "free") + theme_min + xlab("Time (hours)") +
  ylab(expression(paste("Reaction Product (", mu, "mol (MUB) ", g~(DW)^{-1}, ")"))) +
  labs(fill=expression(paste("Substrate\n",mu, "mol ",g~(DW)^{-1}))) +
  theme(legend.position = c(0.85, 0.3))

#Correct outliers
Edata[(Edata$Catchment == "Plesne" & Edata$Horizon == "Litter layer" &
         Edata$Enzyme == "Beta-Glucosidase" & Edata$BenzoicAcid == "False" &
         Edata$SubstrateInitialR == 1 & Edata$InhibitorInitialR == 5 & 
         Edata$ReactionTime < 60), "Product"] <- c(0)

#2. With the benzoic acid
#Weighted 49.488 mg of MUB-G to 100 ml and diluted
#DW is 0.34234
PLO_MUBG_B <- E_calcMUB(dataset = "../HydrolyticRaw/Soil1MUF-G21112022-PLUS.xlsx",
                        MUBconc = c(1, 5, 10, 25, 50, 0),
                        APconc = c(106.83/100, 106.83/20, 106.83/4, 106.83/2, 106.83),#106.83/100, 106.83/20, 106.83/4, 106.83/2
                        Iconc = c(0, 1, 5, 10, 20),
                        Nmeasure = 19, DW = 0.34234, empty = 0)
#Exporting re-calculated data
PLO_MUBG_B$data$Catchment <- c("Plesne")
PLO_MUBG_B$data$Horizon <- c("Litter layer")
PLO_MUBG_B$data$Enzyme <- c("Beta-Glucosidase")
PLO_MUBG_B$data$BenzoicAcid <- c("True")
#Rounding the initial Substrate and Inhibitor concentrations for nicer plotting
PLO_MUBG_B$data$SubstrateInitialR <- round(PLO_MUBG_B$data$SubstrateInitial, 0)
PLO_MUBG_B$data$InhibitorInitialR <- round(PLO_MUBG_B$data$InhibitorInitial, 0)
Edata <- rbind(Edata, PLO_MUBG_B$data)

#Plotting the progress curves
Edata %>% mutate(Tbins = cut(ReactionTime, 
                             breaks = seq(0, max(Edata$ReactionTime), by = 0.1))) %>% 
  group_by(Tbins, Catchment, Horizon, Enzyme, BenzoicAcid, SubstrateInitialR, InhibitorInitialR) %>% 
  summarise(x = mean(ReactionTime/60), y = mean(Product), ySD = sd(Product)) %>% 
  ggplot(aes(x, y)) + geom_point(pch = 21, cex = 6, aes(fill = BenzoicAcid)) + 
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) +
  facet_wrap(SubstrateInitialR~InhibitorInitialR, scales = "free") + theme_min + xlab("Time (hours)") +
  ylab(expression(paste("Reaction Product (", mu, "mol (MUB) ", g~(DW)^{-1}, ")"))) +
  labs(fill=expression(paste("Substrate\n",mu, "mol ",g~(DW)^{-1}))) +
  theme(legend.position = c("top"), legend.direction = "horizontal")

#=============Leucine-Aminopeptidase
#1. Without the benzoic acid
#Weighted 44.075 mg of AMCL to 100 ml and diluted
#DW is 0.34234
PLO_AMCL_0 <- E_calcAMCL(dataset = "../HydrolyticRaw/Soil1AMC-Leu02122022.xlsx",
                        MUBconc = c(0, 1, 5, 10, 25, 50),
                        APconc = c(99.096/100, 99.096/20, 99.096/4, 99.096/2, 99.096),
                        Iconc = c(0, 1, 5, 10, 20),
                        Nmeasure = 15, DW = 0.34234, empty = 0)

#Exporting re-calculated data
PLO_AMCL_0$data$Catchment <- c("Plesne")
PLO_AMCL_0$data$Horizon <- c("Litter layer")
PLO_AMCL_0$data$Enzyme <- c("Leu-Aminopeptidase")
PLO_AMCL_0$data$BenzoicAcid <- c("False")
#Rounding the initial Substrate and Inhibitor concentrations for nicer plotting
PLO_AMCL_0$data$SubstrateInitialR <- round(PLO_AMCL_0$data$SubstrateInitial, 0)
PLO_AMCL_0$data$InhibitorInitialR <- round(PLO_AMCL_0$data$InhibitorInitial, 0)
Edata <- rbind(Edata, PLO_AMCL_0$data)

#Plotting the progress curves
Edata %>% mutate(Tbins = cut(ReactionTime, 
                             breaks = seq(0, max(Edata$ReactionTime), by = 0.1))) %>% 
  group_by(Tbins, Catchment, Horizon, Enzyme, BenzoicAcid, SubstrateInitialR, InhibitorInitialR) %>% 
  filter(Catchment == "Plesne" & Horizon == "Litter layer" & Enzyme == "Leu-Aminopeptidase" & BenzoicAcid == "False") %>% 
  summarise(x = mean(ReactionTime), y = mean(Product), ySD = sd(Product)) %>% 
  ggplot(aes(x, y)) + geom_point(pch = 21, cex = 6, aes(fill = factor(SubstrateInitialR))) + 
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) +
  facet_wrap(~InhibitorInitialR) + theme_min + xlab("Time (hours)") +
  ylab(expression(paste("Reaction Product (", mu, "mol (MUB) ", g~(DW)^{-1}, ")"))) +
  labs(fill=expression(paste("Substrate\n",mu, "mol ",g~(DW)^{-1}))) +
  theme(legend.position = c(0.85, 0.3))

#Correct outliers
Edata[(Edata$Catchment == "Plesne" & Edata$Horizon == "Litter layer" &
        Edata$Enzyme == "Leu-Aminopeptidase" & Edata$BenzoicAcid == "False" &
        #Edata$SubstrateInitialR == 1 & Edata$InhibitorInitialR == 5 & 
        Edata$ReactionTime == 0), "Product"] <- c(0)
Edata[(Edata$Catchment == "Plesne" & Edata$Horizon == "Litter layer" &
         Edata$Enzyme == "Leu-Aminopeptidase" & Edata$BenzoicAcid == "False" &
         Edata$InhibitorInitialR == 0 & 
         Edata$ReactionTime < 10), "Product"] <- c(0)
Edata[(Edata$Catchment == "Plesne" & Edata$Horizon == "Litter layer" &
         Edata$Enzyme == "Leu-Aminopeptidase" & Edata$BenzoicAcid == "False" &
         Edata$InhibitorInitialR == 1 & 
         Edata$ReactionTime < 10), "Product"] <- c(0)

#2. With the benzoic acid
#Weighted 44.075 mg of AMCL to 100 ml and diluted
#DW is 0.34234
PLO_AMCL_B <- E_calcAMCL(dataset = "../HydrolyticRaw/Soil1AMC-Leu02122022-PLUS.xlsx",
                         MUBconc = c(0, 1, 5, 10, 25, 50),
                         APconc = c(99.096/100, 99.096/20, 99.096/4, 99.096/2, 99.096),
                         Iconc = c(0, 1, 5, 10, 20),
                         Nmeasure = 14, DW = 0.34234, empty = 0)

#Exporting re-calculated data
PLO_AMCL_B$data$Catchment <- c("Plesne")
PLO_AMCL_B$data$Horizon <- c("Litter layer")
PLO_AMCL_B$data$Enzyme <- c("Leu-Aminopeptidase")
PLO_AMCL_B$data$BenzoicAcid <- c("True")
#Rounding the initial Substrate and Inhibitor concentrations for nicer plotting
PLO_AMCL_B$data$SubstrateInitialR <- round(PLO_AMCL_B$data$SubstrateInitial, 0)
PLO_AMCL_B$data$InhibitorInitialR <- round(PLO_AMCL_B$data$InhibitorInitial, 0)
Edata <- rbind(Edata, PLO_AMCL_B$data)

#Plotting the progress curves
Edata %>% mutate(Tbins = cut(ReactionTime, 
                             breaks = seq(0, max(Edata$ReactionTime), by = 0.1))) %>% 
  group_by(Tbins, Catchment, Horizon, Enzyme, BenzoicAcid, SubstrateInitialR, InhibitorInitialR) %>% 
  filter(Catchment == "Plesne" & Horizon == "Litter layer" & Enzyme == "Leu-Aminopeptidase" & BenzoicAcid == "True") %>% 
  summarise(x = mean(ReactionTime), y = mean(Product), ySD = sd(Product)) %>% 
  ggplot(aes(x, y)) + geom_point(pch = 21, cex = 6, aes(fill = factor(InhibitorInitialR))) + 
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) +
  facet_wrap(~SubstrateInitialR) + theme_min + xlab("Time (hours)") +
  ylab(expression(paste("Reaction Product (", mu, "mol (MUB) ", g~(DW)^{-1}, ")"))) +
  labs(fill=expression(paste("Substrate\n",mu, "mol ",g~(DW)^{-1}))) +
  theme(legend.position = c(0.85, 0.3))

#========================Plesne catchment - topsoil organic horizon
#=============Beta-Glucosidase
#1. Without the benzoic acid
#Weighted 45.309 mg of MUB-G to 100 ml and diluted
#DW is 0.3423
PLA_MUBG_0 <- E_calcMUB(dataset = "../HydrolyticRaw/Soil1-A-MUF-G03032023.xlsx",
                        MUBconc = c(0, 0, 0, 1, 5, 10),
                        APconc = c(97.805/100, 97.805/20, 97.805/4, 97.805/2, 97.805),
                        Iconc = c(0, 1, 5, 10, 20),
                        Nmeasure = 15, DW = 0.34234, empty = 0)
#Exporting re-calculated data
PLA_MUBG_0$data$Catchment <- c("Plesne")
PLA_MUBG_0$data$Horizon <- c("Organic horizon")
PLA_MUBG_0$data$Enzyme <- c("Beta-Glucosidase")
PLA_MUBG_0$data$BenzoicAcid <- c("False")
#Rounding the initial Substrate and Inhibitor concentrations for nicer PLAtting
PLA_MUBG_0$data$SubstrateInitialR <- round(PLA_MUBG_0$data$SubstrateInitial, 0)
PLA_MUBG_0$data$InhibitorInitialR <- round(PLA_MUBG_0$data$InhibitorInitial, 0)
Edata <- rbind(Edata, PLA_MUBG_0$data)

#Plotting the progress curves
Edata %>% mutate(Tbins = cut(ReactionTime, 
                             breaks = seq(0, max(Edata$ReactionTime), by = 0.1))) %>% 
  group_by(Tbins, Catchment, Horizon, Enzyme, BenzoicAcid, SubstrateInitialR, InhibitorInitialR) %>% 
  filter(Horizon == "Organic horizon" & Enzyme == "Beta-Glucosidase" & BenzoicAcid == "False") %>% 
  summarise(x = mean(ReactionTime/60), y = mean(Product), ySD = sd(Product)) %>% 
  ggplot(aes(x, y)) + geom_point(pch = 21, cex = 6, aes(fill = factor(SubstrateInitialR))) + 
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) +
  facet_wrap(~InhibitorInitialR, scales = "free") + theme_min + xlab("Time (hours)") +
  ylab(expression(paste("Reaction Product (", mu, "mol (MUB) ", g~(DW)^{-1}, ")"))) +
  labs(fill=expression(paste("Substrate\n",mu, "mol ",g~(DW)^{-1}))) +
  theme(legend.position = c(0.85, 0.3))

#2. With the benzoic acid
#Weighted 45.309 mg of MUB-G to 100 ml and diluted
#DW is 0.3423
PLA_MUBG_B <- E_calcMUB(dataset = "../HydrolyticRaw/Soil1-A-MUF-G03032023-PLUS.xlsx",
                        MUBconc = c(0, 0, 0, 1, 5, 10),
                        APconc = c(97.805/100, 97.805/20, 97.805/4, 97.805/2, 97.805),
                        Iconc = c(0, 1, 5, 10, 20),
                        Nmeasure = 14, DW = 0.34234, empty = 0)
#Exporting re-calculated data
PLA_MUBG_B$data$Catchment <- c("Plesne")
PLA_MUBG_B$data$Horizon <- c("Organic horizon")
PLA_MUBG_B$data$Enzyme <- c("Beta-Glucosidase")
PLA_MUBG_B$data$BenzoicAcid <- c("True")
#Rounding the initial Substrate and Inhibitor concentrations for nicer PLAtting
PLA_MUBG_B$data$SubstrateInitialR <- round(PLA_MUBG_B$data$SubstrateInitial, 0)
PLA_MUBG_B$data$InhibitorInitialR <- round(PLA_MUBG_B$data$InhibitorInitial, 0)
Edata <- rbind(Edata, PLA_MUBG_B$data)

#plotting the progress curves
Edata %>% mutate(Tbins = cut(ReactionTime, 
                             breaks = seq(0, max(Edata$ReactionTime), by = 0.1))) %>% 
  group_by(Tbins, Catchment, Horizon, Enzyme, BenzoicAcid, SubstrateInitialR, InhibitorInitialR) %>% 
  summarise(x = mean(ReactionTime/60), y = mean(Product), ySD = sd(Product)) %>% 
  filter(Horizon == "Organic horizon" & Enzyme == "Beta-Glucosidase" & BenzoicAcid == "True") %>% 
  ggplot(aes(x, y)) + geom_point(pch = 21, cex = 6, aes(fill = factor(SubstrateInitialR))) + 
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) +
  facet_wrap(~InhibitorInitialR, scales = "free") + theme_min + xlab("Time (hours)") +
  ylab(expression(paste("Reaction Product (", mu, "mol (MUB) ", g~(DW)^{-1}, ")"))) +
  labs(fill=expression(paste("Substrate\n",mu, "mol ",g~(DW)^{-1}))) +
  theme(legend.position = c("top"), legend.direction = "horizontal")

#=============Leucine-Aminopeptidase
#1. Without the benzoic acid
#Weighted 45.309 mg of AMCL to 100 ml and diluted
#DW is 0.34234
PLA_AMCL_0 <- E_calcAMCL(dataset = "../HydrolyticRaw/Soil1-A-AMC-LEU24032023.xlsx",
                         MUBconc = c(0, 0.25, 0.5, 1, 2, 5),
                         APconc = c(101.87/100, 101.87/20, 101.87/4, 101.87/2, 101.87),
                         Iconc = c(0, 1, 5, 10, 20),
                         Nmeasure = 14, DW = 0.34234, empty = 0)

#Exporting re-calculated data
PLA_AMCL_0$data$Catchment <- c("Plesne")
PLA_AMCL_0$data$Horizon <- c("Organic horizon")
PLA_AMCL_0$data$Enzyme <- c("Leu-Aminopeptidase")
PLA_AMCL_0$data$BenzoicAcid <- c("False")
#Rounding the initial Substrate and Inhibitor concentrations for nicer PLAtting
PLA_AMCL_0$data$SubstrateInitialR <- round(PLA_AMCL_0$data$SubstrateInitial, 0)
PLA_AMCL_0$data$InhibitorInitialR <- round(PLA_AMCL_0$data$InhibitorInitial, 0)
Edata <- rbind(Edata, PLA_AMCL_0$data)

#Plotting the progress curves
Edata %>% mutate(Tbins = cut(ReactionTime, 
                             breaks = seq(0, max(Edata$ReactionTime), by = 0.1))) %>% 
  group_by(Tbins, Catchment, Horizon, Enzyme, BenzoicAcid, SubstrateInitialR, InhibitorInitialR) %>% 
  filter(Horizon == "Organic horizon" & Enzyme == "Leu-Aminopeptidase" & BenzoicAcid == "False") %>% 
  summarise(x = mean(ReactionTime), y = mean(Product), ySD = sd(Product)) %>% 
  ggplot(aes(x, y)) + geom_point(pch = 21, cex = 6, aes(fill = factor(SubstrateInitialR))) + 
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) +
  facet_wrap(~InhibitorInitialR) + theme_min + xlab("Time (hours)") +
  ylab(expression(paste("Reaction Product (", mu, "mol (MUB) ", g~(DW)^{-1}, ")"))) +
  labs(fill=expression(paste("Substrate\n",mu, "mol ",g~(DW)^{-1}))) +
  theme(legend.position = c(0.85, 0.3))

#2. With the benzoic acid
#Weighted 45.309 mg of AMCL to 100 ml and diluted
#DW is 0.34234
PLA_AMCL_B <- E_calcAMCL(dataset = "../HydrolyticRaw/Soil1-A-AMC-LEU24032023-PLUS.xlsx",
                         MUBconc = c(0, 0.25, 0.5, 1, 2, 5),
                         APconc = c(101.87/100, 101.87/20, 101.87/4, 101.87/2, 101.87),
                         Iconc = c(0, 1, 5, 10, 20),
                         Nmeasure = 14, DW = 0.34234, empty = 0)

#Exporting re-calculated data
PLA_AMCL_B$data$Catchment <- c("Plesne")
PLA_AMCL_B$data$Horizon <- c("Organic horizon")
PLA_AMCL_B$data$Enzyme <- c("Leu-Aminopeptidase")
PLA_AMCL_B$data$BenzoicAcid <- c("True")
#Rounding the initial Substrate and Inhibitor concentrations for nicer PLAtting
PLA_AMCL_B$data$SubstrateInitialR <- round(PLA_AMCL_B$data$SubstrateInitial, 0)
PLA_AMCL_B$data$InhibitorInitialR <- round(PLA_AMCL_B$data$InhibitorInitial, 0)
Edata <- rbind(Edata, PLA_AMCL_B$data)

#Plotting the progress curves
Edata %>% mutate(Tbins = cut(ReactionTime, 
                             breaks = seq(0, max(Edata$ReactionTime), by = 0.1))) %>% 
  group_by(Tbins, Catchment, Horizon, Enzyme, BenzoicAcid, SubstrateInitialR, InhibitorInitialR) %>% 
  filter(Horizon == "Organic horizon" & Enzyme == "Leu-Aminopeptidase" & BenzoicAcid == "True") %>% 
  summarise(x = mean(ReactionTime), y = mean(Product), ySD = sd(Product)) %>% 
  ggplot(aes(x, y)) + geom_point(pch = 21, cex = 6, aes(fill = factor(SubstrateInitialR))) + 
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) +
  facet_wrap(~InhibitorInitialR) + theme_min + xlab("Time (hours)") +
  ylab(expression(paste("Reaction Product (", mu, "mol (MUB) ", g~(DW)^{-1}, ")"))) +
  labs(fill=expression(paste("Substrate\n",mu, "mol ",g~(DW)^{-1}))) +
  theme(legend.position = c(0.85, 0.3))

#========================Certovo catchment - litter layer
#=============Beta-Glucosidase
#1. Without the benzoic acid
#Weighted 49.488 mg of MUB-G to 100 ml and diluted
#DW is 0.21875
CTO_MUBG_0 <- E_calcMUB(dataset = "../HydrolyticRaw/Soil2MUF-G25112022.xlsx",
                        MUBconc = c(0, 1, 5, 10, 25, 50),
                        APconc = c(167.18/100, 167.18/20, 167.18/4, 167.18/2, 167.18),
                        Iconc = c(0, 31.29/20, 31.29/4, 31.29/2, 31.29),
                        Nmeasure = 16, DW = 0.21875, empty = 0)
#Exporting re-calculated data
CTO_MUBG_0$data$Catchment <- c("Certovo")
CTO_MUBG_0$data$Horizon <- c("Litter layer")
CTO_MUBG_0$data$Enzyme <- c("Beta-Glucosidase")
CTO_MUBG_0$data$BenzoicAcid <- c("False")
#Rounding the initial Substrate and Inhibitor concentrations for nicer plotting
CTO_MUBG_0$data$SubstrateInitialR <- round(CTO_MUBG_0$data$SubstrateInitial, 0)
CTO_MUBG_0$data$InhibitorInitialR <- round(CTO_MUBG_0$data$InhibitorInitial, 0)
Edata <- rbind(Edata, CTO_MUBG_0$data)

#Plotting the progress curves
Edata %>% mutate(Tbins = cut(ReactionTime, 
                             breaks = seq(0, max(Edata$ReactionTime), by = 0.1))) %>% 
  group_by(Tbins, Catchment, Horizon, Enzyme, BenzoicAcid, SubstrateInitialR, InhibitorInitialR) %>% 
  filter(Catchment == "Certovo" & Horizon == "Litter layer" & Enzyme == "Beta-Glucosidase" & BenzoicAcid == "False") %>% 
  summarise(x = mean(ReactionTime/60), y = mean(Product), ySD = sd(Product)) %>% 
  ggplot(aes(x, y)) + geom_point(pch = 21, cex = 6, aes(fill = factor(SubstrateInitialR))) + 
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) +
  facet_wrap(~InhibitorInitialR) + theme_min + xlab("Time (hours)") +
  ylab(expression(paste("Reaction Product (", mu, "mol (MUB) ", g~(DW)^{-1}, ")"))) +
  labs(fill=expression(paste("Substrate\n",mu, "mol ",g~(DW)^{-1}))) +
  theme(legend.position = c(0.85, 0.3))

#2. With the benzoic acid
#Weighted 49.488 mg of MUB-G to 100 ml and diluted
#DW is 0.21875
CTO_MUBG_B <- E_calcMUB(dataset = "../HydrolyticRaw/Soil2MUF-G25112022-PLUS.xlsx",
                        MUBconc = c(0, 1, 5, 10, 25, 50),
                        APconc = c(167.18/100, 167.18/20, 167.18/4, 167.18/2, 167.18),
                        Iconc = c(0, 31.29/20, 31.29/4, 31.29/2, 31.29),
                        Nmeasure = 15, DW = 0.21875, empty = 0)
#Exporting re-calculated data
CTO_MUBG_B$data$Catchment <- c("Certovo")
CTO_MUBG_B$data$Horizon <- c("Litter layer")
CTO_MUBG_B$data$Enzyme <- c("Beta-Glucosidase")
CTO_MUBG_B$data$BenzoicAcid <- c("True")
#Rounding the initial Substrate and Inhibitor concentrations for nicer plotting
CTO_MUBG_B$data$SubstrateInitialR <- round(CTO_MUBG_B$data$SubstrateInitial, 0)
CTO_MUBG_B$data$InhibitorInitialR <- round(CTO_MUBG_B$data$InhibitorInitial, 0)
Edata <- rbind(Edata, CTO_MUBG_B$data)

#Plotting the progress curves
Edata %>% mutate(Tbins = cut(ReactionTime, 
                             breaks = seq(0, max(Edata$ReactionTime), by = 0.1))) %>% 
  group_by(Tbins, Catchment, Horizon, Enzyme, BenzoicAcid, SubstrateInitialR, InhibitorInitialR) %>% 
  filter(Catchment == "Certovo" & Horizon == "Litter layer" & Enzyme == "Beta-Glucosidase" & BenzoicAcid == "True") %>% 
  summarise(x = mean(ReactionTime/60), y = mean(Product), ySD = sd(Product)) %>% 
  ggplot(aes(x, y)) + geom_point(pch = 21, cex = 6, aes(fill = factor(SubstrateInitialR))) + 
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) +
  facet_wrap(~InhibitorInitialR) + theme_min + xlab("Time (hours)") +
  ylab(expression(paste("Reaction Product (", mu, "mol (MUB) ", g~(DW)^{-1}, ")"))) +
  labs(fill=expression(paste("Substrate\n",mu, "mol ",g~(DW)^{-1}))) +
  theme(legend.position = c("top"), legend.direction = "horizontal")

#=============Leucine-Aminopeptidase
#1. Without the benzoic acid
#Weighted 29.072 mg of AMCL to 100 ml and diluted
#DW is 0.21875
CTO_AMCL_0 <- E_calcAMCL(dataset = "../HydrolyticRaw/CT25OAMCL.xlsx",
                         MUBconc = c(0, 0.25, 0.5, 1, 2, 5),
                         APconc = c(102.29/100, 102.29/20, 102.29/4, 102.29/2, 102.29),
                         Iconc = c(0, 20/20, 20/4, 20/2, 20),
                         Nmeasure = 23, DW = 0.21875, empty = 0)

#Exporting re-calculated data
CTO_AMCL_0$data$Catchment <- c("Certovo")
CTO_AMCL_0$data$Horizon <- c("Litter layer")
CTO_AMCL_0$data$Enzyme <- c("Leu-Aminopeptidase")
CTO_AMCL_0$data$BenzoicAcid <- c("False")
#Rounding the initial Substrate and Inhibitor concentrations for nicer plotting
CTO_AMCL_0$data$SubstrateInitialR <- round(CTO_AMCL_0$data$SubstrateInitial, 0)
CTO_AMCL_0$data$InhibitorInitialR <- round(CTO_AMCL_0$data$InhibitorInitial, 0)
Edata <- rbind(Edata, CTO_AMCL_0$data)

#Plotting the progress curves
Edata %>% mutate(Tbins = cut(ReactionTime, 
                             breaks = seq(0, max(Edata$ReactionTime), by = 0.1))) %>% 
  group_by(Tbins, Catchment, Horizon, Enzyme, BenzoicAcid, SubstrateInitialR, InhibitorInitialR) %>% 
  filter(Catchment == "Certovo" & Horizon == "Litter layer" & Enzyme == "Leu-Aminopeptidase" & BenzoicAcid == "False") %>%
  summarise(x = mean(ReactionTime/60), y = mean(Product), ySD = sd(Product)) %>% 
  ggplot(aes(x, y)) + geom_point(pch = 21, cex = 6, aes(fill = factor(SubstrateInitialR))) + 
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) +
  facet_wrap(~InhibitorInitialR) + theme_min + xlab("Time (hours)") +
  ylab(expression(paste("Reaction Product (", mu, "mol (MUB) ", g~(DW)^{-1}, ")"))) +
  labs(fill=expression(paste("Substrate\n",mu, "mol ",g~(DW)^{-1}))) +
  theme(legend.position = c(0.85, 0.3))

#2. With the benzoic acid
#Weighted 29.072 mg of AMCL to 100 ml and diluted
#DW is 0.21875
CTO_AMCL_B <- E_calcAMCL(dataset = "../HydrolyticRaw/CT25OAMCLBA.xlsx",
                         MUBconc = c(0, 0.25, 0.5, 1, 2, 5),
                         APconc = c(102.29/100, 102.29/20, 102.29/4, 102.29/2, 102.29),
                         Iconc = c(0, 20/20, 20/4, 20/2, 20),
                         Nmeasure = 21, DW = 0.21875, empty = 0)

#Exporting re-calculated data
CTO_AMCL_B$data$Catchment <- c("Certovo")
CTO_AMCL_B$data$Horizon <- c("Litter layer")
CTO_AMCL_B$data$Enzyme <- c("Leu-Aminopeptidase")
CTO_AMCL_B$data$BenzoicAcid <- c("True")
#Rounding the initial Substrate and Inhibitor concentrations for nicer plotting
CTO_AMCL_B$data$SubstrateInitialR <- round(CTO_AMCL_B$data$SubstrateInitial, 0)
CTO_AMCL_B$data$InhibitorInitialR <- round(CTO_AMCL_B$data$InhibitorInitial, 0)
Edata <- rbind(Edata, CTO_AMCL_B$data)

#Plotting the progress curves
Edata %>% mutate(Tbins = cut(ReactionTime, 
                             breaks = seq(0, max(Edata$ReactionTime), by = 0.1))) %>% 
  group_by(Tbins, Catchment, Horizon, Enzyme, BenzoicAcid, SubstrateInitialR, InhibitorInitialR) %>% 
  filter(Catchment == "Certovo" & Horizon == "Litter layer" & Enzyme == "Leu-Aminopeptidase" & BenzoicAcid == "True") %>%
  summarise(x = mean(ReactionTime/60), y = mean(Product), ySD = sd(Product)) %>% 
  ggplot(aes(x, y)) + geom_point(pch = 21, cex = 6, aes(fill = factor(SubstrateInitialR))) + 
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) +
  facet_wrap(~InhibitorInitialR) + theme_min + xlab("Time (hours)") +
  ylab(expression(paste("Reaction Product (", mu, "mol (MUB) ", g~(DW)^{-1}, ")"))) +
  labs(fill=expression(paste("Substrate\n",mu, "mol ",g~(DW)^{-1}))) +
  theme(legend.position = c(0.85, 0.3))

#========================Certovo catchment - topsoil organic horizon
#=============Beta-Glucosidase
#1. Without the benzoic acid
#Weighted 45.309 mg of MUB-G to 100 ml and diluted
#DW is 0.2269
CTA_MUBG_0 <- E_calcMUB(dataset = "../HydrolyticRaw/Soil2-A-MUF-G10032023.xlsx",
                        MUBconc = c(0, 1, 5, 10, 25, 50),
                        APconc = c(147.57/100, 147.57/20, 147.57/4, 147.57/2, 147.57),
                        Iconc = c(0, 30.16/20, 30.16/4, 30.16/2, 30.16),
                        Nmeasure = 14, DW = 0.2269, empty = 0)
#Exporting re-calculated data
CTA_MUBG_0$data$Catchment <- c("Certovo")
CTA_MUBG_0$data$Horizon <- c("Organic horizon")
CTA_MUBG_0$data$Enzyme <- c("Beta-Glucosidase")
CTA_MUBG_0$data$BenzoicAcid <- c("False")
#Rounding the initial Substrate and Inhibitor concentrations for nicer PLAtting
CTA_MUBG_0$data$SubstrateInitialR <- round(CTA_MUBG_0$data$SubstrateInitial, 0)
CTA_MUBG_0$data$InhibitorInitialR <- round(CTA_MUBG_0$data$InhibitorInitial, 0)
Edata <- rbind(Edata, CTA_MUBG_0$data)

#Plotting the progress curves
Edata %>% mutate(Tbins = cut(ReactionTime, 
                             breaks = seq(0, max(Edata$ReactionTime), by = 0.1))) %>% 
  group_by(Tbins, Catchment, Horizon, Enzyme, BenzoicAcid, SubstrateInitialR, InhibitorInitialR) %>% 
  filter(Catchment == "Certovo" & Horizon == "Organic horizon" & Enzyme == "Beta-Glucosidase" & BenzoicAcid == "False") %>% 
  summarise(x = mean(ReactionTime/60), y = mean(Product), ySD = sd(Product)) %>% 
  ggplot(aes(x, y)) + geom_point(pch = 21, cex = 6, aes(fill = factor(SubstrateInitialR))) + 
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) +
  facet_wrap(~InhibitorInitialR, scales = "free") + theme_min + xlab("Time (hours)") +
  ylab(expression(paste("Reaction Product (", mu, "mol (MUB) ", g~(DW)^{-1}, ")"))) +
  labs(fill=expression(paste("Substrate\n",mu, "mol ",g~(DW)^{-1}))) +
  theme(legend.position = c(0.85, 0.3))

#2. With the benzoic acid
#Weighted 45.309 mg of MUB-G to 100 ml and diluted
#DW is 0.2269
CTA_MUBG_B <- E_calcMUB(dataset = "../HydrolyticRaw/Soil2-A-MUF-G10032023-PLUS.xlsx",
                        MUBconc = c(0, 1, 5, 10, 25, 50),
                        APconc = c(147.57/100, 147.57/20, 147.57/4, 147.57/2, 147.57),
                        Iconc = c(0, 30.16/20, 30.16/4, 30.16/2, 30.16),
                        Nmeasure = 14, DW = 0.2269, empty = 0)
#Exporting re-calculated data
CTA_MUBG_B$data$Catchment <- c("Certovo")
CTA_MUBG_B$data$Horizon <- c("Organic horizon")
CTA_MUBG_B$data$Enzyme <- c("Beta-Glucosidase")
CTA_MUBG_B$data$BenzoicAcid <- c("True")
#Rounding the initial Substrate and Inhibitor concentrations for nicer PLAtting
CTA_MUBG_B$data$SubstrateInitialR <- round(CTA_MUBG_B$data$SubstrateInitial, 0)
CTA_MUBG_B$data$InhibitorInitialR <- round(CTA_MUBG_B$data$InhibitorInitial, 0)
Edata <- rbind(Edata, CTA_MUBG_B$data)

#plotting the progress curves
Edata %>% mutate(Tbins = cut(ReactionTime, 
                             breaks = seq(0, max(Edata$ReactionTime), by = 0.1))) %>% 
  group_by(Tbins, Catchment, Horizon, Enzyme, BenzoicAcid, SubstrateInitialR, InhibitorInitialR) %>% 
  summarise(x = mean(ReactionTime/60), y = mean(Product), ySD = sd(Product)) %>% 
  filter(Catchment == "Certovo" & Horizon == "Organic horizon" & Enzyme == "Beta-Glucosidase" & BenzoicAcid == "True") %>% 
  ggplot(aes(x, y)) + geom_point(pch = 21, cex = 6, aes(fill = factor(SubstrateInitialR))) + 
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) +
  facet_wrap(~InhibitorInitialR, scales = "free") + theme_min + xlab("Time (hours)") +
  ylab(expression(paste("Reaction Product (", mu, "mol (MUB) ", g~(DW)^{-1}, ")"))) +
  labs(fill=expression(paste("Substrate\n",mu, "mol ",g~(DW)^{-1}))) +
  theme(legend.position = c("top"), legend.direction = "horizontal")

#Correct unmeasured combinations
Edata[(Edata$Catchment == "Certovo" & Edata$Horizon == "Organic horizon" &
         Edata$Enzyme == "Beta-Glucosidase" &
         Edata$InhibitorInitialR > 7 & Edata$SubstrateInitialR == 142), "Product"] <- c(NA)

#=============Leucine-Aminopeptidase
#1. Without the benzoic acid
#Weighted 45.309 mg of AMCL to 100 ml and diluted
#DW is 0.2269
CTA_AMCL_0 <- E_calcAMCL(dataset = "../HydrolyticRaw/Soil2-A-AMC-LEU28032023.xlsx",
                         MUBconc = c(0, 0.25, 0.5, 1, 2, 5),
                         APconc = c(153.7/100, 153.7/20, 153.7/4, 153.7/2, 153.7),
                         Iconc = c(0, 30.18/20, 30.18/4, 30.18/2, 30.18),
                         Nmeasure = 14, DW = 0.2269, empty = 0)

#Exporting re-calculated data
CTA_AMCL_0$data$Catchment <- c("Certovo")
CTA_AMCL_0$data$Horizon <- c("Organic horizon")
CTA_AMCL_0$data$Enzyme <- c("Leu-Aminopeptidase")
CTA_AMCL_0$data$BenzoicAcid <- c("False")
#Rounding the initial Substrate and Inhibitor concentrations for nicer PLAtting
CTA_AMCL_0$data$SubstrateInitialR <- round(CTA_AMCL_0$data$SubstrateInitial, 0)
CTA_AMCL_0$data$InhibitorInitialR <- round(CTA_AMCL_0$data$InhibitorInitial, 0)
Edata <- rbind(Edata, CTA_AMCL_0$data)

#Plotting the progress curves
Edata %>% mutate(Tbins = cut(ReactionTime, 
                             breaks = seq(0, max(Edata$ReactionTime), by = 0.1))) %>% 
  group_by(Tbins, Catchment, Horizon, Enzyme, BenzoicAcid, SubstrateInitialR, InhibitorInitialR) %>% 
  filter(Catchment == "Certovo" & Horizon == "Organic horizon" & Enzyme == "Leu-Aminopeptidase" & BenzoicAcid == "False") %>% 
  summarise(x = mean(ReactionTime), y = mean(Product), ySD = sd(Product)) %>% 
  ggplot(aes(x, y)) + geom_point(pch = 21, cex = 6, aes(fill = factor(SubstrateInitialR))) + 
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) +
  facet_wrap(~InhibitorInitialR) + theme_min + xlab("Time (hours)") +
  ylab(expression(paste("Reaction Product (", mu, "mol (MUB) ", g~(DW)^{-1}, ")"))) +
  labs(fill=expression(paste("Substrate\n",mu, "mol ",g~(DW)^{-1}))) +
  theme(legend.position = c(0.85, 0.3))

#2. With the benzoic acid
#Weighted 45.309 mg of AMCL to 100 ml and diluted
#DW is 0.2269
CTA_AMCL_B <- E_calcAMCL(dataset = "../HydrolyticRaw/Soil1-A-AMC-LEU24032023-PLUS.xlsx",
                         MUBconc = c(0, 0.25, 0.5, 1, 2, 5),
                         APconc = c(153.7/100, 153.7/20, 153.7/4, 153.7/2, 153.7),
                         Iconc = c(0, 30.18/20, 30.18/4, 30.18/2, 30.18),
                         Nmeasure = 14, DW = 0.2269, empty = 0)

#Exporting re-calculated data
CTA_AMCL_B$data$Catchment <- c("Certovo")
CTA_AMCL_B$data$Horizon <- c("Organic horizon")
CTA_AMCL_B$data$Enzyme <- c("Leu-Aminopeptidase")
CTA_AMCL_B$data$BenzoicAcid <- c("True")
#Rounding the initial Substrate and Inhibitor concentrations for nicer PLAtting
CTA_AMCL_B$data$SubstrateInitialR <- round(CTA_AMCL_B$data$SubstrateInitial, 0)
CTA_AMCL_B$data$InhibitorInitialR <- round(CTA_AMCL_B$data$InhibitorInitial, 0)
Edata <- rbind(Edata, CTA_AMCL_B$data)

#Plotting the progress curves
Edata %>% mutate(Tbins = cut(ReactionTime, 
                             breaks = seq(0, max(Edata$ReactionTime), by = 0.1))) %>% 
  group_by(Tbins, Catchment, Horizon, Enzyme, BenzoicAcid, SubstrateInitialR, InhibitorInitialR) %>% 
  filter(Catchment == "Certovo" & Horizon == "Organic horizon" & Enzyme == "Leu-Aminopeptidase" & BenzoicAcid == "True") %>% 
  summarise(x = mean(ReactionTime), y = mean(Product), ySD = sd(Product)) %>% 
  ggplot(aes(x, y)) + geom_point(pch = 21, cex = 6, aes(fill = factor(SubstrateInitialR))) + 
  geom_errorbar(aes(ymin = y - ySD, ymax = y + ySD)) +
  facet_wrap(~InhibitorInitialR) + theme_min + xlab("Time (hours)") +
  ylab(expression(paste("Reaction Product (", mu, "mol (MUB) ", g~(DW)^{-1}, ")"))) +
  labs(fill=expression(paste("Substrate\n",mu, "mol ",g~(DW)^{-1}))) +
  theme(legend.position = c(0.85, 0.3))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Data evaluation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# #Example of using python odeint function in R - library set is required NOT RUN
# scpy$odeint(WI, y0 = c(100, 0, 0), t = seq(0, 200), args = tuple(pars = c(0.1, 10)))


library(ggplot2)
library(data.table)
library(stringr)
library(dplyr)
library(psych)
library(ggpubr)
library(animation)
library(tidyverse)
library(survival)
library(flexsurv)
library(lubridate)
library(survminer)
library(RCurl)
library(httr)
library(jsonlite)

source("EarlyIncidencePrediction.R")
source("WithinHospModel.R")
source("FigurePlotting.R")
source("IncidenceModel.R")
source("Survival.R")
source("ExpertParameters.R")

readLandKreisNames()
RKIinciList<-getInciForKreiseRKI(c("Emmendi","Freiburg","Hochschwarzwald"),loadData = T)

sum(RKIinciList$all)

singleRtEsti(RKIinciList$all,
                    min(RKIinciList$Date),
                    as.Date("2020-09-28"),
                    pexp(1:100,1/5)-pexp(0:99,1/5),
                    pgamma(1:100,1.87,0.28)-pgamma(0:99,1.87,0.28),
                    runLength=150,
                    underreporting=1)

myEst<-multiRtEsti(RKIinciList$all,
                   min(RKIinciList$Date),
                   as.Date("2020-09-28"),
                   c(1,0),
                   # pexp(1:100,1/3)-pexp(0:99,1/3),
                   pgamma(1:100,1.87,0.28)-pgamma(0:99,1.87,0.28),
                   runLength=150,
                   underreporting=1,
                   numRuns=1000)
warnings()

getReShadedPlot(myEst,ylimits = c(0,6))

source("EarlyIncidencePrediction.R")
source("WithinHospModel.R")
source("FigurePlotting.R")
source("IncidenceModel.R")
source("Survival.R")
source("ExpertParameters.R")

readLandKreisNames()

RKIinciList<-getInciForKreiseRKI(c("Emmendi"),loadData = T)

#Loading the observed data for Emmendingen
obsData<-read.csv("Resources/HospitalData/covid_em.csv",sep=";")
UKFparaList<-readRDS("Resources/UKFsurvivalPara.rds")
warnings()
directICUProp<-0

obsData$Date<-as.Date(obsData$Date)
obsData$Hosp<-obsData$ICU+obsData$Normal
obsData<-obsData[!is.na(obsData$Hosp),]

RKIinciList[1,"thisPop"]

mrRes<-multiInciRun(RKIinciList$all[-length(RKIinciList)],
                    startDate=min(RKIinciList$Date),
                    currentDay=as.Date("2020-11-30"),
                    population=RKIinciList[1,"thisPop"],
                    c(1,0),#pexp(1:100,1/1)-pexp(0:99,1/1),
                    pgamma(1:100,1.87,0.28)-pgamma(0:99,1.87,0.28),
                    numRuns=25,
                    runLength=60,
                    underreporting=1,
                    fixR=T)

admRateEsti<-admRateEstimator(RKIinciList,obsData)

ggplot(admRateEsti)+
  geom_line(aes(x=Date,y=movAvAR))+
  geom_point(aes(x=Date,y=admProps))+
  theme_linedraw()

emWithinHosp<-withinHospitalMultiRunVaryAdmRate(
  numRuns=10,mrRes,
   admRateEsti$movAvAR,
   (1-directICUProp),#pN #Proportion to normal ward
   (directICUProp),#pI #Proportion direct to ICU
   c(1,0),
   getAlphaIDistr(UKFparaList,as.Date("2020-04-28")),
   getAlphaDDistr(UKFparaList,as.Date("2020-04-28")),
   getDeltaHDistr(UKFparaList,as.Date("2020-04-28")),
   getDeltaIDistr(UKFparaList,as.Date("2020-04-28")),
   getDeltaDDistr(UKFparaList,as.Date("2020-04-28"))
  )
  
#This is still needed to plot the points in the BedShadedPlot... Should really rename the variables.
bedObsEm<-data.frame(
  time=obsData$Date,
  timelineICU=obsData$ICU,
  timelineGW=obsData$Normal,
  Hosp=obsData$Hosp,
  Date=obsData$Date
)

getBedShadedPlot(emWithinHosp,
                 ylimits=c(0,50),
                 #bothObs = getObservedBeds(UKFparaList,as.Date("2020-04-28"),allObs=F),
                 bothObs=bedObsEm,
                 startD=min(mrRes$time),
                 cutDate=as.Date("2020-12-12"),
                 xlimits=c(as.Date("2020-02-01"),as.Date("2021-02-01")),
                 plotAll=T
)


#Same example, without the data for October
obsDataNoOct<-obsData[obsData$Date<as.Date("2020-09-29")|obsData$Date>as.Date("2020-11-01"),]
admRateEstiNoOct<-admRateEstimator(RKIinciList,obsDataNoOct)

ggplot(admRateEstiNoOct)+
  geom_line(aes(x=Date,y=movAvAR))+
  geom_line(data=admRateEsti,aes(x=Date,y=movAvAR),linetype=2)+
  geom_point(aes(x=Date,y=admProps))+
  theme_linedraw()

emWithinHospNoOct<-withinHospitalMultiRunVaryAdmRate(
  numRuns=10,mrRes,
  admRateEstiNoOct$movAvAR,
  (1-directICUProp),#pN #Proportion to normal ward
  (directICUProp),#pI #Proportion direct to ICU
  c(1,0),
  getAlphaIDistr(UKFparaList,as.Date("2020-04-28")),
  getAlphaDDistr(UKFparaList,as.Date("2020-04-28")),
  getDeltaHDistr(UKFparaList,as.Date("2020-04-28")),
  getDeltaIDistr(UKFparaList,as.Date("2020-04-28")),
  getDeltaDDistr(UKFparaList,as.Date("2020-04-28"))
)

#This is still needed to plot the points in the BedShadedPlot... Should really rename the variables.
bedObsEmNoOct<-data.frame(
  time=obsDataNoOct$Date,
  timelineICU=obsDataNoOct$ICU,
  timelineGW=obsDataNoOct$Normal,
  Hosp=obsDataNoOct$Hosp,
  Date=obsDataNoOct$Date
)

getBedShadedPlot(emWithinHospNoOct,
                 ylimits=c(0,50),
                 #bothObs = getObservedBeds(UKFparaList,as.Date("2020-04-28"),allObs=F),
                 bothObs=bedObsEmNoOct,
                 startD=min(mrRes$time),
                 cutDate=as.Date("2020-12-12"),
                 xlimits=c(as.Date("2020-02-01"),as.Date("2021-02-01")),
                 plotAll=T
)

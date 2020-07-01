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

source("EarlyIncidencePrediction.R")
source("WithinHospModel.R")
source("FigurePlotting.R")
source("IncidenceModel.R")
source("Survival.R")
source("ExpertParameters.R")

###################################### Forecasts ######################################
############################# Static incidence forecast #######
#Load and combine foreign data
useDataAllDF<-combineRegionsDB(
  list(
    pcmdpcImportRegions(),
    pcmdpcImportProvinces(),
    JHUimportJustCountries()
  ),minCases = 1
)
useData200309DF<-subset(useDataAllDF,Date<as.Date("2020-03-15"))

#Compensate for delay between curves
compTrajectories<-compensateTime(useData200309DF,"Italy")
combiCurve<-getCombinedCurve(compTrajectories,c("Italy","Lombardia","Lodi"),4)
combiCurveMir<-mirrorCurve(combiCurve,52,wind=4)
startDate<-subset(compTrajectories,region=="Germany")[1,"Date"]-subset(compTrajectories,region=="Germany")[1,"Day"]+round(subset(compTrajectories,region=="Germany")[1,"Delay"])

# Run the within-hopsital care path model based on the expert panel to estimate bed demand from the resulting combined curve
inHRuns<-withinHospitalMultiRunvEarly(numRuns=100,data.frame(runNum=1,underlying=290646*combiCurveMir$smoothInciPC),
                                      pL_all,#pL proportion at home
                                      pN_all,#pN #Proportion to normal ward
                                      pI_all,#pI #Proportion direct to ICU
                                      pexp(1:100,1/2)-pexp(0:99,1/2), #alphaH #from confirmation to admission
                                      pexp(1:100,alphaI_all)-pexp(0:99,alphaI_all), #alphaI #admission to ICU from hospita
                                      pexp(1:100,deltaH_all)-pexp(0:99,deltaH_all), #deltaH #Death / Discharge from hospital
                                      pexp(1:100,deltaI_all)-pexp(0:99,deltaI_all), #deltaI #Death / Discharge from ICU
                                      c(1)
)


inHRunExtended<-withinHospitalMultiRunvEarly(numRuns=100,data.frame(runNum=1,underlying=1027351*combiCurveMir$smoothInciPC),
                                      pL_all,#pL proportion at home
                                      pN_all,#pN #Proportion to normal ward
                                      pI_all,#pI #Proportion direct to ICU
                                      pexp(1:100,1/2)-pexp(0:99,1/2), #alphaH #from confirmation to admission
                                      pexp(1:100,alphaI_all)-pexp(0:99,alphaI_all), #alphaI #admission to ICU from hospita
                                      pexp(1:100,deltaH_all)-pexp(0:99,deltaH_all), #deltaH #Death / Discharge from hospital
                                      pexp(1:100,deltaI_all)-pexp(0:99,deltaI_all), #deltaI #Death / Discharge from ICU
                                      c(1)
)
############################# Dynamic incidence forecast ############################################
# load data for all "Kreise" in Baden-Württemberg
casesKreiseBW<-read.csv("Resources/LandkreiseBW_RKI_23062020.csv")

# Select the local Kreise
selectedKreise<-
  casesKreiseBW %>%
  subset(Landkreis %in% c("Breisgau-Hochschwarzwald","Freiburg im Breisgau")) %>%
  group_by(Date) %>%
  summarise(all = sum(Cases)
  ) %>%
  as.data.frame()
selectedKreise$Date<-as.Date(selectedKreise$Date,format="%d/%m/%Y")
selectedKreise<-selectedKreise[order(selectedKreise$Date),]

#Select those numbers reported on or before 5 April
selectedKreiseSub<-subset(selectedKreise,as.Date(Date)<=(as.Date("2020-04-05")))

#Produce the incidence based on the latest data
mrRes<-multiInciRun((getMoveAv(getIncCases(selectedKreise$all),3)),
                    min(selectedKreise$Date),
                    as.Date("2020-04-05"),
                    492000,
                    pexp(1:100,1/7)-pexp(0:99,1/7),
                    pgamma(1:100,1.87,0.28)-pgamma(0:99,1.87,0.28),
                    numRuns=200,
                    runLength=150,
                    underreporting=1)

############################# Dynamic incidence forecast for each day #####################
#Produce the incidence estimates based on incrementally added data (day by day) 
#Needed for SI.
#Prone to crashing... probably a memory issue. Save it step by step.

lapply(1:90,function(dd){
  casessub<-subset(selectedKreise,as.Date(Date)<(as.Date("2020-03-20")+dd))
  saveRDS(multiInciRun(getMoveAv(getIncCases(casessub$all),3),
                       min(casessub$Date),
                       as.Date("2020-03-20")+dd-1,
                       492000,
                       pexp(1:100,1/7)-pexp(0:99,1/7),
                       pgamma(1:100,1.87,0.28)-pgamma(0:99,1.87,0.28),
                       numRuns=200,
                       runLength=150,
                       underreporting=1),paste0("200623-MrRes-",dd,".RDS"))
})

#Combine and save the results:
mrResList<-lapply(1:90,function(x)readRDS(paste0("200623-MrRes-",x,".RDS")))
saveRDS(mrResList,"200626-mrResList.rds")
############################# Dynamic incidence, expert panel care path ############################################

bedPredictionsFixedPara<-lapply(10:90,getBedPredictionFixedPara)
saveRDS(bedPredictionsFixedPara,"200626-allMidBedPredictions.rds")

bedPredictionsFixedPara<-readRDS("200626-allMidBedPredictions.rds")

############################# Survival analysis ############################################
# We are not allowed to share individual patient data
# We therefore provide the results of our survival analysis, used in the following steps

#ourPatientsUKf<-read.csv("/media/Data/CountryData/Germany/Freiburg/COVID/ergebnis.csv",sep=";")
#UKFparaList<-lapply(1:87,function(x){getSurvParas(ourPatientsUKf, as.Date("2020-02-01")+x)})
#saveRDS(UKFparaList,"Resources/UKFsurvivalPara.rds")
# Read the results
UKFparaList<-readRDS("Resources/UKFsurvivalPara.rds")

#Extract the results for the exponential fits of hazards:
expParas<-
  data.frame(
    disGW1=unlist(
      lapply(60:length(UKFparaList),
             function(x)UKFparaList[[x]]$disGW1_exp$res[1])),
    transfGW1=unlist(
      lapply(60:length(UKFparaList),
             function(x)UKFparaList[[x]]$transfGW1_exp$res[1])),
    disICU=unlist(
      lapply(60:length(UKFparaList),
             function(x)UKFparaList[[x]]$disICU_exp$res[1])),
    transfICU=unlist(
      lapply(60:length(UKFparaList),
             function(x)UKFparaList[[x]]$transfICU_exp$res[1])),
    disGW2=unlist(
      lapply(60:length(UKFparaList),
             function(x)UKFparaList[[x]]$disGW2_exp$res[1])),
    
    date=do.call("c",lapply(60:length(UKFparaList),function(x)(UKFparaList[[x]]$lastDate)))
  )

############################# Dynamic forecast for a single date ############################################

admHospPart<-(max(cumsum(getRawAdm(UKFparaList,as.Date("2020-04-05"))$cases))/selectedKreise[selectedKreise$Date==as.Date("2020-04-05"),"all"])
directICUProp<-UKFparaList[[getDateIndex(UKFparaList,as.Date("2020-04-05"))]]$propDirectICU

inHRunWSD<-withinHospitalMultiRunvLate(numRuns=100,mrRes,
                                       (1-admHospPart),#pL proportion at home
                                       (1-directICUProp)*admHospPart,#pN #Proportion to normal ward
                                       (directICUProp)*admHospPart,#pI #Proportion direct to ICU
                                       pexp(1:100,1/2)-pexp(0:99,1/2),
                                       getAlphaIDistr(UKFparaList,as.Date("2020-04-05")),
                                       getAlphaDDistr(UKFparaList,as.Date("2020-04-05")),
                                       getDeltaHDistr(UKFparaList,as.Date("2020-04-05")),
                                       getDeltaIDistr(UKFparaList,as.Date("2020-04-05")),
                                       getDeltaDDistr(UKFparaList,as.Date("2020-04-05"))
)

############################# bed forecasts for all possible dates ############################################

allBedPredictions<-lapply(10:length(mrResList),getBedPrediction)
saveRDS(allBedPredictions,"200629-allBedPredictions.rds")

###################################### Plots and results ############################################

############################# Figure 2 #######
createPerCapitaPlot(useData200309DF,
                    selectRegions=c("Germany","Italy","Lombardia","Lodi"),
                    regionColors=c("Gray","Black","Red","Blue"),
                    xlimits=c(as.Date("2020-02-15"),as.Date("2020-03-15")))
ggsave("Figure2A.svg",width = 6,height = 4)

createPerCapitaPlotComp(compTrajectories,
                        selectRegions=c("Germany","Italy","Lombardia","Lodi"),
                        regionColors=c("Gray","Black","Red","Blue"),
                        xlimits = c(50,95),ylimits=c(10^-7,10^-2)
)
ggsave("Figure2B-aligned.svg",width = 6,height = 4)

ggplot(combiCurveMir)+
  geom_line(aes(x=startDate+Day,y=290646*smoothInciPC),color="#888888")+
  geom_line(data=combiCurve,aes(x=startDate+Day,y=290646*smoothInciPC),color="black")+
  labs(x="",y="Projected confirmed cases UKFreiburg")+
  scale_x_date(limits=c(as.Date("2020-03-14"),as.Date("2020-05-01")))+
  theme_minimal()
ggsave("Figure2C-CombinedIncidence.svg")

getBedShadedPlot(inHRuns,ylimits=c(0,130),startD=startDate+34,xlimits=c(as.Date("2020-03-01"),as.Date("2020-06-01")))
ggsave("Figure2D-Beds-all.svg",width = 6,height=5)

############################# In text results for fig 2#######
icuResOverRuns<-melt(setDT(getICUTimeline(inHRuns)),id="time")%>%
  group_by(variable) %>%
  summarise(max = max(value),
            peakDate = as.Date("2020-01-10")+time[value==max(value)][1]
  ) %>%
  as.data.frame()
quantile(icuResOverRuns$max,probs=c(0.25,0.5,0.75))
median(icuResOverRuns$peakDate)

gwResOverRuns<-melt(setDT(getHospTimeline(inHRuns)),id="time")%>%
  group_by(variable) %>%
  summarise(max = max(value),
            peakDate = as.Date("2020-01-10")+time[value==max(value)][1]
  ) %>%
  as.data.frame()
median(gwResOverRuns$max)
quantile(gwResOverRuns$max,probs=c(0.25,0.5,0.75))
median(gwResOverRuns$peakDate)

max(gwResOverRuns$max)
max(gwResOverRuns$peakDate)

resOverTime<-melt(setDT(getHospTimeline(inHRuns)),id="time")%>%
  group_by(time) %>%
  summarise(med = quantile(value,probs=0.5),
            mea = mean(value),
            Q25 = quantile(value,probs=0.25),
            Q95 = quantile(value,probs=0.75)
  ) %>%
  as.data.frame()

resICUOverTime<-melt(setDT(getICUTimeline(inHRuns)),id="time")%>%
  group_by(time) %>%
  summarise(med = quantile(value,probs=0.5),
            mea = mean(value),
            Q25 = quantile(value,probs=0.25),
            Q95 = quantile(value,probs=0.75)
  ) %>%
  as.data.frame()

max(resOverTime$med)
resOverTime[resOverTime$med==max(resOverTime$med),]
as.Date("2020-01-10")+59

max(resICUOverTime$med)
resICUOverTime[resICUOverTime$med==max(resICUOverTime$med),]
as.Date("2020-01-10")+69

############################# Figure 3 panels ############################

UKFparaList[[64]]$lastDate
plot(UKFparaList[[64]]$disGW1_wb)
plot(UKFparaList[[64]]$disICU_exp)
plot(UKFparaList[[64]]$disGW2_exp)
plot(UKFparaList[[64]]$transfGW1_exp)
plot(UKFparaList[[64]]$transfICU_wb)

ggplot(expParas)+
  geom_line(aes(x=date,y=disGW1,color="Discharge/Death General ward"))+
  geom_line(aes(x=date,y=disICU,color="Discharge/Death Intensive care"))+
  geom_line(aes(x=date,y=disGW2,color="Discharge/Death Step down unit"))+
  geom_line(aes(x=date,y=transfGW1,color="Transfer GW to ICU"))+
  geom_line(aes(x=date,y=transfICU,color="Transfer ICU to SD"))+
  scale_color_manual(limits=c("Discharge/Death General ward",
                              "Discharge/Death Intensive care",
                              "Discharge/Death Step down unit",
                              "Transfer GW to ICU",
                              "Transfer ICU to SD"
  ),
  values=c("black","red","blue","gray","orange"))+
  scale_x_date(limits=c(as.Date("2020-04-01"),as.Date("2020-04-28")))+
  labs(x="Date of estimation",y="Rate")+
  theme_minimal()

############################# Figure 4 panels ############################
#This uses the list of all results
#plotTheseResults<-mrResList[[17]]
#plotBedResults<-allBedPredictions[[17]]$rawResults

#following just on a single date: (April 5)
plotTheseResults<-mrRes
plotBedResults<-inHRunWSD

getIncidencePlot(plotTheseResults,xlimits = c(as.Date("2020-03-01"),as.Date("2020-04-15")))
ggsave("Fig4-A-Incidence.svg",width=6,height=4)

getRePlot(plotTheseResults,ylimits=c(0,7.5),xlimits = c(as.Date("2020-02-14"),as.Date("2020-04-06")))
ggsave("Fig4-B-Rt.svg",width=6,height=5)

getRePlot2(plotTheseResults,
           ylimits=c(0,7.5),
           xlimits = c(as.Date("2020-02-14"),as.Date("2020-04-06"))
      )
      
ggsave("Fig4-B-Rt-2.svg",width=6,height=5)

getSprinklePlotIQR(plotTheseResults,
                   xlimits = c(as.Date("2020-03-01"),as.Date("2020-06-01")),
                   ylimits=c(0,250),
                   cutDate=as.Date("2020-04-05"))

ggsave("Fig4-C-sprinkle.svg",width=6,height=4)

getBedShadedPlot(plotBedResults,
                 ylimits=c(0,175),
                 bothObs = getObservedBeds(UKFparaList,as.Date("2020-04-28"),allObs=F),
                 startD=min(plotTheseResults$time),
                 cutDate=as.Date("2020-04-05"),
                 xlimits=c(as.Date("2020-03-01"),as.Date("2020-07-01"))
                 
                 )
ggsave("Fig4-D-Beds-wObs.svg",width = 6,height = 5)

############################# Results in text for figure 4####################################

UKFparaList[[length(UKFparaList)]]$nAllPats
UKFparaList[[64]]$nAllPats
UKFparaList[[64]]$lastDate
UKFparaList[[64]]$nICUpats

resOverTime<-plotTheseResults %>%
  group_by(time) %>%
  summarise(med = quantile(reportedRaw,probs=0.5),
            mea = mean(reportedRaw),
            Q25 = quantile(reportedEpi,probs=0.25),
            Q95 = quantile(reportedEpi,probs=0.75)
  ) %>%
  as.data.frame()

resOverRun<-plotTheseResults%>%
  group_by(runNum) %>%
  summarise(max = max(reportedRaw),
            peakDate = time[reportedRaw==max(reportedRaw)][1],
            Rtf=mean(recentR0),
            decl=mean(decl)
  ) %>%
  as.data.frame()

median(resOverRun$peakDate)
median(resOverRun$max)
max(resOverRun$max)
min(resOverRun$max)

max(resOverTime$med)

icuResOverRuns<-melt(setDT(getICUTimeline(plotBedResults)),id="time")%>%
  group_by(variable) %>%
  summarise(max = max(value),
            peakDate = min(plotTheseResults$time)+time[value==max(value)][1]
  ) %>%
  as.data.frame()
quantile(icuResOverRuns$max,probs=c(0.25,0.5,0.75))
median(icuResOverRuns$peakDate)

gwResOverRuns<-melt(setDT(getHospTimeline(plotBedResults)),id="time")%>%
  group_by(variable) %>%
  summarise(max = max(value),
            peakDate = min(plotTheseResults$time)+time[value==max(value)][1]
  ) %>%
  as.data.frame()
median(gwResOverRuns$max)
quantile(gwResOverRuns$max,probs=c(0.25,0.5,0.75))
median(gwResOverRuns$peakDate)

max(gwResOverRuns$max)
max(gwResOverRuns$peakDate)

############################# Supplementary figure #######
allForecast4Plots<-lapply(1:(length(mrResList)-9),function(x){produceFore4Plots(x)})
justThePlots<-do.call(c, lapply(1:46,function(x) list(allForecast4Plots[[x]]$theInciPlot,
                                                      allForecast4Plots[[x]]$theRePlot,
                                                      allForecast4Plots[[x]]$theSprinklePlot,
                                                      allForecast4Plots[[x]]$theLateBedPlot)))

allForecastPlots<-lapply(1:length(mrResList),function(x){produceForePlots(x)})
justThePlots<-do.call(c, lapply(1:46,function(x) list(allForecastPlots[[x]]$theInciPlot,
                                                      allForecastPlots[[x]]$theRePlot,
                                                      allForecastPlots[[x]]$theSprinklePlot,
                                                      allForecastPlots[[x]]$theMidBedPlot,
                                                      allForecastPlots[[x]]$theLateBedPlot)))

tryLabs<-unlist(lapply(1:40,
                       function(l)(c(as.character(l+max(mrResList[[1]][!is.na(mrResList[[1]]$ReEsti),"time"])),rep(" ",4)))))
myMultiPage<-lapply(1:13,function(x)ggarrange(plotlist=justThePlots[((x-1)*15)+1:15],ncol = 5,nrow=3,labels = tryLabs[((x-1)*15)+1:15]))
lapply(1:13,function(x) ggsave(paste0("200513-Beds-page",x+10,".svg"),plot=myMultiPage[[x]],height = 8,width=14))

myMultiPage<-lapply(1:13,function(x)ggarrange(plotlist=justThePlots[((x-1)*15)+1:15],ncol = 5,nrow=3,labels = tryLabs[((x-1)*15)+1:15]))
ggsave("200513-Beds-page00.svg",
       plot=ggarrange(plotlist=justThePlots[((2-1)*15)+1:15],
                      ncol = 5,nrow=3,labels = tryLabs[((2-1)*15)+1:15]),
       #myMultiPage[[2]],
       height=8,
       width=14)

lapply(1:13,function(x) print(paste0("200513-Beds-page",str_pad(toString(x),width=2,pad = "0"),".svg")))

####Movie export
videoSlide<-function(x){
  ggsave(paste0("200626-Beds-slide",str_pad(toString(x),width=2,pad = "0"),".png"),
         plot=ggarrange(plotlist=justThePlots[((x-1)*4)+1:4],ncol = 2,nrow=2),
         height = 8,width=8)}

lapply(1:40,videoSlide)

# bed demand forecast based on static model with a population of 1027351 FR-extended:
getBedShadedPlot(inHRunExtended,ylimits=c(0,430),startD=startDate+34,xlimits=c(as.Date("2020-03-01"),as.Date("2020-06-01")))
ggsave("FigureS6-Beds-all.svg",width = 6,height=5)

############################# Export all parameters and AICs for SI table ############################################
allParas<-
  data.frame(
    disGW1=unlist(
      lapply(60:length(UKFparaList),
             function(x)c(
               paste0(round(UKFparaList[[x]]$disGW1_exp$res[1],digits=4)," (",
                      round(UKFparaList[[x]]$disGW1_exp$res[2],digits=4),"-",
                      round(UKFparaList[[x]]$disGW1_exp$res[3],digits=4),") AIC: ",
                      round(UKFparaList[[x]]$disGW1_exp$AIC,2)
               )
               
             ))),
    transfGW1=unlist(
      lapply(60:length(UKFparaList),
             function(x)c(
               paste0(round(UKFparaList[[x]]$transfGW1_exp$res[1],digits=4)," (",
                      round(UKFparaList[[x]]$transfGW1_exp$res[2],digits=4),"-",
                      round(UKFparaList[[x]]$transfGW1_exp$res[3],digits=4),") AIC: ",
                      round(UKFparaList[[x]]$transfGW1_exp$AIC,2)
               )
             ))),
    
    disICU=unlist(
      lapply(60:length(UKFparaList),
             function(x)c(
               paste0(round(UKFparaList[[x]]$disICU_exp$res[1],digits=4)," (",
                      round(UKFparaList[[x]]$disICU_exp$res[2],digits=4),"-",
                      round(UKFparaList[[x]]$disICU_exp$res[3],digits=4),") AIC: ",
                      round(UKFparaList[[x]]$disICU_exp$AIC,2)
               )
             ))),
    
    transfICU=unlist(
      lapply(60:length(UKFparaList),
             function(x)c(
               paste0(round(UKFparaList[[x]]$transfICU_exp$res[1],digits=4)," (",
                      round(UKFparaList[[x]]$transfICU_exp$res[2],digits=4),"-",
                      round(UKFparaList[[x]]$transfICU_exp$res[3],digits=4),") AIC: ",
                      round(UKFparaList[[x]]$transfICU_exp$AIC,2)
               )
             ))),
    
    disGW2=unlist(
      lapply(60:length(UKFparaList),
             function(x)c(
               paste0(round(UKFparaList[[x]]$disGW2_exp$res[1],digits=4)," (",
                      round(UKFparaList[[x]]$disGW2_exp$res[2],digits=4),"-",
                      round(UKFparaList[[x]]$disGW2_exp$res[3],digits=4),") AIC: ",
                      round(UKFparaList[[x]]$disGW2_exp$AIC,2)
               )
               
             ))),
    
    date=do.call("c",lapply(60:length(UKFparaList),function(x)(UKFparaList[[x]]$lastDate)))
  )
write.csv(allParas, "survExpPara.csv")

allParasWB<-
  data.frame(
    disGW1=unlist(
      lapply(60:length(UKFparaList),
             function(x)c(
               paste0("Shape: ",
                      round(UKFparaList[[x]]$disGW1_wb$res[1],digits=4)," (",
                      round(UKFparaList[[x]]$disGW1_wb$res[3],digits=4),"-",
                      round(UKFparaList[[x]]$disGW1_wb$res[5],digits=4),")\n Scale: ",
                      round(UKFparaList[[x]]$disGW1_wb$res[2],digits=4)," (",
                      round(UKFparaList[[x]]$disGW1_wb$res[4],digits=4),"-",
                      round(UKFparaList[[x]]$disGW1_wb$res[6],digits=4),")\n AIC: ",
                      round(UKFparaList[[x]]$disGW1_wb$AIC,2)
               )
               
             ))),
    transfGW1=unlist(
      lapply(60:length(UKFparaList),
             function(x)c(
               paste0("Shape: ",
                      round(UKFparaList[[x]]$transfGW1_wb$res[1],digits=4)," (",
                      round(UKFparaList[[x]]$transfGW1_wb$res[3],digits=4),"-",
                      round(UKFparaList[[x]]$transfGW1_wb$res[5],digits=4),")\n Scale: ",
                      round(UKFparaList[[x]]$transfGW1_wb$res[2],digits=4)," (",
                      round(UKFparaList[[x]]$transfGW1_wb$res[4],digits=4),"-",
                      round(UKFparaList[[x]]$transfGW1_wb$res[6],digits=4),")\n AIC: ",
                      round(UKFparaList[[x]]$transfGW1_wb$AIC,2)
               )
             ))),
    
    disICU=unlist(
      lapply(60:length(UKFparaList),
             function(x)c(
               paste0("Shape: ",round(UKFparaList[[x]]$disICU_wb$res[1],digits=4)," (",
                      round(UKFparaList[[x]]$disICU_wb$res[3],digits=4),"-",
                      round(UKFparaList[[x]]$disICU_wb$res[5],digits=4),")\n Scale: ",
                      round(UKFparaList[[x]]$disICU_wb$res[2],digits=4)," (",
                      round(UKFparaList[[x]]$disICU_wb$res[4],digits=4),"-",
                      round(UKFparaList[[x]]$disICU_wb$res[6],digits=4),")\n AIC: ",
                      round(UKFparaList[[x]]$disICU_wb$AIC,2)
               )
             ))),
    
    transfICU=unlist(
      lapply(60:length(UKFparaList),
             function(x)c(
               paste0("Shape: ",round(UKFparaList[[x]]$disICU_wb$res[1],digits=4)," (",
                      round(UKFparaList[[x]]$transfICU_wb$res[3],digits=4),"-",
                      round(UKFparaList[[x]]$transfICU_wb$res[5],digits=4),")\n Scale: ",
                      round(UKFparaList[[x]]$transfICU_wb$res[2],digits=4)," (",
                      round(UKFparaList[[x]]$transfICU_wb$res[4],digits=4),"-",
                      round(UKFparaList[[x]]$transfICU_wb$res[6],digits=4),")\n AIC: ",
                      round(UKFparaList[[x]]$transfICU_wb$AIC,2)
               )
             ))),
    
    disGW2=unlist(
      lapply(60:length(UKFparaList),
             function(x)c(
               paste0("Shape: ",
                      round(UKFparaList[[x]]$disGW2_wb$res[1],digits=4)," (",
                      round(UKFparaList[[x]]$disGW2_wb$res[3],digits=4),"-",
                      round(UKFparaList[[x]]$disGW2_wb$res[5],digits=4),")\n Scale: ",
                      round(UKFparaList[[x]]$disGW2_wb$res[2],digits=4)," (",
                      round(UKFparaList[[x]]$disGW2_wb$res[4],digits=4),"-",
                      round(UKFparaList[[x]]$disGW2_wb$res[6],digits=4),")\n AIC: ",
                      round(UKFparaList[[x]]$disGW2_wb$AIC,2)
               )
               
             ))),
    
    date=do.call("c",lapply(60:length(UKFparaList),function(x)(UKFparaList[[x]]$lastDate)))
  )
write.csv(allParasWB, "survWBPara.csv")

############################# Extra plots, not in Suppl. Mat.############################################
ggplot(resOverRun)+
  geom_point(aes(x=peakDate,y=max))+
  labs(x="Peak date",y="Peak size")+
  scale_y_log10()+theme_linedraw()
ggsave("FigS-sim-A.svg",width=5,height=5)
ggplot(resOverRun)+
  geom_point(aes(x=Rtf,y=max))+
  labs(x="Most recent Rt",y="Peak size")+
  scale_y_log10()+theme_linedraw()
ggsave("FigS-sim-B.svg",width=5,height=5)
ggplot(resOverRun)+
  geom_point(aes(x=decl,y=max))+
  scale_y_log10()+
  labs(y="Peak size",x="Rt rate of decline")+
  theme_linedraw()
ggsave("FigS-sim-C.svg",width=5,height=5)
ggplot(resOverRun)+
  geom_point(aes(x=Rtf,y=decl))+
  labs(x="Most recent Rt",y="Rt rate of decline")+
  theme_linedraw()
ggsave("FigS-sim-D.svg",width=5,height=5)

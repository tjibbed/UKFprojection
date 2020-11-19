############## plotting functions for early estimate ############## 
createPerCapitaPlot<-function(trajectoriesDB,selectRegions=NULL,regionColors=NULL,xlimits=NULL,ylimits=NULL){
  if(is.null(selectRegions)){
    selectRegions<-unique(trajectoriesDB$region)
  } else {
    trajectoriesDB<-trajectoriesDB[trajectoriesDB$region %in% selectRegions,]
  }
  
  if(is.null(xlimits)){
    xlimits=c(min(trajectoriesDB$Date),max(trajectoriesDB$Date))
  }
  if(is.null(ylimits)){
    withinx<-trajectoriesDB[(trajectoriesDB$Date>=xlimits[1])&(trajectoriesDB$Date<=xlimits[2]),]
    ylimits=c(0.9*min(withinx$CasesPerCapita),1.1*max(withinx$CasesPerCapita))
  }
  ggplot(trajectoriesDB)+
    geom_line(aes(x=Date,y=CasesPerCapita,group=region,colour=region))+
    scale_y_log10(limits=ylimits)+
    (if(is.null(regionColors)){
      (scale_colour_discrete())
    } else {
      (scale_color_manual(values = regionColors,
                          limits = selectRegions))
    })+
    ylab("Cases per capita")+
    xlab("Date")+
    scale_x_date(limits=xlimits)+
    theme_minimal()
}

createPerCapitaPlot<-function(trajectoriesDB,selectRegions=NULL,regionColors=NULL,xlimits=NULL,ylimits=NULL){
  if(is.null(selectRegions)){
    selectRegions<-unique(trajectoriesDB$region)
  } else {
    trajectoriesDB<-trajectoriesDB[trajectoriesDB$region %in% selectRegions,]
  }
  
  if(is.null(xlimits)){
    xlimits=c(min(trajectoriesDB$Date),max(trajectoriesDB$Date))
  }
  if(is.null(ylimits)){
    withinx<-trajectoriesDB[(trajectoriesDB$Date>=xlimits[1])&(trajectoriesDB$Date<=xlimits[2]),]
    ylimits=c(0.9*min(withinx$CasesPerCapita),1.1*max(withinx$CasesPerCapita))
  }
  ggplot(trajectoriesDB)+
    geom_line(aes(x=Date,y=CasesPerCapita,group=region,colour=region))+
    scale_y_log10(limits=ylimits)+
    (if(is.null(regionColors)){
      (scale_colour_discrete())
    } else {
      (scale_color_manual(values = regionColors,
                          limits = selectRegions))
    })+
    ylab("Cases per capita")+
    xlab("Date")+
    scale_x_date(limits=xlimits)+
    theme_minimal()
}

createPerCapitaPlotComp<-function(trajectoriesDB,selectRegions=NULL,regionColors=NULL,xlimits=NULL,ylimits=NULL){
  if(is.null(selectRegions)){
    selectRegions<-unique(trajectoriesDB$region)
  } else {
    trajectoriesDB<-trajectoriesDB[trajectoriesDB$region %in% selectRegions,]
  }
  
  if(is.null(xlimits)){
    xlimits=c(min(trajectoriesDB$CompTime),max(trajectoriesDB$CompTime))
  }
  if(is.null(ylimits)){
    withinx<-trajectoriesDB[(trajectoriesDB$CompTime>=xlimits[1])&(trajectoriesDB$CompTime<=xlimits[2]),]
    ylimits=c(0.9*min(withinx$CasesPerCapita),1.1*max(withinx$CasesPerCapita))
  }
  ggplot(trajectoriesDB)+
    geom_line(aes(x=CompTime,y=CasesPerCapita,group=region,colour=region))+
    scale_y_log10(limits=ylimits)+
    (if(is.null(regionColors)){
      (scale_colour_discrete())
    } else {
      (scale_color_manual(values = regionColors,
                          limits = selectRegions))
    })+
    ylab("Cases per capita")+
    xlab("Days")+
    scale_x_continuous(limits=xlimits)+
    theme_minimal()
}

############# Plots from care-path model#############################
getBedAllLinePlot<-function(runResults,startD=as.Date("2020-01-10"),xlimits=NULL,ylimits=NULL){
  allICUTlFR<-melt(as.data.table(getICUTimeline(runResults)),id.vars="time")
  allhospTlFR<-melt(as.data.table(getHospTimeline(runResults)),id.vars="time")
  if(is.null(xlimits)){
    xlimits=as.Date("2020-01-10")+c(min(runResults$time),max(runResults$time))
  }
  ggplot(allICUTlFR)+
    geom_line(aes(x=startD+time,y=value,group=variable),color="#FF000020")+
    geom_line(data=allhospTlFR,aes(x=startD+time,y=value,group=variable),color="#0000FF20")+
    scale_x_date(limits=xlimits)+
    scale_y_continuous(limits=ylimits)+
    labs(x="time",y="Occupied beds")+
    theme_minimal()
}

getBedShadedPlot<-function(runResults,startD=as.Date("2020-01-10"),plotAll=FALSE,bothObs=NULL,ICUobs=NULL, Hospobs = NULL,xlimits=NULL,ylimits=NULL,cutDate=NULL){
  if(is.null(cutDate)) cutDate<-max(runResults$time)
  if(is.null(xlimits)){
    xlimits=c(min(runResults$time),max(runResults$time))
  }
  thisRunInfo<-extractRunInfo(runResults)
  
  ggplot(thisRunInfo)+
    (if(plotAll) geom_ribbon(aes(x=startD+time-1,ymin=allTotal025,ymax=pmin(allTotal075,ylimits[[2]])),fill="#00000044"))+
    (if(plotAll) geom_ribbon(aes(x=startD+time-1,ymin=allTotal005,ymax=pmin(allTotal095,ylimits[[2]])),fill="#00000044"))+
    geom_ribbon(aes(x=startD+time-1,ymin=hospTotal025,ymax=pmin(hospTotal075,ylimits[[2]])),fill="#0000FF44")+
    geom_ribbon(aes(x=startD+time-1,ymin=hospTotal005,ymax=pmin(hospTotal095,ylimits[[2]])),fill="#0000FF44")+
    geom_ribbon(aes(x=startD+time-1,ymin=icuTotal025,ymax=pmin(icuTotal075,ylimits[[2]])),fill="#FF000044")+
    geom_ribbon(aes(x=startD+time-1,ymin=icuTotal005,ymax=pmin(icuTotal095,ylimits[[2]])),fill="#FF000044")+
    (if(plotAll) geom_line(aes(x=startD+time-1,y=allTotal050),color="black"))+
    geom_line(aes(x=startD+time-1,y=hospTotal050),color="blue")+
    geom_line(aes(x=startD+time-1,y=icuTotal050),color="red")+
    
    (if(!is.null(bothObs)) geom_point(data=subset(bothObs,time<=cutDate),aes(x=time,y=timelineICU),pch=21, size=2,color="black",fill="#FF0000FF"))+
    (if(!is.null(bothObs)) geom_point(data=subset(bothObs,time<=cutDate),aes(x=time,y=timelineGW),pch=21, size=2,color="black",fill="#0000FFFF"))+
    (if((!is.null(bothObs))&plotAll) geom_point(data=subset(bothObs,time<=cutDate),aes(x=time,y=timelineGW+timelineICU),pch=21, size=2,color="black",fill="#000000FF"))+
    
    (if(!is.null(bothObs)) geom_point(data=subset(bothObs,time>cutDate),aes(x=time,y=timelineICU),pch=21, size=2,color="#FF000080",fill="#00000000"))+
    (if(!is.null(bothObs)) geom_point(data=subset(bothObs,time>cutDate),aes(x=time,y=timelineGW),pch=21, size=2,color="#0000FF80",fill="#00000000"))+
    (if((!is.null(bothObs))&plotAll) geom_point(data=subset(bothObs,time>cutDate),aes(x=time,y=timelineGW+timelineICU),pch=21, size=2,color="#00000080",fill="#00000000"))+
    
    (if(!is.null(Hospobs)) geom_point(data=Hospobs,aes(x=time,y=obs),pch=21, size=2,color="black",fill="blue"))+
    (if(!is.null(ICUobs)) geom_point(data=ICUobs,aes(x=time,y=obs),pch=21, size=2,color="black",fill="red"))+
    scale_x_date(limits=xlimits)+
    scale_y_continuous(limits=ylimits)+
    labs(x="time",y="Occupied beds")+
    theme_minimal()
}

produceForePlots<-function(lx){
  startDhere<-mrResList[[lx]][[1,"currentDay"]]
  #startDhere<-max(mrResList[[lx]][!is.na(mrResList[[lx]]$ReEsti),"time"])+1
  
  list(
    
    theInciPlot=getIncidencePlot(mrResList[[lx]],xlimits = c(as.Date("2020-03-01"),as.Date("2020-04-28"))),
    theRePlot=getRePlot2(mrResList[[lx]],ylimits=c(0,7.5),xlimits = c(as.Date("2020-02-01"),startDhere+2)),
    theSprinklePlot=getSprinklePlotIQR(mrResList[[lx]],
                                       xlimits = c(as.Date("2020-03-01"),as.Date("2020-06-01")),
                                       ylimits=c(0,250),
                                       cutDate=startDhere),
    theMidBedPlot=getBedShadedPlot(bedPredictionsFixedPara[[lx]]$rawResults,
                                   ylimits=c(0,175),
                                   bothObs = getObservedBeds(UKFparaList,as.Date("2020-04-28"),allObs=F),
                                   startD=min(mrResList[[lx]]$time),
                                   xlimits=c(as.Date("2020-03-01"),as.Date("2020-07-01")),
                                   cutDate=startDhere),
    theLateBedPlot=getBedShadedPlot(allBedPredictions[[lx]]$rawResults,
                                    ylimits=c(0,175),
                                    bothObs = getObservedBeds(UKFparaList,as.Date("2020-04-28"),allObs=F),
                                    startD=min(mrResList[[lx]]$time),
                                    xlimits=c(as.Date("2020-03-01"),as.Date("2020-07-01")),
                                    cutDate=startDhere)
  )
}

produceFore4Plots<-function(lx){
  startDhere<-mrResList[[lx+9]][[1,"currentDay"]]
  print(startDhere)
  
  list(
    theInciPlot=getIncidencePlot(mrResList[[lx+9]],xlimits = c(as.Date("2020-03-01"),as.Date("2020-06-15"))),
    theRePlot=getRePlot2(mrResList[[lx+9]],
                         ylimits=c(0,7.5),
                         xlimits = c(as.Date("2020-03-01"),as.Date("2020-06-15"))),
    
    theSprinklePlot=getSprinklePlotIQR(mrResList[[lx+9]],
                                       xlimits = c(as.Date("2020-03-01"),as.Date("2020-06-15")),
                                       ylimits=c(0,250),
                                       cutDate=startDhere),
    theLateBedPlot=getBedShadedPlot(allBedPredictions[[lx]]$rawResults,
                                    ylimits=c(0,175),
                                    bothObs = getObservedBeds(UKFparaList,as.Date("2020-04-28"),allObs=F),
                                    startD=min(mrResList[[lx+9]]$time),
                                    xlimits=c(as.Date("2020-03-01"),as.Date("2020-06-15")),
                                    cutDate=startDhere)
  )
}

##################### Plots for dynamic incidence model ###################################
getSprinklePlot<-function(runResults,xlimits=NULL,ylimits=NULL){
  if(is.null(xlimits)){
    xlimits=c(min(runResults$time),max(runResults$time))
  }
  #obsLastDate<-(1+max(runResults[!is.na(runResults$ReEsti),"time"]))
  obsLastDate<-runResults[[1,"currentDay"]]-1
  if(is.null(ylimits)){
    ylimits=c(0,max(runResults$reportedRaw))
  }
  ggplot(runResults)+
    geom_line(aes(x=time,y=reportedRaw,group=runNum),colour="#66666624")+
    geom_step(data=subset(runResults,time<=obsLastDate),
              aes(x=time,y=reportedEpi,group=runNum,color=as.character(runNum)),colour="black",direction="vh")+
    scale_x_date(limits=xlimits)+
    scale_y_continuous(limits=ylimits)+
    labs(x="Date",y="Confirmed cases")+
    theme_minimal()+
    theme(legend.position = "none")
}

getUnderlyingPlot<-function(runResults,xlimits=NULL,ylimits=NULL){
  if(is.null(xlimits)){
    xlimits=c(min(runResults$time),max(runResults$time))
  }
  #obsLastDate<-(1+max(runResults[!is.na(runResults$ReEsti),"time"]))
  obsLastDate<-runResults[[1,"currentDay"]]-1
  if(is.null(ylimits)){
    ylimits=c(0,max(runResults$reportedRaw))
  }
  ggplot(runResults)+
    geom_line(aes(x=time,y=underlying,group=runNum),colour="#66666624")+
    geom_step(data=subset(runResults,time<=obsLastDate),
              aes(x=time,y=reportedEpi,group=runNum,color=as.character(runNum)),colour="black",direction="vh")+
    scale_x_date(limits=xlimits)+
    scale_y_continuous(limits=ylimits)+
    labs(x="Date",y="Confirmed cases")+
    theme_minimal()+
    theme(legend.position = "none")
}

getIncidencePlot<-function(runResults,xlimits=NULL){
  if(is.null(xlimits)){
    xlimits=c(min(runResults$time),max(runResults$time))
  }
  #obsLastDate<-(1+max(runResults[!is.na(runResults$ReEsti),"time"]))
  
  obsLastDate<-runResults[[1,"currentDay"]]-1
  ggplot(runResults)+
    # geom_line(aes(x=time,y=reportedRaw,group=runNum),colour="#66666624")+
    geom_step(data=subset(runResults,time<=obsLastDate),
              aes(x=time,y=reportedEpi,group=runNum,color=as.character(runNum)),colour="black",direction="vh")+
    scale_x_date(limits=xlimits)+
    labs(x="Date",y="Confirmed cases")+
    theme_minimal()+
    theme(legend.position = "none")
}

getRePlot<-function(runResults,xlimits=NULL,ylimits=NULL,showFits=T){
  
  #obsLastDate<-(1+max(runResults[!is.na(runResults$ReEsti),"time"]))
  obsLastDate<-runResults[[1,"currentDay"]]-1
  obsFirstDate<-(min(runResults[!is.na(runResults$ReEsti),"time"]))
  
  if(is.null(xlimits)){
    xlimits=c(obsFirstDate,obsLastDate+7)
  }
  if(is.null(ylimits)){
    ylimits=c(0,1.1*max(runResults$devR0))
  }
  
  ggplot(runResults)+
    geom_point(data=runResults[runResults$time<=(obsLastDate),],aes(x=time,y=ReEsti2),colour="#0000FF07")+
    (if(showFits) geom_line(aes(x=time,y=devR0,group=runNum),colour="#00000012"))+
    geom_hline(yintercept = 1)+
    scale_x_date(limits=xlimits)+
    scale_y_continuous(limits=ylimits)+
    labs(y="Effective Reproduction Number",x="Date")+
    theme_minimal()
}

getRePlot2<-function(runResults,xlimits=NULL,ylimits=NULL){
  
  cutPoint<-runResults$currentDay
  myDF<-runResults[runResults$time<cutPoint,c("time","ReEsti2","ReEsti","reportedRaw")]
  myDF<-myDF[!(myDF$reportedRaw==0),]
  myDF<-myDF[!is.na(myDF$ReEsti),c("time","ReEsti2")]

  myDF$ReEsti2<-round(myDF$ReEsti2,digits = 2)
  
  myUniques<-subset(as.data.frame(table(myDF)),Freq>0)
  myUniques$time<-as.Date(as.character(myUniques$time))
  myUniques$ReEsti2<-as.numeric(as.character(myUniques$ReEsti2))
  
  #obsLastDate<-(1+max(runResults[!is.na(runResults$ReEsti),"time"]))
  obsLastDate<-runResults[[1,"currentDay"]]-1
  obsFirstDate<-(min(runResults[!is.na(runResults$ReEsti),"time"]))
  
  
  if(is.null(xlimits)){
    xlimits=c(obsFirstDate,obsLastDate+7)
  }
  if(is.null(ylimits)){
    ylimits=c(0,1.1*max(runResults$devR0))
  }
  
  ggplot(myUniques)+
    geom_point(aes(x=time,y=ReEsti2,color=Freq))+
    geom_line(data=runResults,aes(x=time,y=devR0,group=runNum),colour="#00000012")+
    geom_hline(yintercept = 1)+
    scale_color_gradient(low="#0000FF07","#0000FFFF")+
    scale_x_date(limits=xlimits)+
    scale_y_continuous(limits=ylimits)+
    labs(y="Effective Reproduction Number",x="Date")+
    theme_minimal()+
    theme(legend.position = "none")
  # theme(legend=none)
}

getReShadedPlot<-function(runResults,xlimits=NULL,ylimits=NULL){
  if(is.null(ylimits)){
    ylimits=c(0,max(runResults$ReEsti))
  }
  ggplot(as.data.frame(
    runResults[(!is.na(runResults$ReEsti))&(runResults$underlying!=0),c("time","ReEsti")] %>% 
    group_by(time)%>%
    summarise(
      meanRe=mean(ReEsti),
      medianRe=median(ReEsti),
      low=quantile(ReEsti,probs=0.05),
      high=quantile(ReEsti,probs=0.95)
    )
  ))+
  geom_ribbon(aes(x=time,ymin=low,ymax=pmin(high,ylimits[[2]])
                  ),fill="#0000FF44")+
  geom_line(aes(x=time,y=meanRe))+
  scale_y_continuous(limits=ylimits)+
  (if(!is.null(xlimits))scale_x_date(limits=xlimits))+
  # geom_line(aes(x=time,y=medianRe),colour="blue")+
  theme_linedraw()
}

getReShadedPlotTwice<-function(runResults1,runResults2,xlimits=NULL,ylimits=NULL){
  if(is.null(ylimits)){
    ylimits=c(0,max(c(runResults$ReEsti,runResults$ReEsti)))
  }
  ggplot(as.data.frame(
    runResults[(!is.na(runResults$ReEsti))&(runResults$underlying!=0),c("time","ReEsti")] %>% 
      group_by(time)%>%
      summarise(
        meanRe=mean(ReEsti),
        medianRe=median(ReEsti),
        low=quantile(ReEsti,probs=0.05),
        high=quantile(ReEsti,probs=0.95)
      )
  ))+
    geom_ribbon(aes(x=time,ymin=low,ymax=pmin(high,ylimits[[2]])
    ),fill="#0000FF44")+
    geom_line(aes(x=time,y=meanRe))+
    scale_y_continuous(limits=ylimits)+
    (if(!is.null(xlimits))scale_x_date(limits=xlimits))+
    # geom_line(aes(x=time,y=medianRe),colour="blue")+
    theme_linedraw()
}

getSprinklePlotIQR<-function(runResults,xlimits=NULL,ylimits=NULL,cutDate=NULL){
  summaryRes<-runResults %>%
    group_by(time) %>%
    summarise(med = quantile(reportedEpi,probs=0.5),
            #  mea = mean(reportedRaw),
              Q05 = quantile(reportedEpi,probs=0.05),
              Q25 = quantile(reportedEpi,probs=0.25),
              Q75 = quantile(reportedEpi,probs=0.75),
              Q95 = quantile(reportedEpi,probs=0.95)         
    ) %>%
    as.data.frame()
  if(is.null(cutDate)) cutDate<-runResults[[1,"currentDay"]]-1
  
  if(is.null(xlimits)){
    xlimits=c(min(runResults$time),max(runResults$time))
  }
  #obsLastDate<-(1+max(runResults[!is.na(runResults$ReEsti),"time"]))
  obsLastDate<-runResults[[1,"currentDay"]]-2
  if(is.null(ylimits)){
    ylimits=c(0,max(runResults$reportedEpi))
  }
  ggplot(runResults)+
    geom_line(data=subset(runResults,time>=cutDate),aes(x=time,y=reportedEpi,group=runNum),colour="#66666624")+
    geom_ribbon(data=summaryRes,
                aes(x=time,ymin=Q25,ymax=pmin(Q75,ylimits[[2]])),fill="#00FF0050")+
    geom_ribbon(data=summaryRes,
                aes(x=time,ymin=Q05,ymax=pmin(Q95,ylimits[[2]])),fill="#00FF0033")+
    geom_step(data=subset(summaryRes,time<=cutDate),
              aes(x=time,y=med),colour="black",direction="vh")+
    geom_step(data=subset(summaryRes,time>=cutDate),
              aes(x=time,y=med),colour="green",direction="vh")+
    scale_x_date(limits=xlimits)+
    scale_y_continuous(limits=ylimits)+
    labs(x="Date",y="Confirmed cases")+
    theme_minimal()+
    theme(legend.position = "none")
}


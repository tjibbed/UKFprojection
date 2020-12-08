################### Early estimate model, expert consensus ###################

createPatientsvEarly<-function(incidence,pL,pN,pI,alphaH,alphaI,deltaH,deltaI,pHosp=NULL){
  incidence<-round(incidence)
  pullPatients<-function(nPat,day){
    tAH<-sample(1:(length(alphaH)),nPat,prob = alphaH,replace = T)
    tAI<-sample(1:(length(alphaI)),nPat,prob = alphaI,replace = T)
    tDH<-sample(1:(length(deltaH)),nPat,prob = deltaH,replace = T)
    tDI<-sample(1:(length(deltaI)),nPat,prob = deltaI,replace = T)
    timeDF<-data.frame(tAH,tAI,tDH,tDI,hospStay<-1,ICUstay<-1,type<-1)
    
    pType<-       # H I H I H I 
      sample(list(c(0,0,0,0,0,0,"L"),                  #Stay at home (low-symptomatic)
                  c(1,1,1,1,1,1,"H"),                  #First to hospital
                  c(1,0,0,1,1,1,"I")                 #Immediate ICU admisison
      ),
      nPat,prob=c(pL, pN, pI),replace = T
      )
    hospChoice<-rep(1,nPat)
    
    typeDF<-as.data.frame(t(matrix(unlist((pType)),nrow=7)))
    typeDF2<-sapply(typeDF[,1:6],function(x)as.numeric(as.character(x)))
    actualTimes<-(timeDF[,1:6]*typeDF2)
    colnames(actualTimes)<-c("tAH","tAI","tDH","tDI","hospStay","ICUstay")
    actualTimes$type<-typeDF[,7]
    actualTimes$hosp<-hospChoice
    actualTimes$day<-day
    return(actualTimes)
  }
  
  allActualTimes<-Reduce(rbind,lapply(1:length(incidence),
                                      function(x){
                                        if(incidence[x]!=0)
                                        {pullPatients(incidence[x],x)}
                                        else {
                                          data.frame(tAH=1,tAI=1,tDH=1,tDI=1,hospStay=1,ICUstay=1,type="L",hosp=-1,day=x)
                                        }
                                      }
  )
  )
  
  allActualTimes<-allActualTimes[allActualTimes$hosp!=-1,]
  hospPats<-allActualTimes[(allActualTimes$type=="H"),]
  #End-of-stay is discharge from hospital:
  hospPats$endStay<-hospPats$tDH
  #Unless admission to ICu comes earlier:
  hospPats[hospPats$tAI<hospPats$tDH,"endStay"]<-hospPats[hospPats$tAI<hospPats$tDH,"tAI"]
  #admission is current day plus time to hospital admission:
  hospadm<-hospPats[,"day"]+hospPats[,"tAH"]-1
  hospdis<-hospadm+hospPats[,"endStay"]
  
  icuPats<-allActualTimes[(allActualTimes$type=="I")|
                            ((allActualTimes$type=="H")&(allActualTimes$tAI<allActualTimes$tDH)),]
  icuPats$endStay<-icuPats$tDI

  ICUadm<-icuPats[,"day"]+icuPats[,"tAH"]+icuPats[,"tAI"]-1
  ICUdis<-ICUadm+icuPats[,"endStay"]
  
  ICUChangeEvents<-data.frame(time=c(ICUadm,ICUdis),
                              event=c(rep(1,length(ICUadm)),rep(-1,length(ICUdis))),
                              hosp=c(icuPats$hosp,icuPats$hosp)
  )
  ICUChangeEvents<-ICUChangeEvents[order(ICUChangeEvents$time),]
  HospChangeEvents<-data.frame(time=c(hospadm,hospdis),
                               event=c(rep(1,length(hospadm)),rep(-1,length(hospdis))),
                               hosp=c(hospPats$hosp,hospPats$hosp))
  HospChangeEvents<-HospChangeEvents[order(HospChangeEvents$time),]
  return(list(hospChangeEvents=HospChangeEvents,ICUChangeEvents=ICUChangeEvents))
}

withinHospitalMultiRunvEarly<-function(numRuns=10,incidenceRuns,pL,pN,pI,alphaH,alphaI,deltaH,deltaI,pHosp=NULL){
  runNames=rep(unique(incidenceRuns$runNum),len=numRuns)
  inHospRuns<-lapply(runNames,
                     function(x){
                       print(x)
                       return(
                         createPatientsvEarly(incidenceRuns[incidenceRuns$runNum==x,"underlying"],pL,pN,pI,alphaH,alphaI,deltaH,deltaI,pHosp=NULL)
                       )
                     }
  )
  return(inHospRuns)
}
################### Late estimate model, local data based ###################

createPatientsvLate<-function(incidence,pL,pN,pI,alphaH,alphaI,alphaD,deltaH,deltaI,deltaD,pHosp=NULL){
  incidence<-round(incidence)
  pullPatients<-function(nPat,day){
    tAH<-sample(1:(length(alphaH)),nPat,prob = alphaH,replace = T)
    tAI<-sample(1:(length(alphaI)),nPat,prob = alphaI,replace = T)
    tAD<-sample(1:(length(alphaD)),nPat,prob = alphaD,replace = T)
    tDH<-sample(1:(length(deltaH)),nPat,prob = deltaH,replace = T)
    tDI<-sample(1:(length(deltaI)),nPat,prob = deltaI,replace = T)
    tDD<-sample(1:(length(deltaD)),nPat,prob = deltaD,replace = T)
    
    timeDF<-data.frame(tAH,tAI,tAD,tDH,tDI,tDD,hospStay<-1,ICUstay<-1,stepStay<-1,type<-1)
    
    pType<-       # H I D H I D H I D 
      sample(list(c(0,0,0,0,0,0,0,0,0,"L"),       #Stay at home (low-symptomatic)
                  c(1,1,1,1,1,1,1,1,1,"H"),       #First to hospital
                  c(1,0,1,0,1,1,1,1,1,"I")        #Immediate ICU admisison
      ),
      nPat,prob=c(pL, pN, pI),replace = T
      )
    #hospChoice<-sample(1:length(pHosp),nPat,prob=pHosp,replace = T)
    hospChoice<-rep(1,nPat)
    
    typeDF<-as.data.frame(t(matrix(unlist((pType)),nrow=10)))
    typeDF2<-sapply(typeDF[,1:9],function(x)as.numeric(as.character(x)))
    actualTimes<-(timeDF[,1:9]*typeDF2)
    colnames(actualTimes)<-c("tAH","tAI","tAD","tDH","tDI","tDD","hospStay","ICUstay","stepStay")
    actualTimes$type<-typeDF[,10]
    actualTimes$hosp<-hospChoice
    actualTimes$day<-day
    return(actualTimes)
  }
  
  allActualTimes<-Reduce(rbind,lapply(1:length(incidence),
                                      function(x){
                                        if(incidence[x]!=0)
                                        {pullPatients(incidence[x],x)}
                                        else {
                                          data.frame(tAH=1,tAI=1,tAD=1,tDH=1,tDI=1,tDD=1,hospStay=1,ICUstay=1,stepStay=1,type="L",hosp=-1,day=x)
                                        }
                                      }
  )
  )
  allActualTimes<-allActualTimes[allActualTimes$hosp!=-1,]

  hospPats<-allActualTimes[(allActualTimes$type=="H"),]
  #End-of-stay is discharge from hospital:
  hospPats$endStay<-hospPats$tDH
  #Unless admission to ICu comes earlier:
  hospPats[hospPats$tAI<hospPats$tDH,"endStay"]<-hospPats[hospPats$tAI<hospPats$tDH,"tAI"]
  #admission is current day plus time to hospital admission:
  hospadm<-hospPats[,"day"]+hospPats[,"tAH"]-1
  #double check: immediate admission
  #hospadm<-hospPats[,"day"]
  hospdis<-hospadm+hospPats[,"endStay"]
  
  icuPats<-allActualTimes[(allActualTimes$type=="I")|
      ((allActualTimes$type=="H")&(allActualTimes$tAI<allActualTimes$tDH)),]
  icuPats$endStay<-icuPats$tDI
  icuPats[icuPats$tAD<icuPats$tDI,"endStay"]<-icuPats[icuPats$tAD<icuPats$tDI,"tAD"]
  #print(icuPats)
  
  ICUadm<-icuPats[,"day"]+icuPats[,"tAH"]+icuPats[,"tAI"]-1
  ICUdis<-ICUadm+icuPats[,"endStay"]
  
  stepdownPats<-allActualTimes[
    ((allActualTimes$type=="I")&(allActualTimes$tAD<allActualTimes$tDI))|
      ((allActualTimes$type=="H")&(allActualTimes$tAD<allActualTimes$tDI)&(allActualTimes$tAI<allActualTimes$tDH)),]

  stepDownAdm<-stepdownPats[,"day"]+stepdownPats[,"tAH"]+stepdownPats[,"tAI"]+stepdownPats[,"tAD"]-1
  stepDownDis<-stepDownAdm+stepdownPats[,"tDD"]
  
  ICUChangeEvents<-data.frame(time=c(ICUadm,ICUdis),
                              event=c(rep(1,length(ICUadm)),rep(-1,length(ICUdis))),
                              hosp=c(icuPats$hosp,icuPats$hosp)
  )
  ICUChangeEvents<-ICUChangeEvents[order(ICUChangeEvents$time),]
  HospChangeEvents<-data.frame(time=c(hospadm,hospdis),
                               event=c(rep(1,length(hospadm)),rep(-1,length(hospdis))),
                               hosp=c(hospPats$hosp,hospPats$hosp))
  HospChangeEventsSD<-data.frame(time=c(stepDownAdm,stepDownDis),
                               event=c(rep(1,length(stepDownAdm)),rep(-1,length(stepDownDis))),
                               hosp=c(stepdownPats$hosp,stepdownPats$hosp))
  HospChangeEvents<-rbind(HospChangeEvents,HospChangeEventsSD)
  HospChangeEvents<-HospChangeEvents[order(HospChangeEvents$time),]
  return(list(hospChangeEvents=HospChangeEvents,ICUChangeEvents=ICUChangeEvents))
}

withinHospitalMultiRunvLate<-function(numRuns=10,incidenceRuns,pL,pN,pI,alphaH,alphaI,alphaD,deltaH,deltaI,deltaD,pHosp=NULL){
  runNames=rep(unique(incidenceRuns$runNum),len=numRuns)
  inHospRuns<-lapply(runNames,
                     function(x){
                       print(x)
                       return(
                         createPatientsvLate(incidenceRuns[incidenceRuns$runNum==x,"underlying"],pL,pN,pI,alphaH,alphaI,alphaD,deltaH,deltaI,deltaD,pHosp=NULL)
                       )
                     }
  )
  return(inHospRuns)
}

###########################Results extraction functions######################################
getHospTimeline<-function(runResults){
  htl<-t(Reduce(rbind,lapply(runResults,function(x){
    timeline=c()
    for(it in 1:350){
      timeline<-c(timeline,sum(c(0,x$hospChangeEvent[(x$hospChangeEvent$time<it),"event"])))
    }
    return(timeline)
  })))
  colnames(htl)<-paste0("R",1:ncol(htl))
  rownames(htl)<-paste0("D",1:nrow(htl))
  
  htl<-as.data.frame(htl)
  htl$time<-1:nrow(htl)
  return(htl)
}


getICUTimeline<-function(runResults){
  itl<-t(Reduce(rbind,lapply(runResults,function(x){
    timeline=c()
    for(it in 1:350){
      timeline<-c(timeline,sum(c(0,x$ICUChangeEvent[(x$ICUChangeEvent$time<it),"event"])))
    }
    return(timeline)
  })))
  colnames(itl)<-paste0("R",1:ncol(itl))
  rownames(itl)<-paste0("D",1:nrow(itl))
  itl<-as.data.frame(itl)
  itl$time<-1:nrow(itl)
  return(itl)
}

extractRunInfo<-function(runResults){
  hospTimeline<-Reduce(rbind,lapply(runResults,function(x){
    timeline=c()
    for(it in 1:550){
      timeline<-c(timeline,sum(c(0,x$hospChangeEvent[(x$hospChangeEvent$time<it),"event"])))
    }
    return(timeline)
  }))
  
  hospInciTimeline<-Reduce(rbind,lapply(runResults,function(x){
    timeline=c()
    for(it in 1:550){
      timeline<-c(timeline,length(c(x$hospChangeEvent[(x$hospChangeEvent$event==1)&(x$hospChangeEvent$time>it)&(x$hospChangeEvent$time<=(it+1)),"event"])))
    }
    return(timeline)
  }))
  
  ICUTimeline<-Reduce(rbind,lapply(runResults,function(x){
    timeline=c()
    for(it in 1:550){
      timeline<-c(timeline,sum(c(0,x$ICUChangeEvent[(x$ICUChangeEvent$time<it),"event"])))
    }
    return(timeline)
  }))
  
  ICUTIncitimeline<-Reduce(rbind,lapply(runResults,function(x){
    timeline=c()
    for(it in 1:550){
      timeline<-c(timeline,length(c(x$ICUChangeEvent[(x$ICUChangeEvent$event==1)&(x$ICUChangeEvent$time>it)&(x$ICUChangeEvent$time<=(it+1)),"event"])))
    }
    return(timeline)
  }))
  rrSelected=ceiling(runif(1)*nrow(hospTimeline))
  resultInfo<-t(rbind(1:550,
                      apply(hospTimeline,2,mean),
                      apply(hospTimeline,2,function(x)quantile(x,probs = c(0.05,0.25,0.5,0.75,0.95))),
                      apply(ICUTimeline,2,mean),
                      apply(ICUTimeline,2,function(x)quantile(x,probs = c(0.05,0.25,0.5,0.75,0.95))),
                      apply(ICUTimeline+hospTimeline,2,mean),
                      apply(ICUTimeline+hospTimeline,2,function(x)quantile(x,probs = c(0.05,0.25,0.5,0.75,0.95))),
                      apply(hospInciTimeline,2,mean),
                      apply(hospInciTimeline,2,function(x)quantile(x,probs = c(0.05,0.25,0.5,0.75,0.95))),
                      apply(ICUTIncitimeline,2,mean),
                      apply(ICUTIncitimeline,2,function(x)quantile(x,probs = c(0.05,0.25,0.5,0.75,0.95))),
                      hospTimeline[rrSelected,],
                      ICUTimeline[rrSelected,],
                      hospInciTimeline[rrSelected,],
                      ICUTIncitimeline[rrSelected,]
  ))
  colnames(resultInfo)<-
    c("time",
      "hospTotalMean","hospTotal005","hospTotal025","hospTotal050","hospTotal075","hospTotal095",
      "icuTotalMean","icuTotal005","icuTotal025","icuTotal050","icuTotal075","icuTotal095",
      "allTotalMean","allTotal005","allTotal025","allTotal050","allTotal075","allTotal095",
      "hospInciMean","hospInci005","hospInci025","hospInci050","hospInci075","hospInci095",
      "icuInciMean","icuInci005","icuInci025","icuInci050","icuInci075","icuInci095",
      "rrHosp","rrICU","rrInciHosp","rrInciICU"
    )
  resultInfo<-as.data.frame(resultInfo)
  return(resultInfo)
}

######################################Middle estimate running############################################
getBedPredictionFixedPara<-function(lx){
  
  startDhere<-max(mrResList[[lx]][!is.na(mrResList[[lx]]$ReEsti),"time"])+1
  print(startDhere)
  
  inHRuns<-withinHospitalMultiRunvEarly(numRuns=100,mrResList[[lx]],
                                        1 - ((1-pL_all)*(290/492)),#pL proportion at home
                                        ((290/492)* pN_all),#pN #Proportion to normal ward
                                        ((290/492)* pI_all),#pI #Proportion direct to ICU
                                        pexp(1:100,1/2)-pexp(0:99,1/2),
                                        pexp(1:100,alphaI_all)-pexp(0:99,alphaI_all), #alphaI #admission to ICU from hospita
                                        pexp(1:100,deltaH_all)-pexp(0:99,deltaH_all), #deltaH #Death / Discharge from hospital
                                        pexp(1:100,deltaI_all)-pexp(0:99,deltaI_all), #deltaI #Death / Discharge from ICU
                                        c(1)
  )
#  thePlot<-getBedShadedPlot(inHRuns,
#                            ylimits=c(0,175),
#                            bothObs = getObservedBeds(UKFparaList,as.Date("2020-04-28"),allObs=F),
#                            startD=min(mrResList[[lx]]$time),
#                            xlimits=c(as.Date("2020-03-01"),as.Date("2020-07-01")),
#                            cutDate=startDhere)
  return(list(
  #  thePlot=thePlot,
    rawResults=inHRuns
  ))
}

######################################Late estimate running############################################

getTimeLineFromCE<-function(ce,tlLength=NULL){
  if(is.null(tlLength)){
    tempDF<-data.frame(time=1:max(ce$time)+1,event=0,hosp=1)
    tempCE<-rbind(ce,tempDF)
  } else {
    tempCE<-rbind(ce[ce$time<=tlLength,],data.frame(time=1:tlLength,event=0,hosp=1))
  }
  tempCE<-tempCE[order(tempCE$time),]
  timelineTemp<-tempCE %>% 
    group_by(time) %>%
    summarise(
      n=n(),
      changes=sum(event)
    ) %<>% mutate(
      absolute=cumsum(changes)
    )
  return(timelineTemp)
}



################### Updated care path model, local data based ###################

getLOSdistr<-function(untilDate,numPatients=10000){
  mypatsTemp<-createPatientsvLate(c(numPatients,0),
                                  0,
                                  (1-directICUProp),
                                  directICUProp,
                                  c(1,0),#pexp(1:100,1/3)-pexp(0:99,1/3),
                                  selectModelOnDate(UKFparaList,untilDate,"transfGW1"),
                                  selectModelOnDate(UKFparaList,untilDate,"transfICU"),
                                  selectModelOnDate(UKFparaList,untilDate,"disGW1"),
                                  selectModelOnDate(UKFparaList,untilDate,"disICU"),
                                  selectModelOnDate(UKFparaList,untilDate,"disGW2")
  )
  
  timelineTemp<-getTimeLineFromCE(rbind(mypatsTemp$hospChangeEvents,mypatsTemp$ICUChangeEvents))
  resLOSdistr<-timelineTemp$absolute/numPatients
  return(resLOSdistr)
}


extraInterPolate<-function(x){
  lenX<-length(x)
  #set firs and last element as first and last known (!is.na) element
  if(is.na(x[1]))x[1]<-x[!is.na(x)][1]
  if(is.na(x[lenX]))x[lenX]<-x[!is.na(x)][sum(!is.na(x))]
  #Then, interpolate
  unknowns<-(which(is.na(x)))
  knowns<-(which(!is.na(x)))
  begins<-unlist(lapply(unknowns,function(y)(max(knowns[knowns<y]))))
  ends<-unlist(lapply(unknowns,function(y)(min(knowns[knowns>y]))))
  distright<-ends-unknowns
  distleft<-unknowns-begins
  interVals<-((x[begins]*distright)+(x[ends]*distleft))/(distright+distleft)
  x[unknowns]<-interVals
  return(x)
}

admRateEstimator<-function(inciDF,bedObservations){
  #inciDF is a data.frame with (at least) following variables:
  # - Date : date of reporting
  # - all  : number of reported cases, in all age classes 
  
  #bedObservations is a data.frame with columns:
  # - Date : date of observation
  # - ICU  : ICU beds occupied
  # - GW   : General ward beds occupied
  # - Hosp : All occupied beds in the hospital (GW+ICU)
  theDistr<-getLOSdistr(max(bedObservations$Date))
  # calculate the proportion admitted among all reported cases in the reversed LOS distribution for each "observed" day.
  pressureCalculator<-function(ddDate){
    thisObsDate=which(inciDF$Date==ddDate)
    calcLen<-min(thisObsDate-1,length(theDistr))
    totalAdmPressure=sum(inciDF[((thisObsDate-calcLen)):(thisObsDate-1),"all"]*theDistr[calcLen:1])
    return(totalAdmPressure)
  }
  #Begin calculating a week after the first case was reported.
  inclRecs<-(inciDF$Date>(min(inciDF$Date)+6))
  resDF<-data.frame(
    Date=inciDF[inclRecs,"Date"],
    admPress=unlist(lapply(inciDF[inclRecs,"Date"],function(x)pressureCalculator(x))),
    cases=inciDF[inclRecs,"all"],
    cLastweek=inciDF[inclRecs,"cLastweek"],
    thisPop=inciDF[inclRecs,"thisPop"],
    RKIcut=inciDF[inclRecs,"RKIcut"]
  )

  resDF=merge(resDF,bedObservations,by="Date",all=T)
  resDF$admProps<-resDF$Hosp/resDF$admPress
  resDF$allAvAR<-(unlist(lapply(resDF$Date,function(x)sum(resDF[resDF$Date<=x,"Hosp"],na.rm = T)/sum(resDF[resDF$Date<=x,"admPress"],na.rm = T))))
  
  calcThese<-resDF[!is.na(resDF$Hosp),]
  resDF[!is.na(resDF$Hosp),"movAvAR"]<-(unlist(lapply(calcThese$Date,function(x)
    sum(calcThese[(calcThese$Date<=(x+15))&(calcThese$Date>=(x-15)),"Hosp"],na.rm=T)/
      sum(calcThese[(calcThese$Date<=(x+15))&(calcThese$Date>=(x-15)),"admPress"],na.rm = T)
    )))
  
  resDF[is.na(resDF$Hosp),"movAvAR"]<-NA
  resDF$movAvAR<-extraInterPolate(resDF$movAvAR)
  
  resDF$movAvWeightAR<-(unlist(lapply(resDF$Date,function(x)
    sum(resDF[(resDF$Date<=(x+15))&(resDF$Date>=(x-15)),"Hosp"],na.rm=F)/
      sum(resDF[(resDF$Date<=(x+15))&(resDF$Date>=(x-15)),"admPress"],na.rm = F)
  )))
  #NOTE: what to do with NAs?
  
  return(resDF)
}

createPatientsVaryAdmRate<-function(incidence,pLv,pN,pI,alphaH,alphaI,alphaD,deltaH,deltaI,deltaD,pHosp=NULL){
  incidence<-round(incidence)
  if((pN+pI)!=1) stop("Normal and ICU parts != 1")
  pullPatients<-function(nPat,pL,day){
    tAH<-sample(1:(length(alphaH)),nPat,prob = alphaH,replace = T)
    tAI<-sample(1:(length(alphaI)),nPat,prob = alphaI,replace = T)
    tAD<-sample(1:(length(alphaD)),nPat,prob = alphaD,replace = T)
    tDH<-sample(1:(length(deltaH)),nPat,prob = deltaH,replace = T)
    tDI<-sample(1:(length(deltaI)),nPat,prob = deltaI,replace = T)
    tDD<-sample(1:(length(deltaD)),nPat,prob = deltaD,replace = T)

    timeDF<-data.frame(tAH,tAI,tAD,tDH,tDI,tDD,hospStay<-1,ICUstay<-1,stepStay<-1,type<-1)
    
    pType<-       # H I D H I D H I D 
      sample(list(c(0,0,0,0,0,0,0,0,0,"L"),       #Stay at home (low-symptomatic)
                  c(1,1,1,1,1,1,1,1,1,"H"),       #First to hospital
                  c(1,0,1,0,1,1,1,1,1,"I")        #Immediate ICU admisison
      ),
      nPat,prob=c(pL, pN*(1-pL), pI*(1-pL)),replace = T
      )
    #hospChoice<-sample(1:length(pHosp),nPat,prob=pHosp,replace = T)
    hospChoice<-rep(1,nPat)
    
    typeDF<-as.data.frame(t(matrix(unlist((pType)),nrow=10)))
    typeDF2<-sapply(typeDF[,1:9],function(x)as.numeric(as.character(x)))
    actualTimes<-(timeDF[,1:9]*typeDF2)
    colnames(actualTimes)<-c("tAH","tAI","tAD","tDH","tDI","tDD","hospStay","ICUstay","stepStay")
    actualTimes$type<-typeDF[,10]
    actualTimes$hosp<-hospChoice
    actualTimes$day<-day
    return(actualTimes)
  }
  
  allActualTimes<-Reduce(rbind,lapply(1:length(incidence),
                                      function(x){
                                        if(incidence[x]!=0)
                                          {
                                          if(x>length(pLv)){
                                            return(pullPatients(incidence[x],1-pLv[length(pLv)],x))
                                          } else {
                                            if(x<=length(pLv)) return(pullPatients(incidence[x],1-pLv[x],x))
                                           
                                          }
                                          }
                                        else {
                                          return(data.frame(tAH=1,tAI=1,tAD=1,tDH=1,tDI=1,tDD=1,hospStay=1,ICUstay=1,stepStay=1,type="L",hosp=-1,day=x))
                                        }
                                      }
  )
  )
  allActualTimes<-allActualTimes[allActualTimes$hosp!=-1,]
  
  hospPats<-allActualTimes[(allActualTimes$type=="H"),]
  #End-of-stay is discharge from hospital:
  hospPats$endStay<-hospPats$tDH
  #Unless admission to ICu comes earlier:
  hospPats[hospPats$tAI<hospPats$tDH,"endStay"]<-hospPats[hospPats$tAI<hospPats$tDH,"tAI"]
  #admission is current day plus time to hospital admission:
  hospadm<-hospPats[,"day"]+hospPats[,"tAH"]-1
  hospdis<-hospadm+hospPats[,"endStay"]
  
  icuPats<-allActualTimes[(allActualTimes$type=="I")|
                            ((allActualTimes$type=="H")&(allActualTimes$tAI<allActualTimes$tDH)),]
  icuPats$endStay<-icuPats$tDI
  icuPats[icuPats$tAD<icuPats$tDI,"endStay"]<-icuPats[icuPats$tAD<icuPats$tDI,"tAD"]

  ICUadm<-icuPats[,"day"]+icuPats[,"tAH"]+icuPats[,"tAI"]-1
  ICUdis<-ICUadm+icuPats[,"endStay"]
  
  stepdownPats<-allActualTimes[
    ((allActualTimes$type=="I")&(allActualTimes$tAD<allActualTimes$tDI))|
      ((allActualTimes$type=="H")&(allActualTimes$tAD<allActualTimes$tDI)&(allActualTimes$tAI<allActualTimes$tDH)),]
  
  stepDownAdm<-stepdownPats[,"day"]+stepdownPats[,"tAH"]+stepdownPats[,"tAI"]+stepdownPats[,"tAD"]-1
  stepDownDis<-stepDownAdm+stepdownPats[,"tDD"]
  
  ICUChangeEvents<-data.frame(time=c(ICUadm,ICUdis),
                              event=c(rep(1,length(ICUadm)),rep(-1,length(ICUdis))),
                              hosp=c(icuPats$hosp,icuPats$hosp)
  )
  ICUChangeEvents<-ICUChangeEvents[order(ICUChangeEvents$time),]
  HospChangeEvents<-data.frame(time=c(hospadm,hospdis),
                               event=c(rep(1,length(hospadm)),rep(-1,length(hospdis))),
                               hosp=c(hospPats$hosp,hospPats$hosp))
  HospChangeEventsSD<-data.frame(time=c(stepDownAdm,stepDownDis),
                                 event=c(rep(1,length(stepDownAdm)),rep(-1,length(stepDownDis))),
                                 hosp=c(stepdownPats$hosp,stepdownPats$hosp))
  HospChangeEvents<-rbind(HospChangeEvents,HospChangeEventsSD)
  HospChangeEvents<-HospChangeEvents[order(HospChangeEvents$time),]
  return(list(hospChangeEvents=HospChangeEvents,ICUChangeEvents=ICUChangeEvents))
}

withinHospitalMultiRunVaryAdmRate<-function(numRuns=10,incidenceRuns,admRateEsti,pN,pI,alphaH,alphaI,alphaD,deltaH,deltaI,deltaD,pHosp=NULL){

  runNames=rep(unique(incidenceRuns$runNum),len=numRuns)
  inHospRuns<-lapply(runNames,
                     function(x){
                       print(x)
                       return(
                         createPatientsVaryAdmRate(incidenceRuns[incidenceRuns$runNum==x,"underlying"],admRateEsti,pN,pI,alphaH,alphaI,alphaD,deltaH,deltaI,deltaD,pHosp=NULL)
                       )
                     }
  )
  return(inHospRuns)
}

getBedPrediction<-function(lx){
  #startDhere<-max(mrResList[[lx]][!is.na(mrResList[[lx]]$ReEsti),"time"])+1
  startDhere<-mrResList[[lx]][[1,"currentDay"]]
  print(startDhere)
  
  #If forecasting further than within-hospital parameters are available:
  if(startDhere>UKFparaList[[length(UKFparaList)]]$lastDate) startDhere<-UKFparaList[[length(UKFparaList)]]$lastDate
  
  admHospPart<-(max(cumsum(getRawAdm(UKFparaList,startDhere)$cases))/selectedKreise[selectedKreise$Date==startDhere,"all"])
  
  directICUProp<-UKFparaList[[getDateIndex(UKFparaList,startDhere)]]$propDirectICU
  
  inHRunWSD<-withinHospitalMultiRunvLate(numRuns=100,mrResList[[lx]],
                                         (1-admHospPart),#pL proportion at home
                                         (1-directICUProp)*admHospPart,#pN #Proportion to normal ward
                                         (directICUProp)*admHospPart,#pI #Proportion direct to ICU
                                         pexp(1:100,1/2)-pexp(0:99,1/2),
                                         getAlphaIDistr(UKFparaList,startDhere),
                                         getAlphaDDistr(UKFparaList,startDhere),
                                         getDeltaHDistr(UKFparaList,startDhere),
                                         getDeltaIDistr(UKFparaList,startDhere),
                                         getDeltaDDistr(UKFparaList,startDhere)
  )
  return(list(
    rawResults=inHRunWSD
  ))
}

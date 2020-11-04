
removeNegatives<-function(cCases,forward=FALSE){
  if(forward) {
    return(unlist(lapply(1:length(cCases),function(x) max(cCases[1:x]))))
  } else {
    return(unlist(lapply(1:length(cCases),function(x) min(cCases[x:length(cCases)]))))
  }
}

getIncCases<-function(cCases){
  cCases<-removeNegatives(cCases)
  cCases-c(0,cCases[-length(cCases)])
}

getMoveAv<-function(l,w){
  starts<-1:(length(l)-w)
  return(c(rep(0,floor(w/2)),unlist(lapply(starts,function(x)mean(l[x:(x+w)])))))
}

produceInfectorInfectee<-function(epicurve,serialDistr,currentDay){
  setToOne<-function(l){l/sum(l)}
  setFirstToOne<-function(l){l/(l[1])}
  
  epicurve<-epicurve[which(epicurve>0)[1]:length(epicurve)] #First entry needs to be a case. 
  
  if(length(serialDistr)<length(epicurve))
    serialDistr<-c(serialDistr,rep(0,length(epicurve)-length(serialDistr)))
  if(length(serialDistr)>length(epicurve))
    serialDistr<-serialDistr[1:length(epicurve)]
  obsProb<-setFirstToOne(rev(cumsum(serialDistr[-length(serialDistr)])))
  
  fGetInf<-function(d) sample(d-(1:(d-1)),epicurve[d],replace=TRUE,prob=serialDistr[1:(d-1)]*rev(epicurve[1:(d-1)]))
  genInfTimes<-unlist(lapply(2:length(epicurve),fGetInf))
  infTimes<-unlist(lapply(1:length(epicurve),
                          function(x)
                            length(genInfTimes[genInfTimes==x])))[-length(epicurve)]
  resInf<-unlist(lapply(1:length(obsProb),
                        function(x)if(infTimes[x]==0){
                          0
                        }else{
                          rnbinom(1,infTimes[x],c(obsProb[x]))
                        }
  ))
  #print(paste0("fPII: resInf ",length(resInf)," epicurve ",length(epicurve) ))
  infInfDF<-data.frame(
             time=currentDay-rev(1:length(resInf)),
             infectors=epicurve[-length(epicurve)],
             infectees=resInf+infTimes
             )

  return(infInfDF)
}

produceRe<-function(epicurve,serialDistr,currentDay){
  infInfDF<-produceInfectorInfectee(epicurve,serialDistr,currentDay)

  infInfDF$ReEsti<-(infInfDF$infectees/infInfDF$infectors)
  if(sum(!(is.na(infInfDF$ReEsti)==(infInfDF$infectors==0)))>0)stop("There is an NA error, please solve")
  infInfDF$ReEsti[is.na(infInfDF$ReEsti)]<-0
  
  return(infInfDF)
}

reconstructBackwards<-function(curve,delayDistr,underreporting=1){
  addedPrevious<-length(delayDistr)
  
  #distribute the reported cases over the previous days.
  actualCurve<-rep(0,addedPrevious+length(curve))
  for(it in 1:length(curve)) {
    nRep<-curve[it]
    if(nRep>0){
      delayCounts<-(table(sample(0:(length(delayDistr)-1),nRep,replace=TRUE,prob=delayDistr)))
      delayDF<-data.frame(delays=(as.numeric(names(delayCounts))),counts=(as.numeric(delayCounts)))
      actualCurve[addedPrevious+it-delayDF$delays]<-actualCurve[addedPrevious+it-delayDF$delays]+delayDF$counts
    }
  }
  #upscale to actual cases.
  myRDDcumsum<-(1-cumsum(rev(c(delayDistr[-1],rep(0,(1+length(actualCurve)-length(delayDistr)))))))/underreporting
  compUR<-rep(0,length(actualCurve))
  for(it in 1:length(actualCurve)) {
    if(actualCurve[it]>0&myRDDcumsum[it]>0){
      compUR[it]<-rnbinom(1,actualCurve[it],myRDDcumsum[it])
    } else {
      compUR[it]<-0
    }
  }
  return(compUR+actualCurve)
}

forwardInci<-function(R0,R0dev,popSize,serialDistr,runLength=150,caseStarter=c(0,1),underreporting=1){
 # print(caseStarter)
  startOffset=length(caseStarter)
  newCurve<-c(caseStarter,rep(0,runLength))
  devBeta<-R0dev/popSize
  
  runDay<-function(d){ #,beta){
    if(d<length(serialDistr)){
      pressure<-rev(serialDistr[1:d])*newCurve[1:d]
    } else {
      pressure<-rev(serialDistr)*newCurve[(d-(length(serialDistr)-1)):(d)]
    }
    if(popSize-sum(newCurve[1:d])<0) stop("Susceptible population under 0, shouldn't happen")
    if(devBeta[d]>1)devBeta[d]<-1
    
    perPersonInfPres<-(1-((1-devBeta[d])^sum(pressure)))
    
    if(perPersonInfPres>1){
      stop(paste0("Infection pressure larger than 1: ", perPersonInfPres))
    }
    newCurve[d+1]<-rbinom(1,popSize-sum(newCurve[1:d]),perPersonInfPres)
    return(newCurve)
  }
  
  for(it in (startOffset):(startOffset+runLength)) {
    newCurve<-runDay(it)#,beta)
  }
  return(newCurve)
}

convertCurveToReported<-function(curve,delayDistr,underreporting=1,exclstarter=c(0,1),replaceStarter=TRUE){
  reportCurve<-rep(0,length(curve))
  for(it in (1:length(curve))) {
    nRep<-rbinom(1,curve[it],1/underreporting)
    if(nRep>0){
      delayCounts<-(table(sample(0:(length(delayDistr)-1),nRep,replace=TRUE,prob=delayDistr)))
      delayDF<-data.frame(delays=(as.numeric(names(delayCounts))),counts=(as.numeric(delayCounts)))
      delayDF<-(delayDF[(delayDF$delays+it)<=length(curve),])
      reportCurve[it+delayDF$delays]<-reportCurve[it+delayDF$delays]+delayDF$counts
    }
  }
  if(replaceStarter)reportCurve[1:length(exclstarter)]<-exclstarter
  return(reportCurve)
}

singleInfInfCounts<-function(epicurve,startDate,currentDay,delayRep,serialinterval,runLength=150,underreporting=1,runNum=1){
  lastday<-((startDate+length(epicurve))-1) # last day in the epi curve
  #print(lastday)
  if(lastday!=(currentDay-1)){ #last day isn't yesterday
    if(lastday>(currentDay-1)){
      epicurve<-epicurve[1:(as.numeric(difftime(currentDay,startDate)))]
      warning("Dates didn't align, fixed by shortening epi curve.") 
    }
    if(lastday<(currentDay-1)) stop("Error: trying to estimate further than the epiCurve. Parameter currentDay is leading, please provide proper epi curve, or other currentDay parameter.")
  }
  lastday<-((startDate+length(epicurve))-1) # last day in the epi curve
  if(lastday!=(currentDay-1)) {
    print(lastday)
    stop("This didn't fix it") 
  }
  
  starterCurve<-reconstructBackwards(epicurve,delayRep,underreporting=underreporting)
  #main reconstructor runs here:
  resDF<-(produceInfectorInfectee(starterCurve,serialinterval,currentDay))
  #resDF$time<-resDF$day
  resDF$currentDay<-currentDay
  resDF$runNum<-runNum
  return(resDF)
}

multiInfInfCounts<-function(epicurve,startDate,currentDay,delayRep,serialinterval,runLength=150,underreporting=1,numRuns=10){
  res<-lapply(1:numRuns,function(x)singleInfInfCounts(epicurve,startDate,currentDay,delayRep,serialinterval,runLength=runLength,underreporting=underreporting,runNum=x))
  return(Reduce(rbind,res))
}


singleRtEsti<-function(epicurve,startDate,currentDay,delayRep,serialinterval,runLength=150,underreporting=1,runNum=1){
  lastday<-((startDate+length(epicurve))-1) # last day in the epi curve
  #print(lastday)
  if(lastday!=(currentDay-1)){ #last day isn't yesterday
    if(lastday>(currentDay-1)){
      epicurve<-epicurve[1:(as.numeric(difftime(currentDay,startDate)))]
      warning("Dates didn't align, fixed by shortening epi curve.") 
    }
    if(lastday<(currentDay-1)) stop("Error: trying to estimate further than the epiCurve. Parameter currentDay is leading, please provide proper epi curve, or other currentDay parameter.")
  }
  lastday<-((startDate+length(epicurve))-1) # last day in the epi curve
  if(lastday!=(currentDay-1)) {
    print(lastday)
    stop("This didn't fix it") 
  }

  starterCurve<-reconstructBackwards(epicurve,delayRep,underreporting=underreporting)
  ReEsti<-produceRe(starterCurve,serialinterval,currentDay)
  
  ReEsti$runNum<-runNum
  ReEsti$currentDay<-currentDay
  ReEsti$underlying=ReEsti$infectors
  return(ReEsti)
}


multiRtEsti<-function(epicurve,startDate,currentDay,delayRep,serialinterval,runLength=150,underreporting=1,numRuns=10){
  res<-lapply(1:numRuns,function(x)singleRtEsti(epicurve,startDate,currentDay,delayRep,serialinterval,runLength=runLength,underreporting=underreporting,runNum=x))
  res<-Reduce(rbind,res)
  
  meanRtTime<-as.data.frame(
    res[res$infectors!=0,]%>% 
      group_by(time)%>%
      summarise(meanRt=mean(ReEsti,na.rm = T))
  )
  rownames(meanRtTime)<-meanRtTime$time
  res$meanRt<-meanRtTime[as.character(res$time),"meanRt"]
  
  return(res)
}

singleRun<-function(epicurve,
                    population,
                    delayRep,
                    serialinterval,
                    currentDay,
                    runLength=150,
                    underreporting=1,
                    fixR=FALSE){
 
  starterCurve<-reconstructBackwards(epicurve,delayRep,underreporting=underreporting)

  ReEsti<-as.numeric(produceRe(starterCurve,serialinterval,currentDay)$ReEsti)
  ReEsti<-c(rep(0,-1+length(starterCurve)-length(ReEsti)),ReEsti,0)
  cStarter<-cumsum(starterCurve)
  ReEsti2<-ReEsti/(1-(cStarter/population))
  interc<-NA
  decl<-NA
  m<-NULL
 
   if(length(ReEsti)<7){#weighted average
    avRe<-sum(ReEsti2*starterCurve)/sum(starterCurve)
  } else {
    avRe<-sum(ReEsti2[(length(ReEsti2)-8):(length(ReEsti2)-2)]*starterCurve[(length(ReEsti2)-8):(length(ReEsti2)-2)])/sum(starterCurve[(length(ReEsti2)-8):(length(ReEsti2)-2)])
  }
  lastR<-ReEsti2[length(ReEsti2)-1]
  
  if(fixR) devR0<-rep(avRe,length(starterCurve)+runLength+1)
  if(!fixR){
    reDF<-subset(data.frame(
      x=1:(length(ReEsti2)),
      y=((ReEsti2)),
      w=starterCurve
    ),y!=0)
  
    lmRes<-lm(reDF$x ~ log(reDF$y))
  
    m<-tryCatch(
      nls(y ~ exp(((x-b)*z)),start=list(b=coef(lmRes)[1],z=1/coef(lmRes)[2]),weights = w,data = reDF),
      error=function(e) return(-1)
    )
    print(lmRes)
    print(m)
    if(!is.numeric(m)){
      lastR0<-as.numeric(exp(((max(reDF$x)-coef(m)["b"])*coef(m)["z"])))
      print(lastR0)
      mNew<-tryCatch(
        nls(y ~ a*exp(((x-max(x))*z)),start=list(a=lastR0,z=as.numeric(coef(m)["z"])),weights = w,data = reDF),
        error=function(e) return(-1)
      )
      print(mNew)
    } else {mNew<-(-1)}
    if(!is.numeric(mNew)){
      interc<-as.numeric(coef(mNew)["a"])
      decl<-as.numeric(coef(mNew)["z"])
    }

  if(!is.numeric(m)){
    if(!fixR){
      if((coef(m)["z"])>0){
        #stop("Increasing Rt detected")
        devR0<-rep(avRe,length(starterCurve)+runLength+1)
      } else {
        devR0<-exp(((1:(length(starterCurve)+runLength+1))-(coef(m)["b"]))*(coef(m)["z"]))
        print(paste0("mean Re last week ",avRe))
      }
    }
   } else {
    devR0<-rep(0,length(starterCurve)+runLength+1)
  }
  }
  print(avRe)
  print(starterCurve)

  newCurve<-forwardInci(avRe,devR0,population,serialinterval,runLength=runLength,caseStarter=starterCurve,underreporting = underreporting)
  newRepCurve<-convertCurveToReported(newCurve,delayRep,underreporting=underreporting,exclstarter=c(rep(0,length(delayRep)),epicurve),replaceStarter=F)
  data.frame(
    ReEsti=c(0,(ReEsti),rep(0,runLength)),
    ReEsti2=c(0,(ReEsti2),rep(0,runLength)),
    devR0=devR0,
    recentR0=avRe,
    lastR0=lastR,
    reportedRaw=newRepCurve,
    underlying=newCurve,
    decl=decl,
    fitR0inter=interc,
    reportedEpi=(c(rep(0,1+(length(starterCurve)-length(epicurve))),epicurve,newRepCurve[(length(starterCurve)+1):(length(newRepCurve)-1)]))
  )
}

singleRunFromHosp<-function(epicurve,population,delayRep,serialinterval,hospDelay,runLength=150,hospitalised=10,underreporting=1,fixR=FALSE){
  confEpiCurve<-reconstructBackwards(epicurve,hospDelay,underreporting=hospitalised)
  return(singleRun(confEpiCurve,population,delayRep,serialinterval,runLength,underreporting,fixR = fixR))
}

multiInciRun<-function(epicurve,startDate,currentDay,population,delayRep,serialinterval,numRuns=125,runLength=150,underreporting=1,fixR=FALSE) {
  if((length(epicurve)+startDate)!=(currentDay-1)){
    if((length(epicurve)+startDate)>(currentDay-1)){
      epicurve<-epicurve[1:(as.numeric(difftime(currentDay,startDate))-1)]
      warning("Dates didn't align, fixed by shortening epi curve.") 
    }
    if((length(epicurve)+startDate)<(currentDay-1)) stop("currentDay parameter is leading, please provide proper epi curve, or other currentDay parameter.")
  }
  if((length(epicurve)+startDate)!=(currentDay-1)) stop("This didn't fix it") 
    
  subset(Reduce(rbind,lapply(1:numRuns,function(x){
    runres<-singleRun(epicurve,
                      population,
                      delayRep,
                      serialinterval,
                      currentDay,
                      runLength,
                      underreporting=underreporting,
                      fixR=fixR)

    runres$runNum<-x
    runres$runPassed<-(sum(runres$devR0)!=0)
    runres$time<-startDate+(1:nrow(runres))-(length(delayRep))
    runres$currentDay<-currentDay
    runres$ReEsti[runres$underlying==0]<-NA
    runres$ReEsti2[runres$underlying==0]<-NA
    runres$startDate<-startDate
    return(runres)
  })),runPassed==TRUE)
}

readLandKreisNames<-function(){
  httr::set_config(httr::config(http_version = 0))
  RKIgotKreise<-GET("https://services7.arcgis.com/mOBPykOjAyBO2ZKk/arcgis/rest/services/RKI_Landkreisdaten/FeatureServer/0/query?where=1%3D1&outFields=*&returnGeometry=false&outSR=4326&f=json")
  tempLKdata<-as.data.frame(fromJSON(content(RKIgotKreise,"text"))$features$attributes)
  LKcodes=unique(tempLKdata[,c("GEN","AGS")])
  LKcodes<<-(LKcodes[!is.na(LKcodes$AGS),])
  return(LKcodes)
}

readRKIdataSingleLK<-function(LK){#Loading a single Landkreis

  if(!(LK %in% LKcodes$AGS)){
    found=FALSE
    if(LK %in% LKcodes$GEN){
      LK<-LKcodes[LKcodes$GEN==LK,"AGS"]
      found=T
    } else {
      detectedPattern<-which(str_detect(LKcodes$GEN,LK))
      if(length(detectedPattern)>1){
        stopstring<-paste0("Landkreis identifier ",LK," ambivalent, multiple matches: ")
        for(ii in 1:length(detectedPattern)) stopstring<-paste0(stopstring,"\n -",LKcodes[detectedPattern[[ii]],"GEN"])
        
        stop(stopstring)
      }
      if(length(detectedPattern)==1){
        LK<-LKcodes[detectedPattern[[1]],"AGS"]
        found=T
      }
    }
    if(!found) stop(paste0("Didn't find identifier: ",LK)) 
  }
  RKIgot<-GET(paste0("https://services7.arcgis.com/mOBPykOjAyBO2ZKk/arcgis/rest/services/RKI_COVID19/FeatureServer/0/query?where=IdLandkreis%20%3D%20'",
                     LK,
                     "'&outFields=*&outSR=4326&f=json"))
  tempRKIdataJSON<-fromJSON(content(RKIgot,"text"))
  tempDF<-as.data.frame(tempRKIdataJSON$features$attributes)
  tempDF$Meldedatum<-as.Date("1970-01-01")+(tempDF$Meldedatum/(3600*24*1000))
  tempDF$Refdatum<-as.Date("1970-01-01")+(tempDF$Refdatum/(3600*24*1000))
  tempDF$Datenstand<-as.Date(substr(tempDF$Datenstand,0,10),format="%d.%m.%Y")
  #print(as.Date("1970-01-01")+max(RKIinci$Meldedatum)/(3600*24*1000))
  return(tempDF)
}

readRKIdatamultiLK<-function(LKs){#Loading multiple Landkreis
  Reduce(rbind,lapply(LKs,readRKIdataSingleLK))
}

readRKIpopLK<-function(){
    RKIgotKreise<-GET("https://services7.arcgis.com/mOBPykOjAyBO2ZKk/arcgis/rest/services/RKI_Landkreisdaten/FeatureServer/0/query?where=1%3D1&outFields=*&returnGeometry=false&outSR=4326&f=json")
    return(fromJSON(RKIgotKreise))
}

getInciForKreiseRKI<-function(l,popLK=NULL, fromDate=NULL,untilDate=NULL,loadData=FALSE){
  if(loadData){
    #https://npgeo-corona-npgeo-de.hub.arcgis.com/datasets/917fc37a709542548cc3be077a786c17_0/data?geometry=-31.470%2C46.269%2C52.378%2C55.886&selectedAttribute=cases7_per_100k
    #if(is.null(popLK)) popLK<-read.csv("https://opendata.arcgis.com/datasets/917fc37a709542548cc3be077a786c17_0.csv")
    if(is.null(popLK)) popLK<-read.csv("Resources/Kreisgrenzen_2017_mit_Einwohnerzahl.csv")
    
  } else {
    if(is.null(popLK))stop("Please provide population data (popLK) and/or case data (dataRKI). Alternatively, set loadData=True to load data from RKI.")
  }
  thisSet<-readRKIdatamultiLK(l)
  totPop<-sum(unlist(lapply(unique(thisSet$IdLandkreis), function(x) ((popLK[!is.na(popLK$AGS),"EWZ"][popLK[!is.na(popLK$AGS),"AGS"]==as.numeric(x)])))))

  RKIinci<-as.data.frame(
    thisSet%>% 
      group_by(Meldedatum) %>%
      summarise(cases=sum(AnzahlFall))
  )
  if(is.null(untilDate)) untilDate<-max(thisSet$Datenstand)
  if(is.null(fromDate)) fromDate<-min(RKIinci$Meldedatum)
  
  RKIinciTempList<-Reduce(rbind,lapply(1:as.numeric(untilDate-fromDate),
                                       function(x) return(data.frame(Date=(fromDate+x),
                                                                     all=sum(RKIinci[as.numeric(as.Date(RKIinci$Meldedatum,origin=as.Date("1900-01-01"))-(fromDate))==x,"cases"])
                                       )
                                       )
  ))
  
  RKIinciTempList$cLastweek<-c(rep(0,6),unlist(lapply(7:nrow(RKIinciTempList),function(x)sum(RKIinciTempList[x-0:6,"all"]))))
  RKIinciTempList$RKIcut<-100000*(RKIinciTempList$cLastweek/totPop)
  RKIinciTempList$thisPop<-totPop
  return(RKIinciTempList)
}


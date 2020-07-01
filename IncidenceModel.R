getIncCases<-function(cumCases){cumCases-c(0,cumCases[-length(cumCases)])}
getMoveAv<-function(l,w){
  starts<-1:(length(l)-w)
  return(c(rep(0,floor(w/2)),unlist(lapply(starts,function(x)mean(l[x:(x+w)])))))
}

produceRe<-function(epicurve,serialDistr){
  setToOne<-function(l){l/sum(l)}
  setFirstToOne<-function(l){l/(l[1])}
  epicurve<-epicurve[which(epicurve>0)[1]:length(epicurve)]
  
  if(length(serialDistr)<length(epicurve))
    serialDistr<-c(serialDistr,rep(0,length(epicurve)-length(serialDistr)))
  if(length(serialDistr)>length(epicurve))
    serialDistr<-serialDistr[1:length(epicurve)]
  
  obsProb<-setFirstToOne(rev(cumsum(serialDistr[-length(serialDistr)])))
  
  fGetInf<-function(d) sample(d-(1:(d-1)),epicurve[d],replace=TRUE,prob=serialDistr[1:(d-1)]*rev(epicurve[1:(d-1)]))
  genInfTimes<-unlist(lapply(2:length(epicurve),fGetInf))
  infTimes<-unlist(lapply(1:length(epicurve),function(x)length(genInfTimes[genInfTimes==x])))[-length(epicurve)]
  resInf<-unlist(lapply(1:length(obsProb),function(x)if(infTimes[x]==0){
    0
  }else{
    rnbinom(1,infTimes[x],obsProb[x])
  }
  ))
  
  esti<-as.data.frame(t((resInf+infTimes)/epicurve[-length(epicurve)]))
  
  if(sum(!(is.na(esti)==(epicurve[-length(epicurve)]==0)))>0)stop("There is an NA error, please solve")
  esti[is.na(esti)]<-0
  
  return(esti)
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
  print(caseStarter)
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

singleRun<-function(epicurve,population,delayRep,serialinterval,runLength=150,underreporting=1){
  starterCurve<-reconstructBackwards(epicurve,delayRep,underreporting=underreporting)
  ReEsti<-as.numeric(produceRe(starterCurve,serialinterval))
  ReEsti<-c(rep(0,-1+length(starterCurve)-length(ReEsti)),ReEsti,0)
  cStarter<-cumsum(starterCurve)
  ReEsti2<-ReEsti/(1-(cStarter/population))
  interc<-NA
  decl<-NA
  
  if(length(ReEsti)<7){#weighted average
    avRe<-sum(ReEsti2*starterCurve)/sum(starterCurve)
  } else {
    avRe<-sum(ReEsti2[(length(ReEsti2)-7):(length(ReEsti2)-1)]*starterCurve[(length(ReEsti2)-7):(length(ReEsti2)-1)])/sum(starterCurve[(length(ReEsti2)-7):(length(ReEsti2)-1)])
  }
  lastR<-ReEsti2[length(ReEsti2)-1]
  
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
    if((coef(m)["z"])>0){
      #stop("Increasing Rt detected")
      devR0<-rep(avRe,length(starterCurve)+runLength+1)
    } else {
      devR0<-exp(((1:(length(starterCurve)+runLength+1))-(coef(m)["b"]))*(coef(m)["z"]))
      print(paste0("mean Re last week ",avRe))
    }
    newCurve<-forwardInci(avRe,devR0,population,serialinterval,runLength=runLength,caseStarter=starterCurve,underreporting = underreporting)
    newRepCurve<-convertCurveToReported(newCurve,delayRep,underreporting=underreporting,exclstarter=c(rep(0,length(delayRep)),epicurve),replaceStarter=F)
  } else {
    devR0<-0
    newCurve<-c(starterCurve,rep(-1,runLength+1))
    newRepCurve<-c(starterCurve,rep(-1,runLength+1))
  }
  
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

singleRunFromHosp<-function(epicurve,population,delayRep,serialinterval,hospDelay,runLength=150,hospitalised=10,underreporting=1){
  confEpiCurve<-reconstructBackwards(epicurve,hospDelay,underreporting=hospitalised)
  return(singleRun(confEpiCurve,population,delayRep,serialinterval,runLength,underreporting))
}

multiInciRun<-function(epicurve,startDate,currentDay,population,delayRep,serialinterval,numRuns=125,runLength=150,underreporting=1) {
  if((length(epicurve)+startDate)!=(currentDay-1)){
    if((length(epicurve)+startDate)>(currentDay-1)){
      epicurve<-epicurve[1:(as.numeric(difftime(currentDay,startDate))-1)]
      warning("Dates didn't align, fixed by shortening epi curve.") 
    }
    if((length(epicurve)+startDate)<(currentDay-1)) stop("currentDay parameter is leading, please provide proper epi curve, or other currentDay parameter.")
  }
  if((length(epicurve)+startDate)!=(currentDay-1)) stop("This didn't fix it") 
    
  subset(Reduce(rbind,lapply(1:numRuns,function(x){
    runres<-singleRun(epicurve,population,delayRep,serialinterval,runLength,underreporting=underreporting)
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

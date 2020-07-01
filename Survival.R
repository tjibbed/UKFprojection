getSurvParas<-function(ourPatients, lastDate){
  
  lastAdm<-ourPatients %>%
    group_by(PIZ) %>%
    summarise(lastFall=max(FALLID))%>%
    as.data.frame
  ourPatients<-ourPatients[ourPatients$FALLID %in% lastAdm$lastFall,]

  ourPatients$adm<-as.POSIXct(ourPatients$AUFENTHALTSBEGINN,format="%Y-%m-%d %H:%M:%S")
  ourPatients$dis<-as.POSIXct(ourPatients$AUFENTHALTSENDE,format="%Y-%m-%d %H:%M:%S")+minutes(1)
  
  ourPatients<-ourPatients[ourPatients$adm<lastDate,]
  ourPatients[ourPatients$dis>lastDate,"ENTLASSCODE"]<-0
  ourPatients[ourPatients$dis>lastDate,"dis"]<-as.POSIXct(paste0(lastDate," 00:00:00"),format="%Y-%m-%d %H:%M:%S")
  ourPatients$los<-difftime(ourPatients$dis,ourPatients$adm,units="days")

  ICUs=c(9283633, 9231625,   9282114, 9211616,9283625, 9211624, 9311637)
  UNZ = 9291806 #British: A&E, USA: ER

  ourPatients$onICU<-0
  ourPatients[ourPatients$KOSTENSTELLENNR %in% ICUs,"onICU"]<-1
  ourPatients[ourPatients$KOSTENSTELLENNR==9291806,"onICU"]<-2

  rawStays<-ourPatients %>%   
    group_by(PIZ) %>%
    summarise(los=sum(los),
            start=min(adm),
            end=max(dis),
            dischargeCodeICU=min(ENTLASSCODE)
    )

  ICUStayFrame<-ourPatients %>%   
    subset(onICU==1) %>%
    group_by(PIZ) %>%
    summarise(losICU=sum(los),
            startICU=min(adm),
            endICU=max(dis),
            dischargeCodeICU=min(ENTLASSCODE)
            )%>%
    as.data.frame()

  GWStayFrame<-ourPatients %>%   
    subset(onICU==0) %>%
    group_by(PIZ) %>%
    summarise(losHosp=sum(los),
            startHosp=min(adm),
            endHosp=max(dis),
          #  losHosp2=max(dis)-min(adm),
            dischargeCodeHosp=min(ENTLASSCODE)
    )%>%
    as.data.frame()

  allStayFrame<-merge(ICUStayFrame,GWStayFrame,by="PIZ",all=T)
  allStayFrame$startICU<-allStayFrame$endICU-allStayFrame$losICU
  allStayFrame[is.na(allStayFrame$startICU),"startICU"]<-allStayFrame[is.na(allStayFrame$startICU),"endHosp"]
  allStayFrame$GW1los<-round(difftime(allStayFrame$startICU,allStayFrame$startHosp,units = "days"),digits = 6)
  allStayFrame$ICUlos<-round(allStayFrame$losICU,digits = 6)
  allStayFrame$GW2los<-round(allStayFrame$losHosp-allStayFrame$GW1los,digits = 6)

  allStayFrame[is.na(allStayFrame$GW1los),"GW1los"]<-0
  allStayFrame[(allStayFrame$GW1los<0),"GW1los"]<-0
  allStayFrame[is.na(allStayFrame$GW2los),"GW2los"]<-0
  allStayFrame[is.na(allStayFrame$losICU),"losICU"]<-0

  #patients in GW1
  allStayFrame$GW1state<-(-1)
  allStayFrame[(allStayFrame$losICU==0)&(allStayFrame$GW1los>0)&(allStayFrame$dischargeCodeHosp==0),"GW1state"]<-0 #still on GW1 (Not discharged. Hasn't been on ICU (yet))
  allStayFrame[(allStayFrame$dischargeCodeHosp!=0)&(allStayFrame$losICU==0),"GW1state"]<-1 #discharged from GW1 (Discharged. Hasn't been on ICU)
  allStayFrame[(allStayFrame$losICU>0)&(allStayFrame$GW1los>0),"GW1state"]<-2 #moved on to ICU (Positive LOS on GW1 and ICU means moved away)

  #patients on ICU
  allStayFrame$ICUstate<-(-1)
  allStayFrame[(allStayFrame$losICU>0)&allStayFrame$dischargeCodeICU==0,"ICUstate"]<-0 #still on ICU (clear, discharge=0 means still here)
  allStayFrame[(allStayFrame$losICU>0)&allStayFrame$dischargeCodeICU>0&allStayFrame$GW2los==0,"ICUstate"]<-1 #discharged from ICU (Discharged, los at GW2 is zero, complete discharge)
  allStayFrame[(allStayFrame$losICU>0)&allStayFrame$dischargeCodeICU>0&allStayFrame$GW2los>0,"ICUstate"]<-2 #moved on to GW2 (Discharged, los at GW2 is positive, moved on)

  #patients in GW2
  allStayFrame$GW2state<-(-1)
  allStayFrame[(allStayFrame$losICU>0)&(allStayFrame$dischargeCodeHosp==0)&(allStayFrame$GW2los>0),"GW2state"]<-0 #still on GW2 (DischargeCode=0, still in hospital, +ve LOS in GW2, so is there)
  allStayFrame[(allStayFrame$losICU>0)&(allStayFrame$dischargeCodeHosp!=0)&(allStayFrame$GW2los>0),"GW2state"]<-1 #Discharged on GW2 (DischargeCode!=0, not in hospital, +ve LOS in GW2, so was there)


  GW1stays<-allStayFrame[allStayFrame$GW1state!=-1,c("GW1los","GW1state")]
  GW2stays<-allStayFrame[allStayFrame$GW2state!=-1,c("GW2los","GW2state")]
  ICUstays<-allStayFrame[allStayFrame$ICUstate!=-1,c("ICUlos","ICUstate")]

  if(nrow(ICUstays)>2){
    allICU_exp<-flexsurvreg(formula = Surv(ICUlos, ICUstate != 0) ~ 1, data = ICUstays, 
                        dist = "exp")
    allICU_wb<-flexsurvreg(formula = Surv(ICUlos, ICUstate != 0) ~ 1, data = ICUstays, 
                        dist = "weibull")
    disICU_exp<-flexsurvreg(formula = Surv(ICUlos, ICUstate == 1) ~ 1, data = ICUstays, 
                        dist = "exp")
    disICU_wb<-flexsurvreg(formula = Surv(ICUlos, ICUstate == 1) ~ 1, data = ICUstays, 
                        dist = "weibull")
    transfICU_exp<-flexsurvreg(formula = Surv(ICUlos, ICUstate == 2) ~ 1, data = ICUstays, 
                        dist = "exp")
    transfICU_wb<-flexsurvreg(formula = Surv(ICUlos, ICUstate == 2) ~ 1, data = ICUstays, 
                        dist = "weibull")
    allGW1_exp<-flexsurvreg(formula = Surv(GW1los, GW1state != 0) ~ 1, data = GW1stays, 
                        dist = "exp")
    disGW1_exp<-flexsurvreg(formula = Surv(GW1los, GW1state == 1) ~ 1, data = GW1stays, 
                        dist = "exp")
    disGW1_wb<-flexsurvreg(formula = Surv(GW1los, GW1state == 1) ~ 1, data = GW1stays, 
                       dist = "weibull")
    transfGW1_exp<-flexsurvreg(formula = Surv(GW1los, GW1state == 2) ~ 1, data = GW1stays, 
                       dist = "exp")
    transfGW1_wb<-flexsurvreg(formula = Surv(GW1los, GW1state == 2) ~ 1, data = GW1stays, 
                       dist = "weibull")
    disGW2_exp<-flexsurvreg(formula = Surv(GW2los, GW2state == 1) ~ 1, data = GW2stays, 
                       dist = "exp")
    disGW2_wb<-flexsurvreg(formula = Surv(GW2los, GW2state == 1) ~ 1, data = GW2stays, 
                       dist = "weibull")
    } else {
      allICU_wb=NULL
      allICU_exp=NULL
      allGW1_exp=NULL
      disGW1_exp=NULL
      disGW1_wb=NULL
      transfGW1_exp=NULL
      transfGW1_wb=NULL
      disICU_exp=NULL
      disICU_wb=NULL
      transfICU_exp=NULL
      transfICU_wb=NULL
      disGW2_exp=NULL
      disGW2_wb=NULL
    }
  nDirectICUPats<-sum(allStayFrame[(allStayFrame$ICUstate!=-1),"GW1state"]==-1)
  nAllPats<-nrow(allStayFrame)
  return(
    list(
      allICU_wb=allICU_wb,
      allGW1_exp=allGW1_exp,
      allICU_exp=allICU_exp,
      GW1stays=GW1stays,
      GW2stays=GW2stays,
      ICUstays=ICUstays,
      rawAdm=nrow(subset(rawStays,(start<=(lastDate))&(start>(lastDate-days(1))))),
      lastDate=lastDate,
      nAllPats=nAllPats,
      propDirectICU=nDirectICUPats/nAllPats,
      nICUpats=sum(subset(ourPatients,onICU==1)$ENTLASSCODE==0),
      nGWpats=sum(subset(ourPatients,onICU==0)$ENTLASSCODE==0),
      disGW1_exp=disGW1_exp,
      disGW1_wb=disGW1_wb,
      transfGW1_exp=transfGW1_exp,
      transfGW1_wb=transfGW1_wb,
      disICU_exp=disICU_exp,
      disICU_wb=disICU_wb,
      transfICU_exp=transfICU_exp,
      transfICU_wb=transfICU_wb,
      disGW2_exp=disGW2_exp,
      disGW2_wb=disGW2_wb
    )
  )
}

selectModel<-function(expMod,wbMod){
  
  if((expMod$AIC)<(wbMod$AIC)){
    print(paste0("Exponential selected. Relative likelihood: ",exp((wbMod$AIC- expMod$AIC)/2)))
    return(
      pexp(1:200,expMod$res[[1]])-pexp(0:199,expMod$res[[1]])
    )
  } else {
    print(paste0("Weibull selected. Relative likelihood: ",exp((expMod$AIC- wbMod$AIC)/2)))
    return(
      pweibull(1:200,shape=wbMod$res[[1]],scale=wbMod$res[[2]])-
        pweibull(0:199,shape=wbMod$res[[1]],scale=wbMod$res[[2]])
    )
  }
}

getDeltaHDistr<-function(paraList,thisDate){
  listEntry<-which(unlist(lapply(1:length(paraList),function(x){paraList[[x]]$lastDate==thisDate})))[[1]]
  disExp<-paraList[[listEntry]]$disGW1_exp
  disWB<-paraList[[listEntry]]$disGW1_wb
  return(selectModel(disExp,disWB))
}
getDeltaIDistr<-function(paraList,thisDate){
  listEntry<-which(unlist(lapply(1:length(paraList),function(x){paraList[[x]]$lastDate==thisDate})))[[1]]
  disExp<-paraList[[listEntry]]$disICU_exp
  disWB<-paraList[[listEntry]]$disICU_wb
  return(selectModel(disExp,disWB))
}
getDeltaDDistr<-function(paraList,thisDate){
  listEntry<-which(unlist(lapply(1:length(paraList),function(x){paraList[[x]]$lastDate==thisDate})))[[1]]
  disExp<-paraList[[listEntry]]$disGW2_exp
  disWB<-paraList[[listEntry]]$disGW2_wb
  return(selectModel(disExp,disWB))
}
getDeltaHDistr<-function(paraList,thisDate){
  listEntry<-which(unlist(lapply(1:length(paraList),function(x){paraList[[x]]$lastDate==thisDate})))[[1]]
  disExp<-paraList[[listEntry]]$disGW1_exp
  disWB<-paraList[[listEntry]]$disGW1_wb
  return(selectModel(disExp,disWB))
}
getDeltaIDistr<-function(paraList,thisDate){
  listEntry<-which(unlist(lapply(1:length(paraList),function(x){paraList[[x]]$lastDate==thisDate})))[[1]]
  disExp<-paraList[[listEntry]]$disICU_exp
  disWB<-paraList[[listEntry]]$disICU_wb
  return(selectModel(disExp,disWB))
}
getAlphaIDistr<-function(paraList,thisDate){
  listEntry<-which(unlist(lapply(1:length(paraList),function(x){paraList[[x]]$lastDate==thisDate})))[[1]]
  disExp<-paraList[[listEntry]]$transfGW1_exp
  disWB<-paraList[[listEntry]]$transfGW1_wb
  return(selectModel(disExp,disWB))
}
getAlphaDDistr<-function(paraList,thisDate){
  listEntry<-which(unlist(lapply(1:length(paraList),function(x){paraList[[x]]$lastDate==thisDate})))[[1]]
  disExp<-paraList[[listEntry]]$transfICU_exp
  disWB<-paraList[[listEntry]]$transfICU_wb
  return(selectModel(disExp,disWB))
}

getDateIndex<-function(paraList,thisDate){
  listEntry<-which(unlist(lapply(1:length(paraList),function(x){paraList[[x]]$lastDate==thisDate})))[[1]]
  return(listEntry)
}

getRawAdm<-function(paraList,thisDate){
  listEntry<-which(unlist(lapply(1:length(paraList),function(x){paraList[[x]]$lastDate==thisDate})))[[1]]
  return(data.frame(
    time=do.call("c",lapply(1:listEntry,function(x)(paraList[[x]]$lastDate))),
    cases=(unlist(lapply(1:listEntry,function(x){(paraList[[x]]$rawAdm)})))
  )
  )
}

getObservedBeds<-function(paraList,thisDate,allObs=FALSE){
  if(allObs){
    listEntry=length(paraList)
  }else{
    listEntry<-which(unlist(lapply(1:length(paraList),function(x){paraList[[x]]$lastDate==thisDate})))[[1]]  
  }
  print(listEntry)
  timelineICU<-unlist(lapply(1:listEntry,function(x)paraList[[x]]$nICUpats))
  timelineGW<-unlist(lapply(1:listEntry,function(x)paraList[[x]]$nGWpats))
  dates<-do.call("c",lapply(1:listEntry,function(x)(paraList[[x]]$lastDate)))
  return(data.frame(
    time=dates,
    timelineICU=timelineICU,
    timelineGW=timelineGW))
}
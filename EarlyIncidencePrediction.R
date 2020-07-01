############## Import functions ############## 

# Italian regions:
pcmdpcImportRegions<-function(){
  ITRegurl <- RCurl::getURL("https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-regioni/dpc-covid19-ita-regioni.csv")
  ITRegdata <- read.csv(text = ITRegurl)
  ItRegionPop<-read.csv("Resources/ItalianRegionPop.tsv",sep = "\t")
  ITRegdata <- merge(ITRegdata,ItRegionPop,by.x="codice_regione",by.y="Region.Code")
  ITRegdata$data<-(as.Date(as.character(ITRegdata$data)))## Making sure date is formatted similarly for all datasets
  ITRegdata$CasesPerCapita<-ITRegdata$totale_casi/ITRegdata$Total.Population
  ITRegdata<-as.data.frame(ITRegdata)
  ITRegdata<-ITRegdata[order(ITRegdata$codice_regione,ITRegdata$data),]
  ITRegdata<-ITRegdata[,c("denominazione_regione","data","totale_casi","Total.Population","CasesPerCapita")]
  colnames(ITRegdata)<-c("region","Date","Cases","Population","CasesPerCapita")
  return(ITRegdata)
}

# Italian provinces:
pcmdpcImportProvinces<-function(){
  ITProvurl <- RCurl::getURL("https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-province/dpc-covid19-ita-province.csv")
  ITProvdata <- read.csv(text = ITProvurl)
  ItProvPop<-read.csv("Resources/ItalianProvincePop.tsv",sep = "\t")
  ITProvdata <- merge(ITProvdata,ItProvPop,by.x="codice_provincia",by.y="Province.Code")
  ITProvdata$data<-(as.Date(as.character(ITProvdata$data)))  ## Making sure date is formatted equally for all datasets
  ITProvdata$CasesPerCapita<-ITProvdata$totale_casi/ITProvdata$Total.Population
  ITProvdata<-as.data.frame(ITProvdata)
  ITProvdata<-ITProvdata[order(ITProvdata$codice_provincia,ITProvdata$data),]
  #Export only useful columns
  ITProvdata<-ITProvdata[,c("denominazione_provincia","data","totale_casi","Total.Population","CasesPerCapita")]
  colnames(ITProvdata)<-c("region","Date","Cases","Population","CasesPerCapita")
  return(ITProvdata)
}

#Countries
JHUimportJustCountries<-function(){
  #Population numbers:
  #https://population.un.org/wpp/Download/Standard/Population/
  JHUurl <- RCurl::getURL("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")
  JHUdata <- read.csv(text = JHUurl)
  justData<-JHUdata[,5:ncol(JHUdata)]
  colnames(JHUdata)<-c(colnames(JHUdata[,1:4]),((substr(colnames(justData),2,100))))
  setDT(JHUdata)
  JHUdata<-melt(JHUdata,id.vars=(1:4))
  JHUdataMelt<-JHUdata[, sum(value),by=list(Country.Region,variable)]
  colnames(JHUdataMelt)<-c("Country.Region","Date","Cases")
  
  JHUdataMelt$Date<-as.Date(JHUdataMelt$Date,format="%m.%d.%y")
  
  countryPops<-read.csv("Resources/worldPopulation.csv")
  colnames(countryPops)<-c("Index","Variant","Country","Notes","Country.code","Type","Parent.code",paste0("Y",1950:2020))
  JHUdataMelt<-merge(JHUdataMelt,countryPops[,c("Country","Y2020")],by.x = "Country.Region",by.y="Country",all.x=TRUE)
  colnames(JHUdataMelt)<-c("region","Date","Cases","Population")
  JHUdataMelt$Population<-1000*as.numeric(str_remove_all(as.character(JHUdataMelt$Population),"[ \t]"))
  JHUdataMelt$CasesPerCapita<-JHUdataMelt$Cases/JHUdataMelt$Population
  return(JHUdataMelt)
}

combineRegionsDB<-function(dataSetList,selectNames=NULL,minCases=0){
  allDB<-Reduce(rbind,dataSetList)
  if(is.null(selectNames)) selectNames<-unique(allDB$region)
  subDB<-allDB[(allDB$region %in% selectNames),]
  subDB<-as.data.frame(subDB[order(subDB$Date),])
  subDB$Day<-as.numeric(difftime(subDB$Date,as.Date("2020-01-01")))
  
  excludedRegions<-unique(subDB[is.na(subDB$CasesPerCapita),"region"])
  if(length(excludedRegions)>0){
    print(paste0("Excluded: ",excludedRegions))
  }
  subDB<-subDB[!is.na(subDB$CasesPerCapita),]
  return(subDB[subDB$Cases>=minCases,])
}

################# functions for calculating the delay ################# 

getLagTimeCurves<-function(inputDF,region1,region2,DEBUG=FALSE){
  R1Curve<-subset(inputDF,region==region1)
  R2Curve<-subset(inputDF,region==region2)
  
  if(DEBUG){
    print(R1Curve)
    print(R2Curve)
  }
  
  if(R1Curve[nrow(R1Curve),"CasesPerCapita"]>R2Curve[nrow(R2Curve),"CasesPerCapita"]){
    testCurve<-R2Curve
    refCurve<-R1Curve
    testRegion=region2
    refRegion=region1
    reverseRef<-FALSE
  }else {
    if(DEBUG)print("REVERSE")
    refCurve<-R2Curve
    testCurve<-R1Curve
    reverseRef<-TRUE
    testRegion=region1
    refRegion=region2
  }
  refPoints<-subset(refCurve,CasesPerCapita>testCurve[nrow(testCurve),"CasesPerCapita"])
  
  if(DEBUG)print(refPoints)
  
  if(nrow(refPoints)==0){
    refPoint<-refCurve[nrow(refCurve),]
  } else{
    refPoint<-refPoints[1,]
  }
  growthRef<-NA
  indexOnRef<-which(refCurve$Day==refPoint$Day)
  if(indexOnRef>7){
    dataYr<-refCurve[(indexOnRef-6):(indexOnRef),"CasesPerCapita"]
    dataXr<-refCurve[(indexOnRef-6):(indexOnRef),"Day"]
    fitRef <- lm(log(dataYr) ~ dataXr)
    growthRef<-coef(fitRef)[2]
    if(DEBUG){
      print(dataXr)
      print(dataYr)
      print(fitRef)
    }
  }
  if(nrow(testCurve)>7){
    if(DEBUG)print(which(refCurve$Day==refPoint$Day))
    
    testPoint<-testCurve[nrow(testCurve),]
    refPointY<-refPoint[1,"CasesPerCapita"]
    refPointX<-refPoint[1,"Day"]
    
    dataY<-testCurve[(nrow(testCurve)-6):(nrow(testCurve)),"CasesPerCapita"]
    dataX<-testCurve[(nrow(testCurve)-6):(nrow(testCurve)),"Day"]
    
    fit <- lm(log(dataY) ~ dataX)
    lag<-((log(refPointY)-coef(fit)[1])/coef(fit)[2]-refPointX)
    
    dt1=log(2)/log(1+coef(fit)[2])
    dt2=log(2)/log(1+growthRef)
    
    if(reverseRef){
      lag<-(-lag)
      dt1=log(2)/log(1+growthRef)
      dt2=log(2)/log(1+coef(fit)[2])
      
    }
    growthTest<-coef(fit)[2]
  } else{
    growthTest<-NA
    dt1<-NA
    dt2<-NA
    lag<-NA
  }
  
  return(data.frame(
    lag=lag,
    refday=refPoint$Day,
    testday=testCurve[nrow(testCurve),"Day"],
    growthRef=growthRef,
    growthTest=growthTest,
    doublingRegion1=dt1,
    doublingRegion2=dt2,
    CasesPerCapita=testCurve[nrow(testCurve),"CasesPerCapita"],
    ref=refRegion,
    test=testRegion,
    region=region1
  ))
}

compensateTime<-function(inputDF,referenceCurve){
  allUseRegions<-unique(as.character(inputDF$region))
  allLags<-data.frame(region=allUseRegions,lags=unlist(lapply(allUseRegions,function(x){getLagTimeCurves(inputDF,referenceCurve,x)[1,"lag"]})))
  for(ii in 1:nrow(allLags)){
    inputDF[as.character(inputDF$region)==allLags[ii,"region"],"Delay"]<-allLags[[ii,"lags"]]
  }
  inputDF$CompTime<-inputDF$Day-inputDF$Delay
  return(inputDF)
}

################# functions for combining the curve ################# 

getIncCases<-function(cumCases){cumCases-c(0,cumCases[-length(cumCases)])}
getMoveAv<-function(l,w){
  starts<-1:(length(l)-w)
  #return(c(rep(0,round(w/2)),unlist(lapply(starts,function(x)mean(l[x:(x+w)])))))
  return(unlist(lapply(starts,function(x)mean(l[x:(x+w)]))))
}

getCombinedCurve<-function(db,inclRegions,wind){
  subdata<-subset(db,region %in% inclRegions)
  start<-(min(floor(subdata$CompTime))+wind)
  end<-(max(floor(subdata$CompTime))-wind)
  getallCasesDay<-function(d){sum(subset(db,CompTime>(d)&CompTime<=(d+1))$Cases)}
  getallPopDay<-function(d){sum(subset(db,CompTime>(d)&CompTime<=(d+1))$Population)}
  popCounts<-unlist(lapply(start:end,getallPopDay))
  casesCounts<-unlist(lapply(start:end,getallCasesDay))
  
  return(
    data.frame(
      Day=start:end,
      combiPopulation=popCounts,
      combiCases=casesCounts,
      combiCasesPC=casesCounts/popCounts,
      combiInciPC=getIncCases(casesCounts/popCounts),
      smoothInciPC=c(rep(0,round(wind/2)),getMoveAv(getIncCases(casesCounts/popCounts),wind))
    )
  )
}

mirrorCurve<-function(cu,mirrorPoint,wind=2){
  curveLength<-nrow(cu)
  newLength<-(mirrorPoint*2)-1
  overlap<-curveLength-mirrorPoint
  
  mirrorInci<-(c(cu[1:mirrorPoint,"combiInciPC"],rev(cu[1:(mirrorPoint-1),"combiInciPC"])))
  
  day<-c(cu[,"Day"],((cu[nrow(cu),"Day"]+1:(mirrorPoint-1-overlap))))
  pop<-c(cu[,"combiPopulation"],rep(0,(mirrorPoint-1-overlap)))
  cas<-c(cu[,"combiCases"],rep(0,(mirrorPoint-1-overlap)))
  
  #combiInciPC #flat-out mirror this
  cpc<-cumsum(mirrorInci)
  smoothCPC=c(rep(0,round(wind/2)),getMoveAv(mirrorInci,wind))
  #smoothCPC=(getMoveAv(mirrorInci,wind))
  
  
  print(day)
  return(
    data.frame(
      Day=day,
      combiPopulation=pop,
      combiCases=cas,
      combiCasesPC=cpc,
      combiInciPC=mirrorInci,
      smoothInciPC=smoothCPC
    )
  )
}

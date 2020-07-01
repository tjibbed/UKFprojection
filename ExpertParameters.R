###########################Early estimate running######################################################

#### Parameters PB ####
losHosp_PB<-11
pSurv_PB<-0.60
losICU_surv_PB<-14
pRapid_PB<-0.30
losICU_rapid_PB<-5
pUltim_PB<-0.1
losICU_ultim_PB<-14
pM_PB<-0.25
pN_PB<-0.25
pI_PB<-0.02
pL_PB<-0.73 #should be written as derived, but is otherwise not needed.
los_ICU_PB<-(pSurv_PB*losICU_surv_PB)+(pRapid_PB*losICU_rapid_PB)+(pUltim_PB*losICU_ultim_PB)
deltaH_PB<-((1-pM_PB)/losHosp_PB)
alphaI_PB<-(pM_PB/losHosp_PB)
deltaI_PB<-1/los_ICU_PB

#### Parameters HB ####
losHosp_HB<-5
pSurv_HB<-0.80
losICU_surv_HB<-11
pRapid_HB<-0.05
losICU_rapid_HB<-5
pUltim_HB<-0.15
losICU_ultim_HB<-21
pM_HB<-0.22
pN_HB<-0.20
pI_HB<-0.03
pL_HB<-0.77 #should be written as derived, but is otherwise not needed.
los_ICU_HB<-(pSurv_HB*losICU_surv_HB)+(pRapid_HB*losICU_rapid_HB)+(pUltim_HB*losICU_ultim_HB)
deltaH_HB<-((1-pM_HB)/losHosp_HB)
alphaI_HB<-(pM_HB/losHosp_HB)
deltaI_HB<-1/los_ICU_HB


#### Parameters WK ####
losHosp_WK<-5
pSurv_WK<-0.80
losICU_surv_WK<-14
pRapid_WK<-0.05
losICU_rapid_WK<-5
pUltim_WK<-0.15
losICU_ultim_WK<-21
pM_WK<-0.25
pN_WK<-0.20
pI_WK<-0.01
pL_WK<-0.79 #should be written as derived, but is otherwise not needed.
los_ICU_WK<-(pSurv_WK*losICU_surv_WK)+(pRapid_WK*losICU_rapid_WK)+(pUltim_WK*losICU_ultim_WK)
deltaH_WK<-((1-pM_WK)/losHosp_WK)
alphaI_WK<-(pM_WK/losHosp_WK)
deltaI_WK<-1/los_ICU_WK


#### Parameters JK ####
losHosp_JK<-5
pSurv_JK<-0.80
losICU_surv_JK<-14
pRapid_JK<-0.05
losICU_rapid_JK<-3
pUltim_JK<-0.15
losICU_ultim_JK<-21
pM_JK<-0.20
pN_JK<-0.20
pI_JK<-0.04
pL_JK<-0.76 #should be written as derived, but is otherwise not needed.
los_ICU_JK<-(pSurv_JK*losICU_surv_JK)+(pRapid_JK*losICU_rapid_JK)+(pUltim_JK*losICU_ultim_JK)
deltaH_JK<-((1-pM_JK)/losHosp_JK)
alphaI_JK<-(pM_JK/losHosp_JK)
deltaI_JK<-1/los_ICU_JK

deltaH_all<-(1-mean(c(pM_PB,pM_HB,pM_WK,pM_JK)))/mean(c(losHosp_PB,losHosp_HB,losHosp_WK,losHosp_JK))
alphaI_all<-(mean(c(pM_PB,pM_HB,pM_WK,pM_JK)))/mean(c(losHosp_PB,losHosp_HB,losHosp_WK,losHosp_JK))
deltaI_all<-(1/mean(c(los_ICU_PB,los_ICU_HB,los_ICU_WK,los_ICU_JK)))
pL_all<-mean(c(pL_PB,pL_HB,pL_WK,pL_JK))
pN_all<-mean(c(pN_PB,pN_HB,pN_WK,pN_JK))
pI_all<-mean(c(pI_PB,pI_HB,pI_WK,pI_JK))

deltaH_all/(deltaH_all+alphaI_all)
alphaI_all/(deltaH_all+alphaI_all)
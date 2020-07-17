# UKFprojection
Projection of ICU bed occupancy in the University Medical Center Freiburg (Uniklinik Freiburg:  UKF)

Provided is the R script we used to produce forcasts of (ICU) bed demand in the University Medical Center Freiburg. 

We will post the manuscript with full description of methods on a preprint server shortly.

The MainAnalysis.R file is of main importance and includes the entire script producing the results and figures included in the manuscript. 

The other scripts contain the functions and variables needed by the main analysis file.
- EarlyIncidencePrediction.r: The static incidence forecast, based on observed reported cases in Italy (on province and region level), compensated for the delay between the trajectoried, and projected to the UKF catchment population. This prediction is only possible at the early stages of the epidemic. 
- IncidenceModel.R: The incidence forecast based on the locally observed reported cases, calculating the local Re(t), and simulating the future cases using a SIR model. These prediction can be made continuously as the epidemic / pandemic progresses

- WithinHospModel.R: The care path model, translating the local incidence to bed occupancy on both the general ward and the ICU.
- ExpertParameters.R: A list of parameters for the care path model, based on expert interviews
- Survival.R: The survival analysis used to estimate the parameters for the care path model from the locally observed admitted patients.

- FigurePlotting.R: functions to produce the various plots in the manuscript based on the above analysis.


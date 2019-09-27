#'@title Two-column coordinates(X,Y) extracted from published Kaplan Meier curves by ScanIt software
#'@description The dataset is extracted from a published Kaplan Meier image by ScanIt software. Locoreginal control events
#' were studied in 424 head and neck cancer patients: 213 in Radiotherapy treatment group and 211 in the Radiotherapy plus cetuximab group.
#' There are 145 pairs of coordinates extracted from the radiation treatment arm, and 136 pairs of coordinates are extracted from the radiation
#' plus arm. For both datasets, the first columns are the times, and the second columns are the survival rates at those time points(in percentage).
#' For each group, there are six patients numbers at risk reported at the months of 0, 10, 20, 30, 40, and 50. Three vectors (i.e., trisk, nrisk.radio and
#' nrisk.radioplus) record those numbers. \cr

#'@name Radiationdata
#'@usage Radiationdata
#'@format List of two dataframes (radio and radioplus) and three vectors
#'@keywords datasets, radiation KM on head and neck cancer
#'@examples
#' radio <- Radiationdata$radio
#' radioplus <- Radiationdata$radioplus
#' trisk <- Radiationdata$trisk
#' nrisk_radio <- Radiationdata$nrisk.radio
#' nrisk_radioplus <- Radiationdata$nrisk.radioplus
#' plot(radio,xlab="time",ylab="survival rates",type="l",
#'     lty=2,col="cyan4",xlim=c(1,70),main="Curves extracted by ScanIt software")
#' lines(radioplus,type="l",col="red4",lty=1)
#' legend("topright", c("Radiotherapy", "Radiotherapy plus cetuximab"),
#'       col = c("cyan4","red4"),lty=c(2,1),text.col = "green4")
#' text(40,80,"Reported Hazard Ratio with 95% CI: 0.69 (0.52,0.89)")
#'
#' pre_radio <- preprocess(dat=radio, trisk=trisk,nrisk=nrisk_radio,maxy=100)
#' est_radio <- getIPD(prep=pre_radio,armind=0,tot.events=NULL)
#' pre_radio_plus <- preprocess(dat=radioplus, trisk=trisk,nrisk=nrisk_radioplus,maxy=100)
#' est_radio_plus <- getIPD(prep=pre_radio_plus,armind=1,tot.events=NULL)
#' library(survival)
#' library(survminer)
#' mixdat <- rbind(est_radio$IPD,est_radio_plus$IPD)
#' coxfit <- coxph(Surv(time,status)~treat,data=mixdat)
#' fit <- survfit(Surv(time, status) ~ treat,data = mixdat)
#' ggsurvplot(fit, data = mixdat, risk.table = TRUE)
#'
#'
#'@references Bonner JA, Harari PM, Giralt J, Azarnia N, Shin DM, Cohen RB, Jones CU, Sur R,
#'Raben D, Jassem J, Ove R, Kies MS, Baselga J, Youssoufian H, Amellal N, Rowinsky EK, Ang KK:
#'Radiotherapy plus Cetuximab for Squamous-Cell Carcinoma of the Head and Neck. N Engl J Med.
#'2006, 354: 567-78. 10.1056/NEJMoa053422.


"Radiationdata"

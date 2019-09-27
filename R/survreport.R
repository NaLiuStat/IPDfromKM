#'@name survreport
#'@title Survival analysis on the reconstructed individual patients data (IPD).
#'@description This function will graph the Kaplan-Meier curves and the cumulative hazard curves from the
#'             reconstructed IPD (the output of \code{\link{getIPD}} function). It will also report
#'             the times when survival rates reach certain percentages, and the survival rates at some
#'             certain time points (for example, at every half-year time points).
#'
#'@usage survreport(ipd1,ipd2=NULL,arms=1,interval=6,s=c(0.75,0.5,0.25))
#'@param ipd1 The first individual patients' lifetable with three columns (time,status, and the treatment indicators), found in the output of \code{\link{getIPD}} function.  \cr
#'@param ipd2 The second indivdual patients' lifetable with three columns (time,status, and the treatment indicators), found in the output of \code{\link{getIPD}} function.  \cr
#'@param arms This indicate the functional analysis of data related to either one or two treatments: can be either 1 or 2.
#'@param interval The time intervals where we want to report the survival rates: the default is at every 6 months time interval.
#'@param s The survival rates want to report the corresponding times
#'@return NULL
#'@export
#'@importFrom survival Surv survfit
#'@importFrom survminer ggsurvplot
#'@importFrom gridExtra grid.arrange

#'@examples
#' ### Get data from the sample dataset=======================
#' radio <- Radiationdata$radio
#' radioplus <- Radiationdata$radioplus
#' trisk <- Radiationdata$trisk
#' nrisk_radio <- Radiationdata$nrisk.radio
#' nrisk_radioplus <- Radiationdata$nrisk.radioplus
#' ### Estimate the IPD for the Radiotherapy treatment group ====================
#' pre_radio <- preprocess(dat=radio, trisk=trisk,nrisk=nrisk_radio,maxy=100)
#' est_radio <- getIPD(prep=pre_radio,armind=0,tot.events=NULL)
#' ### Estimate the IPD for the Radiotherapy plus treatment group ====================
#' pre_radio_plus <- preprocess(dat=radioplus, trisk=trisk,nrisk=nrisk_radioplus,maxy=100)
#' est_radio_plus <- getIPD(prep=pre_radio_plus,armind=1,tot.events=NULL)
#' ### survival report for one arm ===================
#' survreport(ipd1=est_radio$IPD,arms=1,interval=6,s=c(0.8,0.5,0.3))
#' survreport(ipd1=est_radio_plus$IPD,arms=1,interval=10)
#' ### survival report for two arms ===================
#' survreport(ipd1=est_radio$IPD,ipd2=est_radio_plus$IPD,arms=2,interval=8)
#'
#'
#'@references Guyot P, Ades AE, Ouwens MJ, Welton NJ. Enhanced secondary analysis of survival data: reconstructing the data from published Kaplan-Meier survival curves. BMC Med Res Methodol.2012; 1:9.

require("survival")
require("ggplot2")
require("survminer")

survreport<- function(ipd1,ipd2=NULL,arms=1,interval=6, s=c(0.75,0.5,0.25)){

  ## check the input data set
  if (is.null(ipd1)) {stop('Need to input the first IPD dataset:')}
  if ((arms==2) & (is.null(ipd2))) {stop('Need to input the second IPD dataset if want to analyze two arms:')}

  if (sum(colnames(ipd1)==c("time","status","treat"))<3)
  {stop('The dataset should have exactly three columns named as "time","status" and "treat".')}


  if  ((!is.null(ipd2)) & (sum(colnames(ipd1)==c("time","status","treat"))<3))
  {stop('The dataset should have exactly three columns named as "time","status" and "treat".')}

  ### Peform KM analysis and Print the estimated quantities for the first arm
  fit1 <- survfit(Surv(time,status)~1,data=ipd1)
  cat("","\n")
  cat("  The estimated quantitles based on KM analysis of the reconstructed IPD")
  cat("  for the first treatment arm: ","\n")
  qt1 <- data.frame(q = s,
                    qt = quantile(fit1))
  names(qt1) <- c("survival_rate", "time","0.95LCI","O.95UCI")
  rownames(qt1) <- NULL
  print(qt1)

  ### Peform KM analysis and Print the estimated quantities for the second arm
  if (!is.null(ipd2)){
    fit2 <- survfit(Surv(time,status)~1,data=ipd2)
    cat("","\n")
    cat("  The estimated quantitles based on KM analysis of the reconstructed IPD")
    cat("  for the second treatment arm: ","\n")
    qt2 <- data.frame(q = s,
                      qt = quantile(fit2))
    names(qt2) <- c("survival_rate", "time","0.95LCI","O.95UCI")
    rownames(qt2) <- NULL
    print(qt2)
  }
  ### graph KM for one curve
  if (arms==1){
    qt11 <- qt1[!is.na(qt1$time),]
    qt12 <- qt1[is.na(qt1$time),]
    p1 <- ggsurvplot(fit1, data=ipd1, xlab = "Time ", title=" ",
                     censor = T,palette = c("#E7B800"))$plot
    if (nrow(qt11)>0)
      p1 <- p1+geom_segment(data = qt11, aes(x = time, y = survival_rate, xend = time, yend = 0), lty = 2) +
      geom_segment(data = qt11, aes(x = 0, y =survival_rate , xend = time, yend = survival_rate), lty = 2)

    if (nrow(qt12)>0)
      p1 <- p1+geom_segment(data = qt12, aes(x = 0, y =survival_rate , xend = ipd1$time[nrow(ipd1)],
                                             yend = survival_rate), lty = 2)
  }

  if (arms==2){
    qt11 <- qt1[!is.na(qt1$time),]
    qt12 <- qt1[is.na(qt1$time),]
    qt21 <- qt2[!is.na(qt2$time),]
    qt22 <- qt2[is.na(qt2$time),]
    mixdat <- rbind(ipd1,ipd2)
    fitm <- survfit(Surv(time,status)~treat,data=mixdat)
    p1 <- ggsurvplot(fitm, data=mixdat, xlab = "Time ",
                     censor = F,palette = c("#E7B800", "#2E9FDF"))$plot
    if (nrow(qt11)>0)
      p1 <- p1+geom_segment(data = qt11, aes(x = time, y = survival_rate, xend = time, yend = 0), lty = 2) +
      geom_segment(data = qt11, aes(x = 0, y =survival_rate , xend = time, yend = survival_rate), lty = 2)
    if (nrow(qt12)>0)
      p1 <- p1+geom_segment(data = qt12, aes(x = 0, y =survival_rate , xend = ipd1$time[nrow(ipd1)],
                                             yend = survival_rate), lty = 2)
    if (nrow(qt21)>0)
      p1 <- p1+geom_segment(data = qt21, aes(x = time, y = survival_rate, xend = time, yend = 0), lty = 2) +
      geom_segment(data = qt21, aes(x = 0, y =survival_rate , xend = time, yend = survival_rate), lty = 2)
    if (nrow(qt22)>0)
      p1 <- p1+geom_segment(data = qt22, aes(x = 0, y =survival_rate , xend = ipd2$time[nrow(ipd2)],
                                             yend = survival_rate), lty = 2)


  }

  ### graph cumulative hazard
  if (arms==1)
    p2 <- ggsurvplot(fit1,data=ipd1, fun = "cumhaz",palette = c("#E7B800"))$plot
  if (arms==2)
    p2 <- ggsurvplot(fitm, data=mixdat, fun = "cumhaz", palette = c("#E7B800", "#2E9FDF"))$plot

  grid.arrange(p1,p2,ncol=2)


  ### Print survival rates at some certain time points
  nn=floor(ipd1$time[nrow(ipd1)]/interval)
  sur=NULL
  low=NULL
  upp=NULL
  se=NULL
  new.names=NULL
  for (i in 1:nn)
  {
    sur=c(sur,summary(fit1,times=(i*interval))$surv)
    se=c(se,summary(fit1,times=(i*interval))$std.err)
    low=c(low,summary(fit1,times=(i*interval))$lower)
    upp=c(upp,summary(fit1,times=(i*interval))$upper)
    new.names=c(new.names,paste("time = ",i*interval,sep=""))
  }
  survtable=data.frame(survival=sur,se=se,lower=low,upper=upp)
  survtable=round(survtable,4)
  names(survtable)=c("survival_rate","standard_error","0.95LCI","0.95UCI")
  row.names(survtable)=new.names
  if (arms !=2) {
    cat("\n")
    cat(" Survival rates at some certain time points estimated from")
    cat(" reconstructed IPD: ","\n")
    print(survtable)
  }

  if (arms==2){
    nn=floor(ipd2$time[nrow(ipd2)]/interval)
    sur=NULL
    low=NULL
    upp=NULL
    se=NULL
    new.names=NULL
    for (i in 1:nn)
    {
      sur=c(sur,summary(fit2,times=(i*interval))$surv)
      se=c(se,summary(fit2,times=(i*interval))$std.err)
      low=c(low,summary(fit2,times=(i*interval))$lower)
      upp=c(upp,summary(fit2,times=(i*interval))$upper)
      new.names=c(new.names,paste("time = ",i*interval,sep=""))
    }
    survtable2=data.frame(survival=sur,se=se,lower=low,upper=upp)
    survtable2=round(survtable2,4)
    names(survtable2)=c("survival_rate","standard_error","0.95LCI","0.95UCI")
    row.names(survtable2)=new.names
    cat("\n")
    cat(" Survival rates at some certain time points estimated")
    cat("from reconstructed IPD: ","\n")
    ntime <- min(nrow(survtable),nrow(survtable2))
    survtab <- cbind(survtable[1:ntime,],survtable2[1:ntime,])
    survtab[ntime+1,1] <- "Treatment 1"
    survtab[ntime+1,5] <- "Treatment 2"
    row.names(survtab)[ntime+1] <- " "
    survtab[is.na(survtab)] <- ""
    print(survtab[c(ntime+1,1:ntime),])
  }


}




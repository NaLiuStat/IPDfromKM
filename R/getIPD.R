#'@name getIPD
#'@title Reconstruct individual patient data (IPD) from scanned Kaplan-Meier(KM) curves
#'@description This function uses the \code{\link{preprocess}} class as input to map the digitalized curve to the individual patient data.
#'       In addition, the total number of events is an optional input; and the treatment arm can be arbitrarily assigned to label the
#'       patients' group (for example, 0 for control group, and 1 for treatment group). \cr\cr
#'
#'       The output will be the reconstructed IPD in the form of three columns table (time, patient status, and treatment group indicator).
#'
#'@details  The extracted coordinates from scanned KM curve include time points, and the corresponding survival rates. To evaluate the accuracy of our reconstructed process,
#'       the survival rates at each read-in time points will be re-calculated from the estimated IPD, then compared with the original read-in survival rates from the
#'       coordinates. The root mean square error(RMSE) and p values (calculated via a bootstrap two-sampe t test, and Mann-Whitney test) will be reported as measurments
#'       of the discrepancy between estimated and scanned KM curves. \cr\cr
#'
#'       The output includes the reconstructed IPD in the form of lifetable. Using the\code{\link{plotIPD}} function, users can further graph
#'       the survival curve from the reconstructed IPD, and compare it with the read-in dataset from the scanned curve. Simple survival analysis can also be performed
#'       on the reconstructed IPD using the \code{\link{survreport}} function.
#'
#'@usage getIPD(prep,armind=1,tot.events=NULL)
#'@param prep The class object from the \code{\link{preprocess}} function.
#'@param armind The arbitrary lable used as the group indicator for the reconstructed IPD. Typically 0 for the control group and 1 for the treatment group.
#'@param tot.events The total events number. Only available for some published curves, and the default value is NULL.
#'@return An object of class getIPD containing eight components as follow. \cr
#'        IPD:  the estimated individual patient data with three columns (time,status,and treatment group indicator).  \cr
#'        points:   the dataframe with estimations at each read-in time point. \cr
#'        riskmat:   the dataframe of the lower and upper index, and the data estimates for censored patients and events within each time interval.\cr
#'        endpts: the number of patients remaining at the end of trial.\cr
#'        mwtest: the Mann-Whitney test result for comparing distribution of original and estimated curves.\cr
#'        rootp:  the p value of the bootstrap t test for the mean discrepancy of original and estimated curves.\cr
#'        var_surv: the variance of survival rates introduced by the reconstruction procedure.
#'        dt: the dataframe summarized the results.
#'@export
#'@importFrom survival Surv survfit
#'@importFrom boot boot
#'@examples
#'
#' # Radiationdata$radio is a dataset exported from ScanIt software ================
#' radio <- Radiationdata$radio
#'
#' # Load time points when the patients number at risk reported (in month) =========
#' trisk <- Radiationdata$trisk
#'
#' # Load patients number at risk reported at the time points in trisk =============
#' nrisk.radio <- Radiationdata$nrisk.radio
#'
#' # Use the trisk and nrisk as input for preprocess and reconstruction ============
#' pre_radio_1 <- preprocess(dat=Radiationdata$radio, trisk=trisk,
#'              nrisk=nrisk.radio,maxy=100)
#' est_radio_1 <- getIPD(prep=pre_radio_1,armind=0,tot.events=NULL)
#'
#' # Output include reconstructed individual patients data =========================
#' head(est_radio_1$IPD)
#'
#' # When trisk and nrisk were not available, then we must input ====================
#' # total number of initial patients ===============================================
#' pre_radio_2 <- preprocess(dat=Radiationdata$radio, totalpts=213,maxy=100)
#' est_radio_2 <- getIPD(prep=pre_radio_2,armind=0,tot.events=NULL)
#'
#' # Output include reconstructed individual patients data ==========================
#' head(est_radio_2$IPD)

#'@references Guyot P, Ades AE, Ouwens MJ, Welton NJ. Enhanced secondary analysis of survival data: reconstructing the data from published Kaplan-Meier survival curves. BMC Med Res Methodol.2012; 1:9.

require("survival")
require("boot")

getIPD <- function(prep,armind=1,tot.events=NULL){

  ## extract input data from list, and initial vectors----------------------------------------------------------
  dat <- prep[[1]]
  riskmat <- prep[[2]]
  endpts <- prep[[3]]
  ori_dat <- prep[[4]]
  TT <- dat$time
  SS <- dat$sur
  total <- nrow(dat)

  t.risk<-riskmat[,1]
  lower<-riskmat[,2]
  upper<-riskmat[,3]
  n.risk<-riskmat[,4]

  ninterval <- nrow(riskmat)
  ncensor=rep(0,ninterval-1)
  lasti<-rep(1,ninterval)

  cen <- rep(0,total)
  nhat <- rep(n.risk[1]+1,total+1)
  d <- rep(0,total)
  KM.hat<-rep(1,total)

  ## calculate the values for interval 1 to (ninterval-1)------------------------------------
  if (ninterval>1)
  {
    for (i in 1:(ninterval-1)){
      #First approximation of number of censored on interval i
      ncensor[i]<- round(n.risk[i]*SS[lower[i+1]]/SS[lower[i]]- n.risk[i+1])

      #Adjust ncensor until converge to estimation n.hat = n.risk at start of interval (i+1)
      while( ((nhat[lower[i+1]]>n.risk[i+1])&&(ncensor[i]<(n.risk[i]-n.risk[i+1]+1)))
             ||((nhat[lower[i+1]]<n.risk[i+1])&&(ncensor[i]>0)) ){
        if (ncensor[i]<=0){
          cen[lower[i]:upper[i]]<-0
          ncensor[i]<-0
        }

        ## find censor time within interval i:evenly distributed
        if (ncensor[i]>0){
          cen.t<-rep(0,ncensor[i])
          for (j in 1:ncensor[i]){
            cen.t[j]<- TT[lower[i]] +
              j*(TT[lower[(i+1)]]-TT[lower[i]])/(ncensor[i]+1)
          }
          for (k in lower[i]:upper[i])
          {cen[k] <- sum((TT[k]<=cen.t) & (cen.t<TT[k+1]))}
        }


        #Find number of events and patients at risk for each points
        nhat[lower[i]]<-n.risk[i]
        last<-lasti[i]

        # estimate d[k] for interval i
        for (k in lower[i]:upper[i]){
          if (i==1 & k==lower[i]){
            d[k]<-0
            KM.hat[k]<-1
          }
          else {
            d[k]<-round(nhat[k]*(1-(SS[k]/KM.hat[last])))
            KM.hat[k]<-KM.hat[last]*(1-(d[k]/nhat[k]))
          }
          nhat[k+1]<-nhat[k]-d[k]-cen[k]
          if (d[k] != 0) last<-k
        }  # end for k loop

        # update ncensor for interval i
        ncensor[i]<- ncensor[i]+(nhat[lower[i+1]]-n.risk[i+1])
      } # end while loop

      ## repare for next interval
      n.risk[i+1]<-nhat[lower[i+1]]
      lasti[(i+1)]<-last
    }  # end for loop
  }  # end if ninterval>1



  ## calculate the value for the last interval-----------------------------------------------------------

  ## for the case have more than one intervals, assume the last interval has the same censor rate as before
  if (ninterval>1)
  { ## the number of events happend in the last interval
    if (is.null(tot.events))
       {leftd <- 0}
    else
       {temp <- sum(d[1:(lower[ninterval]-1)])
       leftd <- tot.events-temp
       leftd <- ifelse(leftd<=0,0,leftd)
    }
    ## patients number at the end of trial
    mm <- ifelse(is.null(endpts),0,endpts)
    ## the number of censoring happened in the last interval
    ncensor[ninterval] <- min(mean(ncensor[1:(ninterval-1)])*(TT[total]-t.risk[ninterval])
                              /(t.risk[ninterval]-t.risk[ninterval-1]),(n.risk[ninterval]-mm-leftd))
  }

  ## for the case only have one interval (no nrisk and trisk information), assume ncensor=0
  if (ninterval==1)
  {if (is.null(tot.events))
  {ncensor[ninterval] <- 0}  ## since no censor information available, assume ncensor=0
    else {ncensor[ninterval] <- n.risk[ninterval]-tot.events}
  }


  ## estimate the censoring for each point
  if (ncensor[ninterval]<=0){
    cen[lower[ninterval]:upper[ninterval]]<-0
    ncensor[ninterval]<-0
  }
  if ((ncensor[ninterval]>0) && (upper[ninterval]>lower[ninterval])){
    cen.t<-rep(0,ncensor[ninterval])
    for (j in 1:ncensor[ninterval]){
      cen.t[j]<- TT[lower[ninterval]] +
        j*(TT[upper[(ninterval)]]-TT[lower[ninterval]])/(ncensor[ninterval]+1)
    }
    for (k in lower[ninterval]:upper[ninterval])
    {cen[k] <- sum((TT[k]<=cen.t) & (cen.t<TT[k+1]))}
    cen[is.na(cen)] <- 0
  }
  if ((ncensor[ninterval]>0) && (upper[ninterval]==lower[ninterval]))
  {cen[upper[ninterval]] <- ncensor[ninterval]}

  ## estimate number of events and patients at risk for each point
  nhat[lower[ninterval]]<-n.risk[ninterval]
  last<-lasti[ninterval]
  for (k in lower[ninterval]:upper[ninterval]){
    if (k==1){
      d[k]<-0
      KM.hat[k]<-1
    }
    else {
      if (KM.hat[last]!=0) {d[k]<-round(nhat[k]*(1-(SS[k]/KM.hat[last])))}
      else {d[k] <- 0}
      KM.hat[k]<-KM.hat[last]*(1-(d[k]/nhat[k]))
    }
    nhat[k+1]<-nhat[k]-d[k]-cen[k]
    if (nhat[k+1]<0) {
      nhat[k+1] <- 0
      cen[k] <- nhat[k]-d[k]
    }

    if (d[k] != 0) last<-k
  }  # end for k loop



  ### summarize the estimation of number of censor numbers, event numbers, and patients numbers -------------------
  ### at risk for each interval -----------------------------------------------------------------------------------
  event.hat <- rep(0,ninterval)
  n.risk.hat <- rep(0,ninterval)
  censor.hat <- rep(0,ninterval)
  for (i in 1: ninterval){
    censor.hat[i] <- sum(cen[lower[i]:upper[i]])
    n.risk.hat[i] <- nhat[lower[i]]
    event.hat[i] <- sum(d[lower[i]:upper[i]])
  }
  riskmat <- cbind (riskmat,n.risk.hat,censor.hat,event.hat)
  Points=cbind(time=TT,surv=SS,risk=nhat[1:total],censor=cen,event=d)


  ## reconstruct individual patients records--------------------------------------------------------------------
  ipd=data.frame(time=numeric(0),status=integer(0),treat=integer(0))
  for (i in 1:total){
    if (d[i]>0) {
      time <- rep(dat$time[i],d[i])
      status <- rep(1,d[i])
      treat <- rep(armind,d[i])
      ipd <- rbind(ipd,cbind(time,status,treat))
    }
    if (cen[i]>0) {
      if (i < total) {time <- rep((dat$time[i]+dat$time[i+1])/2,cen[i])}
      if (i == total) {time <- rep(dat$time[i],cen[i])}
      status <- rep(0,cen[i])
      treat <- rep(armind,cen[i])
      ipd <- rbind(ipd,cbind(time,status,treat))
    }
  }  # end for loop

  if (nhat[total+1]>0)
  {
    time <- rep(dat$time[total],nhat[total+1])
    status <- rep(0,nhat[total+1])
    treat <- rep(armind,nhat[total+1])
    ipd <- rbind(ipd,cbind(time,status,treat))
  }



  ## comparing the survival rates calculated from the estimated IPD and the read-in data -----------------------------------------------------
  ## using

  fit <- survfit(Surv(ipd$time,ipd$status)~1,data=ipd)
  ## find survival estimation from ipd for each time points
  ## ss2 vector has the survival rates from reconstruct, ss1 has the survival rates from the original read in,
  ## and tt has all the non-repeated time points
  ## diff_s vector records the differences between estimated and original survival rates
  ss1 <- c()
  tt <- c()
  for (i in 1:(total-1))
  {if (TT[i]!=TT[i+1])
  {ss1 <- c(ss1,SS[i])
  tt <- c(tt,TT[i])}
  }
  ll <- length(tt)
  ss2 <- rep(0,ll)
  diff_s <- rep(0,ll)
  index <- 1
  for (i in 1:ll){
    if (tt[i]<fit$time[1]) {ss2[i] <- 1}
    else{
      while ((fit$time[index]<=tt[i]) && (index<length(fit$time))) {index <- (index+1)}
      if (index>1) {ss2[i] <- fit$surv[index-1]}
    }
    diff_s[i] <- ss2[i]-ss1[i]
  }

  ## find rmse, standard deviation of unexplained variance
  rmse <- (sum(diff_s^2)/(ll-1))^0.5

  ## variance introduced from the reconstruction procesure= mse + variance by read in error
  del_s <- 1/500
  var_surv <- rmse^2+del_s^2

  ## Mann-Whitney test
  mwtest <- wilcox.test(ss1,ss2)

  ## bootstrap
  nBoot <- 1000
  bootdat <- cbind(ss1,ss2)
  tstat <- function(bootdat) {
    tobj <- t.test(bootdat[,1],bootdat[,2])
    t <- tobj$statistic
    return(t)
  }
  bootstat <- function(bootdat, indices) {
    d <- bootdat[indices,]
    t <- tstat(d)
    return(t)
  }
  b <- boot(data=bootdat,statistic=bootstat,R=nBoot)
  bootp <- length(which(b$t > b$t0))/nBoot

  ## print out results----------------------------------------------------------
  cat("","\n")
  cat("              Total number of patients is ", nrow(ipd), "\n")
  cat("  The summary of the estimation is as follow", "\n")
  print(riskmat)
  cat("  The root mean square error(RMSE) of the estimations is", rmse, "\n")
  cat("  The max absolute error of the estimation is ", max(abs(diff_s)), "\n")
  cat("  The variance of survival rates introduced by the reconstruction procedure is  ", var_surv, "\n")

  cat("","\n")
  cat("        Using bootstrap t test to compare the discrepancy of original and estimated curves","\n")
  cat("date:  original survival rates, estimated survival rates  ", "\n")
  cat("nBoot= ", nBoot, ",  p-value= ",bootp, "\n")
  cat("Alternative hypothesis: the mean discrepancy is not equal to 0","\n")

  cat("","\n")
  cat("        Using Mann-Whitney test to compare the distributions of original and estimated curves","\n")
  cat("date:  original survival rates, estimated survival rates  ", "\n")
  cat("W= ", mwtest$statistic, ",  p-value= ",mwtest$p.value, "\n")
  cat("Alternative hypothesis: true location shift is not equal to 0","\n")

  ## make a DT table form output summary
  type=c("Total Number of Patients", "Root Mean Square Error", "Variance of Survival Rates from the Reconstruction",
         "Mann-Whiteney Test P-value","Bootstrap T Test P-value")
  value=c(nrow(ipd),rmse,var_surv,mwtest$p.value,bootp)
  d=data.frame(ResultSummary=type,value=round(value,3))

  ## return
  re <- list(IPD=ipd, Points=Points,riskmat=riskmat,endpts=endpts,mwtest=mwtest,bootp=bootp,var_surv=var_surv,dt=d)
  return(re)
}   # function getIPD


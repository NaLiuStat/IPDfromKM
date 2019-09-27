#'@name preprocess
#'@title Preprocess the X,Y coordinates dataset extracted from Kaplan-Meier(KM) curves
#'@description Before using the \code{\link{getIPD}} function to estimate the individual patients records (IPD), the coordinates of the survival
#'             curves exported from the digitized software need to be pretreated. \cr\cr
#'             The preprocess function will do three things.
#'             First,sort the data by time progression and discard any unnecessary points to avoid duplications.
#'             Second, adjust the monotonicity to make sure that the survival rates are non-increasing with time.
#'             Third, split the curve into intervals. For each starting point of an interval, we find the number of patients at risk,
#'             the time at which the number at risk was provided, and the first and the last index of the extracted points coordinates within each time interval.
#'@details We recommend using \href{https://www.digitizeit.de}{DigitizeIt} or \href{https://www.amsterchem.com/scanit.html}{ScanIt} software to extract
#'         the X,Y coordinates from the scanned Kaplan-Meier images. The coordinates datasets exported from the digitize softwares are either .csv
#'         or .txt files with two columns: the first column is the time, and the second column is the survival rates. \cr\cr
#'         For most of Kaplan-Meier curves, we can also find the numbers of patients at risk at normally 5-10 or more time points under the x-axis. Those time points
#'         and the numbers of patients should be manually input as the form of vectors. However, when those information is not availabe, we can leave the "trisk" and
#'         "nrisk" parameter as "NULL". But the total initial patients number "totalpts" can not be "NULL" in this case.  \cr\cr
#'         Sample data can be found in \code{\link{Radiationdata}}.
#'@usage preprocess(dat,trisk=NULL,nrisk=NULL,totalpts=NULL,maxy=100)
#'@param dat The two columns dataset read by some digitized software programs. The first column should be the times, and the second column should be the survival rates.
#'@param trisk A list of times points when the numbers of patients at risk were reported. This often can be found under the x-axis of the KM curve. The default value is NULL.
#'@param nrisk A list of the numbers of patients at risk reported at the time points in the trisk vector. This often can be found under the x-axis of the KM curve. The default value is NULL.
#'@param totalpts The initial number of patients, with a default value of NULL. However, when both trisk and nrisk are NULL, we need this number for the estimation.
#'@param maxy The scale of survival rates. This is 100 when the percentages are provided; 1 when the decimal numbers are provided.
#'@return An object of class preprocess containing four components as below. \cr
#'        originaldat: the read-in dataset.\cr
#'        newdat: the two column data after preprocess in the form of (time,survival). \cr
#'        intervalIndex: the table reporting the lower and upper index of the extracted coordinates within each time interval.\cr
#'        endpts: the number of patients remaining at the end of the trial.\cr
#'


#'@export
#'@examples
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
#'
#'@references Guyot P, Ades AE, Ouwens MJ, Welton NJ. Enhanced secondary analysis of survival data:
#'reconstructing the data from published Kaplan-Meier survival curves. BMC Med Res Methodol.2012; 1:9.


preprocess <- function(dat,trisk=NULL,nrisk=NULL,totalpts=NULL,maxy=100){

  ## Formalize the dataset to appropriate values and orders-----------------------------------------------------------
  indat <- dat
  ## check the input data set
  if (ncol(dat)!=2) {stop('The dataset should have exactly two columns.')}

  ## name the columns with "time" and "sur"
  names(dat) <- c("time","sur")

  ## order the data with time
  dat <- dat[order(dat$time),]

  ## rescale to survival rate to within [0,1] in case using percentage scale
  if (maxy==100) dat$sur <- dat$sur/100

  ## formalize dataset to make sure that the start point was (0,1), all survival
  ## rates within [0,1], and time always positive.
  dat <- dat[which(dat$sur<=1),]
  dat <- dat[which(dat$sur>=0),]
  dat <- dat[which(dat$time>=0),]
  if (dat[1,1]!=0) dat <- rbind(c(0,1),dat)
  if (dat[1,1]==0 & dat[1,2]<1) dat <- rbind(c(0,1),dat)

  ## make sure there is no repeat reads
  dat<-unique(dat)

  ## remove outlier points
  nl <- nrow(dat)
  d <- rep(0,nl-1)
  for (i in 1:(nl-1)) {d[i] <- abs(dat$sur[i+1]-dat$sur[i])}
  ndx <- order(d) [(nl-1):min(ceiling(0.99*nl),(nl-1))]
  rdx <- NULL
  for (j in 1:length(ndx)){
    if (is.element(ndx[j]+1,ndx)) rdx <- c(rdx,ndx[j]+1)
  }
  if (!is.null(rdx)) {dat <- dat[-rdx,]}

  ## make sure one time point has at most two survival reads
  uni.t <-  unique(dat[,1])
  newdat  <- data.frame(time=numeric(),sur=numeric())
  for(i in 1:length(uni.t)){
    temp.time <-  uni.t[i]
    subdat <-  dat[which(dat[,1] == temp.time),]
    if(nrow(subdat) > 1){
      minrow <- c(temp.time,min(subdat$sur))
      maxrow <- c(temp.time,max(subdat$sur))
      newdat[nrow(newdat)+1,] <- maxrow
      newdat[nrow(newdat)+1,] <- minrow
    } else newdat[nrow(newdat)+1,] <-  subdat
  }
  rownames(newdat) <- NULL

  ## As time increasing, make sure survival rates were non-increasing
  m <- nrow(newdat)
  for(i in 2:m)
    if(newdat$sur[i] > newdat$sur[i-1])
      newdat$sur[i] <-  newdat$sur[i-1]
  newdat <- unique(newdat)

  ## find out the interval index table (riskmat)-------------------------------------------------------------

  ## for the case have information of trisk and nrisk, find riskmat table
  if (!(is.null(nrisk))){

    ## trim nrisk to have at most one zero at the end, then match the length of nrisk and trisk
    n.int <- length(nrisk)
    while ((nrisk[n.int]==0) && (nrisk[n.int-1]==0)) n.int <- n.int-1
    nrisk <- nrisk[1:n.int]
    trisk <- trisk[1:n.int]

    ## extract upper and lower index for each interval
    lower <- NULL
    upper <- NULL
    temp.t <-  newdat[,1]
    event.t <-  NULL

    for(i in 1:n.int){
      if (i<n.int) {my.id <- which(temp.t >= trisk[i] & temp.t < trisk[(i+1)])}
      ## right open last interval
      else {my.id <-  which(temp.t >= trisk[i])}

      if(length(my.id) > 0){
        event.t <-  c(event.t, trisk[i])
        lower <-  c(lower, my.id[1])
        upper <-  c(upper, my.id[length(my.id)])
      }
    }
  } ## finish the !isnull(nrisk) case

  ## when trisk and nrisk not available, assign riskmat
  if (is.null(nrisk)){
    event.t <- 0
    lower <- 1
    upper <- dim(newdat)[1]
    if (is.null(totalpts)) { stop ("At least need to know total initial number of patients.")}
    else {nrisk <- totalpts}
  }

  ## record the endpts: when the curve end at the horizontal segment and there are patients alive at the end of the trial,
  ##                    record the number of patients left in the endpts (useful for IPD estimaiton); otherwise, null
  if (!is.null(trisk))
  {if (trisk[length(trisk)] != event.t[length(event.t)]) {endpts=nrisk[length(nrisk)]}
    else {endpts=NULL}
  }
  if (is.null(trisk)) {endpts=NULL}

  ## output ---------------------------------------------------------------------------------

  riskmat = data.frame(t.risk=event.t, lower=lower, upper=upper,
                       n.risk=nrisk[1:length(event.t)])
  re=list(newdat=newdat,intervalIndex=riskmat,endpts=endpts,originaldat=indat)
  cat("Total points read from Kaplan-Meier curve=", nrow(indat), "\n")
  cat("Points left after preprocess =", nrow(newdat), "\n")
  cat("The indexes for each reported interval", "\n")
  print(riskmat)
  return (re)
}


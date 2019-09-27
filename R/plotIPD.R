#'@name plotIPD
#'@title Plot graphs comparing the survival estimates from the reconstructed individual patient data(IPD) with the data from the read-in coordinates
#'@description This function plots the survival curve from the reconstructed individual patient data (IPD),
#'              and compares it with the scanned curve extracted from published image. The output includes three graphs:
#'              the graph comparing the two Kaplan-Meier curves from the read-in coordinates, and from the reconstructed IPD;
#'              the graph comparing the reported patients number at risk and the estimated patient number at risk;
#'              and the graph plotting the discrepancy between the estimated survival rates and the read-in survival rates.
#'
#'@usage plotIPD(est,ori_dat,maxy=100)
#'@param est The list returned from the \code{\link{getIPD}} function.
#'@param ori_dat The two-column read-in coordinates data used as the input of the\code{\link{preprocess}} function. The first column should be the times, and the second column should be the survival rates.
#'@param maxy The scale of survival rates. 100 when the percentages are provided; and 1 when the decimal numbers are provided.
#'@return NULL


#'@export
#'@importFrom survival Surv survfit
#'@importFrom ggplot2 ggplot
#'@importFrom gridExtra grid.arrange


#'@examples
#' # Radiationdata$radio is a dataset exported from ScanIt software
#' radio <- Radiationdata$radio
#' # Load time points when the patients number at risk reported (in month)
#' trisk <- Radiationdata$trisk
#'
#' ##### Use the trisk and nrisk as input ==========
#' pre_radio_1 <- preprocess(dat=Radiationdata$radio, trisk=trisk,nrisk=nrisk.radio,maxy=100)
#' est_radio_1 <- getIPD(prep=pre_radio_1,armind=0,tot.events=NULL)
#' # Output include reconstructed individual patients data
#' head(est_radio_1$IPD)
#' # Plot the estimations
#' plotIPD (est=est_radio_1,ori_dat=Radiationdata$radio,maxy=100)
#'
#' ##### When trisk and nrisk were not available, then we must input total ========
#' ##### number of initial patients ===============================================
#' pre_radio_2 <- preprocess(dat=Radiationdata$radio, totalpts=213,maxy=100)
#' est_radio_2 <- getIPD(prep=pre_radio_2,armind=0,tot.events=NULL)
#' # Output include reconstructed individual patients data
#' head(est_radio_2$IPD)
#' # Plot the estimations
#' plotIPD (est=est_radio_2,ori_dat=Radiationdata$radio,maxy=100)
#'
#'@references Guyot P, Ades AE, Ouwens MJ, Welton NJ. Enhanced secondary analysis of survival data: reconstructing the data from published Kaplan-Meier survival curves. BMC Med Res Methodol.2012; 1:9.

require("survival")
require("ggplot2")
require("survminer")
require("gridExtra")

plotIPD <- function (est,ori_dat,maxy=100)

{
  ## unpack data
  ipd <- est$IPD
  poi <- est$Points
  riskmat <- est$riskmat
  endpts <- est$endpts
  total <- nrow(poi)
  TT <- est$Points[,1]
  SS <- est$Points[,2]
  t.risk<-riskmat[,1]
  lower<-riskmat[,2]
  upper<-riskmat[,3]
  n.risk<-riskmat[,4]
  ninterval <- nrow(riskmat)
  nhat <- poi[,3]

  ## find survival estimation from ipd for each time points-------------------------------------------------------
  ## ss2 vector has the survival rates from reconstruct, ss1 has the survival rates from the original read in,
  ## and tt has all the non-repeated time points
  ## diff_s vector records the differences between estimated and original survival rates
  fit <- survfit(Surv(ipd$time,ipd$status)~1,data=ipd)
  ## ss1
  ss1 <- c()
  tt <- c()
  for (i in 1:(total-1))
  {if (TT[i]!=TT[i+1])
  {ss1 <- c(ss1,SS[i])
  tt <- c(tt,TT[i])}
  }
  ll <- length(tt)
  ## ss2
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

  ## begin graph --------------------------------------------------------------------------------------

  ### graph 1
  gtt <- rep(tt,2)
  gss <- c(ss1,ss2)
  tp <- c(rep("read-in co-ordinates",length(tt)),rep("estimated IPD",length(tt)))
  gtab <- data.frame("time"=gtt,"survival"=gss,"From"=tp)
  p1 <- ggplot() +
    geom_line(data = gtab, aes(x = time, y = survival, color = From), size = 1)+
    labs(subtitle="Compare KM curves from the read-in \n Coordinates with the one from the estimated IPD",
         x="time",
         y="survival rates")


  ## graph2: comparing estimated risk and true risk
  if (length(t.risk)>1)
  {  gtab1 <- data.frame(time=TT,n=nhat[1:total],type=rep("Estimated",total))
  gtab2 <- data.frame(time=t.risk,n=n.risk,type=rep("Reported",length(n.risk)))
  zz <- rbind(gtab1,gtab2)
  p2 <- ggplot(zz, aes(x=time, y=n, color=type)) +
    geom_point(data=zz[zz$type=="Reported",]) +
    geom_line(data=zz[zz$type=="Estimated", ]) +
    scale_color_manual("Patients numbers at risk",
                       values = c("Reported" = "red4", "Estimated" = "blue4"))+
    labs(subtitle="Compare the estimated and reported \n patients numbers at risk",
         x="time",
         y="Patients number at risk")
  }


  ## graph3: calculate discprency suvival rates and graph
  gtab3 <- data.frame("time"=tt,"Difference"=diff_s)
  p3 <- ggplot(gtab3, aes(x=time, y=Difference, color="blue4")) +
    geom_point() + ylim(-0.15,0.15)+
    geom_hline(yintercept=0, color = "red4")+
    labs(subtitle="Discrepancy of estimated and read in \n survival rate",
         x="time",
         y=expression(paste(Delta," suvival rate")))+
    theme(legend.position='none')

  if (length(t.risk)>1) {
    grid.arrange(p1,p2,p3,nrow=2)
  }
  if (length(t.risk)==1) {grid.arrange(p1,p3,nrow=2)}


}

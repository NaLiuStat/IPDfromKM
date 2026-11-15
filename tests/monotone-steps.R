library(IPDfromKM)

# Five deaths at 1,5,5,10,10. The two largest adjacent drops are genuine.
# 200 coordinates are the minimum that selects two differences in this filter.
dat <- data.frame(time=c(seq(0,.99,length.out=99),
                         seq(1,4.9,length.out=99),5,10),
                  surv=c(rep(1,99),rep(.8,99),.4,0))
p <- preprocess(dat,totalpts=5,maxy=1)
stopifnot(identical(p$newdat$time, dat$time),
          identical(p$newdat$sur, dat$surv),
          identical(p$originaldat, dat))
set.seed(725)
ipd <- suppressWarnings(getIPD(p))$IPD
stopifnot(nrow(ipd) == 5, all(ipd$status == 1),
          identical(sort(ipd$time), c(1,5,5,10,10)))
fit <- survfit(Surv(time,status) ~ 1, ipd)
stopifnot(summary(fit)$table['median'] == 5,
          abs(summary(fit,times=5)$surv - .4) < 1e-12)

# Percentage input must retain the same event times and step heights.
p <- preprocess(transform(dat,surv=100*surv),totalpts=5,maxy=100)
stopifnot(identical(p$newdat$time, dat$time),
          isTRUE(all.equal(p$newdat$sur, dat$surv)))

# Still remove an isolated downward click followed by a return to the plateau.
spike <- dat
spike$surv[199:200] <- c(.2,.8)
p <- preprocess(spike,totalpts=5,maxy=1)
stopifnot(!5 %in% p$newdat$time, tail(p$newdat$sur,1) == .8)

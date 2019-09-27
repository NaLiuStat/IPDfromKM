# IPDfromKM R package

This package is an implementation to reconstruct individual patient data (IPD) from coordinates extracted from 
Kaplan-Meier (KM) survival curves. Digitize software, such as DigitizeIt(for MAC or windows) or ScanIt(for windows) 
can be used to get the coordinates from the published KM curves. 
To increasing accuracy, calculations also can include optional information, like reported at risk patient numbers 
(reported at the beginning of time intervals), total number of initial patients, and total number of events.

## Installing

The package can be intalled in R studio with version (>= 3.5.0)

```
install.packages("devtools")
devtools::install_github("NaLiuStat/IPDfromKM")
```



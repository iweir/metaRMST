# metaRMST
An R package for meta-analysis of RMSTD

#### R Installation Instructions
```
install.packages("devtools")
library(devtools)
devtools::install_github("iweir/metaRMST")
library(metaRMST)
```

#### Example 
```
# read in built-in dataset 
dat <- AorticStenosisTrials

# demonstration of meta-analysis to obtain combined effect by multivariate meta-analysis model (method="mvma")
result <- metaRMSTD(dat, time_horizons=c(12,24,36), MA_method="mvma")

# generate figure:
obj <- RMSTcurves(dat, time_horizons=c(12,24,36))
RMSTplot(obj)
```

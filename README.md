# metaRMST
An R package for meta-analysis of RMSTD (available on CRAN at: https://CRAN.R-project.org/package=metaRMST)

Reference: Weir, IR., Tian, L., and Trinquart, L. (2019). Multivariate meta-analysis model for the difference in restricted mean survival times. Biostatistics. 

#### R Installation Instructions
```
# from CRAN:
install.packages("metaRMST")
library(metaRMST)

# or you may install directly from github:
install.packages("devtools")
library(devtools)
devtools::install_github("iweir/metaRMST")
library(metaRMST)
```

#### Example 
```
# load library
library(metaRMST)

# read in built-in dataset 
data(AorticStenosisTrials)

# demonstration of meta-analysis to obtain combined effect by multivariate meta-analysis model (method="mvma")
mvma_res <- metaRMSTD(AorticStenosisTrials, time_horizons=c(12,24,36), MA_method="mvma")
mvma_res$REresult

# generate figure:
obj <- RMSTcurves(AorticStenosisTrials, time_horizons=c(12,24,36), tmax=40, nboot=50)
RMSTplot(obj, xlim=c(0,40), ylim=c(-0.25,2.75), yby=0.5, ylab="RMSTD (mos)", xlab="Time (mos)")
```

![](https://cranlogs.r-pkg.org/badges/grand-total/metaRMST?color=yellow)

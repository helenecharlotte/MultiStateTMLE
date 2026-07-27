### master.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Feb  5 2026 (10:09) 
## Version: 
## Last-Updated: May 13 2026 (20:04) 
##           By: Helene
##     Update #: 17
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

run.master <- function(install.simevent = FALSE, file.path = "") {
    
    library(data.table)
    library(survival)
    library(ggplot2)
    library(xtable)
    library(gridExtra)
    library(zoo)
    library(nleqslv)
    library(foreach)
    library(doParallel)
    library(parallel)
    library(devtools)
    if (install.simevent) remotes::install_github("miclukacova/simevent")
    library(simevent)

    source(paste0("./", file.path, "R/prepare.initial.R"))
    source(paste0("./", file.path, "R/compute.Q.clever.per.id.R"))
    source(paste0("./", file.path, "R/tmle.alpha.fun.R"))
    source(paste0("./", file.path, "R/calibration.fun.R"))    
    source(paste0("./", file.path, "R/estimate.alpha.fun.R"))
    source(paste0("./", file.path, "R/estimate.derivative.R"))
}

######################################################################
### master.R ends here

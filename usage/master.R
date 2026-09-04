### master.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Feb  5 2026 (10:09) 
## Version: 
## Last-Updated: Sep  4 2026 (13:58) 
##           By: Helene
##     Update #: 24
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

    source(paste0("./", file.path, "R/sim.from.data.R"))
    source(paste0("./", file.path, "R/sim.generic.R"))
    source(paste0("./", file.path, "R/prepare.initial.R"))
    source(paste0("./", file.path, "R/compute.Q.clever.per.id.R"))
    source(paste0("./", file.path, "R/tmle.alpha.fun.R"))
    source(paste0("./", file.path, "R/calibration.curve.fun.R"))    
    source(paste0("./", file.path, "R/calibration.fun.R"))
    source(paste0("./", file.path, "R/make.calibrated.contrasts.R"))        
    source(paste0("./", file.path, "R/estimate.alpha.fun.R"))
    source(paste0("./", file.path, "R/estimate.derivative.R"))
}

######################################################################
### master.R ends here

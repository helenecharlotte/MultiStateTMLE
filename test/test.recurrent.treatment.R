### test.recurrent.treatment.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Feb  8 2026 (21:00) 
## Version: 
## Last-Updated: Feb  9 2026 (14:29) 
##           By: Helene
##     Update #: 41
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:



if (system("echo $USER",intern=TRUE)%in%c("jhl781")){
    setwd("//projects/biostat01/people/jhl781/research/TMLE-from-2020june/targeted-max-likelihood/")
    .libPaths( c( .libPaths(), "//projects/biostat01/people/jhl781/tmp"))
    local <- FALSE
} else {
    setwd("~/research/TMLE-from-2020june/targeted-max-likelihood/")
    local <- TRUE
}

library(devtools)
load_all('./simevent/')
library(simevent)

if (system("echo $USER",intern=TRUE)%in%c("jhl781")){
    setwd("//projects/biostat01/people/jhl781/research/TMLE-from-2020june/targeted-max-likelihood/MultiStateTMLE/")
    .libPaths( c( .libPaths(), "//projects/biostat01/people/jhl781/tmp"))
} else {
    setwd("~/research/TMLE-from-2020june/targeted-max-likelihood/MultiStateTMLE/")
}


source("./usage/master.R")
source("./extra/simRecurrentTreatment.R")

run.master()


alpha <- 0.5
use.cores <- 3
N <- 500
m <- 1


set.seed(m+33339)

N <- 500


eta <- c(0.1, 0.025, 0.085, 0.1, 0.1, 0.025)
nu <- rep(1.1, 6)

eta[3] <- 2.5*eta[3]
nu[3] <- nu[3]

set.seed(2020115)
    
plotEventData(dt <- simRecurrentTreatment(N = N, 
                                          eta = eta, 
                                          nu = nu,
                                          beta_L_Z = 1.5, beta_L_D = 1.25,
                                          beta_W_Z = 1.5, beta_W_D = 1.25,
                                          #beta_Z_D = -0.5, beta_Z_L = -1.25, beta_Z_W = -1.25, beta_Z_D2 = -0.5,
                                          beta_Z_D = 0, beta_Z_L = 0, beta_Z_W = 0, beta_Z_D2 = 0,
                                          beta_L0_L = 1, beta_L0_W = 1,
                                          beta_W_L = 1.5, beta_L_W = 0.5,
                                          beta_L0_Z = 1,
                                          lower = 10^(-150),      
                                          upper = 1e3
                                          ))

setnames(dt, "ID", "id")
setnames(dt, "Time", "time")
setnames(dt, "Delta", "delta")


dt[time <= 3, table(delta)]

table(dt[time <= 3, sum(delta == 1), by = "id"][[2]])
table(dt[time <= 3, sum(delta == 2), by = "id"][[2]])

run.master()

initial.fit <- prepare.initial(dt,
                               tau = 3,
                               fit.types = list(
                                   Z = list(model = "Surv(tstart, tstop, delta == 2)~L0+L+W+Z",
                                            fit = "cox"),
                                   L = list(model = "Surv(tstart, tstop, delta == 3)~L0+Z+W",
                                            fit = "cox",
                                            at.risk = function(dt) (dt[["L"]] == 0)),
                                   W = list(model = "Surv(tstart, tstop, delta == 4)~L0+Z+L",
                                            fit = "cox",
                                            at.risk = function(dt) (dt[["W"]] == 0)),
                                   outcome1 = list(model = "Surv(tstart, tstop, delta == 1)~L0+Z+L+W",
                                                   fit = "cox"),
                                   cr1 = list(model = "Surv(tstart, tstop, delta == 5)~L0+Z+L+W",
                                              fit = "cox"),
                                   censoring = list(model = "Surv(tstart, tstop, delta == 0)~L0+Z+L+W",
                                                    fit = "cox")
                               ),
                               max.count = 1,
                               browse = FALSE,
                               verbose = TRUE)


initial.fit$process.types

(test1 <- tmle.alpha.fun(initial.fit = initial.fit,
                         target = "outcome1",
                         tau = 3,
                         alpha = 0.5,
                         #browse = TRUE,
                         verbose = TRUE))

(testz <- tmle.alpha.fun(initial.fit = initial.fit,
                         target = "z",
                         tau = 3,
                         alpha = 0.5,
                        # browse = TRUE,
                         verbose = TRUE))

(test2 <- tmle.alpha.fun(initial.fit = initial.fit,
                         target = "cr1",
                         tau = 3,
                         alpha = 0.5,
                         #browse = TRUE,
                         verbose = TRUE))



#-------------------------------------------------------------------------------------------#
## compute true values 
#-------------------------------------------------------------------------------------------#

N2 <- 1e6
tau <- 3

which.alphas <- 1

for (alpha in which.alphas) {

    int_effect <- lapply(1:10, function(rep) {
        data_alpha <-  simRecurrentTreatment(N = N2,
                                             cens = 0,
                                             eta = c(eta[1:2], eta[3]*alpha, eta[4:length(eta)]), 
                                             nu = nu,
                                             beta_L_Z = 1.5, beta_L_D = 1.25,
                                             beta_W_Z = 1.5, beta_W_D = 1.25,
                                             #beta_Z_D = -0.5, beta_Z_L = -1.25, beta_Z_W = -1.25, beta_Z_D2 = -0.5,
                                             beta_Z_D = 0, beta_Z_L = 0, beta_Z_W = 0, beta_Z_D2 = 0,
                                             beta_L0_L = 1, beta_L0_W = 1,
                                             beta_W_L = 1.5, beta_L_W = 0.5,
                                             beta_L0_Z = 1,
                                             lower = 10^(-150),      
                                             upper = 1e3
                                             )

        data_alpha[, idN := 1:.N, by = "ID"]
        data_alpha[, N := .N, by = "ID"]
        
        #Proportion of subjects dying before some time $\tau$ in treatment group
        prop_alpha <- mean(data_alpha[, sum(Delta == 1 & Time < tau), by = "ID"][[2]]) # with intervention
        prop_alpha_death <- mean(data_alpha[, any(Delta == 5 & Time < tau), by = "ID"][[2]]) # with intervention
        prop_alpha_Z <- mean(data_alpha[, sum(Delta == 2 & Time < tau), by = "ID"][[2]]) # with intervention

        tmp.out <- list(effect_z = prop_alpha_Z,
                        effect_1 = prop_alpha,
                        effect_2 = prop_alpha_death)
         
        if (any(class(tmp.out) == "try-error")) {
            print(paste0("alpha = ", alpha, ", rep = ", rep))
            return(list(effect_z = NA,
                        effect_1 = NA,
                        effect_2 = NA))
        }
         
         return(tmp.out)
     })

     print(truth.out <- c(truth_1 = mean(sapply(int_effect, function(rep) rep$effect_1), na.rm = TRUE),
                          truth_2 = mean(sapply(int_effect, function(rep) rep$effect_2), na.rm = TRUE),
                          truth_z = mean(sapply(int_effect, function(rep) rep$effect_z), na.rm = TRUE)))

    saveRDS(truth.out,
            file = paste0("./output/",
                          "sim-study-recurrent-treatment-",
                          "true-estimands-alpha-", alpha,
                          ".rds"))
     
}



######################################################################
### test.recurrent.treatment.R ends here

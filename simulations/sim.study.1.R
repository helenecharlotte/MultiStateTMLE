### sim.study.1.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Feb  5 2026 (10:12) 
## Version: 
## Last-Updated: Feb  6 2026 (20:21) 
##           By: Helene
##     Update #: 66
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

source("./usage/master.R")
source("./extra/simTreatmentCR.R")

run.master()

#-------------------------------------------------------------------------------------------#
## simulate 1 dataset and apply function 
#-------------------------------------------------------------------------------------------#

if (FALSE) {
    
    N <- 500
    eta <- c(0.1, 0.025, 0.085, 0.1, 0.1, 0.025)

    set.seed(20205)
    
    plotEventData(dt <- simTreatmentCR(N = N, 
                                       eta = eta, 
                                       nu = rep(1.1, 6),
                                       beta_L_Z = 1.5, beta_L_D = 1.25,
                                       beta_W_Z = 1.5, beta_W_D = 1.25,
                                       beta_Z_D = -0.5, beta_Z_L = -1.25, beta_Z_W = -1.25,
                                       beta_L0_L = 1, beta_L0_W = 1,
                                       beta_W_L = 0.5, beta_L_W = 0.5,
                                       beta_L0_Z = 1,
                                       lower = 10^(-150),      
                                       upper = 1e3
                                       ))

    setnames(dt, "ID", "id")
    setnames(dt, "Time", "time")
    setnames(dt, "Delta", "delta")

    initial.fit <- prepare.initial(dt,
                                   tau = 3,
                                   fit.types = list(
                                       Z = list(model = "Surv(tstart, tstop, delta == 2)~L0+L+W",
                                                fit = "cox",
                                                at.risk = function(dt) (dt[["Z"]] == 0)),
                                       L = list(model = "Surv(tstart, tstop, delta == 3)~L0+Z+W",
                                                fit = "cox",
                                                at.risk = function(dt) (dt[["L"]] == 0)),
                                       W = list(model = "Surv(tstart, tstop, delta == 4)~L0+Z+L",
                                                fit = "cox",
                                                at.risk = function(dt) (dt[["W"]] == 0)),
                                       outcome1 = list(model = "Surv(tstart, tstop, delta == 1)~L0+Z+L",
                                                       fit = "cox"),
                                       cr1 = list(model = "Surv(tstart, tstop, delta == 5)~L0+Z+L",
                                                  fit = "cox"),
                                       censoring = list(model = "Surv(tstart, tstop, delta == 0)~L0+Z+L",
                                                        fit = "cox")
                                   ),
                                   browse = FALSE,
                                   verbose = TRUE)

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
                             #browse = TRUE,
                             verbose = TRUE))
}

#-------------------------------------------------------------------------------------------#
## simulation study 
#-------------------------------------------------------------------------------------------#    

N <- 500
M <- 500

use.cores <- 50

eta <- c(0.1, 0.025, 0.085, 0.1, 0.1, 0.025)

which.alphas <- c(0.5, 1, 1.5)

output.weights <- c(1, 0.99, 0.975) 

target.parameters <- c("target", "auxiliary")

which.misspecify <- list(c(), c("1", "L", "W")) #, c("1", "L", "W"), c("1", "W"), c("1", "L")

for (misspecify in list(c(), which.misspecify)) { 

    est.lists <- lapply(which.alphas, function(alpha) list())
    est.aux.lists <- lapply(which.alphas, function(alpha) list())

    names(est.lists) <- which.alphas
    names(est.aux.lists) <- which.alphas
    
    model.1 <- "Surv(tstart, tstop, delta == 1)~L0+Z+W+L"

    if (("1" %in% misspecify) & ("W"  %in% misspecify)) {
        model.1 <- gsub("\\+W", "", model.1)
    }

    if (("1" %in% misspecify) & ("L"  %in% misspecify)) {
        model.1 <- gsub("\\+L", "", model.1)
    }

    model.2 <- "Surv(tstart, tstop, delta == 5)~L0+Z+W+L"

    if (("1" %in% misspecify) & ("W"  %in% misspecify)) {
        model.2 <- gsub("\\+W", "", model.2)
    }

    if (("1" %in% misspecify) & ("L"  %in% misspecify)) {
        model.2 <- gsub("\\+L", "", model.2)
    }

    model.Z <- "Surv(tstart, tstop, delta == 2)~L0+L+W"

    if (("Z" %in% misspecify) & ("W"  %in% misspecify)) {
        model.Z <- gsub("\\+W", "", model.Z)
    }

    if (("Z" %in% misspecify) & ("L"  %in% misspecify)) {
        model.Z <- gsub("\\+L", "", model.Z)
    }

    if ("L" %in% misspecify) {
        model.L <- "Surv(tstart, tstop, delta == 3)~L0+W"
    } else {
        model.L <- "Surv(tstart, tstop, delta == 3)~L0+Z+W"
    }

    if ("W" %in% misspecify) {
        model.W <- "Surv(tstart, tstop, delta == 4)~L0+L"
    } else {
        model.W <- "Surv(tstart, tstop, delta == 4)~L0+Z+L"
    }

    for (m in 1:M) {

        set.seed(m+33339)
        
        dt <- simTreatmentCR(N = N, 
                             eta = eta, 
                             nu = rep(1.1, 6),
                             beta_L_Z = 1.5, beta_L_D = 1.25,
                             beta_W_Z = 1.5, beta_W_D = 1.25,
                             beta_Z_D = -0.5, beta_Z_L = -1.25, beta_Z_W = -1.25,
                             beta_L0_L = 1, beta_L0_W = 1,
                             beta_W_L = 0.5, beta_L_W = 0.5,
                             beta_L0_Z = 1,
                             lower = 10^(-150),      
                             upper = 1e3
                             )

        setnames(dt, "ID", "id")
        setnames(dt, "Time", "time")
        setnames(dt, "Delta", "delta")

        initial.fit <- prepare.initial(dt,
                                       tau = 3,
                                       fit.types = list(
                                           Z = list(model = model.Z,
                                                    fit = "cox",
                                                    at.risk = function(dt) (dt[["Z"]] == 0)),
                                           L = list(model = model.L,
                                                    fit = "cox",
                                                    at.risk = function(dt) (dt[["L"]] == 0)),
                                           W = list(model = model.W,
                                                    fit = "cox",
                                                    at.risk = function(dt) (dt[["W"]] == 0)),
                                           outcome1 = list(model = model.1,
                                                           fit = "cox"),
                                           cr1 = list(model = model.2,
                                                      fit = "cox"),
                                           censoring = list(model = "Surv(tstart, tstop, delta == 0)~L0+Z+L",
                                                            fit = "cox")
                                       ),
                                       browse = FALSE,
                                       verbose = TRUE)

        for (alpha in which.alphas) {

            if ("target" %in% target.parameters) {
                est.lists[[as.character(alpha)]][[m]] <- 
                    tmle.alpha.fun(initial.fit = initial.fit,
                                   target = "outcome1",
                                   tau = 3,
                                   alpha = alpha,
                                   use.cores = use.cores,
                                   output.weights = output.weights,
                                   output.eic = TRUE,
                                   verbose = TRUE)
            
                saveRDS(est.lists[[as.character(alpha)]],
                        file = paste0("./output/",
                                      "sim-study-1-",
                                      "estimation-target-alpha", alpha,
                                      ifelse("1" %in% misspecify, "-misspecify-1", ""),
                                      ifelse("Z" %in% misspecify, "-misspecify-Z", ""),
                                      ifelse("L" %in% misspecify, "-misspecify-L", ""),
                                      ifelse("W" %in% misspecify, "-misspecify-W", ""),
                                      "-M", M, "-N", N, 
                                      ".rds"))

            }

            if ("auxiliary" %in% target.parameters) {
                est.aux.lists[[as.character(alpha)]][[m]] <- 
                    tmle.alpha.fun(initial.fit = initial.fit,
                                   target = "z",
                                   tau = 3,
                                   alpha = alpha,
                                   use.cores = use.cores,
                                   output.weights = output.weights,
                                   output.eic = TRUE,
                                   verbose = TRUE)

                saveRDS(est.aux.lists[[as.character(alpha)]],
                        file = paste0("./output/",
                                      "sim-study-1-",
                                      "estimation-auxiliary-alpha", alpha,
                                      ifelse("1" %in% misspecify, "-misspecify-1", ""),
                                      ifelse("Z" %in% misspecify, "-misspecify-Z", ""),
                                      ifelse("L" %in% misspecify, "-misspecify-L", ""),
                                      ifelse("W" %in% misspecify, "-misspecify-W", ""),
                                      "-M", M, "-N", N, 
                                      ".rds"))
            }
        }
    }
}


#-------------------------------------------------------------------------------------------#
## compute true values 
#-------------------------------------------------------------------------------------------#

N2 <- 1e6
tau <- 3

for (alpha in which.alphas) {

    int_effect <- lapply(1:100, function(rep) {
        data_alpha <- simTreatmentCR(N = N2,
                                     cens = 0,
                                     eta = c(eta[1:2], eta[3]*alpha, eta[4:length(eta)]), 
                                     nu = rep(1.1, 6),
                                     beta_L_Z = 1.5, beta_L_D = 1.25,
                                     beta_W_Z = 1.5, beta_W_D = 1.25,
                                     beta_Z_D = -0.5, beta_Z_L = -1.25, beta_Z_W = -1.25,
                                     beta_L0_L = 1, beta_L0_W = 1,
                                     beta_W_L = 0.5, beta_L_W = 0.5,
                                     beta_L0_Z = 1,
                                     lower = 10^(-150),      
                                     upper = 1e3
                                     )

        data_alpha[, idN := 1:.N, by = "ID"]
        data_alpha[, N := .N, by = "ID"]
        
        #Proportion of subjects dying before some time $\tau$ in treatment group
        prop_alpha <- data_alpha[idN == N, mean(Delta == 1 & Time < tau)] # with intervention
        prop_alpha_Z <- mean(data_alpha[, any(Delta == 2 & Time < tau)[1], by = "ID"][[2]]) # with intervention

         tmp.out <- list(effect_z = prop_alpha_Z,
                         effect_1 = prop_alpha)
         
        if (any(class(tmp.out) == "try-error")) {
            print(paste0("alpha = ", alpha, ", rep = ", rep))
            return(list(effect_z = NA,
                        effect_1 = NA))
        }
         
         return(tmp.out)
     })

     print(truth.out <- c(truth_1 = mean(sapply(int_effect, function(rep) rep$effect_1), na.rm = TRUE),
                          truth_z = mean(sapply(int_effect, function(rep) rep$effect_z), na.rm = TRUE)))

    saveRDS(truth.out,
            file = paste0("./output/",
                          "sim-study-1-",
                          "true-estimands-alpha-", alpha,
                          ".rds"))
     
}

######################################################################
### sim.study.1.R ends here

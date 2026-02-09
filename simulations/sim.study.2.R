### sim.study.2.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Feb  5 2026 (11:25) 
## Version: 
## Last-Updated: Feb  6 2026 (08:08) 
##           By: Helene
##     Update #: 71
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

source("./usage/master.R")
source("./extra/simDropin.R")

run.master()


#-------------------------------------------------------------------------------------------#
## simulate 1 dataset and apply function 
#-------------------------------------------------------------------------------------------#

if (FALSE) {
    
    N <- 500

    set.seed(20205)
    
    plotEventData((dt <- simDropin(N = 500,
                                   beta_L_D = 1.1,
                                   beta_L_Z = 5,
                                   beta_Z_D = -3,
                                   beta_Z_L = -3,
                                   beta_A0_L = -2.75,
                                   eta = c(0.25, 0.5*0.05, 0.05, 4*0.25),
                                   nu = c(1.1, 0.8*1.0, 0.2*0.2, 1.1),
                                   beta_A0_D = -0.5))[Time <= 4])

    setnames(dt, "ID", "id")
    setnames(dt, "Time", "time")
    setnames(dt, "Delta", "delta")

    initial.fit <- prepare.initial(dt,
                                   tau = 3,
                                   fit.types = list(
                                       Z = list(model = "Surv(tstart, tstop, delta == 2)~L0+A0+L",
                                                fit = "cox",
                                                at.risk = function(dt) (dt[["Z"]] == 0)),
                                       L = list(model = "Surv(tstart, tstop, delta == 3)~L0+A0+Z",
                                                fit = "cox",
                                                at.risk = function(dt) (dt[["L"]] == 0)),
                                       outcome1 = list(model = "Surv(tstart, tstop, delta == 1)~L0+A0+Z+L",
                                                       fit = "cox"),
                                       censoring = list(model = "Surv(tstart, tstop, delta == 0)~L0+A0+Z+L",
                                                        fit = "cox")
                                   ),
                                   fit.treatment = "A0~L0",
                                   a = 0:1,
                                   browse = FALSE,
                                   verbose = TRUE)

    (test1.a1 <- tmle.alpha.fun(initial.fit = initial.fit,
                                target = "outcome1",
                                tau = 3,
                                a = 1,
                                alpha = 0.5,
                                #browse = TRUE,
                                verbose = TRUE))

    (testz.a1 <- tmle.alpha.fun(initial.fit = initial.fit,
                                target = "z",
                                tau = 3,
                                a = 1,
                                alpha = 0.5,
                                #browse = TRUE,
                                verbose = TRUE))

    (test1.a0 <- tmle.alpha.fun(initial.fit = initial.fit,
                                target = "outcome1",
                                tau = 3,
                                a = 0,
                                alpha = 0.5,
                                #browse = TRUE,
                                verbose = TRUE))

    (testz.a0 <- tmle.alpha.fun(initial.fit = initial.fit,
                                target = "z",
                                tau = 3,
                                a = 0,
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
which.as <- 1:0

output.weights <- c(1, 0.99, 0.975) 

target.parameters <- c("target", "auxiliary")

which.misspecify <- list(c(), c("1", "L"))

for (a in which.as) {
    
    for (misspecify in which.misspecify) {

        est.lists <- lapply(which.alphas, function(alpha) list())
        est.aux.lists <- lapply(which.alphas, function(alpha) list())

        names(est.lists) <- which.alphas
        names(est.aux.lists) <- which.alphas
    
        model.1 <- "Surv(tstart, tstop, delta == 1)~L0+A0+Z+L"

        if (("1" %in% misspecify) & ("L"  %in% misspecify)) {
            model.1 <- gsub("\\+A0", "", gsub("\\+L", "", model.1))
        }

        model.Z <- "Surv(tstart, tstop, delta == 2)~L0+A0+L"

        if (("Z" %in% misspecify) & ("L"  %in% misspecify)) {
            model.Z <- gsub("\\+L", "", model.Z)
        }

        model.L <- "Surv(tstart, tstop, delta == 3)~L0+A0+Z"

        if ("L" %in% misspecify) {
            model.L <- "Surv(tstart, tstop, delta == 3)~L0"
        }

        for (m in 1:M) {

            set.seed(m+33339)
        
            dt <- simDropin(N = N,
                            beta_L_D = 1.1,
                            beta_L_Z = 5,
                            beta_Z_D = -3,
                            beta_Z_L = -3,
                            beta_A0_L = -2.75,
                            eta = c(0.25, 0.5*0.05, 0.05, 4*0.25),
                            nu = c(1.1, 0.8*1.0, 0.2*0.2, 1.1),
                            beta_A0_D = -0.5)

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
                                               outcome1 = list(model = model.1,
                                                               fit = "cox"),
                                               censoring = list(model = "Surv(tstart, tstop, delta == 0)~L0+A0+Z+L",
                                                                fit = "cox")
                                           ),
                                           fit.treatment = "A0~L0",
                                           a = 0:1,
                                           browse = FALSE,
                                           verbose = TRUE)

            for (alpha in which.alphas) {

                if ("target" %in% target.parameters) {
                    est.lists[[as.character(alpha)]][[m]] <- 
                        tmle.alpha.fun(initial.fit = initial.fit,
                                       target = "outcome1",
                                       tau = 3,
                                       alpha = alpha,
                                       a = a, 
                                       use.cores = use.cores,
                                       output.weights = output.weights,
                                       output.eic = TRUE,
                                       verbose = TRUE)
            
                    saveRDS(est.lists[[as.character(alpha)]],
                            file = paste0("./output/",
                                          "sim-study-2-",
                                          "estimation-target-alpha", alpha, "-a", a,
                                          ifelse("1" %in% misspecify, "-misspecify-1", ""),
                                          ifelse("Z" %in% misspecify, "-misspecify-Z", ""),
                                          ifelse("L" %in% misspecify, "-misspecify-L", ""),
                                          "-M", M, "-N", N, 
                                          ".rds"))

                }

                if ("auxiliary" %in% target.parameters) {
                    est.aux.lists[[as.character(alpha)]][[m]] <- 
                        tmle.alpha.fun(initial.fit = initial.fit,
                                       target = "z",
                                       tau = 3,
                                       alpha = alpha,
                                       a = a, 
                                       use.cores = use.cores,
                                       output.weights = output.weights,
                                       output.eic = TRUE,
                                       verbose = TRUE)

                    saveRDS(est.aux.lists[[as.character(alpha)]],
                            file = paste0("./output/",
                                          "sim-study-2-",
                                          "estimation-auxiliary-alpha", alpha, "-a", a, 
                                          ifelse("1" %in% misspecify, "-misspecify-1", ""),
                                          ifelse("Z" %in% misspecify, "-misspecify-Z", ""),
                                          ifelse("L" %in% misspecify, "-misspecify-L", ""),
                                          "-M", M, "-N", N, 
                                          ".rds"))
                }
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

    for (a in 0:1) {

        int_effect <- lapply(1:10, function(rep) {
            data_alpha <- simDropin(N = N2, 
                                    cens = 0,
                                    generate.A0 = function(N,L0) rep(a,N),
                                    beta_L_D = 1.1,
                                    beta_L_Z = 5,
                                    beta_Z_D = -3,
                                    beta_Z_L = -3,
                                    beta_A0_L = -2.75,
                                    eta = c(0.25, 0.5*0.05, alpha*0.05, 4*0.25),
                                    nu = c(1.1, 0.8*1.0, 0.2*0.2, 1.1),
                                    beta_A0_D = -0.5,
                                    max_iter = 1e3)

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
                              "sim-study-2-",
                              "true-estimands-alpha-", alpha, "-a", a,
                              ".rds"))

    }
}



######################################################################
### sim.study.2.R ends here

### tmle.alpha.fun.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Feb  4 2026 (12:51) 
## Version: 
## Last-Updated: Sep  1 2026 (06:21) 
##           By: Helene
##     Update #: 499
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

tmle.alpha.fun <- function(target = "z",
                           tau = 1.2,
                           alpha = 1, z.name = "z", ## old options (but new one not made yet)
                           alpha.list = NULL, ## new option (not made yet)
                           ## should be on the form: alpha.list(z = function(time, covariates))
                           a = NULL,
                           initial.fit = NULL,
                           dt = NULL,
                           years.lost = NULL,
                           only.first = NULL,
                           target.only.in.state = NULL, # function(states) (states[["L"]]==states[["Z"]]),
                           target.by.state = FALSE,
                           conv.const = 1,
                           one.step = FALSE,
                           verbose = FALSE,
                           max.iter = 10,
                           min.iter = 0,
                           use.cores = 50,
                           truncate.weights = 0,
                           output.convergence = FALSE,
                           output.eic = FALSE,
                           output.weights = NULL,
                           output.a.weights = NULL,
                           browse = FALSE,
                           verbose.exponential = FALSE,
                           ...) {

    if (length(initial.fit) == 0) {
        initial.fit <-
            prepare.initial(tau = tau,
                            a = a,
                            verbose = verbose,
                            ...)
    }

    tmp.long <- data.table::copy(initial.fit$tmp.long)
    depend.matrix <- data.table::copy(initial.fit$depend.matrix)

    n <- length(unique(tmp.long[["id"]]))

    process.names <- initial.fit$process.names
    process.deltas <- initial.fit$process.deltas
    process.types <- initial.fit$process.types

    cens.process.id <- initial.fit$cens.process.id

    at.risks <- initial.fit$at.risks
      
    z.process.id <- (1:length(process.names))[tolower(process.names) == tolower(z.name)]

    if (length(alpha.list)>0) {
        message("using new alpha option: check")
        for (alpha.kk in 1:length(alpha.list)) {
            tmp.long[, paste0("alpha.", names(alpha.list)[alpha.kk]) := #
                                do.call(alpha.list[[alpha.kk]], .SD), .SDcols = names(formals(alpha.list[[alpha.kk]]))]
        }
        z.present <- TRUE
    } 

    if (length(alpha.list) == 0) {
        if (length(z.process.id) == 0) {
            z.present <- FALSE
            z.process.id <- Inf
            z.name <- "not.applied"
            if (!identical(alpha, 1)) warning("No 'z' process found in initial.fit; alpha-scaling will be ignored (treated as 1).")
        } else {
            z.present <- TRUE
            z.name <- process.names[z.process.id]
        }
    }

    #--------------------------------    
    #-- initialization step:
    
    states <- copy(depend.matrix)
    state.names <- names(states)
    setkeyv(states, c("state", state.names[state.names %in% process.names]))

    target.id <- (1:length(process.names))[tolower(target) == tolower(process.names)]
    target.name <- process.names[target.id]
    
    if (length(only.first)>0) {
        if (only.first) { ## & z.present
            ## process.types[[z.name]] <- "one.jump"
            ## tmp.long <- tmp.long[tmp.long[[z.name]] <= 1]
            process.types[[target.name]] <- "one.jump"
            tmp.long <- tmp.long[tmp.long[[target.name]] <= 1]
        }
    }

    tmp.long[, at.risk := (time <= final.time)]

    if (target.name %in% names(at.risks)) {
        # states <- states[at.risks[[target.name]](states)]
        if (process.types[target.name] == "one.jump") { # <- fixme
            tmp.long[, at.risk := at.risk*(at.risks[[target.name]](tmp.long))]
        }
        states[, at.risk := at.risks[[target.name]](depend.matrix)]
    }

    if (length(only.first)>0 & z.present) {
        if (only.first) {
            tmp.long[, at.risk := at.risk*(tmp.long[[z.name]] <= 1)]
            states[, at.risk := at.risk & (states[[z.name]] <= 1)]
        }
    }

    if (length(target.only.in.state)>0) {
        states[, target.only.in.state := target.only.in.state(states)]
    }

    for (process.name in process.names[!(process.names %in% names(at.risks))]) {
        tmp.long[, (paste0("at.risk.", process.name)) := at.risk]
    }
    
    P.prefix <- ifelse(length(a)>0, paste0("P.a", a, "."), paste0("P."))

    if (length(alpha.list)>0) {
        for (alpha.kk in 1:length(alpha.list)) {
            tmp.long[, paste0("alpha.", names(alpha.list)[alpha.kk]) := #
                                do.call(alpha.list[[alpha.kk]], .SD), .SDcols = names(formals(alpha.list[[alpha.kk]]))]
        }
    } 
    
    #--------------------------------    
    #-- compute clever weights:

    if (z.present) {
        if (length(alpha.list)>0) {
            tmp.long[, clever.weight.alpha := 1]
            for (alpha.kk in 1:length(alpha.list)) {
                tmp.long[, (paste0("cum.hazard.", names(alpha.list)[alpha.kk])) := cumsum(get(paste0("at.risk.", names(alpha.list)[alpha.kk]))*
                                                                                          get(paste0("P.", names(alpha.list)[alpha.kk]))), by = "id"]
                tmp.long[, (paste0("cum.hazard.", names(alpha.list)[alpha.kk], ".1")) :=
                               c(0, get(paste0("cum.hazard.", names(alpha.list)[alpha.kk]))[-.N]),
                         by = "id"]
                tmp.long[, clever.weight.alpha := clever.weight.alpha*
                               get(paste0("alpha.", names(alpha.list)[alpha.kk]))^get(names(alpha.list)[alpha.kk])*
                               exp(-(get(paste0("alpha.", names(alpha.list)[alpha.kk]))-1)*
                                   get(paste0("cum.hazard.", names(alpha.list)[alpha.kk], ".1")))]
            }
        } else {
            tmp.long[, cum.hazard.z := cumsum(get(paste0("at.risk.", process.names[z.process.id]))*get(paste0("P.", process.names[z.process.id]))), by = "id"]
            tmp.long[, cum.hazard.z.1 := c(0, cum.hazard.z[-.N]), by = "id"]
            tmp.long[, clever.weight.alpha := alpha^get(z.name)*exp(-(alpha-1)*cum.hazard.z.1)]
        }
        if (truncate.weights>0) {
            no.truncated <- tmp.long[time <= final.time & clever.weight.alpha > truncate.weights, length(unique(id))]
            tmp.long[time <= final.time & clever.weight.alpha > truncate.weights, clever.weight.alpha := truncate.weights]
        }
    } else {
        tmp.long[, clever.weight.alpha := 1]
    }
    
    if (length(a)>0) {
        tmp.long[, clever.weight := tmp.long[[paste0("clever.weight.a", a)]]]
    }

    if (length(output.weights)>0) {
        ### alpha weights
        tmp.weight <- tmp.long[time <= final.time, max(clever.weight.alpha), by = "id"]
        out.weights <- sapply(output.weights, function(ppp) {
            as.numeric(quantile(tmp.weight[[2]], p = ppp))
        })
        names(out.weights) <- c(paste0("q", output.weights*100))
        if (truncate.weights) out.weights <- c(out.weights, no.truncated = no.truncated)
        ### censoring*a weights
        tmp.cens.weight <- tmp.long[time <= final.time, max(clever.weight), by = "id"]
        out.cens.weights <- sapply(output.weights, function(ppp) {
            as.numeric(quantile(tmp.cens.weight[[2]], p = ppp))
        })
        names(out.cens.weights) <- c(paste0("q", output.weights*100))
    }

    if (length(output.a.weights)>0) {
        tmp.weight <- tmp.long[time <= final.time, max(clever.weight), by = "id"]
        out.a.weights <- sapply(output.a.weights, function(ppp) {
            as.numeric(quantile(tmp.weight[[2]], p = ppp))
        })
        names(out.a.weights) <- c(paste0("q", output.a.weights*100))
    }
           
    #--------------------------------    
    #-- counterfactual P.Z:

    if (length(alpha.list)>0) {
        for (alpha.kk in 1:length(alpha.list)) {
            for (varname in grep(paste0(P.prefix, names(alpha.list)[alpha.kk]), names(tmp.long), value = TRUE)) {
                if (varname != paste0(P.prefix, names(alpha.list)[alpha.kk]))
                    tmp.long[[varname]] <- tmp.long[[paste0("alpha.", names(alpha.list)[alpha.kk])]]*tmp.long[[varname]]
            }
        }
    } else {
        for (varname in grep(paste0(P.prefix, z.name), names(tmp.long), value = TRUE)) {
            if (varname != paste0(P.prefix, z.name))
                tmp.long[[varname]] <- alpha*tmp.long[[varname]]
        }
    }

    #--------------------------------    
    #-- tmle iterations:

    converged <- FALSE
    clever.Q.label <- ifelse(length(years.lost)>0, "clever.Q.years.lost.", "clever.Q.")
    
    if (browse) browser()
    
    for (iter in 1:max.iter) {

        print(paste0("iter = ", iter))
        
        setkey(tmp.long, id, time)

        tmp.long[["gamma.mapping.warning"]] <- (tmp.long[["id"]] == 1)*(iter == 1)*(verbose) 

        dt_list <- split(tmp.long[, !(names(tmp.long) %in% c("Q",
                                                             grep(clever.Q.label, names(tmp.long), value = TRUE))),
                                  with = FALSE], by = "id", keep.by = TRUE)
       
        t2 <- system.time({
            dt_list <- mclapply(
                dt_list,
                compute.Q.clever.per.id,
                states = states,
                process.types = process.types,
                P.prefix = P.prefix,
                parameter = target.name,
                get.years.lost = (length(years.lost)>0),
                years.lost.block.size = years.lost,
                clever.by.state = target.by.state,
                mc.cores = min(detectCores()-1, use.cores)
            )
        })

        if (verbose) print(t2)

        tmp.long <- rbindlist(dt_list)

        #-- current estimator for target parameter:  
        target.est <- mean(tmp.long[, Q[1], by = "id"][[2]])

        if (iter == 1) {
            g.est <- target.est
        }

        eic <- tmp.long[at.risk == 1,
                        Q[1]-target.est,
                        by = "id"][[2]]

        clever.ids <- (1:length(process.names))[sapply(process.names, function(process.name) paste0(clever.Q.label, process.name) %in% names(tmp.long))]

        for (process.jj in (1:length(process.names))[clever.ids]) {
            name.jj <- process.names[process.jj]
            if (length(alpha.list)>0) {
                if (name.jj %in% names(alpha.list)) {
                    eic <- eic + tmp.long[at.risk == 1, sum(get(paste0("alpha.", name.jj))*clever.weight*clever.weight.alpha*(
                        get(paste0(clever.Q.label, process.names[process.jj])))*((delta == process.deltas[process.jj]) - get(paste0("at.risk.", name.jj))*get(paste0("P.", name.jj)))),
                        by = "id"][[2]]
                } else {
                    eic <- eic + tmp.long[at.risk == 1, sum(clever.weight*clever.weight.alpha*(
                        get(paste0(clever.Q.label, process.names[process.jj])))*((delta == process.deltas[process.jj]) - get(paste0("at.risk.", name.jj))*get(paste0("P.", name.jj)))),
                        by = "id"][[2]]
                }
            } else {
                eic <- eic + tmp.long[at.risk == 1, sum(alpha^(process.jj == z.process.id)*clever.weight*clever.weight.alpha*(
                    get(paste0(clever.Q.label, process.names[process.jj])))*((delta == process.deltas[process.jj]) - get(paste0("at.risk.", name.jj))*get(paste0("P.", name.jj)))),
                    by = "id"][[2]]
            }
        }
        
        if (iter == 1) {
            eic.init <- copy(eic)
            target.se <- sqrt(mean(eic^2/n))
        }

        if (one.step) {
            one.step.est <- mean(eic + target.est)
            break()
        }

        ## if (iter == 2) browser()

        print(paste0("eic equation solved at = ", abs(mean(eic))))

        if (iter>min.iter) {
            if (abs(mean(eic)) <= conv.const*target.se/(log(n))) {
                converged <- TRUE
                print(paste0("converged after ", iter, " iterations"))
                break()
            }
        }

        if (length(alpha.list)>0) {
            if (target.by.state) {
                target.fun <- function(eps, process.jj) {
                    name.jj <- process.names[process.jj]
                    if (name.jj %in% names(alpha.list)) {
                        mean(tmp.long[at.risk == 1, sum(get(paste0("alpha.", name.jj))*clever.weight*clever.weight.alpha*(
                            get(paste0(clever.Q.label, process.names[process.jj])))*((delta == process.deltas[process.jj]) - get(paste0("at.risk.", name.jj))*get(paste0("P.", name.jj))*exp(eps*(
                                get(paste0(clever.Q.label, process.names[process.jj])))))), by = "id"][[2]])
                    } else {
                        mean(tmp.long[at.risk == 1, sum(clever.weight*clever.weight.alpha*(
                            get(paste0(clever.Q.label, process.names[process.jj])))*((delta == process.deltas[process.jj]) - get(paste0("at.risk.", name.jj))*get(paste0("P.", name.jj))*exp(eps*(
                                get(paste0(clever.Q.label, process.names[process.jj])))))), by = "id"][[2]])
                    }
                }
            } else {
                target.fun <- function(eps, process.jj) {
                    name.jj <- process.names[process.jj]
                    if (name.jj %in% names(alpha.list)) {
                        mean(tmp.long[at.risk == 1, sum(get(paste0("alpha.", name.jj))*clever.weight*clever.weight.alpha*(
                            get(paste0(clever.Q.label, process.names[process.jj])))*((delta == process.deltas[process.jj]) - get(paste0("at.risk.", name.jj))*get(paste0("P.", name.jj))*exp(eps))), by = "id"][[2]])
                    } else {
                        mean(tmp.long[at.risk == 1, sum(clever.weight*clever.weight.alpha*(
                            get(paste0(clever.Q.label, process.names[process.jj])))*((delta == process.deltas[process.jj]) - get(paste0("at.risk.", name.jj))*get(paste0("P.", name.jj))*exp(eps))), by = "id"][[2]])
                    }
                }
            }
        } else { 
            if (target.by.state) {
                target.fun <- function(eps, process.jj) {
                    name.jj <- process.names[process.jj]
                    mean(tmp.long[at.risk == 1, sum((alpha)^(process.jj == z.process.id)*clever.weight*clever.weight.alpha*(
                        get(paste0(clever.Q.label, process.names[process.jj])))*((delta == process.deltas[process.jj]) - get(paste0("at.risk.", name.jj))*get(paste0("P.", name.jj))*exp(eps*(
                            get(paste0(clever.Q.label, process.names[process.jj])))))), by = "id"][[2]])
                }
            } else {
                target.fun <- function(eps, process.jj) {
                    name.jj <- process.names[process.jj]
                    mean(tmp.long[at.risk == 1, sum((alpha)^(process.jj == z.process.id)*clever.weight*clever.weight.alpha*(
                        get(paste0(clever.Q.label, process.names[process.jj])))*((delta == process.deltas[process.jj]) - get(paste0("at.risk.", name.jj))*get(paste0("P.", name.jj))*exp(eps))), by = "id"][[2]])
                }
            }
        }

        ## browser()
        
        for (process.jj in (1:length(process.names))[clever.ids]) {
            name.jj <- process.names[process.jj]
            eps.jj <- nleqslv(0.00, function(eps) target.fun(eps, process.jj))$x
            if (verbose) print(paste0("eps.", process.names[process.jj], " = ", eps.jj))
            if (FALSE & abs(eps.jj) > 10) {
                eps.jj <- eps.jj/10
                ## message(print(paste0("eps.", process.names[process.jj], " ill-defined")))
                ## if (verbose) print(paste0("eps.", process.names[process.jj], " = ", eps.jj))
            }
            if (target.by.state) {
                tmp.long[[(paste0(P.prefix, name.jj))]] <-
                    tmp.long[[(paste0(P.prefix, name.jj))]]*exp(eps.jj*tmp.long[[paste0(clever.Q.label, process.names[process.jj])]])
            } else {
                tmp.long[[(paste0(P.prefix, name.jj))]] <-
                    tmp.long[[(paste0(P.prefix, name.jj))]]*exp(eps.jj)
            }
            if (P.prefix != "P.") {
                if (target.by.state) {
                    tmp.long[[paste0("P.", name.jj)]] <- tmp.long[[paste0("P.", name.jj)]] *
                        exp(eps.jj * tmp.long[[paste0(clever.Q.label, name.jj)]])
                } else {
                    tmp.long[[(paste0("P.", name.jj))]] <- tmp.long[[(paste0("P.", name.jj))]]*exp(eps.jj)
                }
            }
            for (state.jj in depend.matrix[, unique(state)]) {
                if (target.by.state) {
                    tmp.long[[(paste0(P.prefix, name.jj, ".", state.jj))]] <-
                        tmp.long[[(paste0(P.prefix, name.jj, ".", state.jj))]]*exp(eps.jj*tmp.long[[paste0(clever.Q.label, process.names[process.jj], ".", state.jj)]])
                } else {
                    tmp.long[[(paste0(P.prefix, name.jj, ".", state.jj))]] <-
                        tmp.long[[(paste0(P.prefix, name.jj, ".", state.jj))]]*exp(eps.jj)
                }
            }
        }

        if (length(alpha.list)>0) {
            tmp.long[, clever.weight.alpha := 1]
            for (alpha.kk in 1:length(alpha.list)) {
                tmp.long[, (paste0("cum.hazard.", names(alpha.list)[alpha.kk])) := cumsum(get(paste0("at.risk.", names(alpha.list)[alpha.kk]))*
                                                                                          get(paste0("P.", names(alpha.list)[alpha.kk]))), by = "id"]
                tmp.long[, (paste0("cum.hazard.", names(alpha.list)[alpha.kk], ".1")) :=
                               c(0, get(paste0("cum.hazard.", names(alpha.list)[alpha.kk]))[-.N]),
                         by = "id"]
                tmp.long[, clever.weight.alpha := clever.weight.alpha*
                               get(paste0("alpha.", names(alpha.list)[alpha.kk]))^get(names(alpha.list)[alpha.kk])*
                               exp(-(get(paste0("alpha.", names(alpha.list)[alpha.kk]))-1)*
                                   get(paste0("cum.hazard.", names(alpha.list)[alpha.kk], ".1")))]
            }
        } else if (z.present) {
            tmp.long[, cum.hazard.z := cumsum(get(paste0("at.risk.", process.names[z.process.id]))*get(paste0("P.", process.names[z.process.id]))), by = "id"]
            tmp.long[, cum.hazard.z.1 := c(0, cum.hazard.z[-.N]), by = "id"]
            tmp.long[, clever.weight.alpha := (alpha)^get(z.name)*
                           exp(-(alpha-1)*cum.hazard.z.1)]
        }

        names.P <- names(tmp.long)[substr(names(tmp.long), 1, 2) == "P."]

        for (name.P in names.P) {
            if (any(tmp.long[[name.P]]>1)) {
                if (verbose.exponential) print(paste0("transform ", name.P, " with 1-exp to avoid values >1"))
                tmp.long[, (name.P) := 1-exp(-tmp.long[[name.P]])]
            }
        }

    }

    if (one.step) {
        out <- list(estimate = c(one.step.est = one.step.est, se = target.se,
                                 g.est = g.est))
    } else {
        out <- list(estimate = c(tmle.est = target.est, se = target.se,
                                 g.est = g.est))
    }

    if (output.convergence) {
        out[[length(out)+1]] <- c(iter = iter, 
                                  eic.solved.at = abs(mean(eic)),
                                  converged = converged)
        names(out)[length(out)] <- "convergence"
    }

    if (length(output.weights)>0) {
        out[[length(out)+1]] <- out.weights
        names(out)[length(out)] <- "weights"
        out[[length(out)+1]] <- out.cens.weights
        names(out)[length(out)] <- "cens.weights"
    }

    if (length(output.a.weights)>0) {
        out[[length(out)+1]] <- out.a.weights
        names(out)[length(out)] <- "a.weights"
    }

    if (output.eic) {
        out[[length(out)+1]] <- eic.init
        names(out)[length(out)] <- "eic"
    }

    return(out)
}

######################################################################
### tmle.alpha.fun.R ends here

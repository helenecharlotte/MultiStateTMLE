### tmle.alpha.fun.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Feb  4 2026 (12:51) 
## Version: 
## Last-Updated: Feb 10 2026 (12:20) 
##           By: Helene
##     Update #: 209
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
                           alpha = 1, z.name = "z",
                           a = NULL,
                           initial.fit = NULL,
                           verbose = FALSE,
                           max.iter = 10,
                           use.cores = 50,
                           truncate.weights = 0,
                           output.convergence = FALSE,
                           output.eic = FALSE,
                           output.weights = NULL,
                           browse = FALSE,
                           ...) {

    if (length(initial.fit) == 0) {
        initial.fit <-
            prepare.initial(tau = tau,
                            a = a,
                            verbose = verbose,
                            ...)
    }

    tmp.long <- initial.fit$tmp.long

    n <- length(unique(tmp.long[["id"]]))

    depend.matrix <- initial.fit$depend.matrix

    process.names <- initial.fit$process.names
    process.deltas <- initial.fit$process.deltas
    process.types <- initial.fit$process.types

    cens.process.id <- initial.fit$cens.process.id

    at.risks <- initial.fit$at.risks
      
    z.process.id <- (1:length(process.names))[tolower(process.names) == z.name]

    if (length(z.process.id) == 0) {
        z.present <- FALSE
        z.process.id <- Inf
        z.name <- "not.applied"
        if (!identical(alpha, 1)) warning("No 'z' process found in initial.fit; alpha-scaling will be ignored (treated as 1).")
    } else {
        z.present <- TRUE
        z.name <- process.names[z.process.id]
    }

    #--------------------------------    
    #-- initialization step:
    
    states <- copy(depend.matrix)
    state.names <- names(states)
    setkeyv(states, c("state", state.names[state.names %in% process.names]))

    target.id <- (1:length(process.names))[tolower(target) == tolower(process.names)]
    target.name <- process.names[target.id]

    tmp.long[, at.risk := (time <= final.time)]

    if (target.name %in% names(at.risks)) {
        # states <- states[at.risks[[target.name]](states)]
        tmp.long[, at.risk := at.risk*(at.risks[[target.name]](tmp.long))]
        states[, at.risk := at.risks[[target.name]](depend.matrix)]
    }

    for (process.name in process.names[!(process.names %in% names(at.risks))]) {
        tmp.long[, (paste0("at.risk.", process.name)) := at.risk]
    }
    
    P.prefix <- ifelse(length(a)>0, paste0("P.a", a, "."), paste0("P."))

    #--------------------------------    
    #-- compute clever weights:

    if (z.present) {
        tmp.long[, cum.hazard.z := cumsum(get(paste0("at.risk.", process.names[z.process.id]))*get(paste0("P.", process.names[z.process.id]))), by = "id"]
        tmp.long[, cum.hazard.z.1 := c(0, cum.hazard.z[-.N]), by = "id"]
        tmp.long[, clever.weight.alpha := alpha^get(z.name)*exp(-(alpha-1)*cum.hazard.z.1)]
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
        tmp.weight <- tmp.long[time <= final.time, max(clever.weight.alpha), by = "id"]
        out.weights <- sapply(output.weights, function(ppp) {
            as.numeric(quantile(tmp.weight[[2]], p = ppp))
        })
        names(out.weights) <- c(paste0("q", output.weights*100))
        if (truncate.weights) out.weights <- c(out.weights, no.truncated = no.truncated)
    }
           
    #--------------------------------    
    #-- counterfactual P.Z:

    for (varname in grep(paste0(P.prefix, z.name), names(tmp.long), value = TRUE)) {
        if (varname != paste0(P.prefix, z.name))
            tmp.long[[varname]] <- alpha*tmp.long[[varname]]
    }

    #--------------------------------    
    #-- tmle iterations:

    converged <- FALSE
    if (browse) browser()
    
    for (iter in 1:max.iter) {

        print(paste0("iter = ", iter))
        
        setkey(tmp.long, id, time)

        dt_list <- split(tmp.long[, !(names(tmp.long) %in% c("Q",
                                                             grep("clever.Q", names(tmp.long), value = TRUE))),
                                  with = FALSE], by = "id", keep.by = TRUE)

        if (FALSE) {
            dt_list <- tmp.long[, !(names(tmp.long) %in% c("Q",
                                                           grep("clever.Q", names(tmp.long), value = TRUE))),
                                with = FALSE]

            profvis(compute.Q.clever.per.id(dt_list[id <= 40], states = states,
                                            process.types = process.types,
                                            P.prefix = P.prefix,
                                            #browse = TRUE,
                                            parameter = target.name))

            compute.Q.clever.per.id(dt_list[id <= 40], states = states,
                                    process.types = process.types,
                                    P.prefix = P.prefix,
                                    browse = TRUE,
                                    parameter = target.name)
        }

        
        t2 <- system.time({
            dt_list <- mclapply(
                dt_list,
                compute.Q.clever.per.id,
                states = states,
                process.types = process.types,
                P.prefix = P.prefix,
                parameter = target.name,
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

        clever.ids <- (1:length(process.names))[sapply(process.names, function(process.name) paste0("clever.Q.", process.name, "0") %in% names(tmp.long))]

        for (process.jj in (1:length(process.names))[clever.ids]) {
            name.jj <- process.names[process.jj]
            eic <- eic + tmp.long[at.risk == 1, sum(alpha^(process.jj == z.process.id)*clever.weight*clever.weight.alpha*(
                get(paste0("clever.Q.", process.names[process.jj], "1")) - get(paste0("clever.Q.", process.names[process.jj], "0")))*((delta == process.deltas[process.jj]) - get(paste0("at.risk.", name.jj))*get(paste0("P.", name.jj)))),
                by = "id"][[2]]
        }
        
        if (iter == 1) {
            eic.init <- copy(eic)
            target.se <- sqrt(mean(eic^2/n))
        }

        print(paste0("eic equation solved at = ", abs(mean(eic))))

        if (abs(mean(eic)) <= target.se/(log(n))) {
            converged <- TRUE
            break(print(paste0("finished after ", iter, " iterations")))
        }

        target.fun <- function(eps, process.jj) {
            name.jj <- process.names[process.jj]
            mean(tmp.long[at.risk == 1, sum((alpha)^(process.jj == z.process.id)*clever.weight*clever.weight.alpha*(
                get(paste0("clever.Q.", process.names[process.jj], "1")) - get(paste0("clever.Q.", process.names[process.jj], "0")))*((delta == process.deltas[process.jj]) - get(paste0("at.risk.", name.jj))*get(paste0("P.", name.jj))*exp(eps))), by = "id"][[2]])
        }

        for (process.jj in (1:length(process.names))[clever.ids]) {
            name.jj <- process.names[process.jj]
            eps.jj <- nleqslv(0.00, function(eps) target.fun(eps, process.jj))$x
            if (verbose) print(paste0("eps.", process.names[process.jj], " = ", eps.jj))
            tmp.long[[(paste0(P.prefix, name.jj))]] <- tmp.long[[(paste0(P.prefix, name.jj))]]*exp(eps.jj)
            if (P.prefix != "P.") tmp.long[[(paste0("P.", name.jj))]] <- tmp.long[[(paste0("P.", name.jj))]]*exp(eps.jj)
            for (state.jj in depend.matrix[, unique(state)]) {
                tmp.long[[(paste0(P.prefix, name.jj, ".", state.jj))]] <- tmp.long[[(paste0(P.prefix, name.jj, ".", state.jj))]]*exp(eps.jj)
            }
        }

        if (z.present) {
            tmp.long[, cum.hazard.z := cumsum(get(paste0("at.risk.", process.names[z.process.id]))*get(paste0("P.", process.names[z.process.id]))), by = "id"]
            tmp.long[, cum.hazard.z.1 := c(0, cum.hazard.z[-.N]), by = "id"]
            tmp.long[, clever.weight.alpha := alpha^get(z.name)*exp(-(alpha-1)*cum.hazard.z.1)]
        }

    }

    out <- list(estimate = c(tmle.est = target.est, se = target.se,
                             g.est = g.est))

    if (output.convergence) {
        out[[length(out)+1]] <- c(iter = iter, 
                                  eic.solved.at = abs(mean(eic)),
                                  converged = converged)
        names(out)[length(out)] <- "convergence"
    }

    if (length(output.weights)>0) {
        out[[length(out)+1]] <- out.weights
        names(out)[length(out)] <- "weights"
    }

    if (output.eic) {
        out[[length(out)+1]] <- eic.init
        names(out)[length(out)] <- "eic"
    }

    return(out)
}

######################################################################
### tmle.alpha.fun.R ends here

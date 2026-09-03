### sim.generic.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Aug 28 2026 (09:46) 
## Version: 
## Last-Updated: Aug 30 2026 (19:50) 
##           By: Helene
##     Update #: 290
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

sim.generic <- function(baseline = list(),
                        processes = list(),
                        effects = list(),
                        sim.object = list(),
                        cens = 1,
                        alpha.intervention = list(),
                        baseline.intervention = list(),
                        n = 500,
                        browse = FALSE) {

    if (length(baseline) == 0 | length(processes) == 0 | length(effects) == 0) {
        baseline <- sim.object$baseline
        processes <- sim.object$processes
        effects <- sim.object$effects
    }
    
    baseline.vars <- names(baseline)
    
    add_cov <- copy(baseline)

    process.names <- names(processes)

    which.cens <- process.names[sapply(processes, function(process) process[["type"]] == "censoring")]
    which.terminal <- process.names[sapply(processes, function(process) process[["type"]] == "terminal")]
    which.one.jump <- setdiff(process.names[sapply(processes, function(process) process[["type"]] == "one.jump")],
                              c(which.terminal, which.cens))

    process.order <- c(which.cens, which.terminal,
                       setdiff(process.names,
                               c(which.cens, which.terminal)))

    eta <- sapply(processes[process.order], function(process) process[["eta"]])
    nu <- sapply(processes[process.order], function(process) process[["nu"]])

    if (length(baseline.intervention)>0) {
        for (bname in names(baseline.intervention)) {
            add_cov[[bname]] <- function(N) rep(baseline.intervention[[bname]], N)
        }
    }

    if (length(alpha.intervention)>0) {
        for (alphaname in names(alpha.intervention)) {
            eta[names(processes[process.order]) == alphaname] <-
                alpha.intervention[[alphaname]]*eta[names(processes[process.order]) == alphaname]
        }
    }
    
    at_risk <- function(events) {

        out <- numeric(length(process.order))

        names(out) <- process.order

        ## censoring
        out[process.order %in% which.cens] <- cens

        ## terminal events
        out[process.order %in% which.terminal] <- 1

        ## one jump
        for (one.jump in which.one.jump) {
            idx <- which(process.order == one.jump)
            out[idx] <- as.numeric(events[idx] == 0)
        }

        ## recurrent
        out[setdiff(process.order,
                    c(which.cens,
                      which.terminal,
                      which.one.jump))] <- 1

        return(out)
    }

    if (!("A0" %in% baseline.vars)) {
        add_A0 <- 1
    } else {
        add_A0 <- 0
    }

    if (!("L0" %in% baseline.vars)) {
        add_L0 <- 1
    } else {
        add_L0 <- 0
    }
    
    beta <- matrix(
        0,
        nrow = length(process.order)+length(baseline.vars)+add_A0+add_L0,
        ncol = length(process.order)
    )

    rownames(beta) <- c(if (add_L0) "L0", if (add_A0) "A0", baseline.vars, process.order)
    colnames(beta) <- process.order
   
    for (effect in effects) {
        beta[effect[1], effect[2]] <- as.numeric(effect[3])
    }

    override_beta <- NULL

    if (browse) browser()

    data <- simEventData(N = n, beta = beta, eta = eta, nu = nu,
                         max_cens = Inf,
                         max_events = 50,
                         at_risk = at_risk,
                         lower = 1e-25, upper = 1e8,
                         term_deltas = 0:length(which.terminal),
                         #gen_L0 = add_cov[["L0"]],
                         #gen_A0 = {if ("A0" %in% names(add_cov)) function(N, L0) add_cov[["A0"]](N) else NULL},
                         add_cov = add_cov[!(names(add_cov) %in% c("A0", "L0"))],
                         override_beta = override_beta
                         )

    if (add_L0) {
        data[["L0"]] <- NULL
    }

    if (add_A0) {
        data[["A0"]] <- NULL
    }

    for (jj in 0:length(which.terminal)) {
        data[[paste0("N", jj)]] <- NULL
    }

    setnames(data, paste0("N", (length(which.terminal)+1):(length(process.order)-1)),
             setdiff(process.order, c(which.cens, which.terminal)))

    setnames(data, c("Delta", "Time", "ID"), c("delta", "time", "id"))

    return(data)
 

}

######################################################################
### sim.generic.R ends here

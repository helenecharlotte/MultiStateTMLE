### sim.from.data.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Aug 29 2026 (10:18) 
## Version: 
## Last-Updated: Aug 31 2026 (14:21) 
##           By: Helene
##     Update #: 42
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

sim.from.data <- function(n = 500,
                          sim.parameters,
                          cens = 1,
                          alpha.intervention = list(),
                          baseline.intervention = list(),
                          browse = FALSE,
                          verbose = FALSE) {

    processes <- names(sim.parameters)[
        names(sim.parameters) != "baseline.summary"
    ]

    baseline.vars <- unique(names(sim.parameters$baseline.summary))

    add_cov <- lapply(baseline.vars, function(x) {
        svar <- sim.parameters$baseline.summary[[x]]
        type <- svar$type
        if (type == "numeric") {
            if ((svar$type == "numeric" &&
                 svar$min %in% c(0, -1) &&
                 svar$max %in% c(1, 0) &&
                 length(unique(c(svar$min, svar$max))) == 2)) {
                p <- svar$mean
                return(function(N) stats::rbinom(N, 1, p))
            } else {
                mu <- svar$mean
                sd <- svar$sd
                return(function(N) stats::rnorm(N, mean = mu, sd = sd))
            }
        } else {
            probs <- unlist(svar$proportions)
            vals <- seq_along(names(probs))
            return(function(N) sample(vals, N, replace = TRUE, prob = probs))
        }
    })

    names(add_cov) <- baseline.vars

    process.names <- sim.parameters$model.structure$process.names
    process.deltas <- sim.parameters$model.structure$process.deltas
    
    which.cens <- process.names[sim.parameters$model.structure$cens.process.id]
    which.terminal <- setdiff(process.names[sim.parameters$model.structure$process.types == "terminal"],
                              which.cens)
    which.one.jump <- setdiff(process.names[sim.parameters$model.structure$process.types == "one.jump"],
                              c(which.terminal, which.cens))

    process.order <- process.names[order(process.deltas)]
        #c(which.cens, which.terminal, setdiff(process.names,
        #                                     c(which.cens, which.terminal)))

    eta <- as.numeric(sapply(sim.parameters[process.order], function(sim.param) sim.param[["weibull.parameters"]][substr(names(sim.param[["weibull.parameters"]]), 1, 3) == "eta"]))
    nu <- as.numeric(sapply(sim.parameters[process.order], function(sim.param) sim.param[["weibull.parameters"]][substr(names(sim.param[["weibull.parameters"]]), 1, 2) == "nu"]))

    if (length(baseline.intervention)>0) {
        for (bname in names(baseline.intervention)) {
            add_cov[[bname]] <- function(N) rep(baseline.intervention[[bname]], N)
        }
    }
    
    if (length(alpha.intervention)>0) {
        for (alphaname in names(alpha.intervention)) {
            eta[process.order == alphaname] <-
                alpha.intervention[[alphaname]]*eta[process.order == alphaname]
        }
    }
    
    if (verbose) {
        print(paste0("eta = ", eta))
        print(paste0("nu = ", nu))
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

    override_beta <- NULL

    for (proc in process.order) {

        cp <- sim.parameters[[proc]]$cox.parameters

        for (v in names(cp)){

            if (v %in% c(baseline.vars, process.order)) {
                beta[v, proc] <- cp[v]
            } else {
                bvar <- baseline.vars[sapply(baseline.vars, function(baseline.var)
                    length(grep(baseline.var, v, value = TRUE))>0)]
                bvalue <- gsub(bvar, "", v)
                
                svar <- sim.parameters$baseline.summary[[bvar]]
                probs <- unlist(svar$proportions)
                vals <- seq_along(names(probs))

                bval <- vals[names(probs) == bvalue][1]
                
                out_vec <- cp[v]
                #names(out_vec) <- paste0("N", (0:(length(process.order)-1))[process.order == proc])
                names(out_vec) <- paste0(
                    "N", match(proc, process.order)-1)
                #out_list <- list(out_vec)
                #names(out_list) <- paste0("(", bvar, "==", bvalue, ")")
                expr <- paste0("(", bvar, "==", bval, ")")
                if (is.null(override_beta[[expr]])) {
                    override_beta[[expr]] <- out_vec
                } else {
                    override_beta[[expr]] <- c(override_beta[[expr]], out_vec)
                }
                ##override_beta[[length(override_beta)+1]] <- out_vec
                ##names(override_beta)[length(override_beta)] <- expr
                #override_beta[[length(override_beta)+1]] <-
                #    out_list
            }

        }
    }

    if (verbose) print(beta)

    if (browse) browser()

    term.processes <- c(which.cens, which.terminal)
    term.deltas <- match(term.processes, process.order) - 1L

    non.term.processes <- setdiff(process.order, term.processes)
    non.term.deltas <- match(non.term.processes,
                             process.order) - 1L
        
    data <- simEventData(N = n, beta = beta, eta = eta, nu = nu,
                         max_cens = Inf,
                         max_events = 50,
                         at_risk = at_risk,
                         lower = 1e-45, upper = 1e10,
                         term_deltas = term.deltas,#0:length(which.terminal),
                         gen_L0 = add_cov[["L0"]],
                         gen_A0 = {if ("A0" %in% names(add_cov)) function(N, L0) add_cov[["A0"]](N) else NULL},
                         add_cov = add_cov[!(names(add_cov) %in% c("A0", "L0"))],
                         override_beta = override_beta
                         )

    if (add_L0) {
        data[["L0"]] <- NULL
    }

    if (add_A0) {
        data[["A0"]] <- NULL
    }

    for (jj in term.deltas) {
        data[[paste0("N", jj)]] <- NULL
    }

    setnames(data, paste0("N", non.term.deltas),
             non.term.processes)

    setnames(data, c("Delta", "Time", "ID"), c("delta", "time", "id"))

    return(data)

}

######################################################################
### sim.from.data.R ends here

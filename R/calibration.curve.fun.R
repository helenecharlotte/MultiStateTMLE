### calibration.curve.fun.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Sep  3 2026 (20:04) 
## Version: 
## Last-Updated: Sep  4 2026 (08:40) 
##           By: Helene
##     Update #: 94
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

calibration.curve.fun <- function(initial.fit = NULL,
                                  a = NULL,
                                  alpha.grid = seq(0, 5, length = 10),
                                  browse = FALSE,
                                  verbose = TRUE,
                                  output.eic = FALSE,
                                  tau = 1.2,
                                  tau.z = tau,
                                  use.cores = 50,
                                  target = "outcome", # c("outcome", "cr")
                                  z.name = "z",
                                  min.iter = 1, 
                                  target.by.state = FALSE,
                                  block.size.z = NULL,
                                  block.size.target = NULL,
                                  ...) {

    a.fixed <- a
    use.cores.fixed <- use.cores
    verbose.fixed <- verbose
    tau.fixed <- tau
    z.name.fixed <- z.name
    block.size.z.fixed <- block.size.z
    block.size.target.fixed <- block.size.target
    target.by.state.fixed <- target.by.state
    min.iter.fixed <- min.iter

    fit.local <- data.table::copy(initial.fit)
    
    n <- length(unique(initial.fit[["tmp.long"]][["id"]]))

    tmle.alpha.fun.fixed <- function(tau = tau.z,
                                     alpha = 1,
                                     parameter = z.name,
                                     output.eic = TRUE,
                                     a = a.fixed,
                                     block.size = NULL,
                                     output.weights = FALSE) {
        tmle.alpha.fun(initial.fit = data.table::copy(fit.local),
                       tau = tau,
                       alpha = alpha,
                       a = a,
                       target = parameter,
                       z.name = z.name.fixed,
                       use.cores = use.cores.fixed,
                       verbose = verbose.fixed,
                       output.eic = output.eic,
                       years.lost = block.size,
                       min.iter = min.iter.fixed,
                       target.by.state = target.by.state.fixed,
                       output.weights = output.weights,
                       ...)
    }

    if (browse) browser()

    alpha.out <- lapply(alpha.grid, function(alpha) {

        if (alpha > 0) {
            calibration.est <- tmle.alpha.fun.fixed(parameter = z.name,
                                                    tau = tau.z,
                                                    alpha = alpha,
                                                    block.size = block.size.z,
                                                    output.weights = TRUE,
                                                    output.eic = TRUE)

            est.deriv.auxiliary <- estimate.derivative(alpha,
                                                       parameter = z.name,
                                                       fun = tmle.alpha.fun.fixed,
                                                       block.size = block.size.z,
                                                       h = 0.3*n^{-1/6}*sqrt(mean(calibration.est$eic^2)))

            alpha.eic <- (1/est.deriv.auxiliary)*(0 - calibration.est$eic)

        } else {
            calibration.est <- list(estimate = c(tmle.est = 0, se = 0), eic = 0)
            alpha.eic <- 0
            est.deriv.auxiliary <- Inf
        }

        out.target <- lapply(target, function(target1) {
            
            target.est <- tmle.alpha.fun.fixed(tau = tau,
                                               parameter = target1,
                                               alpha = alpha,
                                               block.size = block.size.target,
                                               output.eic = TRUE)

            target.se.crude <- target.est[["estimate"]]["se"]

            if (alpha > 0) {

                est.deriv.target <- estimate.derivative(alpha,
                                                        parameter = target1,
                                                        tau = tau,
                                                        fun = tmle.alpha.fun.fixed,
                                                        block.size = block.size.target,
                                                        h = 0.3*n^{-1/6}*sqrt(mean(calibration.est$eic^2)))

                target.se <- sqrt(mean((target.est[["eic"]] + est.deriv.target*alpha.eic)^2)/n)

            } else {
                est.deriv.target <- Inf
                target.se <- sqrt(mean((target.est[["eic"]])^2)/n)
            }

            return(list(target.est = target.est, target.se.crude = target.se.crude,
                        target.se = target.se, est.deriv.target = est.deriv.target))
        })

        if (length(target)>1) {
            names(out.target) <- target
        } else {
            names(out.target) <- ""
        }

        out.estimate <- c(
            alpha = alpha, 
            exposure.est = calibration.est$estimate[[grep("tmle.est|one.step.est", names(calibration.est$estimate), value = TRUE)]],
            exposure.se = calibration.est$estimate[["se"]],
            target.est = sapply(out.target, function(out.target.1) as.numeric(out.target.1[["target.est"]][["estimate"]][grep("tmle.est|one.step.est", names(out.target.1[["target.est"]][["estimate"]]), value = TRUE)])), #["tmle.est"]
            target.se = sapply(out.target, function(out.target.1) out.target.1[["target.se"]])## ,
            ## target.se.crude = sapply(out.target, function(out.target.1) as.numeric(out.target.1[["target.se.crude"]]))
        )
    
        out.checks <- c(est.deriv.auxiliary = est.deriv.auxiliary,
                        est.deriv.target = sapply(out.target, function(out.target.1) out.target.1[["est.deriv.target"]]))

        if (output.eic) {
            target.eic <- lapply(out.target, function(out.target.1) out.target.1[["target.est"]][["eic"]])
            eic <- lapply(out.target, function(out.target.1) out.target.1[["target.est"]][["eic"]] + out.target.1[["est.deriv.target"]]*alpha.eic)
            if (length(target) == 1) {
                target.eic <- target.eic[[1]]
                eic <- eic[[1]]
            }
            return(list(estimate = out.estimate,
                        checks = out.checks,
                        alpha.eic = alpha.eic,
                        calibration.eic = calibration.est$eic,
                        target.eic = target.eic,
                        eic = eic))
        } else {
            return(list(estimate = out.estimate,
                        checks = out.checks))
        }
    })

    names(alpha.out) <- paste0("alpha = ", alpha.grid)

    estimate.table <- data.table::rbindlist(
                                      lapply(
                                          alpha.out,
                                          function(x) {
                                              as.list(x$estimate)
                                          }
                                      ),
                                      fill = TRUE
                                  )

    checks.list <- lapply(
        alpha.out,
        function(x) x$checks
    )

    names(checks.list) <- paste0(
        "alpha = ",
        alpha.grid
    )

    if (output.eic) {
        eic.list <- lapply(
            alpha.out,
            function(x) list(alpha.eic = x$alpha.eic, 
                             calibration.eic = x$calibration.eic,
                             target.eic = x$target.eic,
                             eic = x$eic
                             ))
        names(eic.list) <- paste0(
            "alpha = ",
            alpha.grid
        )
    } else {
        eic.list <- NULL
    }

    return(
        list(
            estimate = estimate.table,
            checks = checks.list,
            eics = eic.list
        ))
}

######################################################################
### calibration.curve.fun.R ends here

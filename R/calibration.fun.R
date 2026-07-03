### calibration.fun.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: May 13 2026 (19:31) 
## Version: 
## Last-Updated: Jul  3 2026 (09:56) 
##           By: Helene
##     Update #: 138
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

calibration.fun <- function(initial.fit = NULL,
                            a =  NULL, 
                            theta = 0.5, rho = NULL, delta = NULL,
                            rho.1a = NULL,
                            browse = FALSE,
                            verbose = TRUE,
                            output.eic = FALSE,
                            tau = 1.2,
                            tau.z = tau,
                            use.cores = 50,
                            target = "outcome",
                            z.name = "z",
                            ...) {

    a.fixed <- a
    use.cores.fixed <- use.cores
    verbose.fixed <- verbose
    tau.fixed <- tau
    z.name.fixed <- z.name

    fit.local <- data.table::copy(initial.fit)
    
    n <- length(unique(initial.fit[["tmp.long"]][["id"]]))

    tmle.alpha.fun.fixed <- function(tau = tau.z,
                                     alpha = 1,
                                     parameter = z.name,
                                     output.eic = TRUE,
                                     a = a.fixed) {
        tmle.alpha.fun(initial.fit = data.table::copy(fit.local),
                       tau = tau,
                       alpha = alpha,
                       a = a,
                       target = parameter,
                       z.name = z.name.fixed,
                       use.cores = use.cores.fixed,
                       verbose = verbose.fixed,
                       output.eic = output.eic,
                       ...)
    }

    if (length(rho)>0) {
        est.aux.1 <- tmle.alpha.fun.fixed(alpha = 1,
                                          a = a,
                                          output.eic = TRUE)
        theta <- rho*est.aux.1[["estimate"]][grep("tmle.est|one.step.est", names(est.aux.1[["estimate"]]), value = TRUE)]
        theta.eic <- rho*est.aux.1$eic
    } else if (length(rho.1a)>0) {
        est.aux.1a <- tmle.alpha.fun.fixed(alpha = 1,
                                           a = 1-a,
                                           output.eic = TRUE)
        theta <- rho.1a*est.aux.1a[["estimate"]][grep("tmle.est|one.step.est", names(est.aux.1[["estimate"]]), value = TRUE)]
        theta.eic <- rho.1a*est.aux.1a$eic
    } else if (length(delta)>0) {
        est.aux.1 <- tmle.alpha.fun.fixed(alpha = 1,
                                          a = a,
                                          output.eic = TRUE)
        theta <- delta + est.aux.1[["estimate"]][grep("tmle.est|one.step.est", names(est.aux.1[["estimate"]]), value = TRUE)]
        theta.eic <- est.aux.1$eic
    }  else {
        theta.eic <- 0
    }

    theta.se <- sqrt(mean(theta.eic^2)/n)

    if (theta<0) {
        theta <- 0
    } else if (theta>1) {
        theta <- 1
    }   

    if (browse) browser()

    if (theta>1 | theta<0) {
        stop(paste0("inadmissible target: theta = ", theta))
    }

    if (theta > 0) {
        est.alpha <- estimate.alpha.fun(fun = tmle.alpha.fun.fixed,
                                        c_n = n^{-1/2}/log(n),
                                        verbose = verbose,
                                        theta = theta)

        est.deriv.auxiliary <- estimate.derivative(est.alpha$alpha.hat,
                                                   parameter = z.name,
                                                   fun = tmle.alpha.fun.fixed, 
                                                   h = 0.3*n^{-1/6}*sqrt(mean(est.alpha$eic^2)))

        alpha.eic <- (1/est.deriv.auxiliary)*(theta.eic - est.alpha$eic)
    } else {
        est.alpha <- list(alpha.hat = 0,
                          eic = 0,
                          converged = TRUE,
                          dist = 0,
                          c.n = 0,
                          grid = data.frame(alpha = 0,
                                            psi   = 0))
        alpha.eic <- 0
        est.deriv.auxiliary <- Inf
    }

    alpha.se <- sqrt(mean(alpha.eic^2)/n)

    target.est <- tmle.alpha.fun.fixed(tau = tau,
                                       parameter = target,
                                       alpha = est.alpha$alpha.hat,
                                       output.eic = TRUE)

    target.se.crude <- target.est[["estimate"]]["se"]

    if (theta > 0) {

        est.deriv.target <- estimate.derivative(est.alpha$alpha.hat,
                                                parameter = target,
                                                tau = tau,
                                                fun = tmle.alpha.fun.fixed, 
                                                h = 0.3*n^{-1/6}*sqrt(mean(est.alpha$eic^2)))

        target.se <- sqrt(mean((target.est[["eic"]] + est.deriv.target*alpha.eic)^2)/n)

    } else {
        est.deriv.target <- Inf
        target.se <- sqrt(mean((target.est[["eic"]])^2)/n)
    }

    out.estimate <- c(alpha.est = est.alpha$alpha.hat,
                      alpha.se = alpha.se,
                      target.est = as.numeric(target.est[["estimate"]][grep("tmle.est|one.step.est", names(target.est[["estimate"]]), value = TRUE)]), #["tmle.est"]
                      target.se = target.se,
                      target.se.crude = as.numeric(target.se.crude),
                      theta = theta, theta.se = theta.se)
    
    out.checks <- c(est.deriv.auxiliary = est.deriv.auxiliary,
                    est.alpha.converged = est.alpha$converged,
                    est.alpha.dist = est.alpha$dist,
                    est.alpha.cn = est.alpha$c.n,
                    est.deriv.target = est.deriv.target)

    if (output.eic) {
        return(list(estimate = out.estimate,
                    checks = out.checks,
                    alpha.eic = alpha.eic,
                    target.eic = target.est[["eic"]],
                    eic = target.est[["eic"]] + est.deriv.target*alpha.eic))
    } else {
        return(list(estimate = out.estimate,
                    checks = out.checks))
    }

}

######################################################################
### calibration.fun.R ends here

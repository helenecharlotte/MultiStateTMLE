### make.calibrated.contrasts.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Sep  4 2026 (10:04) 
## Version: 
## Last-Updated: Sep  4 2026 (10:21) 
##           By: Helene
##     Update #: 18
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

make.calibrated.contrasts <- function(
                                      target = NULL,
                                      treatment.fit,
                                      placebo.fit,
                                      calibrated.fit,
                                      conf.level = 0.95
                                      ) {

    treatment.est <-
        treatment.fit$estimate[["tmle.est"]]

    placebo.est <-
        placebo.fit$estimate[["tmle.est"]]

    calibrated.est <-
        calibrated.fit$estimate[[ifelse(length(target)>0, paste0("target.est.", target), "target.est")]]

    treatment.eic <- treatment.fit$eic
    placebo.eic <- placebo.fit$eic
    calibrated.eic <- if (length(target)>0) calibrated.fit$eic[[target]] else calibrated.fit$eic

    n <- length(treatment.eic)

    if (length(placebo.eic) != n ||
        length(calibrated.eic) != n) {
        stop("The EIC vectors must have the same length.")
    }

    contrast <- function(
        estimate.1,
        estimate.0,
        eic.1,
        eic.0
    ) {

        estimate <- estimate.1 - estimate.0
        eic <- eic.1 - eic.0

        se <- stats::sd(eic) / sqrt(n)

        critical.value <-
            stats::qnorm(1 - (1 - conf.level) / 2)

        c(
            estimate = estimate,
            se = se,
            lower = estimate - critical.value * se,
            upper = estimate + critical.value * se
        )
    }

    total <- contrast(
        estimate.1 = treatment.est,
        estimate.0 = placebo.est,
        eic.1 = treatment.eic,
        eic.0 = placebo.eic
    )

    mediated <- contrast(
        estimate.1 = treatment.est,
        estimate.0 = calibrated.est,
        eic.1 = treatment.eic,
        eic.0 = calibrated.eic
    )

    calibrated <- contrast(
        estimate.1 = calibrated.est,
        estimate.0 = placebo.est,
        eic.1 = calibrated.eic,
        eic.0 = placebo.eic
    )

    out <- rbind(total, mediated, calibrated)
    out <- data.table(target = target, contrast = rownames(out), out)
    
    decomposition.error <-
        total[["estimate"]] -
        mediated[["estimate"]] -
        calibrated[["estimate"]]

    eic.decomposition.error <-
        max(abs(
            (treatment.eic - placebo.eic) -
            (treatment.eic - calibrated.eic) -
            (calibrated.eic - placebo.eic)
        ))

    attr(out, "decomposition.error") <-
        decomposition.error

    attr(out, "eic.decomposition.error") <-
        eic.decomposition.error

    return(out)
}

######################################################################
### make.calibrated.contrasts.R ends here

### simDropin.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Feb  5 2026 (11:26) 
## Version: 
## Last-Updated: Feb  5 2026 (13:28) 
##           By: Helene
##     Update #: 2
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

simDropin <- function(N,
                      beta_L_A = 0, beta_L_Z = 3, beta_L_D = 0.5, beta_L_C = 0,
                      beta_A_L = -0.5,  beta_A_Z = -0.5, beta_A_D = -1, beta_A_C = 0,
                      beta_Z_L = -2, beta_Z_A = 0, beta_Z_D = -3, beta_Z_C = 0,
                      beta_L0_L = 1, beta_L0_A = 1, beta_L0_Z = 0.1, beta_L0_D = 1, beta_L0_C = 0,
                      beta_A0_L = -2.5, beta_A0_A = 0, beta_A0_Z = 0, beta_A0_D = -0.5, beta_A0_C = 0,
                      eta = c(0.25, 0.05, 0.05, 0.25),
                      nu = c(1.1, 1.0, 0.2, 1.1),
                      adherence = FALSE,
                      followup = Inf, cens = 1,
                      max_iter = 500,
                      generate.A0 = function(N, L0) stats::rbinom(N, 1, 0.5)) {

    Time <- A0 <- N0 <- N1 <- ID <- NULL

    if (adherence) { #-- should add one more process
        if (length(eta) == 4) eta <- c(eta,eta[4])
        if (length(nu) == 4) nu <- c(nu,nu[4])
        at_risk <- function(events) {
            return(c(
                cens, # You might be at risk for censoring
                1, # If you have not died yet you are at risk for dying
                as.numeric(events[3] == 0), # You are at risk for drop-in initiation if you have not initiated yet
                as.numeric(events[4] == 0), # You are only at risk for a change in the covariate process if you have not experienced a change yet
                as.numeric(events[5] == 0))) # You are only at risk for a change in the treatment process if you have not experienced a change yet
        }
    } else {
        at_risk <- function(events) {
            return(c(
                cens, # You might be at risk for censoring
                1, # If you have not died yet you are at risk for dying
                as.numeric(events[3] == 0), # You are at risk for drop-in initiation if you have not initiated yet
                as.numeric(events[4] == 0))) # You are only at risk for a change in the covariate process if you have not experienced a change yet
        }
    }

    beta <- matrix(ncol = length(eta), nrow = 2+length(eta))

    if (length(eta) == 4) {
        # The effect of L0 on the processes C, D, Z, L
        beta[1,] <- c(beta_L0_C, beta_L0_D, beta_L0_Z, beta_L0_L)
        # A0 is 0
        beta[2,] <- c(beta_A0_C, beta_A0_D, beta_A0_Z, beta_A0_L)
        # The processes C and D are terminal
        beta[c(3,4),] <- 0
        # The effect of Z on the processes C, D, Z, L
        beta[5,] <- c(beta_Z_C, beta_Z_D, 0, beta_Z_L)
        # The effect of L on the processes C, D, Z, L
        beta[6,] <- c(beta_L_C, beta_L_D, beta_L_Z, 0)
    } else {
        # The effect of L0 on the processes C, D, Z, L, A
        beta[1,] <- c(beta_L0_C, beta_L0_D, beta_L0_Z, beta_L0_L, beta_L0_A)
        # A0 is 0
        beta[2,] <- c(beta_A0_C, beta_A0_D, beta_A0_Z, beta_A0_L, beta_A0_A)
        # The processes C and D are terminal
        beta[c(3,4),] <- 0
        # The effect of Z on the processes C, D, Z, L, A
        beta[5,] <- c(beta_Z_C, beta_Z_D, 0, beta_Z_L, beta_Z_A)
        # The effect of L on the processes C, D, Z, L, A
        beta[6,] <- c(beta_L_C, beta_L_D, beta_L_Z, 0, beta_L_A)
        # The effect of A on the processes C, D, Z, L, A
        beta[7,] <- c(beta_A_C, beta_A_D, beta_A_Z, beta_A_L, 0)
    }

    data <- simEventData(N,
                         beta = beta,
                         eta = eta,
                         nu = nu,
                         max_cens = followup,
                         at_risk = at_risk,
                         gen_A0 = generate.A0,
                         max_iter = max_iter,
                         lower = 1e-300,
                         upper = 1e20)

    setnames(data, c("N2", "N3"), c("Z", "L"))

    if (length(eta)>4) setnames(data, c("N4"), c("A"))

    data[, N0 := NULL]
    data[, N1 := NULL]    

    return(data)
}

######################################################################
### simDropin.R ends here

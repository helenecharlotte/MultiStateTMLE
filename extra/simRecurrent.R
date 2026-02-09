### simRecurrent.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Feb  6 2026 (19:07) 
## Version: 
## Last-Updated: Feb  7 2026 (20:53) 
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

simRecurrent <- function(N, beta_L_Z = 1, beta_L_D = 1, beta_Z_D = -1,
                         beta_W_Z = 1, beta_W_D = 1,
                         beta_Z_L = -0.5, beta_Z_W = -0.5, beta_L0_Z = 1,
                         beta_L0_L = 1, beta_L0_W = 0,
                         beta_W_L = 0.5, beta_L_W = 0.5,
                         beta_L0_D = 1, beta_L0_C = 0, beta_L_C = 0, beta_W_C = 0,
                         beta_Z_C = 0, eta = rep(0.1,5), nu = rep(1.1,5),
                         ##
                         beta_L0_D2 = 1, beta_Z_D2 = -0.5, beta_L_D2 = 0, beta_W_D2 = 1, 
                         ##
                         followup = Inf, cens = 1, lower = 10^(-15), upper = 200,
                         beta_L_Z_prime = 0, beta_L_D_prime = 0, beta_Z_D_prime = 0,
                         beta_Z_L_prime = 0, beta_L0_Z_prime = 0, beta_L0_L_prime = 0,
                         beta_L0_D_prime = 0, beta_L0_C_prime = 0, beta_L_C_prime = 0,
                         beta_Z_C_prime = 0, t_prime = NULL, at_risk_cov = NULL,
                         override_beta = list("N1 >= 1" = c("N1" = 2),
                                              "N1 >= 2" = c("N1" = 1))){

  Time <- A0 <- N0 <- N1 <- ID <- NULL

   
    at_risk <- function(events) {
        return(c(
            cens, # You might be at risk for censoring
            1, # If you have not died yet you are at risk for dying
            as.numeric(events[3] == 0), # You are at risk for an treatment if you have not had one yet
            as.numeric(events[4] == 0), # You are only at risk for a change in the covariate process if you have not experienced a change yet
            1)) # If you have not died yet you are at risk for dying
    }

    beta <- matrix(ncol = 4+1, nrow = 6+1)
    # The effect of L0 on the processes C, D, Z, L
    beta[1,] <- c(beta_L0_C, beta_L0_D, beta_L0_Z, beta_L0_L, beta_L0_D2)
    # A0 is 0
    beta[2,] <- 0
    # The processes C and D are terminal
    beta[c(3,4),] <- 0
    # The effect of Z on the processes C, D, Z, L
    beta[5,] <- c(beta_Z_C, beta_Z_D, 0, beta_Z_L, beta_Z_D2)
    # The effect of L on the processes C, D, Z, L
    beta[6,] <- c(beta_L_C, beta_L_D, beta_L_Z, 0, beta_L_D2)
    # The processes D2 is terminal
    beta[7,] <- 0

    if(!is.null(t_prime)){
        tv_eff <- matrix(ncol = 4, nrow = 6)
        # The change in effect of L0 on the processes C, D, Z, L
        tv_eff[1,] <- c(beta_L0_C_prime, beta_L0_D_prime, beta_L0_Z_prime, beta_L0_L_prime)
        # A0 is 0
        tv_eff[2,] <- 0
        # The processes C and D are terminal
        tv_eff[c(3,4),] <- 0
        # The change in effect of Z on the processes C, D, Z, L
        tv_eff[5,] <- c(beta_Z_C_prime, beta_Z_D_prime, 0, beta_Z_L_prime)
        # The change in effect of L on the processes C, D, Z, L
        tv_eff[6,] <- c(beta_L_C_prime, beta_L_D_prime, beta_L_Z_prime, 0)

        data <- simEventTV(N, beta = beta, eta = eta, nu = nu, max_cens = followup,
                           at_risk = at_risk, t_prime = t_prime, tv_eff = tv_eff,
                           at_risk_cov = at_risk_cov, upper = upper, lower = lower)
    }
    else{
        data <- simEventData(N, beta = beta, eta = eta, nu = nu, max_cens = followup,
                             at_risk = at_risk, at_risk_cov = at_risk_cov,
                             term_deltas = c(0,4),
                             override_beta = override_beta,
                             lower = lower, upper = upper)
    }

    data[, c("N0", "N1", "N4", "A0") := NULL]
    colnames(data)[c(5,6)] <- c("Z", "L")

  return(data)
}

######################################################################
### simRecurrent.R ends here

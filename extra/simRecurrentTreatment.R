#' Simulate Event History Data with Treatment and Time-Dependent Covariate
#'
#' Simulates event history data with four types of events representing censoring (0), death (1), treatment (2), and covariate change (3).
#' Death and censoring are terminal events; treatment and covariate events can occur only once.
#'
#' Event intensities are modeled using Weibull hazards with parameters \eqn{\nu} (scale) and \eqn{\eta} (shape),
#' and covariate effects controlled by specified \code{beta} parameters. For example, \code{beta_L_A} quantifies
#' the effect of covariate \code{L = 1} on the hazard of treatment.
#'
#' @param N Integer. Number of individuals to simulate.
#' @param beta_L_A Numeric. Effect of covariate \code{L = 1} on treatment hazard. Default 1.
#' @param beta_L_D Numeric. Effect of covariate \code{L = 1} on death hazard. Default 1.
#' @param beta_A_D Numeric. Effect of treatment \code{A = 1} on death hazard. Default -1.
#' @param beta_L0_A Numeric. Effect of baseline covariate \code{L0} on treatment hazard. Default 1.
#' @param beta_A_L Numeric. Effect of treatment \code{A = 1} on covariate hazard. Default -0.5.
#' @param beta_L0_L Numeric. Effect of baseline covariate \code{L0} on covariate hazard. Default 1.
#' @param beta_L0_D Numeric. Effect of baseline covariate \code{L0} on death hazard. Default 1.
#' @param beta_L0_C Numeric. Effect of baseline covariate \code{L0} on censoring hazard. Default 0.
#' @param beta_L_C Numeric. Effect of covariate \code{L = 1} on censoring hazard. Default 0.
#' @param beta_A_C Numeric. Effect of treatment \code{A = 1} on censoring hazard. Default 0.
#' @param eta Numeric vector of length 4. Shape parameters for Weibull intensities, parameterized as \eqn{\eta \nu t^{\nu - 1}}. Default is \code{rep(0.1, 4)}.
#' @param nu Numeric vector of length 4. Scale parameters for Weibull hazards. Default is \code{rep(1.1, 4)}.
#' @param followup Numeric. Maximum censoring time. Defaults to \code{Inf} (no censoring).
#' @param cens Integer (0 or 1). Indicates if censoring is possible. Default 1.
#' @param op Integer (0 or 1). Indicates if treatment (operation) is possible. Default 1.
#' @param lower Numeric. Lower bound for root finding (inverse cumulative hazard). Default \code{1e-15}.
#' @param upper Numeric. Upper bound for root finding (inverse cumulative hazard). Default 200.
#' @param beta_L_A_prime Numeric. Additional effect of covariate \code{L = 1} on treatment hazard. Default 0.
#' @param beta_L_D_prime Numeric. Additionalffect of covariate \code{L = 1} on death hazard. Default 0.
#' @param beta_A_D_prime Numeric. Effect of treatment \code{A = 1} on death hazard. Default 0.
#' @param beta_L0_A_prime Numeric. Effect of baseline covariate \code{L0} on treatment hazard. Default 0.
#' @param beta_A_L_prime Numeric. Effect of treatment \code{A = 1} on covariate hazard. Default 0.
#' @param beta_L0_L_prime Numeric. Effect of baseline covariate \code{L0} on covariate hazard. Default 0.
#' @param beta_L0_D_prime Numeric. Effect of baseline covariate \code{L0} on death hazard. Default 0.
#' @param beta_L0_C_prime Numeric. Effect of baseline covariate \code{L0} on censoring hazard. Default 0.
#' @param beta_L_C_prime Numeric. Effect of covariate \code{L = 1} on censoring hazard. Default 0.
#' @param beta_A_C_prime Numeric. Effect of treatment \code{A = 1} on censoring hazard. Default 0.
#' @param t_prime Numeric scalar or NULL. Time point where effects change (optional).
#' @param at_risk_cov Function. Function determining if an individual is at risk for each event type, given their covariates. Takes a numeric vector covariates and returns a binary vector. Default returns 1 for all events.
#'
#' @return A \code{data.frame} with columns:
#' \itemize{
#'   \item \code{ID} - Individual identifier.
#'   \item \code{Time} - Event time.
#'   \item \code{Delta} - Event type (0=censoring, 1=death, 2=treatment, 3=covariate change).
#'   \item \code{L0} - Baseline covariate.
#'   \item \code{L} - Time-dependent covariate.
#'   \item \code{A} - Treatment status.
#' }
#'
#' @export
#'
#' @examples
#' simTreatment(10)

simRecurrentTreatment <- function(N, beta_L_Z = 1, beta_L_D = 1, beta_Z_D = -1,
                                  beta_W_Z = 1, beta_W_D = 1,
                                  beta_Z_L = -0.5, beta_Z_W = -0.5, beta_L0_Z = 1,
                                  beta_L0_L = 1, beta_L0_W = 0,
                                  beta_W_L = 0.5, beta_L_W = 0.5,
                                  beta_L0_D = 1, beta_L0_C = 0, beta_L_C = 0, beta_W_C = 0,
                                  beta_Z_C = 0, eta = rep(0.1,6), nu = rep(1.1,6),
                                  ##
                                  beta_L0_D2 = 1, beta_Z_D2 = -0.5, beta_L_D2 = 0, beta_W_D2 = 1, 
                                  ##
                                  followup = Inf, cens = 1, lower = 10^(-15), upper = 200,
                                  beta_L_Z_prime = 0, beta_L_D_prime = 0, beta_Z_D_prime = 0,
                                  beta_Z_L_prime = 0, beta_L0_Z_prime = 0, beta_L0_L_prime = 0,
                                  beta_L0_D_prime = 0, beta_L0_C_prime = 0, beta_L_C_prime = 0,
                                  beta_Z_C_prime = 0, t_prime = NULL, at_risk_cov = NULL){

    Time <- A0 <- N0 <- N1 <- ID <- NULL
   
    at_risk <- function(events) {
        return(c(
            cens, # You might be at risk for censoring
            1, # If you have not died yet you are at risk for dying
            1, #as.numeric(events[3] == 0), # You are at risk for an treatment if you have not had one yet
            as.numeric(events[4] == 0), # You are only at risk for a change in the covariate process if you have not experienced a change yet
            as.numeric(events[5] == 0),  # You are only at risk for a change in the other covariate process if you have not experienced a change yet
            1)) # If you have not died yet you are at risk for dying
    }

    beta <- matrix(ncol = 4+2, nrow = 6+2)
    # The effect of L0 on the processes C, D, Z, L, W
    beta[1,] <- c(beta_L0_C, beta_L0_D, beta_L0_Z, beta_L0_L, beta_L0_W, beta_L0_D2)
    # A0 is 0
    beta[2,] <- 0
    # The processes C and D are terminal
    beta[c(3,4),] <- 0
    # The effect of Z on the processes C, D, Z, L, W
    beta[5,] <- c(beta_Z_C, beta_Z_D, 0, beta_Z_L, beta_Z_W, beta_Z_D2)
    # The effect of L on the processes C, D, Z, L, W
    beta[6,] <- c(beta_L_C, beta_L_D, beta_L_Z, 0, beta_L_W, beta_L_D2)
    # The effect of W on the processes C, D, Z, L, W
    beta[7,] <- c(beta_W_C, beta_W_D, beta_W_Z, beta_W_L, 0, beta_W_D2)
    # The process D2 is terminal
    beta[8,] <- 0

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
                           at_risk = at_risk, at_risk_cov = at_risk_cov, upper = upper,
                           term_deltas = c(0,1,5),
                           override_beta = list("N2 >= 1" = c("N2" = 0.75),
                                                "N2 >= 1" = c("N1" = -0.5),
                                                "N2 >= 1" = c("N3" = -1.25),
                                                "N2 >= 1" = c("N4" = -1.25),
                                                "N2 >= 1" = c("N5" = -0.5)),
                           lower = lower)
  }

    data[, c("N0", "N1", "A0", "N5") := NULL]

    setnames(data, c("N2", "N3", "N4"), c("Z", "L", "W")) 

  return(data)
}

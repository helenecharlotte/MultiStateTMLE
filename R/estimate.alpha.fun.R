### estimate.alpha.fun.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: May 13 2026 (19:38) 
## Version: 
## Last-Updated: Jul  2 2026 (17:58) 
##           By: Helene
##     Update #: 52
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

estimate.alpha.fun <- function(theta,
                               fun,
                               c_n,
                               alpha_init = 1,
                               expand_up = 2,#1.25,
                               expand_down = 0.25,#0.5,#0.8,
                               max_iter = 100,
                               alpha_min = 1e-3,
                               alpha_max = 100,
                               use.cores = 50,
                               verbose = FALSE,
                               trace_every = 1L) {

  if (missing(c_n) || is.null(c_n) || !is.finite(c_n)) {
    stop("Please supply a finite c_n explicitly.")
  }

  if (!is.finite(alpha_min) || !is.finite(alpha_max) || alpha_min <= 0 || alpha_max <= alpha_min) {
    stop("Require 0 < alpha_min < alpha_max.")
  }

    if (!is.finite(alpha_init)) stop("alpha_init must be finite.")

    # A little extra room for expansion + bracketing
    max_evals <- 2L * max_iter + 5L

    alpha_hist <- numeric(max_evals)
    psi_hist   <- numeric(max_evals)
    se_hist    <- numeric(max_evals)
    eic_hist   <- vector("list", max_evals)

    n_eval <- 0L
    best_idx <- NA_integer_
    best_dist <- Inf

    cache <- new.env(parent = emptyenv())

    record_eval <- function(rec) {
        n_eval <<- n_eval + 1L
        alpha_hist[n_eval] <<- rec$alpha
        psi_hist[n_eval]   <<- rec$psi
        se_hist[n_eval]    <<- rec$se
        eic_hist[[n_eval]]  <<- rec$eic

        d <- abs(rec$psi - theta)
        if (d < best_dist) {
            best_dist <<- d
            best_idx <<- n_eval
        }
        invisible(NULL)
    }

    eval_psi <- function(alpha) {
        alpha <- max(min(alpha, alpha_max), alpha_min)
        key <- sprintf("%.15g", alpha)

        if (exists(key, envir = cache, inherits = FALSE)) {
            out <- get(key, envir = cache, inherits = FALSE)
            out$cached <- TRUE
            return(out)
        }
    
        psi <- fun(
            alpha = alpha
        )

    out <- list(
      alpha = alpha,
      psi   = as.numeric(psi[["estimate"]][[grep("tmle.est|one.step.est", names(psi[["estimate"]]), value = TRUE)]]), #[["tmle.est"]]
      se    = as.numeric(psi[["estimate"]][["se"]]),
      eic   = psi[["eic"]],
      cached = FALSE
    )

    assign(key, out, envir = cache)
    out
  }

  progress <- function(stage, iter, alpha, psi, lo = NA_real_, hi = NA_real_, cached = FALSE) {
    if (!verbose) return(invisible(NULL))

    msg <- sprintf(
      "%s iter=%03d alpha=%.8g psi=%.8g |psi-theta|=%.3g",
      stage, iter, alpha, psi, abs(psi - theta)
    )

    if (is.finite(lo) && is.finite(hi)) {
      msg <- paste0(msg, sprintf(" bracket=[%.8g, %.8g] width=%.3g", lo, hi, hi - lo))
    }

    if (cached) {
      msg <- paste0(msg, " [cache hit]")
    }

    cat(msg, "\n")
    if (interactive()) flush.console()
    invisible(NULL)
  }

  # Initial evaluation
  rec0 <- eval_psi(alpha_init)
  record_eval(rec0)
  progress("init", 0L, rec0$alpha, rec0$psi, cached = rec0$cached)

  if (best_dist <= c_n) {
    ord <- order(alpha_hist[seq_len(n_eval)])
    return(list(
      alpha.hat = alpha_hist[best_idx],
      converged = TRUE,
      dist = best_dist,
      c.n = c_n,
      grid = data.frame(
        alpha = alpha_hist[seq_len(n_eval)][ord],
        psi   = psi_hist[seq_len(n_eval)][ord],
        se    = se_hist[seq_len(n_eval)][ord]
      ),
      eic = eic_hist[[best_idx]]
    ))
  }

  # This assumes psi(alpha) increases with alpha, like your current code.
  # If your function is decreasing instead, flip the inequalities below.
  if (rec0$psi < theta) {
    # Search upward to find an upper bracket
    lo <- rec0$alpha
    f_lo <- rec0$psi
    hi <- NA_real_
    f_hi <- NA_real_

    alpha <- rec0$alpha
    for (iter in seq_len(max_iter)) {
      alpha_new <- min(alpha_max, max(alpha_min, alpha * expand_up))
      if (alpha_new <= alpha * (1 + 1e-15)) break

      rec <- try(eval_psi(alpha_new))

      if (inherits(rec, "try-error")) {

          expand_up <- expand_up*0.9
          
      } else {
      
          if (!isTRUE(rec$cached)) record_eval(rec)
          progress("expand↑", iter, rec$alpha, rec$psi, cached = rec$cached)

          if (rec$psi >= theta) {
              hi <- rec$alpha
              f_hi <- rec$psi
              break
          }
          
          lo <- rec$alpha
          f_lo <- rec$psi
          alpha <- rec$alpha
      }
    }
  } else {
    # Search downward to find a lower bracket
    hi <- rec0$alpha
    f_hi <- rec0$psi
    lo <- NA_real_
    f_lo <- NA_real_

    alpha <- rec0$alpha
    for (iter in seq_len(max_iter)) {
      alpha_new <- max(alpha_min, min(alpha_max, alpha * expand_down))
      if (alpha_new >= alpha * (1 - 1e-15)) break

      rec <- eval_psi(alpha_new)
      if (!isTRUE(rec$cached)) record_eval(rec)
      progress("expand↓", iter, rec$alpha, rec$psi, cached = rec$cached)

      if (rec$psi <= theta) {
        lo <- rec$alpha
        f_lo <- rec$psi
        break
      }

      hi <- rec$alpha
      f_hi <- rec$psi
      alpha <- rec$alpha
    }
  }

  bracketed <- is.finite(lo) && is.finite(hi) && lo < hi &&
               is.finite(f_lo) && is.finite(f_hi) &&
               f_lo <= theta && theta <= f_hi

  if (bracketed) {
    for (iter in seq_len(max_iter)) {
      if (abs(f_hi - f_lo) < .Machine$double.eps) break

      # Secant step, clipped to the bracket; fallback to midpoint if needed
      alpha_new <- lo + (theta - f_lo) * (hi - lo) / (f_hi - f_lo)
      if (!is.finite(alpha_new) || alpha_new <= lo || alpha_new >= hi) {
        alpha_new <- 0.5 * (lo + hi)
      }

      if (alpha_new <= alpha_min) alpha_new <- alpha_min
      if (alpha_new >= alpha_max) alpha_new <- alpha_max

      rec <- eval_psi(alpha_new)
      if (!isTRUE(rec$cached)) record_eval(rec)
      progress("search", iter, rec$alpha, rec$psi, lo, hi, cached = rec$cached)

      if (abs(rec$psi - theta) <= c_n) break

      if (rec$psi < theta) {
        lo <- rec$alpha
        f_lo <- rec$psi
      } else {
        hi <- rec$alpha
        f_hi <- rec$psi
      }

      if (abs(hi - lo) <= .Machine$double.eps * max(1, abs(lo), abs(hi))) break
    }
  } else if (verbose) {
    cat("No bracket found; returning the best point seen.\n")
    if (interactive()) flush.console()
  }

  # Best point seen
  alpha_best <- alpha_hist[best_idx]
  converged <- best_dist <= c_n

  ord <- order(alpha_hist[seq_len(n_eval)])

    list(
        alpha.hat = alpha_best,
        converged = converged,
        dist = best_dist,
        c.n = c_n,
        grid = data.frame(
            alpha = alpha_hist[seq_len(n_eval)][ord],
            psi   = psi_hist[seq_len(n_eval)][ord],
            se    = se_hist[seq_len(n_eval)][ord]
        ),
        eic = eic_hist[[best_idx]]
    )
}


if (FALSE) {
estimate.alpha.fun <- function(theta,
                               fun,
                               parameter = "z",
                               tau = 3,
                               c_n = n^{-1/2}/log(n),
                               alpha_init = 1,
                               expand_up = 1.25,
                               expand_down = 0.8,
                               max_iter = 100,
                               alpha_min = 1e-3,
                               alpha_max = 100,
                               use.cores = 50,
                               ...) {



    # storage
    alpha_hist   <- numeric()
    psi_hist     <- numeric()
    se_hist      <- numeric()
    eic          <- numeric()

    # helper with caching
    eval_psi <- function(alpha) {
        idx <- which(abs(alpha_hist - alpha) < 1e-12)
        if (length(idx) > 0) return(psi_hist[idx[1]])
        psi <- fun(tau = tau, alpha = alpha, parameter = parameter,
                   use.cores = use.cores, output.eic = TRUE, ...)
        alpha_hist <<- c(alpha_hist, alpha)
        psi_hist   <<- c(psi_hist, psi[["estimate"]][grep("tmle.est|one.step.est", names(psi[["estimate"]]), value = TRUE)])
        se_hist    <<- c(se_hist, psi[["estimate"]]["se"])
        eic        <<- psi[["eic"]]
        psi[["estimate"]]["est"]
    }

    # step 0
    alpha_prev <- NA
    alpha_curr <- alpha_init
    psi_curr   <- eval_psi(alpha_curr)

    if (abs(psi_curr - theta) <= c_n) {
        return(list(
            alpha.hat = alpha_curr,
            converged = 1,
            dist = abs(psi_curr - theta),
            c.n = c_n,
            grid = data.frame(alpha = alpha_hist,
                              psi   = psi_hist,
                              se    = se_hist),
            eic = eic
        ))
    }

    # step 1
    if (psi_curr < theta - c_n) {
        alpha_next <- max(0, expand_up * alpha_curr)
    } else {
        alpha_next <- max(0, expand_down * alpha_curr)
    }

    # main loop
    for (m in 1:max_iter) {

        alpha_next <- max(min(alpha_next, alpha_max), alpha_min)

        psi_next   <- eval_psi(alpha_next)

        if (abs(psi_next - theta) <= c_n) break

        # update rule
        if (psi_next < theta - c_n) {
            alpha_new <- if (!is.na(alpha_prev) && alpha_next < alpha_curr)
                             (alpha_next + alpha_curr) / 2
                         else
                             expand_up * alpha_next
        } else {
            alpha_new <- if (!is.na(alpha_prev) && alpha_next > alpha_curr)
                             (alpha_next + alpha_curr) / 2
                         else
                             expand_down * alpha_next
        }

        if (alpha_new < 0) alpha_new <- 0

        alpha_prev <- alpha_curr
        alpha_curr <- alpha_next
        alpha_next <- alpha_new
    }

    # ordered grid
    ord <- order(alpha_hist)

    dist <- abs(psi_hist - theta)
    best_idx <- which.min(dist)
    alpha_best <- alpha_hist[best_idx]

    if (min(dist) <= c_n) {
        converged <- TRUE
    } else {
        converged <- FALSE
        eval_psi(alpha_best)
    }
   
    list(
        alpha.hat = alpha_best,
        converged = converged,
        dist = min(dist),
        c.n = c_n,
        grid = data.frame(alpha = alpha_hist[ord],
                          psi   = psi_hist[ord])
    )
}

}

######################################################################
### estimate.alpha.fun.R ends here

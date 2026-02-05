### compute.Q.clever.per.id.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Feb  4 2026 (08:39) 
## Version: 
## Last-Updated: Feb  5 2026 (09:39) 
##           By: Helene
##     Update #: 25
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

#' compute_Q_and_clever_per_id: unified recursion for finite product state spaces
#'
#' @param dt_id data.table with time-ordered rows for a single id. Column names must include hazard columns with prefix P.prefix
#'        of the form: P.<name>.<s> where <name> is a process/outcome name and <s> is the integer state index (1..S).
#' @param states data.table enumerating the product state space. Must have S rows and columns named "<proc>.1", "<proc2>.1", ...
#'        and a column `state` giving the state index (1..S) matching hazard column second index.
#' @param P.prefix prefix used for hazard columns in dt_id, default "P." (expects columns like "P.Z.1", "P.outcome1.3", ...)
#' @param state.idx.col column name in dt_id that holds the state index (value in 1..S)
#' @param parameter string: the name of the target process/outcome. If "target" (default) we use the first discovered terminal/outcome name.
#' @param process.deltas optional numeric vector (length = number of stateful processes) - NOT used to determine caps; only validated if provided.
#' @param compute.clever logical, whether to compute clever covariates (default TRUE)
#' @return data.table: original dt_id augmented with Q and clever.Q.<name>0 / clever.Q.<name>1 for each discovered name.

compute.Q.clever.per.id <- function(dt_id,
                                    states,
                                    P.prefix = "P.",
                                    state.idx.col = "state",
                                    parameter = "target",
                                    compute.clever = TRUE) {

    # --- checks / sizes
    S  <- nrow(states)
    if (is.null(S) || S < 1) stop("'states' must be a data.table with S >= 1 rows")
    Tn <- nrow(dt_id)
    if (Tn < 1) stop("dt_id must have at least one row")

    # --- stateful process names (columns in states except 'state' and optional 'at.risk')
    state_cols_all <- setdiff(names(states), c("at.risk", state.idx.col))
    state_processes <- state_cols_all
    M_proc <- length(state_processes)

    # --- parse hazard columns in dt_id with pattern P.<name>.<s>
    pcols <- grep(paste0("^", P.prefix), names(dt_id), value = TRUE)
    if (length(pcols) == 0) {
        warning("No hazard columns found with prefix '", P.prefix, "'. Returning Q=0")
        out <- copy(dt_id)
        out[, Q := 0]
        return(out[])
    }
    re <- paste0("^", P.prefix, "([^\\.]+)\\.(\\d+)$")
    regs <- regmatches(pcols, regexec(re, pcols))
    found_names <- unique(vapply(regs, function(x) if(length(x) >= 3) x[2] else NA_character_, ""))
    found_names <- found_names[!is.na(found_names)]
    if (length(found_names) == 0) stop("No hazard columns matched the expected pattern '", P.prefix, "<name>.<s>'")

    # build columns_by_name and hazard matrices (Tn x S) for each name
    columns_by_name <- setNames(lapply(found_names, function(x) rep(NA_character_, S)), found_names)
    names(columns_by_name) <- found_names
    for (i in seq_along(pcols)) {
        r <- regs[[i]]
        if (length(r) < 3) next
        nm <- r[2]; sidx <- as.integer(r[3])
        if (is.na(sidx) || sidx < 1 || sidx > S) next
        columns_by_name[[nm]][sidx] <- pcols[i]
    }
    hazard_mats <- list()
    for (nm in found_names) {
        mat <- matrix(0, nrow = Tn, ncol = S)
        for (s in seq_len(S)) {
            colname <- columns_by_name[[nm]][s]
            if (!is.na(colname) && colname %in% names(dt_id)) mat[, s] <- as.numeric(dt_id[[colname]])
        }
        hazard_mats[[nm]] <- mat
    }

    # partition names into stateful (present in states) and terminal (others)
    stateful_names <- intersect(found_names, state_processes)
    terminal_names  <- setdiff(found_names, stateful_names)
    stateful_names <- sort(stateful_names)
    terminal_names <- sort(terminal_names)

    # ---- proc caps from states (maximum observed count for each process)
    proc_caps <- integer(M_proc)
    for (j in seq_len(M_proc)) {
        vals <- as.integer(states[[ state_processes[j] ]])
        if (all(is.na(vals))) stop("state column ", state_processes[j], " appears entirely NA")
        proc_caps[j] <- max(vals, na.rm = TRUE)
    }

    # ---- build gamma: stay and jump successor indices (capped)
    key_vec <- do.call(paste, c(lapply(states[, ..state_processes], as.character), sep = ","))
    state_lookup <- setNames(states[[state.idx.col]], key_vec)
    gamma_stay_idx <- states[[state.idx.col]]
    gamma_jump_idx <- vector("list", M_proc)
    for (j in seq_len(M_proc)) {
        gamma_jump_idx[[j]] <- integer(S)
        for (s in seq_len(S)) {
            parts <- as.integer(strsplit(key_vec[s], ",")[[1]])
            parts[j] <- pmin(parts[j] + 1L, proc_caps[j])
            nk <- paste(parts, collapse = ",")
            if (is.null(state_lookup[[nk]])) stop("gamma mapping error: successor key not found: ", nk)
            gamma_jump_idx[[j]][s] <- as.integer(state_lookup[[nk]])
        }
    }

    # ---- classify target and special cases
    target_name <- parameter
    if (! (target_name %in% c(stateful_names, terminal_names)) ) {
        stop("parameter must be one of discovered names (stateful or terminal). Found: ", paste(c(stateful_names, terminal_names), collapse = ", "))
    }
    target_in_states <- target_name %in% stateful_names
    has_at_risk <- "at.risk" %in% names(states)

    # target types:
    # - "terminal": not among states
    # - "recurrent": among states and no at.risk column
    # - "state-with-atrisk": among states and states$at.risk provided (applies only to target)
    if (!target_in_states) {
        target_type <- "terminal"
    } else if (has_at_risk) {
        target_type <- "state-with-atrisk"
    } else {
        target_type <- "recurrent"
    }

  target_idx_in_statecols <- if (target_in_states) which(state_processes == target_name) else NA_integer_
  target_one_jump <- FALSE
  if (target_in_states) target_one_jump <- (proc_caps[target_idx_in_statecols] == 1)

    target_can_jump_from_state <- function(s) {
        if (target_type == "state-with-atrisk") {
            isTRUE(states$at.risk[s])
        } else {
            TRUE
        }
    }

    # ---- backward recursion: compute Q_t (vector of length S) for each time t
    Q_by_time <- vector("list", Tn)
    # Terminal condition at T: use immediate reward (hazard) for terminal targets, or last hazard for stateful targets
    if (target_type == "terminal") {
        Q_T <- hazard_mats[[target_name]][Tn, ]
        # if target also happens to be recorded in states (indicator), force post-jump states to 1
        if (target_name %in% names(states)) {
            Q_T[ as.integer(states[[target_name]]) >= 1 ] <- 1
        }
    } else {
        # stateful target: initial immediate reward uses target process hazard at time T
        if (target_name %in% names(hazard_mats)) {
            Q_T <- hazard_mats[[target_name]][Tn, ]
        } else {
            Q_T <- rep(0, S)
        }
        # if one-jump target, any state with target already >=1 should have Q_T = 1
        if (target_one_jump) Q_T[ as.integer(states[[target_name]]) >= 1 ] <- 1
    }
    Q_by_time[[Tn]] <- Q_T

    # recursion backwards
    for (tt in (Tn - 1):1) {
        Qn <- Q_by_time[[tt + 1]]

        # sum hazards of terminal outcomes and stateful processes
        Pout_sum <- if (length(terminal_names) > 0) Reduce(`+`, lapply(terminal_names, function(nm) hazard_mats[[nm]][tt, ]), init = rep(0, S)) else rep(0, S)
        Pproc_sum <- if (length(stateful_names) > 0) Reduce(`+`, lapply(stateful_names, function(nm) hazard_mats[[nm]][tt, ]), init = rep(0, S)) else rep(0, S)
        Pstay <- 1 - Pout_sum - Pproc_sum
        Pstay[Pstay < 0 & Pstay > -1e-12] <- 0

        Qt <- Pstay * Qn[ gamma_stay_idx ]

        # add contributions from each stateful process
        for (j in seq_along(state_processes)) {
            procj <- state_processes[j]
            Pj <- if (procj %in% names(hazard_mats)) hazard_mats[[procj]][tt, ] else rep(0, S)

            if (procj == target_name && target_type == "state-with-atrisk") {
                # only target has special at.risk logic
                can_jump <- vapply(seq_len(S), function(ss) target_can_jump_from_state(ss), logical(1))
                # where jump allowed: usual jump to gamma_jump
                Qt[can_jump] <- Qt[can_jump] + Pj[can_jump] * Qn[ gamma_jump_idx[[j]][can_jump] ]
                # where jump not allowed: redirect mass to 'stay' (i.e. successor is stay)
                Qt[!can_jump] <- Qt[!can_jump] + Pj[!can_jump] * Qn[ gamma_stay_idx[!can_jump] ]
            } else {
                # usual processing
                Qt <- Qt + Pj * Qn[ gamma_jump_idx[[j]] ]
            }
        }

        # add immediate reward if target is terminal
        if (target_type == "terminal") {
            Qt <- Qt + hazard_mats[[target_name]][tt, ]
        } else {
            # if target is one-jump and state already has target >=1, enforce Q=1
            if (target_one_jump) {
                already_idx <- which(as.integer(states[[target_name]]) >= 1)
                if (length(already_idx) > 0) Qt[already_idx] <- 1
            }
        }

        Q_by_time[[tt]] <- Qt
    }

    # --- extract Q per-row and clever covariates
    out <- copy(dt_id)
    # Q per-row:
    Q_row <- numeric(Tn)
    for (tt in seq_len(Tn)) {
        s <- dt_id[[state.idx.col]][tt]
        Q_row[tt] <- Q_by_time[[tt]][s]
    }
    out[, Q := Q_row]

    if (!compute.clever) return(out[])

    # prepare clever matrices: columns for every discovered name (stateful + terminal)
    all_names <- sort(unique(c(stateful_names, terminal_names)))
    clever0 <- matrix(NA_real_, nrow = Tn, ncol = length(all_names), dimnames = list(NULL, all_names))
    clever1 <- matrix(NA_real_, nrow = Tn, ncol = length(all_names), dimnames = list(NULL, all_names))

    for (tt in seq_len(Tn)) {
        s <- dt_id[[state.idx.col]][tt]
        Qn <- if (tt < Tn) Q_by_time[[tt + 1]] else Q_by_time[[tt]]

        # terminal names:
        for (nm in terminal_names) {
            clever0[tt, nm] <- Qn[s]
            if (nm == target_name && target_type == "terminal") {
                # target terminal -> clever1 = 1
                clever1[tt, nm] <- 1
            } else {
                # for non-target terminal: clever1 should reflect whether target count already achieved at state s
                if (target_in_states) {
                    clever1[tt, nm] <- as.numeric(as.integer(states[[target_name]][s]) >= 1)
                } else {
                    # fallback if target not among states
                    clever1[tt, nm] <- 0 
                }
            }
        }

        # stateful processes:
        for (j in seq_along(state_processes)) {
            procj <- state_processes[j]
            s_no <- gamma_stay_idx[s]
            s_yes <- gamma_jump_idx[[j]][s]

            clever0[tt, procj] <- Qn[s_no]

            if (procj == target_name) {
                # target stateful process special cases:
                if (target_one_jump) {
                    # one-jump target -> clever1 = 1 for all states (by your rule)
                    clever1[tt, procj] <- 1
                } else if (target_type == "state-with-atrisk" && !target_can_jump_from_state(s)) {
                    # not at risk at this state -> behaves like terminal at this state
                    clever1[tt, procj] <- 1
                } else {
                    # usual recurrent case
                    clever1[tt, procj] <- Qn[s_yes]
                }
            } else {
                # non-target stateful -> usual post-jump expected future
                clever1[tt, procj] <- Qn[s_yes]
            }
        }
    }

    # attach clever columns to out
    for (nm in all_names) {
        out[[ paste0("clever.Q.", nm, "0") ]] <- clever0[, nm]
        out[[ paste0("clever.Q.", nm, "1") ]] <- clever1[, nm]
    }

    return(out[])
}



######################################################################
### compute.Q.clever.per.id.R ends here

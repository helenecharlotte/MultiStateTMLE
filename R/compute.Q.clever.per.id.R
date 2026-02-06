### compute.Q.clever.per.id.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Feb  4 2026 (08:39) 
## Version: 
## Last-Updated: Feb  6 2026 (21:19) 
##           By: Helene
##     Update #: 67
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
                                    process.types,
                                    P.prefix = "P.",
                                    state.idx.col = "state",
                                    parameter = "target",
                                    compute.clever = TRUE) {
    ## Safe minimal dependencies on data.table already present in your project
    # normalize/handle empty states table
    if (is.null(states) || nrow(states) == 0) {
        states <- data.table::data.table(state = 1L)
        empty_states_flag <- TRUE
    } else {
        empty_states_flag <- FALSE
    }

    S  <- nrow(states)
    if (is.null(S) || S < 1) stop("'states' must be a data.table with S >= 1 rows")
    Tn <- nrow(dt_id)
    if (Tn < 1) stop("dt_id must have at least one row")

    # --- stateful process names that are actually columns in 'states' (exclude at.risk and state index)
    state_cols_all <- setdiff(names(states), c("at.risk", state.idx.col))
    state_processes <- state_cols_all    # initial list (may be augmented below)
    # If process.types says the target is "recurrent" but target not in states,
    # force the target into state_processes so the stateful loop executes.
    if (!missing(process.types) && !is.null(process.types) &&
        !is.null(parameter) && parameter %in% names(process.types) &&
        identical(process.types[[parameter]], "recurrent") &&
        !(parameter %in% state_processes)) {
        state_processes <- c(parameter, state_processes)
    }

    # keep track of which of the state_processes are actually represented as columns in 'states'
    stateful_in_states <- intersect(state_processes, state_cols_all)

    # --- parse hazard columns that look like P.<name>.<s> (we purposely restrict to columns ending with .<digits>)
    pcols <- grep(paste0("^", P.prefix, "[^\\.]+\\.[0-9]+$"), names(dt_id), value = TRUE)
    if (length(pcols) == 0) {
        warning("No hazard columns found with prefix '", P.prefix, "'. Returning dt_id with Q=0")
        out <- data.table::copy(dt_id)
        out[, Q := 0]
        return(out[])
    }
    re <- paste0("^", P.prefix, "([^\\.]+)\\.(\\d+)$")
    regs <- regmatches(pcols, regexec(re, pcols))
    found_names <- unique(vapply(regs, function(x) if (length(x) >= 3) x[2] else NA_character_, FUN.VALUE = ""))
    found_names <- found_names[!is.na(found_names)]
    if (length(found_names) == 0) stop("No hazard columns matched the expected pattern '", P.prefix, "<name>.<s>'")

    # construct columns_by_name and hazard matrices (Tn x S)
    columns_by_name <- setNames(lapply(found_names, function(x) rep(NA_character_, S)), found_names)
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

    # partition discovered names
    stateful_names <- intersect(found_names, state_processes)
    terminal_names  <- setdiff(found_names, stateful_names)

    # if process.types provided, restrict terminal_names to those declared terminal
    if (!missing(process.types) && !is.null(process.types)) {
        declared_terminal <- names(process.types)[sapply(process.types, function(x) identical(x, "terminal"))]
        terminal_names <- intersect(terminal_names, declared_terminal)
    }

    stateful_names <- sort(stateful_names)
    terminal_names  <- sort(terminal_names)

    # ---- proc_caps and gamma mapping for the columns that *are actually in states*
    if (length(stateful_in_states) > 0) {
        proc_caps <- integer(length(stateful_in_states))
        for (j in seq_along(stateful_in_states)) {
            vals <- as.integer(states[[ stateful_in_states[j] ]])
            if (all(is.na(vals))) stop("state column ", stateful_in_states[j], " appears entirely NA")
            proc_caps[j] <- max(vals, na.rm = TRUE)
        }
        # build the key vector for the *actual* state-columns
        key_vec_real <- do.call(paste, c(lapply(states[, ..stateful_in_states], as.character), sep = ","))
        state_lookup_real <- setNames(states[[ state.idx.col ]], key_vec_real)
    } else {
        proc_caps <- integer(0)
        key_vec_real <- as.character(states[[ state.idx.col ]])
        state_lookup_real <- setNames(states[[ state.idx.col ]], key_vec_real)
    }

    # global stay mapping (maps row -> state index)
    gamma_stay_idx <- states[[ state.idx.col ]]

    # build gamma_jump_idx aligned to 'state_processes' order.
    # For processes actually represented in states we compute the real successor;
    # for processes not present in states we use identity mapping (so Q_next[ gamma_jump_idx[[j]] ] == Q_next).
    M_proc_all <- length(state_processes)
    gamma_jump_idx <- vector("list", M_proc_all)
    for (j_all in seq_len(M_proc_all)) {
        procj_name <- state_processes[j_all]
        if (procj_name %in% stateful_in_states) {
            # index within the actual-state-columns vector
            j_in <- which(stateful_in_states == procj_name)
            gamma_jump_idx[[j_all]] <- integer(S)
            for (s in seq_len(S)) {
                parts <- as.integer(strsplit(key_vec_real[s], ",")[[1]])
                parts[j_in] <- pmin(parts[j_in] + 1L, proc_caps[j_in])
                nk <- paste(parts, collapse = ",")
                if (is.null(state_lookup_real[[nk]])) stop("gamma mapping error: successor key not found: ", nk)
                gamma_jump_idx[[j_all]][s] <- as.integer(state_lookup_real[[nk]])
            }
        } else {
            # process not represented in 'states' -> identity successor mapping
            gamma_jump_idx[[j_all]] <- seq_len(S)
        }
    }

    # ---- target classification (we expect process.types to tell us terminal/recurrent/one.jump/state-with-atrisk)
    target_name <- parameter
    if (!(target_name %in% c(stateful_names, terminal_names) || target_name %in% state_processes)) {
        stop("parameter must be one of discovered names (stateful or terminal) or a recurrent target declared in process.types. Found: ",
             paste(c(stateful_names, terminal_names, state_processes), collapse = ", "))
    }
    target_in_states <- target_name %in% stateful_in_states
    target_type <- NULL
    if (!missing(process.types) && !is.null(process.types) && target_name %in% names(process.types)) {
        target_type <- process.types[[target_name]]
    } else {
        # fallback heuristic: if target is in terminal_names -> terminal, else if in stateful_in_states -> recurrent, else one.jump
        target_type <- if (target_name %in% terminal_names) "terminal" else if (target_name %in% stateful_in_states) "recurrent" else "one.jump"
    }
    target_one_jump <- identical(target_type, "one.jump")

    # helper: at.risk for state-with-atrisk (only meaningful if states has at.risk column)
    target_can_jump_from_state <- function(s) {
        if (identical(target_type, "state-with-atrisk") && ("at.risk" %in% names(states))) {
            isTRUE(states$at.risk[s])
        } else {
            TRUE
        }
    }

    # ---- backward recursion: compute Q_t (vector length S) for each time t
    Q_by_time <- vector("list", Tn)

    # Terminal condition (time T)
    if (target_name %in% names(hazard_mats)) {
        Q_T <- hazard_mats[[ target_name ]][Tn, ]
    } else {
        Q_T <- rep(0, S)
    }
    # if one-jump and states contains an indicator column for the target, those states have Q_T = 1
    if (target_one_jump && target_name %in% names(states)) {
        Q_T[ as.integer(states[[ target_name ]]) >= 1 ] <- 1
    }
    Q_by_time[[Tn]] <- Q_T

    # recursion backwards
    for (tt in (Tn - 1):1) {
        Qn <- Q_by_time[[tt + 1]]

        Pout_sum  <- if (length(terminal_names) > 0) Reduce(`+`, lapply(terminal_names, function(nm) hazard_mats[[nm]][tt, ]), init = rep(0, S)) else rep(0, S)
        Pproc_sum <- if (length(stateful_names)  > 0) Reduce(`+`, lapply(stateful_names, function(nm) hazard_mats[[nm]][tt, ]), init = rep(0, S)) else rep(0, S)
        Pstay <- 1 - Pout_sum - Pproc_sum
        Pstay[Pstay < 0 & Pstay > -1e-12] <- 0

        Qt <- Pstay * Qn[ gamma_stay_idx ]

        # contributions from each process in the stateful loop (note: state_processes may include target even if not in states)
        for (j in seq_along(state_processes)) {
            procj <- state_processes[j]
            Pj <- if (procj %in% names(hazard_mats)) hazard_mats[[procj]][tt, ] else rep(0, S)

            # special handling for "state-with-atrisk" target: redirect mass where jump not allowed
            if (procj == target_name && identical(target_type, "state-with-atrisk")) {
                can_jump <- vapply(seq_len(S), function(ss) target_can_jump_from_state(ss), logical(1))
                if (any(can_jump)) Qt[can_jump] <- Qt[can_jump] + Pj[can_jump] * Qn[ gamma_jump_idx[[j]][can_jump] ]
                if (any(!can_jump)) Qt[!can_jump] <- Qt[!can_jump] + Pj[!can_jump] * Qn[ gamma_stay_idx[!can_jump] ]
            } else {
                # general case: use jump successor mapping (identity if process not represented in states)
                Qt <- Qt + Pj * Qn[ gamma_jump_idx[[j]] ]
            }
        }

        # add immediate reward (1 * hazard) for the target if hazard available (keeps same convention as original code)
        if (identical(target_type, "terminal") | identical(target_type, "recurrent")) {
            if (target_name %in% names(hazard_mats)) {
                Qt <- Qt + hazard_mats[[ target_name ]][tt, ]
            }
        }

        # enforce one-jump terminal states -> Q = 1 where appropriate
        if (target_one_jump && target_name %in% names(states)) {
            already_idx <- which(as.integer(states[[ target_name ]]) >= 1)
            if (length(already_idx) > 0) Qt[already_idx] <- 1
        }

        Q_by_time[[tt]] <- Qt
    }

    # --- extract Q per-row
    out <- data.table::copy(dt_id)
    Q_row <- numeric(Tn)
    for (tt in seq_len(Tn)) {
        s <- as.integer(dt_id[[ state.idx.col ]][tt])
        Q_row[tt] <- Q_by_time[[tt]][ s ]
    }
    out[, Q := Q_row]

    if (!compute.clever) return(out[])

    # ---------- clever covariates
    # produce clever columns for all processes we care about: include any state_processes and terminal_names
    all_names <- sort(unique(c(state_processes, terminal_names)))
    clever0 <- matrix(NA_real_, nrow = Tn, ncol = length(all_names), dimnames = list(NULL, all_names))
    clever1 <- matrix(NA_real_, nrow = Tn, ncol = length(all_names), dimnames = list(NULL, all_names))

    for (tt in seq_len(Tn)) {
        s <- as.integer(dt_id[[ state.idx.col ]][tt])
        Qn <- if (tt < Tn) Q_by_time[[tt + 1]] else Q_by_time[[tt]]

        # terminal processes
        for (nm in terminal_names) {
            clever0[tt, nm] <- Qn[s]
            if (nm == target_name && identical(target_type, "terminal")) {
                clever1[tt, nm] <- 1
            } else {
                # non-target terminal: indicate whether target count already achieved (if target encoded in states)
                if (target_name %in% names(states)) {
                    clever1[tt, nm] <- as.numeric(as.integer(states[[ target_name ]][ s ]) >= 1)
                } else {
                    clever1[tt, nm] <- 0
                }
            }
        }

        # stateful processes
        for (j in seq_along(state_processes)) {
            procj <- state_processes[j]
            s_no <- gamma_stay_idx[s]
            s_yes <- gamma_jump_idx[[j]][s]   # identity if proc not in states

            clever0[tt, procj] <- Qn[s_no]

            if (procj == target_name) {
                # target-specific rules
                if (identical(target_type, "one.jump")) {
                    clever1[tt, procj] <- 1
                } else if (identical(target_type, "state-with-atrisk") && !target_can_jump_from_state(s)) {
                    clever1[tt, procj] <- 1
                } else if (identical(target_type, "recurrent")) {
                    # recurrent target: clever1 = 1 + expected future after jump
                    if (procj %in% stateful_in_states) {
                        clever1[tt, procj] <- 1 + Qn[s_yes]
                    } else {
                        # target not represented in states -> successor is same state (identity)
                        clever1[tt, procj] <- 1 + Qn[s]
                    }
                } else {
                    # default: post-jump expected future
                    clever1[tt, procj] <- Qn[s_yes]
                }
            } else {
                # non-target stateful process: usual post-jump expected future
                clever1[tt, procj] <- Qn[s_yes]
            }
        }
    }

    # attach clever cols
    for (nm in all_names) {
        out[[ paste0("clever.Q.", nm, "0") ]] <- clever0[, nm]
        out[[ paste0("clever.Q.", nm, "1") ]] <- clever1[, nm]
    }

    return(out[])
}





######################################################################
### compute.Q.clever.per.id.R ends here

#' Compute Q and clever covariates per id using a discrete-state backward recursion
#'
#' The function expects a time-ordered `dt_id` (rows for a single id) that
#' contains hazard columns named like `P.<name>.<s>` for each discovered
#' process `<name>` and state index `s` in `1..S`. It uses the product-space
#' `states` table to map state indices and to build successor mappings.
#'
#' @param dt_id data.table with time-ordered rows for a single id. Must contain columns for hazards of the form P.<name>.<s>.
#' @param states data.table enumerating the product state space. Must have S rows and columns for each stateful process (named by process) plus a column `state` giving indices 1..S.
#' @param process.types named character vector or list describing each discovered process type. Allowed values include: "terminal", "recurrent", "one.jump", "state-with-atrisk". If omitted, a small heuristic is used.
#' @param P.prefix prefix used for hazard columns in dt_id (default "P."). The function looks for columns matching `^P.prefix<name>.<s>$`.
#' @param state.idx.col name of the column in dt_id that holds the state index (value in 1..S). Default "state".
#' @param parameter character: name of the target process/outcome. If "target" (default), the first discovered terminal/outcome name is used.
#' @param process.deltas optional numeric vector (length = number of discovered processes). Not used internally for mapping, only validated if provided.
#' @param compute.clever logical; whether to compute and append clever.Q.<name>0 / clever.Q.<name>1 (default TRUE).
#' @return data.table: `dt_id` augmented with column `Q` and columns `clever.Q.<name>0` and `clever.Q.<name>1` for each discovered name.
#' @examples
#' # compute.Q.clever.per.id(dt_id = some_dt_for_one_id, states = depend.matrix, process.types = process.types)
#' @export
#' 
compute.Q.clever.per.id <- function(dt_id,
                                    states,
                                    process.types = NULL,
                                    P.prefix = "P.",
                                    state.idx.col = "state",
                                    parameter = "target",
                                    process.deltas = NULL,
                                    compute.clever = TRUE,
                                    browse = FALSE, browse2 = FALSE,
                                    get.years.lost = FALSE,
                                    years.lost.block.size = 10,
                                    clever.by.state = FALSE) {

    requireNamespace("data.table")
    # Defensive checks & normalization
    if (is.null(states) || nrow(states) == 0) {
        states <- data.table::data.table(state = 1L)
    }
    S <- length(unique(states[[state.idx.col]]))#nrow(states)
    if (S < 1) stop("'states' must have at least one row")

    Tn <- nrow(dt_id)
    if (Tn < 1) stop("dt_id must have at least one row (time-ordered)")

    # state columns (exclude 'at.risk' and the state index)
    state_cols_all <- setdiff(names(states), c("at.risk", state.idx.col, "target.only.in.state",
                                               sstates <- names(states)[substr(names(states), 1, 2) == "s."]))
    state_processes <- state_cols_all

    time_var_match <- names(states)[substr(names(states), 1, 5) == "time."]

    # canonical full ordering of state columns (same order as columns in 'states')
    full_state_cols <- state_cols_all

    derived_processes <- state_processes[!(state_processes %in% names(process.types))]
    j_derived <- which(state_processes == derived_processes)
    state_processes <- setdiff(state_processes, derived_processes)

    if ("target.only.in.state" %in% names(states)) {
        is.target.only.in.state <- TRUE
        target.only.in.state <- states[[ "target.only.in.state" ]]
    } else {
        is.target.only.in.state <- FALSE
    }

    # If the declared target is recurrent but not in states, ensure loops still consider it
    if (!is.null(process.types) && !is.null(parameter) &&
        parameter %in% names(process.types) &&
        identical(process.types[[parameter]], "recurrent") &&
        !(parameter %in% state_processes)) {
        state_processes <- c(parameter, state_processes)
    }
    # Track which processes actually have columns in 'states'
    stateful_in_states <- intersect(state_processes, state_cols_all)

    # Parse hazard columns in dt_id
    pcols <- grep(paste0("^", P.prefix, "[^\\.]+\\.[0-9]+$"), names(dt_id), value = TRUE)
    if (length(pcols) == 0L) {
        warning("No hazard columns found with prefix '", P.prefix, "'. Returning Q=0")
        out <- data.table::copy(dt_id)
        out[, Q := 0]
        return(out[])
    }
    re <- paste0("^", P.prefix, "([^\\.]+)\\.(\\d+)$")
    regs <- regmatches(pcols, regexec(re, pcols))
    found_names <- unique(vapply(regs, function(x) if (length(x) >= 3) x[2] else NA_character_, FUN.VALUE=""))
    found_names <- found_names[!is.na(found_names)]
    if (length(found_names) == 0L) stop("No hazard columns matched the expected pattern '", P.prefix, "<name>.<s>'")

    # Build columns_by_name and hazard matrices (Tn x S)
    columns_by_name <- setNames(lapply(found_names, function(x) rep(NA_character_, S)), found_names)
    for (i in seq_along(pcols)) {
        r <- regs[[i]]
        if (length(r) < 3) next
        nm <- r[2]; sidx <- as.integer(r[3])
        if (is.na(sidx) || sidx < 1L || sidx > S) next
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

    # Partition discovered names into stateful (those in state_processes) and terminal (the rest)
    stateful_names <- intersect(found_names, state_processes)
    terminal_names  <- setdiff(found_names, stateful_names)
    # Optionally restrict terminal_names to those declared terminal in process.types
    if (!is.null(process.types)) {
        declared_terminal <- names(process.types)[sapply(process.types, function(x) identical(x, "terminal"))]
        if (length(declared_terminal) > 0L) terminal_names <- intersect(terminal_names, declared_terminal)
    }
    stateful_names <- sort(stateful_names)
    terminal_names <- sort(terminal_names)

    # ---- build successor mappings
    # For the actual state columns (stateful_in_states) compute proc_caps and key vectors
    if (length(stateful_in_states) > 0L) {
        proc_caps <- integer(length(stateful_in_states))
        for (j in seq_along(stateful_in_states)) {
            vals <- as.integer(states[[ stateful_in_states[j] ]])
            proc_caps[j] <- max(vals, na.rm = TRUE)
        }
        if (length(derived_processes)>0) {
            stateful_derived <- c(stateful_in_states, derived_processes)
            key_vec_real <- do.call(paste, c(lapply(states[, ..stateful_derived], as.character), sep = ","))
        } else {
            key_vec_real <- do.call(paste, c(lapply(states[, ..stateful_in_states], as.character), sep = ","))
        }
        state_lookup_real <- setNames(states[[ state.idx.col ]], key_vec_real)
    } else {
        proc_caps <- integer(0)
        key_vec_real <- as.character(states[[ state.idx.col ]])
        state_lookup_real <- setNames(states[[ state.idx.col ]], key_vec_real)
    }

    # Stay mapping: maps row -> state index
    gamma_stay_idx <- unique(states[[ state.idx.col ]])
    # Build gamma_jump_idx aligned with state_processes
    M_proc_all <- length(state_processes)
    gamma_jump_idx <- vector("list", M_proc_all)
    ##browser()
    if (TRUE) {
        for (j_all in seq_len(M_proc_all)) {
            procj_name <- state_processes[j_all]
            if (procj_name %in% stateful_in_states & length(procj_name_s <- grep(procj_name, sstates, value = TRUE))) {
                gamma_jump_idx[[j_all]] <- states[[procj_name_s]]
            } else {
                gamma_jump_idx[[j_all]] <- seq_len(S)
            }
        }
    } else {
        for (j_all in seq_len(M_proc_all)) {
            procj_name <- state_processes[j_all]
            is_derived <- length(grep(procj_name, derived_processes))>0
            if (procj_name %in% stateful_in_states) {
                j_in <- which(stateful_in_states == procj_name)
                gamma_jump_idx[[j_all]] <- integer(S)
                for (s in seq_len(S)) {
                    parts <- as.integer(strsplit(key_vec_real[s], ",")[[1]])
                    parts[j_in] <- pmin(parts[j_in] + 1L, proc_caps[j_in])
                    if (is_derived) { # & parts[j_in] == proc_caps[j_in] ## <- fix: adapt to multiple!
                        for (derived_var in grep(procj_name, derived_processes, value = TRUE)) {
                            j_derived <- which(c(state_processes, derived_processes) == derived_var)
                            if (length(grep("T\\.", derived_var))>0) {
                                parts[j_derived] <- parts[j_derived]
                            } else {                   
                                parts[j_derived] <- 1-parts[j_derived] # <- this need be fixed for more general derived.vars
                            }
                        }
                    }
                    nk <- paste(parts, collapse = ",")
                    if (!(nk %in% names(state_lookup_real))) {
                        ## this was added to handle the case that some jumps are not admissible/possible
                        ##(is.null(state_lookup_real[[nk]])) {
                        ## stop("gamma mapping error: successor key not found: ", nk)
                        if (any(dt_id[["gamma.mapping.warning"]] == 1)) {
                            message("gamma mapping warning: successor key not found: ", nk)
                        }
                        parts <- as.integer(strsplit(key_vec_real[s], ",")[[1]])
                        parts[j_in] <- pmin(parts[j_in], proc_caps[j_in])
                        nk <- paste(parts, collapse = ",")
                    } 
                    gamma_jump_idx[[j_all]][s] <- as.integer(state_lookup_real[[nk]])

                }
            } else {
                # process not represented in states -> identity successor mapping
                gamma_jump_idx[[j_all]] <- seq_len(S)
            }
        }
    }
    
    # ---- classify target
    target_name <- parameter
    if (identical(parameter, "target")) {
        # choose default: first terminal found, else first found name
        if (length(terminal_names) > 0L) target_name <- terminal_names[1L] else target_name <- found_names[1L]
    }
    if (!(target_name %in% c(found_names, state_processes))) {
        stop("parameter must be one of discovered names or a recurrent target declared. Found: ", paste(c(found_names, state_processes), collapse = ", "))
    }
    target_in_states <- target_name %in% stateful_in_states
    if (!is.null(process.types) && target_name %in% names(process.types)) {
        target_type <- process.types[[target_name]]
    } else {
        target_type <- if (target_name %in% terminal_names) "terminal" else if (target_name %in% stateful_in_states) "recurrent" else "one.jump"
    }
    target_one_jump <- identical(target_type, "one.jump")
    target_terminal <- identical(target_type, "terminal")
    target_recurrent <- identical(target_type, "recurrent")
    target_with_atrisk <- identical(target_type, "state-with-atrisk")
    target_terminal_or_recurrent <- target_terminal | target_recurrent
    target_name_in_states <- target_name %in% names(states)

    is_recurrent_state_processes <-
        target_recurrent*(state_processes == target_name)

    ## target_can_jump_from_state <- function(s) {
    ##     if (target_with_atrisk && ("at.risk" %in% names(states))) {
    ##         isTRUE(states$at.risk[s])
    ##     } else TRUE
    ## }

    # ---- backward recursion

    ## ---------- Precompute hazard sums (one-time, outside tt-loop)
    if (length(terminal_names) > 0L) {
        Pout_sum_mat <- Reduce(`+`, lapply(terminal_names, function(nm) hazard_mats[[nm]]),
                               init = matrix(0, nrow = Tn, ncol = S))
    } else {
        Pout_sum_mat <- matrix(0, nrow = Tn, ncol = S)
    }
    if (length(stateful_names) > 0L) {
        Pproc_sum_mat <- Reduce(`+`, lapply(stateful_names, function(nm) hazard_mats[[nm]]),
                                init = matrix(0, nrow = Tn, ncol = S))
    } else {
        Pproc_sum_mat <- matrix(0, nrow = Tn, ncol = S)
    }

    zeroS <- rep(0, S)

    M <- length(state_processes)
    hazard_arr <- array(0, dim = c(Tn, S, M))
    
    for (j in seq_along(state_processes)) {
        nm <- state_processes[j]
        if (nm %in% names(hazard_mats)) {
            hazard_arr[,,j] <- hazard_mats[[nm]]   # matrix (Tn x S)
        } else {
            hazard_arr[,,j] <- matrix(0, nrow = Tn, ncol = S) # cheap once
        }
    }

    S1 <- nrow(states)
    gamma_jump_idx_mat <- matrix(NA_integer_, nrow = S1, ncol = M)
    for (j in seq_len(M)) gamma_jump_idx_mat[, j] <- gamma_jump_idx[[j]]

    states1 <- unique(states, by = "state")

    # --- attach rowwise Q
    out <- data.table::copy(dt_id)
    s_vec <- as.integer(out[[ state.idx.col ]])
    s_vec_row_index <- as.integer(out[[ paste0(state.idx.col, ".row.index") ]])

    all_names <- sort(unique(c(state_processes, terminal_names)))

    no_all <- length(all_names)

    if (target_name_in_states) {
        target_observed_by_state <- as.integer(states[[ target_name ]] >= 1)  # length S
    } else {
        target_observed_by_state <- integer(S)  # zeros
    }

    ##if (browse) browser()    

    if (length(time_var_match)>0) { ## fix: if there is more than one such variable

        if (TRUE) { ## THIS, is where I am currenty at (April, 2)

            state_times <- states[, time_var_match, with = FALSE]
            inf_rows <- lapply(time_var_match, function(time_var_match_jj) is.infinite(state_times[[(time_var_match_jj)]]))
            names(inf_rows) <- time_var_match
            
            uniq_time <- unique(dt_id[, time_var_match, with = FALSE])[, row := 1:.N]
            
            ### this one below should give row of uniq_time!
            tt_time <- uniq_time[dt_id, on = time_var_match][["row"]]
            ## tt_time <- lapply(time_var_match, function(time_var_match_jj)
            ##     dt_id[[time_var_match_jj]]
            ##     )

            rows_by_time <- vector("list", nrow(uniq_time))
            names(rows_by_time) <- uniq_time[["row"]]

            for (tm in 1:nrow(uniq_time)) {
                rows_by_time[[tm]] <- TRUE
                for (time_var_match_jj in time_var_match) {
                    rows_by_time[[tm]] <- rows_by_time[[tm]] &
                        (state_times[[time_var_match_jj]] == uniq_time[[time_var_match_jj]][tm] |
                         inf_rows[[time_var_match_jj]])
                }
            }

        } else {
        
            tt_time <- lapply(time_var_match, function(t_rep) c())#matrix(0, nrow = nrow(dt_id), ncol = length(time_var_match))
            rows_by_time <- lapply(time_var_match, function(t_rep) c())
            names(tt_time) <- names(rows_by_time) <- time_var_match
            for (time_var_match_jj in time_var_match) {
                state_time <- states[[time_var_match_jj]]
                inf_rows   <- is.infinite(state_time)
                tt_time[[time_var_match_jj]] <- dt_id[[time_var_match_jj]]
                uniq_time <- unique(tt_time[[time_var_match_jj]])
                rows_by_time[[time_var_match_jj]] <- setNames(vector("list", length(uniq_time)),
                                                              as.character(uniq_time))
                for (tm in uniq_time) {
                    rows_by_time[[time_var_match_jj]][[as.character(tm)]] <-
                        (state_time == tm | inf_rows)#c(inf_rows, which(state_time == tm))
                }
            }
        }
    }

    if (FALSE) {
        st_time_mat <- as.matrix(states[, ..time_var_match])
        st_inf_mat  <- is.infinite(st_time_mat)
        cache <- new.env(parent = emptyenv())

        get_rows_gamma <- function(tt) {
            key <- paste(vapply(time_var_match, function(v) dt_id[[v]][tt], numeric(1)), collapse = "|")
            ans <- cache[[key]]
            if (!is.null(ans)) return(ans)

            cur <- vapply(time_var_match, function(v) dt_id[[v]][tt], numeric(1))
            cmp <- sweep(st_time_mat, 2, cur, FUN = "==") | st_inf_mat
            ans <- rowSums(cmp) == length(time_var_match)

            cache[[key]] <- ans
            ans
        }
    }
    
    compute.Q.up.to.Tn <- function(Tn) {

        # terminal condition at last row: use hazard if present or zero
        if (target_name %in% names(hazard_mats)) {
            if (is.target.only.in.state) {
                Q_T <- target.only.in.state*hazard_mats[[ target_name ]][Tn, ]
            } else {
                Q_T <- hazard_mats[[ target_name ]][Tn, ]
            }
        } else {
            Q_T <- rep(0, S)
        }
        if (target_one_jump && target_name_in_states) {
            Q_T[ as.integer(states1[[ target_name ]]) >= 1 ] <- 1
        }
        
        Q_mat <- matrix(0, nrow = Tn, ncol = S)
        Q_mat[Tn, ] <- Q_T
        Q_row <- numeric(Tn)

        if (compute.clever) {
            # ---------- clever covariates
            clever0 <- matrix(NA_real_, nrow = Tn, ncol = length(all_names), dimnames = list(NULL, all_names))
            clever1 <- matrix(NA_real_, nrow = Tn, ncol = length(all_names), dimnames = list(NULL, all_names))
            if (!target_recurrent && target_name_in_states) { ## <- 20/3: this needs checking: !target_recurrent && 
                target_observed_by_state <- as.integer(states1[[ target_name ]] >= 1)  # length S
            } else {
                target_observed_by_state <- integer(S)  # zeros
            }
        }
        ##if (browse) browser()

        if (Tn>1) {
            for (tt in (Tn - 1):1) {

                Qn <- Q_mat[tt+1,]

                Pout_sum <- Pout_sum_mat[tt, ]
                Pproc_sum <- Pproc_sum_mat[tt, ]

                Pstay <- 1 - Pout_sum - Pproc_sum

                Qt <- Pstay * Qn[ gamma_stay_idx ]

                P_states <- hazard_arr[tt, ,]

                ## if (browse) browser()
                if (FALSE) {
                    if (length(time_var_match)>0) { ## fix: if there is more than one such variable 
                        rows_gamma <- lapply(time_var_match, function(t_rep) c())
                        names(rows_gamma) <- time_var_match
                        for (time_var_match_jj in time_var_match) {
                            rows_gamma[[time_var_match_jj]] <-
                                rows_by_time[[time_var_match_jj]][[as.character(tt_time[[time_var_match_jj]][tt])]]
                        }
                        rows_gamma <- apply(do.call("cbind", rows_gamma), 1, function(row) all(row)) 
                        Qn_jump_mat <- matrix(Qn[gamma_jump_idx_mat[rows_gamma,]], nrow = S, ncol = M)
                    } else {
                        Qn_jump_mat <- matrix(Qn[gamma_jump_idx_mat], nrow = S, ncol = M)
                    }
                } else {
                    if (length(time_var_match) > 0L) {
                        if (TRUE) {

                            if (TRUE) { ## THIS, is where I am currenty at (April, 2)

                                rows_gamma <- rows_by_time[[as.character(tt_time[[tt]])]]
                                Qn_jump_mat <- matrix(Qn[gamma_jump_idx_mat[rows_gamma,]], nrow = S, ncol = M)
                                
                            } else {
                                rows_gamma_list <- lapply(time_var_match, function(v) {
                                    rows_by_time[[v]][[as.character(tt_time[[v]][tt])]]
                                })

                                rows_gamma <- Reduce(`&`, rows_gamma_list)

                                idx_gamma <- which(rows_gamma)

                                ##Qn_jump_mat <- matrix(0, nrow = S, ncol = M)
                                if (length(idx_gamma) > 0L) {
                                    Qn_jump_mat <- matrix(
                                        Qn[gamma_jump_idx_mat[idx_gamma, , drop = FALSE]],
                                        nrow = length(idx_gamma),
                                        ncol = M
                                    )
                                }
                            }
                        } else {
                            rows_gamma <- get_rows_gamma(tt)
                            idx_gamma <- which(rows_gamma)
                            Qn_jump_mat <- matrix(Qn[gamma_jump_idx_mat[idx_gamma, , drop = FALSE]],
                                                  nrow = length(idx_gamma), ncol = M)
                        }
                    } else {
                        Qn_jump_mat <- matrix(Qn[gamma_jump_idx_mat], nrow = S, ncol = M)
                    }
                }
                
                Qt <- Qt + rowSums(P_states * Qn_jump_mat)

                # add immediate reward for terminal and recurrent targets
                if (target_terminal_or_recurrent) {
                    if (is.target.only.in.state) {
                        Qt <- Qt + target.only.in.state*hazard_mats[[ target_name ]][tt, ]
                    } else {
                        Qt <- Qt + hazard_mats[[ target_name ]][tt, ]
                    }
                }

                # enforce one-jump terminal states -> Q=1 where already observed
                if (target_one_jump && target_name_in_states) {
                    already_idx <- which(as.integer(states1[[ target_name ]]) >= 1)
                    if (length(already_idx) > 0L) Qt[already_idx] <- 1
                }
                
                Q_mat[tt, ] <- Qt

            }
        } 

        for (tt in seq_len(Tn)) {
            s <- s_vec[tt]
            Q_row[tt] <- Q_mat[tt,s]
        }
        ##Q_row <- Q_mat[cbind(seq_len(Tn), s_vec)]

        clever.fun <- function(s_vec) {

            # --- vectorized clever computation (after Q_mat is computed) ---
            clever0 <- matrix(NA_real_, nrow = Tn, ncol = no_all, dimnames = list(NULL, all_names))
            clever1 <- matrix(NA_real_, nrow = Tn, ncol = no_all, dimnames = list(NULL, all_names))

            tt_idx <- seq_len(max(1, Tn-1))   # we compute for tt = 1..Tn-1 (use Q_mat[tt+1, ...])
            m <- (Tn>1)*length(tt_idx)

            # 1) clever0 rows for tt=1..Tn-1:
            if (m > 0) {
                s_no_vec <- gamma_stay_idx[s_vec[1:m]]               # length m
                # value for each tt: Q_mat[tt+1, s_no_vec[tt]]
                clever0_vals <- Q_mat[cbind(tt_idx + 1L, s_no_vec)]
                # fill rows 1..m with repeated scalar across columns
                clever0[1:m, ] <- matrix(rep(clever0_vals, times = no_all), nrow = m, ncol = no_all)
            }
            # last row (tt = Tn): use Q_mat[Tn, ..] as in original code
            clever0[Tn, ] <- Q_mat[Tn, gamma_stay_idx[s_vec[Tn]]]    # scalar repeated

            # 2) terminal processes:
            if (length(terminal_names) > 0L) {
                # default for each tt: target_observed_by_state[s_vec[tt]]
                term_base <- target_observed_by_state[s_vec]          # length Tn
                # make m x p_term matrix
                if (m > 0) clever1[1:m, terminal_names] <- matrix(rep(term_base[1:m],
                                                                      length(terminal_names)),
                                                                  nrow = m,
                                                                  ncol = length(terminal_names))
                clever1[Tn, terminal_names] <- term_base[Tn]

                # if the target is a terminal and should be set to 1:
                if (target_terminal && target_name %in% terminal_names) {
                    if (is.target.only.in.state) {
                        clever1[, target_name] <- 1*target.only.in.state[s_vec[1:(m+1)]]
                    } else {
                        clever1[, target_name] <- 1L
                    }
                }
            }

            ##if (browse) browser()
            # 3) stateful processes (vectorized)
            if (length(state_processes) > 0L) {
                # succ_idx_mat: Tn x M where each row tt lists successor state indices for s_vec[tt]
                if (FALSE) {
                    if (length(time_var_match) > 0L) {
                        if (TRUE) {
                            st <- data.table::copy(states)[, row := .I]
                            data.table::setkeyv(st, c(state.idx.col, time_var_match))
                            dt_key <- dt_id[, .(state = s_vec, time = get(time_var_match))]
                            s_exact <- st[dt_key, row]
                            s_inf   <- st[data.table::data.table(state = s_vec, time = Inf), row]
                            s_lookup <- data.table::fcoalesce(s_exact, s_inf)
                        } else if (FALSE) {
                            s_lookup <- lapply(time_var_match, function(t_rep) c())
                            names(s_lookup) <- time_var_match
                            for (time_var_match_jj in time_var_match) {
                                st <- data.table::copy(states)[, row := .I]
                                data.table::setkeyv(st, c(state.idx.col, time_var_match))
                                dt_key <- dt_id[, .(state = s_vec, time = get(time_var_match_jj))]
                                s_exact <- st[dt_key, row]
                                s_inf   <- st[data.table::data.table(state = s_vec, time = Inf), row]
                                s_lookup[[time_var_match_jj]] <- data.table::fcoalesce(s_exact, s_inf)
                            }
                        } else {                       
                            st <- data.table::copy(states)[, row := .I]
                            data.table::setkeyv(st, c(state.idx.col, time_var_match))
                            dt_key <- dt_id[, .(state = s_vec)][, time := Inf]
                            for (time_var_match_jj in time_var_match) {
                                s_jj <- s_vec %in% st[st[[time_var_match_jj]]<Inf][[state.idx.col]]
                                dt_key[["time"]][s_jj] <- dt_id[[time_var_match_jj]][s_jj]
                                ##dt_key <- dt_id[, .(state = s_vec, time = get(time_var_match_jj))]
                            }
                            ##dt_key <- dt_id[, .(state = s_vec, time = get(time_var_match_jj))]
                            ##dt_key <- dt_id[, .(state = s_vec, time = get(time_var_match_jj))]
                            ##dt_key <- dt_id[, c(list(state = s_vec), .SD), .SDcols = time_var_match]
                            s_exact <- st[dt_key, row]
                            s_inf   <- st[data.table::data.table(state = s_vec, time = Inf), row]
                            s_lookup <- data.table::fcoalesce(s_exact, s_inf)
                        }
                    } else {
                        s_lookup <- s_vec
                    }
                }

                ##succ_idx_mat <- gamma_jump_idx_mat[s_lookup, , drop = FALSE]   # Tn x M
                succ_idx_mat <- gamma_jump_idx_mat[s_vec_row_index, , drop = FALSE]   # Tn x M
                
                if (m > 0) {
                    # build long index to extract Q_mat[tt+1, succ_idx_mat[tt, j]] for all tt,j
                    rows_long  <- rep(tt_idx + 1L, times = ncol(succ_idx_mat))  # length m*M
                    cols_long  <- as.integer(as.vector(succ_idx_mat[1:m, , drop = FALSE])) # column-major flatten
                    vals_long  <- Q_mat[cbind(rows_long, cols_long)]
                    # reshape to m x M
                    clever1_stateful <- matrix(vals_long, nrow = m, ncol = ncol(succ_idx_mat))
                    # add recurrence offsets if present
                    if (any(is_recurrent_state_processes != 0L)) {
                        clever1_stateful <- clever1_stateful + matrix(rep(is_recurrent_state_processes, each = m),
                                                                      nrow = m, ncol = length(is_recurrent_state_processes), byrow = FALSE)
                    }
                    # assign into clever1 rows 1..m for columns named state_processes
                    clever1[1:m, state_processes] <- clever1_stateful
                }
                # last row: tt = Tn
                clever1[Tn, state_processes] <- Q_mat[Tn, as.integer(succ_idx_mat[Tn, ])] +
                    (is_recurrent_state_processes * 1) # if recurrence offset (careful with types)
            }

            return(list(clever1 = clever1,
                        clever0 = clever0))
        }
            
        # At this point clever0/c1 are filled; they match the row-wise scalar logic in your loop.

        if (compute.clever) {
            clever.s.obs <- clever.fun(s_vec)
            out.list <- list(Q_row = Q_row,
                             clever1 = clever.s.obs$clever1,
                             clever0 = clever.s.obs$clever0)
            if (clever.by.state) {
                for (s in seq_len(S)) {
                    clever.s <- clever.fun(rep(s, length = length(s_vec)))
                    out.list[[paste0("clever1.", s)]] <- clever.s$clever1
                    out.list[[paste0("clever0.", s)]] <- clever.s$clever0
                    out.list[[paste0("Q.", s)]] <- Q_mat[, s]
                }
                return(out.list)
            } else{
                for (s in seq_len(S)) {
                    out.list[[paste0("Q.", s)]] <- Q_mat[, s]
                }
                return(out.list)
            }
        } else {
            return(list(Q_row = Q_row))
        }
    }

    if (browse) browser()

    Q.out <- compute.Q.up.to.Tn(Tn)
  
    out[, Q := Q.out$Q_row]

    if (compute.clever) {
        # attach clever columns
        for (nm in all_names) {
            out[[ paste0("clever.Q.", nm, "0") ]] <- Q.out$clever0[, nm] # <- fix: can remove again
            out[[ paste0("clever.Q.", nm, "1") ]] <- Q.out$clever1[, nm] # <- fix: can remove again
            out[[ paste0("clever.Q.", nm) ]] <- Q.out$clever1[, nm] - Q.out$clever0[, nm]
            if (clever.by.state) {
                for (s in seq_len(S)) {
                    out[[ paste0("clever.Q.", nm, ".", s) ]] <-
                        Q.out[[paste0("clever1.", s)]][, nm] - Q.out[[paste0("clever0.", s)]][, nm]
                }
            }
        }
        for (s in seq_len(S)) {
            out[[ paste0("Q.", s) ]] <-
                Q.out[[paste0("Q.", s)]]
        }
    }

    if (get.years.lost) {

        # k block size for approximation
        k <- years.lost.block.size
        starts <- unique(c(seq(1L, Tn, by = k), Tn))
        ends   <- pmin(starts + k - 1L, Tn)
        # safe representative: use block right endpoint (>= every tt in block)
        ## rep_idx <- ends
        rep_idx <- floor((starts + ends) / 2)

        # map every time index to block index 1..length(starts)
        block_of <- findInterval(1:Tn, starts)  # block index for each tt (starts sorted)

        Q.years.lost <- 0
        clever.years.lost.matrix <- matrix(0, nrow = Tn, ncol = length(all_names),
                                           dimnames = list(NULL, all_names))

        if (clever.by.state) {
            clever.years.lost.list <- list()
            for (s in seq_len(S)) {
                clever.years.lost.list[[s]] <- matrix(0, nrow = Tn, ncol = length(all_names),
                                                           dimnames = list(NULL, all_names))
            }
        }
        
        dt.vec <- diff(c(0, dt_id[["time"]]))    # length Tn
        last_block <- NA_integer_
        Q.out.block <- NULL

        # loop descending as you had it
        for (tt in Tn:1) {

            b <- block_of[tt]   # which block tt belongs to
            if (is.na(last_block) || b != last_block) {
                # new block encountered -> compute once for block representative
                m_rep <- rep_idx[b]     # this is >= tt for any tt in that block
                Q.out.block <- compute.Q.up.to.Tn(m_rep)   # returns list(Q_row, clever1, clever0)
                last_block <- b
            }

            # reuse Q.out.block for all tt in the same block
            Q.years.lost <- Q.years.lost + dt.vec[tt] * Q.out.block$Q_row[1]

            # update rows 1:tt with today's contribution (this is heavy; see next paragraph)

            clever_diff <- Q.out.block$clever1 - Q.out.block$clever0    # matrix m_rep x P
            
            nrep <- nrow(clever_diff)      # = m_rep
            imax <- min(tt, nrep)
            
            clever.years.lost.matrix[1:imax, ] <- clever.years.lost.matrix[1:imax, ] +
                dt.vec[tt] * out[[ "at.risk" ]][1:imax] * clever_diff[1:imax, ]

            if (clever.by.state) {
                for (s in seq_len(S)) {
                    clever_diff <- Q.out.block[[paste0("clever1.", s)]] - Q.out.block[[paste0("clever0.", s)]]    # matrix m_rep x P
                    clever.years.lost.list[[s]][1:imax, ] <- clever.years.lost.list[[s]][1:imax, ] +
                        dt.vec[tt] * out[[ "at.risk" ]][1:imax] * clever_diff[1:imax, ]
                }
            }
        }

        out [[ "Q" ]] <- Q.years.lost

        for (nm in all_names) {
            out[[ paste0("clever.Q.years.lost.", nm) ]] <- clever.years.lost.matrix[, nm]
            if (clever.by.state) {
                for (s in seq_len(S)) {
                    out[[ paste0("clever.Q.years.lost.", nm, ".", s) ]] <-
                        clever.years.lost.list[[s]][, nm]
                }
            }
        }
    }

    if (browse2) browser()
       
    return(out[])
}

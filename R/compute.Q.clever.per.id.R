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
                                    years.lost.block.size = 10) {

    requireNamespace("data.table")
    # Defensive checks & normalization
    if (is.null(states) || nrow(states) == 0) {
        states <- data.table::data.table(state = 1L)
    }
    S <- nrow(states)
    if (S < 1) stop("'states' must have at least one row")

    Tn <- nrow(dt_id)
    if (Tn < 1) stop("dt_id must have at least one row (time-ordered)")

    # state columns (exclude 'at.risk' and the state index)
    state_cols_all <- setdiff(names(states), c("at.risk", state.idx.col))
    state_processes <- state_cols_all

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
        key_vec_real <- do.call(paste, c(lapply(states[, ..stateful_in_states], as.character), sep = ","))
        state_lookup_real <- setNames(states[[ state.idx.col ]], key_vec_real)
    } else {
        proc_caps <- integer(0)
        key_vec_real <- as.character(states[[ state.idx.col ]])
        state_lookup_real <- setNames(states[[ state.idx.col ]], key_vec_real)
    }

    # Stay mapping: maps row -> state index
    gamma_stay_idx <- states[[ state.idx.col ]]

    # Build gamma_jump_idx aligned with state_processes
    M_proc_all <- length(state_processes)
    gamma_jump_idx <- vector("list", M_proc_all)
    for (j_all in seq_len(M_proc_all)) {
        procj_name <- state_processes[j_all]
        if (procj_name %in% stateful_in_states) {
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
            # process not represented in states -> identity successor mapping
            gamma_jump_idx[[j_all]] <- seq_len(S)
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

    target_can_jump_from_state <- function(s) {
        if (target_with_atrisk && ("at.risk" %in% names(states))) {
            isTRUE(states$at.risk[s])
        } else TRUE
    }

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

    gamma_jump_idx_mat <- matrix(NA_integer_, nrow = S, ncol = M)
    for (j in seq_len(M)) gamma_jump_idx_mat[, j] <- gamma_jump_idx[[j]]

    # --- attach rowwise Q
    out <- data.table::copy(dt_id)
    s_vec <- as.integer(out[[ state.idx.col ]])

    all_names <- sort(unique(c(state_processes, terminal_names)))

    no_all <- length(all_names)

    if (target_name_in_states) {
        target_observed_by_state <- as.integer(states[[ target_name ]] >= 1)  # length S
    } else {
        target_observed_by_state <- integer(S)  # zeros
    }

    compute.Q.up.to.Tn <- function(Tn) {

        # terminal condition at last row: use hazard if present or zero
        if (target_name %in% names(hazard_mats)) {
            Q_T <- hazard_mats[[ target_name ]][Tn, ]
        } else {
            Q_T <- rep(0, S)
        }
        if (target_one_jump && target_name_in_states) {
            Q_T[ as.integer(states[[ target_name ]]) >= 1 ] <- 1
        }
        
        Q_mat <- matrix(0, nrow = Tn, ncol = S)
        Q_mat[Tn, ] <- Q_T
        Q_row <- numeric(Tn)

        if (compute.clever) {
            # ---------- clever covariates
            clever0 <- matrix(NA_real_, nrow = Tn, ncol = length(all_names), dimnames = list(NULL, all_names))
            clever1 <- matrix(NA_real_, nrow = Tn, ncol = length(all_names), dimnames = list(NULL, all_names))
            if (target_name_in_states) {
                target_observed_by_state <- as.integer(states[[ target_name ]] >= 1)  # length S
            } else {
                target_observed_by_state <- integer(S)  # zeros
            }
        }

        if (Tn>1) {
            for (tt in (Tn - 1):1) {

                Qn <- Q_mat[tt+1,]

                Pout_sum <- Pout_sum_mat[tt, ]
                Pproc_sum <- Pproc_sum_mat[tt, ]

                Pstay <- 1 - Pout_sum - Pproc_sum

                Qt <- Pstay * Qn[ gamma_stay_idx ]

                P_states <- hazard_arr[tt, ,]

                Qn_jump_mat <- matrix(Qn[gamma_jump_idx_mat], nrow = S, ncol = M)
                
                Qt <- Qt + rowSums(P_states * Qn_jump_mat)
        
                # add immediate reward for terminal and recurrent targets
                if (target_terminal_or_recurrent) {
                    Qt <- Qt + hazard_mats[[ target_name ]][tt, ]
                }

                # enforce one-jump terminal states -> Q=1 where already observed
                if (target_one_jump && target_name_in_states) {
                    already_idx <- which(as.integer(states[[ target_name ]]) >= 1)
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
            if (m > 0) clever1[1:m, terminal_names] <- matrix(rep(term_base[1:m], length(terminal_names)),
                                                              nrow = m, ncol = length(terminal_names))
            clever1[Tn, terminal_names] <- term_base[Tn]

            # if the target is a terminal and should be set to 1:
            if (target_terminal && target_name %in% terminal_names) {
                clever1[, target_name] <- 1L
            }
        }

        # 3) stateful processes (vectorized)
        if (length(state_processes) > 0L) {
            # succ_idx_mat: Tn x M where each row tt lists successor state indices for s_vec[tt]
            succ_idx_mat <- gamma_jump_idx_mat[s_vec, , drop = FALSE]   # Tn x M

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

        # At this point clever0/c1 are filled; they match the row-wise scalar logic in your loop.

        if (compute.clever) {
            return(list(Q_row = Q_row,
                        clever1 = clever1,
                        clever0 = clever0))
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
            out[[ paste0("clever.Q.", nm, "0") ]] <- Q.out$clever0[, nm]
            out[[ paste0("clever.Q.", nm, "1") ]] <- Q.out$clever1[, nm]
        }
    }

    weight.vec <- out[["clever.weight"]]*out[["clever.weight.alpha"]]

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
                dt.vec[tt] * weight.vec[1:imax] * out[[ "at.risk" ]][1:imax] * clever_diff[1:imax, ]
        }

        out [[ "Q" ]] <- Q.years.lost

        for (nm in all_names) {
            out[[ paste0("clever.Q.years.lost.", nm) ]] <- clever.years.lost.matrix[, nm]
        }
    }

    
    if (browse2) browser()

    if (FALSE & get.years.lost) {

            Q.years.lost <- 0
            clever.years.lost.matrix <- matrix(0, nrow = Tn, ncol = length(all_names),
                                               dimnames = list(NULL, all_names))
            dt.vec <- diff(c(0, dt_id[["time"]]))
            for (tt in Tn:1) {

                Q.out.tt <- compute.Q.up.to.Tn(tt) # gets E[N^x(t)|S_k=S]

            
                Q.years.lost <- Q.years.lost + dt.vec[tt]*Q.out.tt$Q_row[1] # contribution to integral: E[N^x(t)|S_0=S]dt
                clever.years.lost.matrix[1:tt,] <- clever.years.lost.matrix[1:tt,] +
                    dt.vec[tt]*weight.vec[tt]*(Q.out.tt$clever1[, ] - Q.out.tt$clever0[, ]) # contribution to clever covariates: (clever1_k_t - clever0_k_t)dt
            
            }

    }

    if (FALSE & get.years.lost) {
        
        source("./try/compute.Q.up.to.m.R")

        ################################

        save.Q <- numeric(Tn)
        save.k.matrix <- matrix(0, nrow = Tn, ncol = length(all_names),
                                dimnames = list(NULL, all_names))
        dt_vec <- c(diff(dt_id[["time"]]), 0)

        for (m in 1:(Tn - 1)) {
            # Q_up_to_m: rows 1..m, cols 1..S
            Q_up_to_m <- compute.Q.up.to.m(m, Q_mat, hazard_arr, gamma_jump_idx_mat,
                                           Pout_sum_mat, Pproc_sum_mat, gamma_stay_idx)

            dt_m <- dt_vec[m]

            # for each start k <= m: contribution for the time slice (T_m, T_{m+1}]
            for (k in seq_len(nrow(Q_up_to_m))) {
                s_k <- s_vec[k]                # observed state at time k
                tmp.k <- Q_up_to_m[k, ]        # cumulative up to m when starting (k, *)
                s_no <- gamma_stay_idx[s_k]    # successor index if no jump

                # subtracted value (expected cumulative up to m under no-jump)
                sub_val <- tmp.k[s_no]

                # terminal processes: scalar per terminal
                clever1_term_vals <- rep(as.integer(target_observed_by_state[s_k]), length(terminal_names))
                if (target_terminal && target_name %in% terminal_names) {
                    ix <- which(terminal_names == target_name)
                    clever1_term_vals[ix] <- 1L
                }

                # update terminal columns (vectorized over terminal_names)
                save.k.matrix[k, terminal_names] <- save.k.matrix[k, terminal_names] +
                    dt_m * (clever1_term_vals - sub_val)

                # stateful processes
                succ_idx_vec <- gamma_jump_idx_mat[s_k, ]   # length M
                clever1_vals_stateful <- tmp.k[succ_idx_vec] + is_recurrent_state_processes
                save.k.matrix[k, state_processes] <- save.k.matrix[k, state_processes] +
                    dt_m * (clever1_vals_stateful - sub_val)

                # accumulate integrateted expectation for start k (if you want it)
                # Q_up_to_m[k, s_k] is E[N(t)] at times in the m-th interval for start k
                save.Q[k] <- save.Q[k] + dt_m * tmp.k[s_k]
            }
        }

        out <- cbind(out, save.k.matrix)

        out[, Q.years.lost := save.Q]

    }
       
    return(out[])
}

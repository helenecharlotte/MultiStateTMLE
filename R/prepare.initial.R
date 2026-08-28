### prepare.initial.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Feb  4 2026 (08:47) 
## Version: 
## Last-Updated: Aug 28 2026 (09:55) 
##           By: Helene
##     Update #: 920
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

prepare.initial <- function(dt,
                            tau = 1.2, 
                            fit.types = list(
                                Z = list(model = "Surv(tstart, tstop, delta == 2)~L0+L+A0",
                                         fit = "cox",
                                         at.risk = function(dt) (dt[["Z"]] == 0)),
                                L = list(model = "Surv(tstart, tstop, delta == 3)~L0+Z+A0",
                                         fit = "cox",
                                         at.risk = function(dt) (dt[["L"]] == 0)),
                                #W = list(model = "Surv(tstart, tstop, delta == 4)~L0+Z",
                                #         fit = "cox",
                                #         at.risk = function(dt) (dt[["W"]] == 0)),
                                outcome1 = list(model = "Surv(tstart, tstop, delta == 1)~L0+Z+L+A0",
                                                fit = "cox"),
                                #death = list(model = "Surv(tstart, tstop, delta == 5)~L0+Z+L",
                                #             fit = "cox"),
                                censoring = list(model = "Surv(tstart, tstop, delta == 0)~L0+Z+L+A0",
                                                 fit = "cox")
                            ),
                            use.exponential = FALSE,
                            verbose.exponential = FALSE,
                            browse = FALSE,
                            max.count = 1,
                            a = NULL,
                            derived.vars = list(),
                            depend.time = list(),
                            fit.treatment = NULL,
                            prune.states = FALSE,
                            ##### hal parameters:
                            cut.time = 35,
                            cut.one.way = 15,
                            cut.Tk = 3, cut.Tk.values = NULL,
                            max.Tk = max.count,
                            two.way = NULL,
                            reduce.NK = 0.025,
                            hal.sl = list(c(cut.one.way = cut.one.way,
                                            cut.Tk = cut.Tk,
                                            cut.Tk.values = cut.Tk.values,
                                            two.way = two.way)),
                            screen.two.way = FALSE,
                            lambda.cvs = c((9:1)/10,(9:2)/10^2, seq(1/10^2, 1/10^5, length = 100)),#c(sapply(1:5, function(jjj) (9:1)/(10^jjj))),
                            event.dependent.cv = FALSE,
                            npenalize.vars = NULL,
                            V = 10,
                            seed.hal = NULL,
                            reduce.seed.dependence = FALSE,
                            penalize.time = FALSE,
                            use.cores = 50,
                            use.cores.prediction = 1,#5,
                            verbose.hal = FALSE,
                            browse.hal = FALSE,
                            cv.glmnet = FALSE, #currently not supported
                            #####################
                            verbose = FALSE,
                            return.parameters.for.simulation = FALSE) {

    dt <- copy(dt)

    n <- length(unique(dt[["id"]]))

    #-- if no "fit" is specified for types, then specify it as cox:
    for (nm in names(fit.types)) {
        if (!is.list(fit.types[[nm]])) fit.types[[nm]] <- list(model = fit.types[[nm]], fit = "cox")
    }

    any.hal <- any(which.hal <- sapply(fit.types, function(fit.type) fit.type$fit == "hal"))

    #-- if no "fit" is specified for treatment model, then specify it as glm: 
    if (!is.list(fit.treatment)) {
        fit.treatment <- list(model = fit.treatment, fit = "glm")
    }

    at.risk.ids <- (1:length(fit.types))[sapply(fit.types, function(fit.type) ("at.risk" %in% names(fit.type)))]

    extract_dt_cols <- function(at_risk_fun) {
        expr <- body(at_risk_fun)

        cols <- character()

        recurse <- function(e) {
            if (is.call(e) && identical(e[[1]], as.name("[["))) {
                # dt[["Z"]] → extract "Z"
                if (length(e) == 3 && is.character(e[[3]])) {
                    cols <<- c(cols, e[[3]])
                }
            }
            if (is.call(e) || is.expression(e)) {
                lapply(as.list(e), recurse)
            }
        }

        recurse(expr)
        unique(cols)
    }

    #-- covariates/predictors extracted from models: 
    varnames <- unique(c(
        unlist(lapply(fit.types, function(fit.type) {
            if (length(fit.type[["model"]])>0) return(strsplit(strsplit(fit.type[["model"]], "~")[[1]][2], "\\+|\\*")[[1]])
        })),
        #unlist(as.character(lapply(fit.types[at.risk.ids], function(fit.type) {
        #    setdiff(as.character(body(fit.type[["at.risk"]])[[2]][[2]]), c("[[", "dt"))
        #})))
        unlist(lapply(fit.types[at.risk.ids], function(fit.type) {
            extract_dt_cols(fit.type[["at.risk"]])
        }))
    ))

    varnames <- sapply(varnames, function(varname) gsub(" ", "", varname))

    N.vars <- NULL

    #-- remove time vars for now
    if (length(depend.time)>0) {
        for (time.jj in 1:length(depend.time)) {
            varnames <- setdiff(varnames, paste0("T.", names(depend.time)[time.jj]))
        }
    }

    #-- in case of interactions: 
    if (length(interaction.names <- grep(":", varnames, value = TRUE))>0) {
        varnames <- unique(c(varnames[!varnames %in% interaction.names],
                             str_split(interaction.names, ":")[[1]]))
    }

    #-- corresponding delta values: 
    process.deltas <- c(as.numeric(unlist(lapply(fit.types, function(fit.type) {
        tmp.delta <- strsplit(fit.type[["model"]], "~")[[1]][1]
        as.numeric(substr(tmp.delta, nchar(tmp.delta)-1, nchar(tmp.delta)-1))
    }))))
    
    #-- names of types
    process.names <- names(fit.types)
    state.names <- process.names[process.names %in% varnames[varnames %in% process.names]]
    state.deltas <- process.deltas[process.names %in% state.names]

    #-- process types
    process.types <- lapply(1:length(process.names), function(process.jj) {
        delta.jj <- process.deltas[process.jj]
        is.terminal <- all(dt[, (delta == delta.jj)*(time < max(time)), by = "id"][[2]] == 0)
        is.recurrent <- any(dt[delta == delta.jj, .N, by = "id"][["N"]]>1)
        return(ifelse(is.terminal, "terminal", ifelse(is.recurrent, "recurrent", "one.jump")))
    })
    names(process.types) <- process.names

    at.risks <- lapply(fit.types[at.risk.ids], function(fit.type) fit.type[["at.risk"]])
    names(at.risks) <- names(fit.types)[at.risk.ids]

    which.recurrent <- intersect(state.names, process.names[sapply(process.types, function(process.type) process.type == "recurrent")])

    #-- name of baseline treatment variable:
    if (length(fit.treatment[["model"]])>0) {
        A0.name <- strsplit(fit.treatment[["model"]], "~")[[1]][1]
        varnames <- unique(c(varnames, A0.name))
    } else {
        A0.name <- NULL 
    }

    if (typeof(dt[["time"]]) != "double") {
        warning("NB: the time variable is not numeric - will be converted")
        dt[, time := as.numeric(time)]
    }

    ##dt[diff(c(0, time)) == 0, time := time+1e-6, by = "id"]
    
    dt[, tstart := c(0, time[-.N]), by = "id"]
    dt[, tstop := time]

    if (length(which.recurrent) > 0) {
        for (var.recurrent in which.recurrent) {
            dt[, (var.recurrent) := cumsum(delta == process.deltas[process.names == var.recurrent]), by = "id"]
        }
    }

    ### this added June 4
    if (length(state.names)>0) {
        for (process.jj in 1:length(state.names[state.names %in% names(fit.types)])) {
            dt[, (state.names[process.jj]) := cumsum(1*(delta == state.deltas[process.jj])), by = "id"]
        }
    }

    for (varname in state.names) {
        dt[, (varname) := c(0, get(varname)[-.N]), by = "id"]
    }

    if (length(which.recurrent) > 0) {
        for (var.recurrent in which.recurrent) {
            for (count.jj in 1:max.count) {
                dt[, (paste0(var.recurrent, ".", count.jj)) := (dt[[var.recurrent]] >= count.jj)]
            }
        }
    }

    if (length(derived.vars)>0) {
        for (derived.jj in 1:length(derived.vars)) {
            derived.var <- names(derived.vars)[derived.jj]
            derived.fun <- derived.vars[[derived.jj]]
            dt[, (derived.var) := derived.fun(dt)]
            state.names <- unique(c(state.names, derived.var))
        }
    }

    if (length(depend.time)>0) {
        for (time.jj in 1:length(depend.time)) {
            time.var <- names(depend.time)[time.jj]
            time.fun <- depend.time[[time.jj]]
            suppressWarnings(
                dt[, (paste0("T.", time.var)) := min(time[c(get(time.var)[-1],get(time.var)[.N]) == 1]), by = "id"]
            )
            dt[tstart < get(paste0("T.", time.var)), (paste0("T.", time.var)) := 0] #time.fun(tstart)
            dt[tstart >= get(paste0("T.", time.var)), (paste0("T.", time.var)) := as.numeric(time.fun(get(paste0("T.", time.var))))]
            varnames <- c(varnames, paste0("T.", time.var))
            state.names <- c(state.names, paste0("T.", time.var))
        }
    }

    #--------------------------------
    #-- intervention part; for clever weight estimation:
    # (we start with cox models, if HAL is specified this is fitted later)

    cens.process.id <- (1:length(process.deltas))[process.deltas == 0]

    dt[, idN := 1:.N, by = "id"]
    if (fit.treatment[["fit"]] == "glm" & length(fit.treatment[["model"]])>0) {
        fit.A0 <- glm(as.formula(fit.treatment[["model"]]), data=dt[idN == 1], family=binomial)
        if (verbose) print(summary(fit.A0))
    } else if (length(fit.treatment[["model"]])>0)  {
        print("NB: need to incorporate other estimations methods than glm for treatment")
    }

    #--------------------------------
    #-- outcome / clever covariate part:
    # (we start with cox models, if HAL is specified this is fitted later)
    ##browser()
    fit.cox.types <- lapply(1:length(fit.types), function(fit.type.jj) {
        model.jj <- fit.types[[fit.type.jj]][["model"]]
        if (length(depend.time)>0) {
            for (time.jj in 1:length(depend.time)) {
                model.jj <-
                    gsub(paste0("T.", names(depend.time)[time.jj]),
                         paste0("as.factor(", paste0("T.", names(depend.time)[time.jj]), ")"),
                         model.jj)
            }
        }
        if (length(which.recurrent) > 0) {
            for (var.recurrent in which.recurrent) {
                model.jj <- gsub(paste0("\\+", var.recurrent),
                                 paste0(paste0("+", var.recurrent, ".", 1:max.count), collapse = ""),
                                 model.jj)
            }
        }
        if (fit.type.jj %in% at.risk.ids) {
            tmp.cox <- coxph(as.formula(model.jj),
                             data = dt[fit.types[[fit.type.jj]][["at.risk"]](dt)], 
                             control = coxph.control(timefix = FALSE))
        } else {
            tmp.cox <- coxph(as.formula(model.jj),
                             data = dt,
                             control = coxph.control(timefix = FALSE))
        }
        if (verbose) message("------------------------------")
        if (verbose) message(paste0("fit.type = ", names(fit.types)[fit.type.jj]))
        if (verbose) print(tmp.cox)
        tmp.type <- suppressWarnings(setDT(basehaz(tmp.cox, centered=TRUE)))[, (paste0("dhazard.", names(fit.types)[fit.type.jj])) := c(hazard[1],diff(hazard))][, (paste0("hazard.", names(fit.types)[fit.type.jj])) := hazard][, -"hazard", with = FALSE]
        return(list(fit.cox = tmp.cox, tmp.type = tmp.type[tmp.type[[paste0("dhazard.", names(fit.types)[fit.type.jj])]]>0]))
    })

    if (return.parameters.for.simulation) {

        # baseline dataset: one row per subject
        baseline.dt <- dt[idN == 1]

        # baseline covariates
        baseline.vars <- setdiff(
            varnames,
            c(process.names, "time", "tstart", "tstop", "delta")
        )

        baseline.summary <- lapply(baseline.vars, function(v) {

            x <- baseline.dt[[v]]

            if (is.numeric(x)) {
                list(
                    type = "numeric",
                    mean = mean(x, na.rm = TRUE),
                    sd = sd(x, na.rm = TRUE),
                    median = median(x, na.rm = TRUE),
                    min = min(x, na.rm = TRUE),
                    max = max(x, na.rm = TRUE)
                )
            } else {
                tab <- table(x, useNA = "ifany")
                list(
                    type = class(x)[1],
                    frequencies = as.list(tab),
                    proportions = as.list(prop.table(tab))
                )
            }

        })
        
        names(baseline.summary) <- baseline.vars

        parameters.for.simulation <- lapply(fit.cox.types, function(fit.cox.type) {
            
            tmp.type <- fit.cox.type[["tmp.type"]]
            tmp.type[, Haz := cumsum(tmp.type[[2]])]
    
            weib_fit <- lm(log(Haz) ~ log(time), data = tmp.type)
            coef_weib <- coef(weib_fit)

            return(list(cox.parameters = coef(fit.cox.type[["fit.cox"]]),
                        weibull.parameters = c(eta = exp(coef_weib[1]),
                                               nu = coef_weib[2])))
        })

        names(parameters.for.simulation) <- names(fit.types)

        parameters.for.simulation$model.structure <- list(
            process.names = process.names,
            process.deltas = process.deltas,
            process.types = process.types,
            cens.process.id = cens.process.id,
            at.risks = at.risks)

        parameters.for.simulation[["baseline.summary"]] <- baseline.summary

        parameters.for.simulation[["n.subjects"]] <- n
       
        return(parameters.for.simulation)
    }

    #-- get all unique times; 
    unique.times <- sort(unique(dt[["time"]]))
    unique.times <- unique.times[unique.times <= tau]
    all.times <- data.table(expand.grid(time = unique.times,
                                        id = unique(dt[["id"]])))

    #-- collect data with all time-points;
    tmp.inner <- merge(dt[, time.obs := time], all.times, by = c("id", "time"), all = TRUE)[order(id, time.obs)]
    
    for (varname in varnames) {
        tmp.inner[, (varname) := na.locf(get(varname)), by = "id"]
    }

    #--------------------------------
    #-- expanded dataset to work with for TMLE:

    tmp.long <- tmp.inner[order(id, time)][is.na(delta), delta := 0][, -c("tstop"), with = FALSE]

    for (fit.type.jj in 1:length(fit.types)) {
        tmp.long <- merge(tmp.long, fit.cox.types[[fit.type.jj]]["tmp.type"][[1]], by = "time", all = TRUE)
    }

    tmp.long <- tmp.long[order(id,time)][!is.na(id)]

    tmp.long[, time.obs := nafill(time.obs, "nocb"), by = "id"]
    tmp.long[is.na(time.obs), time.obs := -Inf]

    tmp.long[, final.time := max(time.obs), by = "id"]

    for (fit.type.jj in 1:length(fit.types)) {
        tmp.long[is.na(get(paste0("dhazard.", names(fit.types)[fit.type.jj]))), (paste0("dhazard.", names(fit.types)[fit.type.jj])) := 0]
    }

    if (length(state.names)>0) {
        for (process.jj in 1:length(state.names[state.names %in% names(fit.types)])) {
            tmp.long[, (state.names[process.jj]) := cumsum(1*(delta == state.deltas[process.jj])), by = "id"]
        }
    }
    
    if (length(state.names)>0) {
        for (varname in state.names[state.names %in% names(fit.types)]) {
            tmp.long[, (varname) := c(0, get(varname)[-.N]), by = "id"]
        }
    }

    if (length(which.recurrent) > 0) {
        for (var.recurrent in which.recurrent) {
            for (count.jj in 1:max.count) {
                tmp.long[, (paste0(var.recurrent, ".", count.jj)) := (tmp.long[[var.recurrent]] >= count.jj)]
            }
        }
    }

    if (length(derived.vars)>0) {
        for (derived.jj in 1:length(derived.vars)) {
            derived.var <- names(derived.vars)[derived.jj]
            derived.fun <- derived.vars[[derived.jj]]
            tmp.long[, (derived.var) := derived.fun(tmp.long)]
        }
    }

    tmp.long[, tstart := c(0, time[-.N]), by = "id"]

    if (length(depend.time)>0) {
        for (time.jj in 1:length(depend.time)) {
            time.var <- names(depend.time)[time.jj]
            time.value <- as.numeric(strsplit(time.var, "\\.")[[1]][2])
            if (is.na(time.value)) time.value <- 1
            time.fun <- depend.time[[time.jj]]
            suppressWarnings(
                tmp.long[, (paste0("T.", time.var)) := min(time[c(get(time.var)[-1],get(time.var)[.N]) == time.value]), by = "id"]
            )
            tmp.long[tstart < get(paste0("T.", time.var)), (paste0("T.", time.var)) := 0] #time.fun(tstart)
            tmp.long[tstart >= get(paste0("T.", time.var)), (paste0("T.", time.var)) := as.numeric(time.fun(get(paste0("T.", time.var))))]
            tmp.long[, (paste0("time.", time.var)) := as.numeric(time.fun(time))]
        }
    }
    
    if (length(a) > 0) {
        tmp.long[["pi.A0.1"]] <- predict(fit.A0, newdata = tmp.long, type = "response")
    }

    #--------------------------------    
    #-- hal

    if (browse.hal) browser()
    
    if (any.hal) {

        deltas.hal <- as.numeric(sapply(fit.types[which.hal], function(fit.type) {
            strsplit(strsplit(gsub(" ", "", fit.type$model), "delta==")[[1]][2], ")")[[1]][1]
        }))
            
        tmp.long[, idN := 1:.N, by = "id"]

        N.vars <- state.names
        T.vars <- grep("T.", varnames, value = TRUE)

        for (N.var in N.vars[N.vars %in% process.names]) {
            N.var.count <- tmp.long[, (get(N.var))[.N], by = "id"]
            NK <- ceiling(quantile(N.var.count[[2]], probs = 1-reduce.NK))
            for (jnk in 1:NK) {
                T.vars <- unique(c(T.vars, paste0("T.", N.var, ".", jnk)))
            }
        }

        N.vars <- unique(c(N.vars, T.vars))
        
        baseline.vars <- setdiff(varnames, N.vars)

        ##T.vars2 <- unique(c(T.vars, unlist(lapply(hal.sl, function(hsl) hsl[grep("T.vars", names(hsl))]))))
        ##T.vars2 <- T.vars2[T.vars2 != ""]

        for (T.var in T.vars) {
            time.value <- as.numeric(strsplit(T.var, "\\.")[[1]][3])
            if (is.na(time.value)) time.value <- 1
            time.var <- strsplit(T.var, "\\.")[[1]][2]
            suppressWarnings(
                tmp.long[, (T.var) := min(time[c(get(time.var)[-1],get(time.var)[.N]) == time.value]),
                         by = "id"]
            )
            tmp.long[tstart < get(T.var), (T.var) := 0] 
        }

        tmp.long[, (varnames) := lapply(.SD, function(x) {
            if (is.character(x)) {
                return(as.numeric(as.factor(x)))
            } else if (!is.numeric(x)) {
                return(as.numeric(x))
            } else return(x)
        }), .SDcols = varnames]

        grid.times <- unique.times[floor(seq(1, length(unique.times), length=cut.time+2))[-c(1, cut.time+2)]]
        #grid.times <- grid.times[grid.times != 0.00001]

        tmp.long[, grid.period := as.numeric(cut(tstart, unique(c(0, grid.times, Inf)), include.lowest=TRUE, right=FALSE))]
        tmp.long[, grid.time := c(0, grid.times, Inf)[grid.period]]

        by.vars <- c("id", "grid.period")
        if (length(N.vars)>0) by.vars <- c(by.vars, N.vars)

        for (delta.value in deltas.hal) {
           
            name.hal <- names(which.hal)[deltas.hal == delta.value]
            is.name.hal.at.risk <- (1:length(fit.types))[names(fit.types) == name.hal] %in% at.risk.ids

            if (is.name.hal.at.risk) {
                tmp.long[, (paste0("at.risk.", delta.value)) :=
                               fit.types[[name.hal]][["at.risk"]](tmp.long)]
            } else {
                tmp.long[, (paste0("at.risk.", delta.value)) := (time<=final.time)]
            }
            
            if (delta.value == 0) {
                tmp.long[, (paste0("ddd.", delta.value)) := sum(get(paste0("at.risk.", delta.value)) & time==final.time & delta==delta.value), by = by.vars]
            } else {
                tmp.long[, (paste0("ddd.", delta.value)) := sum(get(paste0("at.risk.", delta.value)) & (time<=final.time)*(delta==delta.value)), by = by.vars]
            }

            #tmp.long[, (paste0("risk.time.", delta.value)) := 0]
            tmp.long[time <= final.time & get(paste0("at.risk.", delta.value)), (paste0("risk.time.", delta.value)) := sum(time - tstart), by = by.vars]
        }

        tmp.long.reduced <- unique(tmp.long[time<=final.time], by=by.vars)

        make.hal.formula <- function(cut.one.way = cut.one.way,
                                     cut.Tk = cut.Tk,
                                     cut.Tk.values = cut.Tk.values,
                                     two.way = two.way,
                                     N.vars = N.vars,
                                     T.vars = T.vars,
                                     max.count = max.count) {

            grid.baseline.vars <- lapply(baseline.vars, function(var) {
                sort(unique(tmp.long[idN == 1, var, with = FALSE][floor(seq(1, n, length=cut.one.way))])[[1]])
            })
            names(grid.baseline.vars) <- baseline.vars

            if (!all(T.vars == "")) {
                if (length(cut.Tk.values)>0 & length(cut.Tk.values[[1]])>0) {
                    grid.T.vars <- lapply(T.vars, function(T.var) {
                        cut.Tk.values
                    })
                } else {
                    grid.T.vars <- lapply(T.vars, function(T.var) {
                        tmp.var <- strsplit(gsub("T.", "", T.var), "\\.")[[1]]
                        N.var <- tmp.var[1]
                        N.count <- tmp.var[2]
                        unique.T.var <- sort(unique(
                            tmp.long[time <= tau & get(N.var) >= N.count][[T.var]][c(0,diff(tmp.long[time <= tau & get(N.var) >= N.count][[T.var]]))>0]))
                        return(unique.T.var[floor(seq(1, length(unique.T.var), length = min(cut.Tk, length(unique.T.var)/3)))])
                    })
                }
                names(grid.T.vars) <- T.vars
            } else {
                T.vars <- NULL
            }

            if (all(N.vars == "")) {
                N.vars <- NULL
            }
            
            N.vars <- unique(c(setdiff(N.vars, grep("T.", N.vars, value = TRUE)), T.vars))

            hal.formula.main <- paste0("delta ~ -1 + ", ifelse(length(A0.name)>0, paste0(A0.name," + "), ""),
                                       paste0("(grid.time >=", grid.times, ")", collapse = "+"),
                                       " + ", paste0(sapply(setdiff(baseline.vars, A0.name), function(var)
                                           paste0("(", var, ">=", grid.baseline.vars[[var]], ")", collapse = "+")), collapse = "+"),
                                       ifelse(length(setdiff(N.vars, T.vars))>0,
                                              paste0(" + ", paste0(sapply(setdiff(N.vars, T.vars), function(var)
                                                  paste0("(", var, ">=", 1:max.count, ")", collapse = "+")), collapse = "+")), ""),
                                       ifelse(length(T.vars)>0, paste0(" + ", paste0(sapply(T.vars, function(var)
                                           paste0("(", var, ">=", grid.T.vars[[var]], ")", collapse = "+")), collapse = "+")), ""))
           
            if (length(two.way) > 0 & length(two.way[[1]]) > 0) {
                two.way.matrix <- do.call("rbind", lapply(two.way, function(two.way.jj) {
                    two.way.inner <- unique(
                        do.call("rbind",
                                combn(unlist(lapply(seq_along(two.way.jj[1:2]), function(ii) {
                                    two.way.var <- two.way.jj[ii]
                                    k1 <- as.numeric(two.way.jj[3])
                                    k2 <- ifelse(length(two.way.jj) >= 4, as.numeric(two.way.jj[4]), k1)
                                    kk <- if (ii == 1) k1 else k2
                                    if (two.way.var == "time") {
                                        return(paste0("(grid.time >=",
                                                      grid.times[seq(1, length(grid.times), length = kk)],
                                                      ")"))
                                    } else if (two.way.var %in% baseline.vars) {
                                        return(paste0("(", two.way.var, ">=",
                                                      grid.baseline.vars[[two.way.var]][seq(1, length(grid.baseline.vars[[two.way.var]]),
                                                                                            length = min(length(grid.baseline.vars[[two.way.var]]), kk))],
                                                      ")"))
                                    } else if (two.way.var %in% setdiff(N.vars, T.vars)) {
                                        return(paste0("(", two.way.var, ">=", (1:max.count)[seq(1, max.count,
                                                                                                length = kk)],
                                                      ")"))
                                    } else if (two.way.var %in% T.vars) {
                                        return(paste0("(", two.way.var, ">=", grid.T.vars[[two.way.var]][seq(1, length(grid.T.vars[[two.way.var]]),
                                                                                                             length = min(length(grid.T.vars[[two.way.var]]), kk))], ")"))
                                    } else if (two.way.var == A0.name) {
                                        return(A0.name)
                                    }
                                })), 2, simplify = FALSE)))
                    two.way.inner <-
                        two.way.inner[
                            sapply(str_split(two.way.inner[,1], ">="), function(twi) gsub("\\(", "", twi[1])) !=
                            sapply(str_split(two.way.inner[,2], ">="), function(twi) gsub("\\(", "", twi[1])), ]
                    if (substr(two.way.jj[1], 1, 2) == "T." &
                        substr(two.way.jj[2], 1, 2) == "T.") {
                        two.way.inner <-
                            two.way.inner[
                                sapply(str_split(two.way.inner[,1], ">="), function(twi) as.numeric(gsub(")", "", twi[2]))) < 
                                sapply(str_split(two.way.inner[,2], ">="), function(twi) as.numeric(gsub(")", "", twi[2]))), ]
                    }
                    if (substr(two.way.jj[1], 1, 9) == "time" &
                        substr(two.way.jj[2], 1, 2) == "T.") {
                        two.way.inner <-
                            two.way.inner[
                                sapply(str_split(two.way.inner[,1], ">="), function(twi) as.numeric(gsub(")", "", twi[2]))) > 
                                sapply(str_split(two.way.inner[,2], ">="), function(twi) as.numeric(gsub(")", "", twi[2]))), ]
                    }
                    if (substr(two.way.jj[2], 1, 9) == "time" &
                        substr(two.way.jj[1], 1, 2) == "T.") {
                        two.way.inner <-
                            two.way.inner[
                                sapply(str_split(two.way.inner[,1], ">="), function(twi) as.numeric(gsub(")", "", twi[2]))) <
                                sapply(str_split(two.way.inner[,2], ">="), function(twi) as.numeric(gsub(")", "", twi[2]))), ]
                    }
                    return(two.way.inner[apply(two.way.inner, 1, function(two.way.inner.row) {
                        !((length(grep(two.way.jj[1], two.way.inner.row[1]))>0 &
                           length(grep(two.way.jj[1], two.way.inner.row[2]))>0) |
                          (length(grep(two.way.jj[2], two.way.inner.row[1]))>0 &
                           length(grep(two.way.jj[2], two.way.inner.row[2]))>0))}),])
                }))
                hal.formula.main <- paste0(hal.formula.main,
                                           " + ", paste0(two.way.matrix[,1], ":", two.way.matrix[,2], collapse = "+"))
            }
            return(hal.formula.main)
        }

        make.hal <- function(hal.jj) {

            if (sum(names(hal.jj) %in% c("cut.one.way", "max.count", "T.vars", "N.vars",
                                         "cut.Tk", "cut.Tk.values", "two.way"))>0) {
                if (length(grep("max.count", hal.jj))>0) {
                    hal.jj["max.count"] <- min(hal.jj["max.count"], NK)
                }
                return(hal.jj) 
            }

            if (sum(names(hal.jj) %in% c("components", "interactions")) == 0) {

                split.hal.jj <- strsplit(hal.jj, "_")
                if (length(split.hal.jj[[1]])>1) {
                    use.times <- gsub("time", "", split.hal.jj[[1]][2])
                } else {
                    use.times <- NULL
                }

                hal.jj <- split.hal.jj[[1]][1]
                
                if (hal.jj == "default") {
                    hal.jj <- list(
                        components = list(
                            baseline = list(complexity = "medium"),
                            history  = list(complexity = "medium", time.mode = "auto")
                        ),
                        interactions = list(
                            time.with = c("history", A0.name),
                            history.with = c("history"),
                            baseline.with = c("baseline", A0.name)
                        )
                    )
                } else if (hal.jj == "history.interactions") {
                    hal.jj <- list(
                        components = list(
                            baseline = list(complexity = "medium"),
                            history  = list(complexity = "medium", time.mode = "auto")
                        ),
                        interactions = list(
                            time.with = c("history"),
                            history.with = c("history")
                        )
                    )
                } else if (hal.jj == "dense") {
                    hal.jj <- list(
                        components = list(
                            baseline = list(complexity = "high"),
                            history  = list(complexity = "high", time.mode = "auto")
                        ),
                        interactions = list(
                            time.with = c("history", A0.name),
                            baseline.with = c("baseline", A0.name)
                        )
                    )
                } else if (hal.jj == "full2") {
                    hal.jj <- list(
                        components = list(
                            baseline = list(complexity = "high"),
                            history  = list(complexity = "high", time.mode = "auto")
                        ),
                        interactions = list(
                            time.with = c("history", A0.name, "baseline"),
                            history.with = c("history", "baseline", A0.name), 
                            baseline.with = c("baseline", A0.name)
                        )
                    )
                } else if (hal.jj == "full") {
                    hal.jj <- list(
                        components = list(
                            baseline = list(complexity = "high"),
                            history  = list(complexity = "high", time.mode = "auto")
                        ),
                        interactions = list(
                            time.with = c("history", A0.name),
                            history.with = c("history"), 
                            baseline.with = c("baseline", A0.name)
                        )
                    )
                } else if (hal.jj == "no.interactions") {
                    hal.jj <- list(
                        components = list(
                            baseline = list(complexity = "medium"),
                            history  = list(complexity = "medium", time.mode = "auto")
                        )
                    )
                } else if (hal.jj == "baseline") {
                    hal.jj <- list(
                        components = list(
                            baseline = list(complexity = "medium"),
                            history  = list(complexity = "none", time.mode = "none")
                        )
                    )
                } else if (hal.jj == "sparse") {
                    hal.jj <- list(
                        components = list(
                            baseline = list(complexity = "medium"),
                            history  = list(complexity = "medium", time.mode = "none")
                        )
                    )
                } else {
                    message("NB: the choice of hal not supported - sparse will be fitted")
                    hal.jj <- list(
                        components = list(
                            baseline = list(complexity = "medium"),
                            history  = list(complexity = "medium", time.mode = "none")
                        )
                    )
                }

                if (length(use.times)>0) {
                    hal.jj[["components"]][length(hal.jj[["components"]])+1] <-
                        list(list(use.times = use.times))
                    names(hal.jj[["components"]])[length(hal.jj[["components"]])] <-
                        "time"
                }
            }

            if (sum(names(hal.jj) %in% c("components", "interactions"))>0) {
                hal.out <- c()
                if (length(baseline <- hal.jj[["components"]][["baseline"]])>0) {
                    if (baseline$complexity == "medium") {
                        hal.out <- c(hal.out, cut.one.way = 20)
                    } else if (baseline$complexity == "high") {
                        hal.out <- c(hal.out, cut.one.way = 35)
                    } else {
                        hal.out <- c(hal.out, cut.one.way = 10)
                    }
                }
                if (length(history <- hal.jj[["components"]][["history"]])>0) {
                    if (history$complexity == "medium") {
                        hal.out <- c(hal.out, max.count = min(3, NK),
                                     cut.Tk = 5)
                    } else if (history$complexity == "high") {
                        hal.out <- c(hal.out, max.count = min(10, NK),
                                     cut.Tk = 10)
                    } else if (history$complexity == "low") {
                        hal.out <- c(hal.out, max.count = min(2, NK),
                                     cut.Tk = 2)
                    } else {
                        hal.out <- c(hal.out, max.count = 0,
                                     cut.Tk = 0)
                    }
                    if (!("time.mode" %in% names(history))) {
                        history$time.mode <- "auto"
                    }
                    if (history$time.mode == "auto") {
                        hal.out <- c(hal.out, T.vars = grep(paste0(".", 1:hal.out["max.count"], collapse = "|"),
                                                            T.vars, value = TRUE))                             
                    } else if (history$time.mode == "first" & hal.out["max.count"] >= 1) {
                        hal.out <- c(hal.out, T.vars = grep(grep("\\.1$", T.vars, value = TRUE),
                                                            T.vars, value = TRUE))
                    } else {
                        hal.out <- c(hal.out, T.vars = "")
                    }
                }
                if (TRUE) { ## this is test, not recommended
                    if (length(time <- hal.jj[["components"]][["time"]])>0) {
                        if ("complexity" %in% names(time)>0) {
                            if (time$complexity == "medium") {
                                hal.out <- c(hal.out, penalize.time = 0.4)
                            } else if (time$complexity == "high") {
                                hal.out <- c(hal.out, penalize.time = 0)
                            } else if (time$complexity == "none") {
                                hal.out <- c(hal.out, penalize.time = 1)
                            } else {
                                hal.out <- c(hal.out, penalize.time = 0.7)
                            }
                        }
                        if ("use.times" %in% names(time)>0) {
                            hal.out <- c(hal.out, use.times = time$use.times)
                        }
                    }
                }
                
                if (sum(names(hal.jj) %in% c("interactions"))>0) {

                    two.way.out <- list()

                    if (length(time.with <- hal.jj[["interactions"]][["time.with"]])>0) {
                        if ("history" %in% time.with) {
                            for (T.var in hal.out[grep("T.vars", names(hal.out))]) {
                                two.way.out[[length(two.way.out)+1]] <-
                                    c("time", T.var, 30, 10)
                            }
                            for (N.var in setdiff(N.vars, grep("T.", N.vars, value = TRUE))) {
                                two.way.out[[length(two.way.out)+1]] <-
                                    c("time", N.var, 10, 10)
                            }
                        }
                        if (length(A0.name)>0) {
                            if (A0.name %in% time.with) {
                                two.way.out[[length(two.way.out)+1]] <-
                                    c("time", A0.name, 10)
                            }
                        }
                        if ("baseline" %in% time.with) {
                            for (baseline.var in setdiff(baseline.vars, A0.name)) {
                                two.way.out[[length(two.way.out)+1]] <-
                                    c("time", baseline.var, 10)
                            }
                        }
                    }

                    ### fix: need to make this: 
                    if (length(history.with <- hal.jj[["interactions"]][["history.with"]])>0) {
                        if ("history" %in% history.with) {
                            if ((length.T.vars <- length(T.vars <- hal.out[grep("T.vars", names(hal.out))])) > 1) {
                                for (tt1 in 1:(length.T.vars-1)) {
                                    for (tt2 in (tt1+1):length(T.vars)) {
                                        two.way.out[[length(two.way.out)+1]] <-
                                            c(T.vars[tt1], T.vars[tt2], min(10, hal.out["cut.Tk"]))
                                    }
                                }
                            }
                        }
                        if (length(A0.name)>0) {
                            if (A0.name %in% history.with) {
                                for (T.var in hal.out[grep("T.vars", names(hal.out))]) {
                                    two.way.out[[length(two.way.out)+1]] <-
                                        c(A0.name, T.var, 10)
                                }
                                for (N.var in setdiff(N.vars, grep("T.", N.vars, value = TRUE))) {
                                    two.way.out[[length(two.way.out)+1]] <-
                                        c(A0.name, N.var, 10)
                                }
                            }
                        }
                        if ("baseline" %in% history.with) {
                            for (baseline.var in setdiff(baseline.vars, A0.name)) {
                                for (T.var in hal.out[grep("T.vars", names(hal.out))]) {
                                    two.way.out[[length(two.way.out)+1]] <-
                                        c(baseline.var, T.var, 10)
                                }
                                for (N.var in setdiff(N.vars, grep("T.", N.vars, value = TRUE))) {
                                    two.way.out[[length(two.way.out)+1]] <-
                                        c(baseline.var, N.var, 10)
                                }
                            }
                        }
                    }

                    if (length(baseline.with <- hal.jj[["interactions"]][["baseline.with"]])>0) {
                        if ("baseline" %in% baseline.with) {
                            if ((length.baseline <- length(baseline.covars <- setdiff(baseline.vars, A0.name)))>1) {
                                for (baseline.j1 in 1:(length.baseline-1)) {
                                    for (baseline.j2 in (baseline.j1+1):length(baseline.covars)) {
                                        two.way.out[[length(two.way.out)+1]] <-
                                            c(baseline.covars[baseline.j1], baseline.covars[baseline.j2], 10)
                                    }
                                }
                            }
                        }
                        if (length(A0.name)>0) {
                            if (A0.name %in% baseline.with) {
                                for (baseline.j1 in 1:(length(baseline.covars <- setdiff(baseline.vars, A0.name)))) {
                                    two.way.out[[length(two.way.out)+1]] <-
                                        c(baseline.covars[baseline.j1],  A0.name, 10)
                                }
                            }
                        }
                    }
                } else {

                    two.way.out <- list()

                }

                hal.out <- c(hal.out, two.way = unique(two.way.out))
            }

            return(hal.out)
        }
        
        ##if (length(hal.sl)>1) browser()

        ## FIX! need to make this work
        ### browser()
        fit.hals.list <- lapply(1:length(hal.sl), function(hal.jj) {

            hal.parameters <- make.hal(hal.jj = hal.sl[[hal.jj]])

            if (length(grep("penalize.time", names(hal.parameters)))>0) {
                penalize.time <- as.numeric(hal.parameters[grep("penalize.time", names(hal.parameters))])
            }

            if (length(grep("use.times", names(hal.parameters)))>0) {
                use.times <- as.numeric(hal.parameters[grep("use.times", names(hal.parameters))])
            } else {
                use.times <- 1
            }

            if (verbose) print(paste0(hal.jj, "/", length(hal.sl), " (hal-sl)"))

            hal.formula.main <- make.hal.formula(cut.one.way = ifelse(
                                                     "cut.one.way" %in% names(hal.parameters),
                                                     as.numeric(hal.parameters[["cut.one.way"]]),
                                                     cut.one.way),
                                                 cut.Tk = if (
                                                     length(grep("cut.Tk", names(hal.parameters)))>0) {
                                                              as.numeric(hal.parameters[grep("cut.Tk", names(hal.parameters))])
                                                          } else cut.Tk,
                                                 cut.Tk.values = if (
                                                     length(grep("cut.Tk.values", names(hal.parameters)))>0) {
                                                                     hal.parameters[grep("cut.Tk.values", names(hal.parameters))]
                                                                 } else cut.Tk.values,
                                                 two.way = if (
                                                     length(which.two.way <- grep("two.way", names(hal.parameters)))>0) {
                                                               lapply(which.two.way, function(tw) hal.parameters[[tw]])
                                                           } else two.way,
                                                 N.vars = if (
                                                     length(grep("N.vars", names(hal.parameters)))>0) {
                                                              hal.parameters[grep("N.vars", names(hal.parameters))]
                                                          } else N.vars,
                                                 T.vars = if (
                                                     length(grep("T.vars", names(hal.parameters)))>0) {
                                                              hal.parameters[grep("T.vars", names(hal.parameters))]
                                                          } else T.vars,
                                                 max.count = ifelse(
                                                     "max.count" %in% names(hal.parameters),
                                                     as.numeric(hal.parameters[["max.count"]]),
                                                     max.count))
        
            X <- Matrix(
                model.matrix(formula(hal.formula.main),
                             data = tmp.long.reduced), sparse = TRUE)

            X <- X[, !(colnames(X) %in% grep("FALSE", colnames(X), value = TRUE))]

            ##col_ones <- Matrix::colMeans(X != 0)
            ##X <- X[, col_ones >= min.no.of.ones]
        
            x.vector <- hash_sparse_rows_dgC(X)
            ### browser()
            tmp.long.reduced[, x:=x.vector]

            if (FALSE) {
                
                mf <- model.frame(
                    formula(hal.formula.main),
                    data = tmp.long.reduced,
                    na.action = na.omit
                )

                dim(mf)

                attr(mf, "na.action")
                length(attr(mf, "na.action"))

                mf <-  model.frame(
                    formula(hal.formula.main),
                    data = tmp.long.reduced,
                    na.action = na.pass
                )

                dim(mf)


                
                vars <- all.vars(formula(hal.formula.main))
                vars
                colSums(is.na(tmp.long.reduced[, vars, drop = FALSE]))

            }

            if (length(seed.hal) == 0) seed.hal <- sample(1e8, 1)

            fit.hals <- lapply(deltas.hal, function(delta.value) {

                ##message("-------")
                if (verbose) message(paste0("delta = ", delta.value))
                ##message("-------")

                ##name.hal <- names(which.hal)[deltas.hal == delta.value]
                ##is.name.hal.at.risk <- (1:length(fit.types))[names(fit.types) == name.hal] %in% at.risk.ids

                tmp.long.reduced[, risk.time := get(paste0("risk.time.", delta.value))]

                first.hal <- fit.hal(
                    hal.pseudo.dt = tmp.long.reduced[tmp.long.reduced[[paste0("at.risk.", delta.value)]]],
                    X.hal = X[tmp.long.reduced[[paste0("at.risk.", delta.value)]],], 
                    delta.var = "delta",
                    delta.value = delta.value,
                    time.var = "time",
                    treatment = A0.name,
                    lambda.cvs = lambda.cvs,
                    penalize.time = if (penalize.time<1) {penalize.time} else {screen.two.way|penalize.time},#FALSE,
                    penalize.treatment = screen.two.way,#FALSE,
                    npenalize.vars = npenalize.vars,
                    event.dependent.cv = event.dependent.cv,
                    reduce.seed.dependence = reduce.seed.dependence,
                    return.train.fit = FALSE,
                    use.times = use.times,
                    check.support.for.basis.functions = 0.01,
                    V = V,
                    parallelize.cve = min(detectCores()-1, use.cores),
                    cv.glmnet = cv.glmnet,
                    verbose = verbose.hal,
                    return.cve = TRUE,
                    seed = seed.hal)

                if (screen.two.way) {
                    get.coef <- coef(first.hal$hal.fit, s = first.hal$lambda)
                    coef.dt <- data.table(non.zero = get.coef@Dimnames[[1]][get.coef@i+1],
                                          coef = get.coef@x)#[non.zero %in% c("Aobs", "Y.dummy >= 1TRUE")]
                    first.hal.coefs <- coef.dt[["coef"]]
                    which.two.way <- names(first.hal.coefs) <- gsub("TRUE|FALSE", "", paste0(coef.dt[["non.zero"]]))
                    two.way <- do.call(
                        "rbind",
                        combn(which.two.way[which.two.way != "(Intercept)"], 2, simplify = FALSE))
                    two.way <- two.way[apply(two.way, 1, function(two.way.row) {
                        (strsplit(two.way.row[1], ">=")[[1]][1] != strsplit(two.way.row[2], ">=")[[1]][1])
                    }),]
                    hal.formula.main.2 <-
                        paste0(hal.formula.main,
                               " + ", paste0("(", two.way[,1], "):(", two.way[,2], ")", collapse = "+"))
                    X.2 <- Matrix(
                        model.matrix(formula(hal.formula.main.2),
                                     data = tmp.long.reduced), sparse = TRUE)
                    return(second.hal <- fit.hal(
                               hal.pseudo.dt = tmp.long.reduced[tmp.long.reduced[[paste0("at.risk.", delta.value)]]], 
                               X.hal = X.2[tmp.long.reduced[[paste0("at.risk.", delta.value)]]], 
                               delta.var = "delta",
                               delta.value = delta.value,
                               time.var = "time",
                               treatment = A0.name,
                               lambda.cvs = lambda.cvs,
                               penalize.time = FALSE|penalize.time,
                               penalize.treatment = FALSE,
                               event.dependent.cv = event.dependent.cv,
                               reduce.seed.dependence = reduce.seed.dependence,
                               return.train.fit = FALSE,
                               use.times = use.times,
                               check.support.for.basis.functions = TRUE,
                               V = V,
                               parallelize.cve = min(detectCores()-1, use.cores),
                               cv.glmnet = cv.glmnet,
                               verbose = verbose.hal,
                               return.cve = TRUE,
                               seed = seed.hal))
                } else {
                    return(first.hal)
                }
            })

        })

        if (length(fit.hals.list) == 1) {
            fit.hals <- fit.hals.list[[1]]
            if ((!verbose.hal) & verbose) {
                lapply(fit.hals, function(fh) {
                    message("--------------------------------------------")
                    message(paste0("delta = ", fh$delta.value))
                    message("--------------------------------------------")
                    print(coef(fh$hal.fit, s=fh$lambda.cv))
                })
            }
        } else if (length(fit.hals.list) > 1) {
            fit.hals <- lapply(1:length(fit.hals.list[[1]]), function(fh.process.index) {
                if (verbose) {
                    message("--------------------------------------------")
                    message(paste0("delta = ", fit.hals.list[[1]][[fh.process.index]]$delta.value))
                    message("--------------------------------------------")
                }
                # here:
                inner.cve <- sapply(fit.hals.list, function(fh) fh[[fh.process.index]]$cve)
                picked <- (1:length(inner.cve))[inner.cve == min(inner.cve)]
                if (verbose) print(inner.cve)
                if (verbose) message(paste0("picked hal: ", picked))
                fh <- fit.hals.list[[picked]][[fh.process.index]]
                fh.coef <- try(coef(fh$hal.fit, s=fh$lambda.cv))
                if (FALSE & inherits( fh.coef, "try-error")) {
                    inner.cve <- inner.cve[-picked]
                    picked <- (1:length(inner.cve))[inner.cve == min(inner.cve)]
                    fh <- fit.hals.list[[picked]][[fh.process.index]]
                    fh.coef <- try(coef(fh$hal.fit, s=fh$lambda.cv))
                    if (inherits( fh.coef, "try-error")) {
                        inner.cve <- inner.cve[-picked]
                        picked <- (1:length(inner.cve))[inner.cve == min(inner.cve)]
                        fh <- fit.hals.list[[picked]][[fh.process.index]]
                        fh.coef <- try(coef(fh$hal.fit, s=fh$lambda.cv))
                    }
                }
                if (verbose) {
                    print(coef(fh$hal.fit, s=fh$lambda.cv))
                }
                return(fh)
            })
        }

        tmp.long[, risk.time := time-tstart]

        ##browser()

        hal.coefs <- lapply(fit.hals, function(fh) {
            delta <- fh$delta.value
            get.coef <- coef(fh$hal.fit, s = fh$lambda)
            coef.dt <- data.table(non.zero = get.coef@Dimnames[[1]][get.coef@i+1],
                                  coef = get.coef@x)#[non.zero %in% c("Aobs", "Y.dummy >= 1TRUE")]
            coefs <- coef.dt[["coef"]]
            names(coefs) <- paste0(coef.dt[["non.zero"]])
            return(coefs)
        })

        non.zero.vars <- unique(unlist(sapply(hal.coefs[deltas.hal != 0],
                                              function(hal.coef) names(hal.coef))))

        non.zero.state.vars <- lapply(N.vars, function(N.var) {
            tmp.vars <- gsub("TRUE", "", grep(N.var, non.zero.vars, value = TRUE))
            if (length((tmp.interaction <- grep(":", tmp.vars, value = TRUE)) > 0)) {
                tmp.interaction.vars <- unlist(lapply(strsplit(tmp.interaction, ":"), function(interaction.split)
                    grep(N.var, interaction.split, value = TRUE)))
                tmp.vars <- unique(c(setdiff(tmp.vars, tmp.interaction), tmp.interaction.vars))
            }
            return(sort(unique(tmp.vars[!tmp.vars %in% grep(paste0("T.", N.var), tmp.vars, value = TRUE)])))
        })
        names(non.zero.state.vars) <- N.vars
        ## the following is to make sure that N.x >= y is always included when T.x.y is
        for (T.var in T.vars) {
            if (length(non.zero.state.vars[[T.var]])>0) {
                time.value <- as.numeric(strsplit(T.var, "\\.")[[1]][3])
                if (is.na(time.value)) time.value <- 1
                time.var <- strsplit(T.var, "\\.")[[1]][2]
                if (all(non.zero.state.vars[[time.var]] != paste0(time.var, " >= ", time.value))) {
                    non.zero.state.vars[[time.var]] <-
                        sort(c(non.zero.state.vars[[time.var]], paste0(time.var, " >= ", time.value)))
                }
            }
        }
        ##message("FIX: CHECK this below")
        non.zero.state.vars <- non.zero.state.vars[(1:length(non.zero.state.vars))[sapply(non.zero.state.vars, function(non.zero) {
            length(non.zero)>0
        })]]
    }
    
    #--------------------------------    
    #-- training part is over
    
    tmp.long <- tmp.long[time <= tau]

    ### browser()


    if (any.hal) {

        hal.vars.list <- lapply(1:length(fit.hals), function(kk) {
            test.get.beta <- try(out.beta <- fit.hals[[kk]][["hal.fit"]]$beta@Dimnames[[1]],
                                 silent = TRUE)
            if (any(class(test.get.beta) == "try-error")) {
                out.beta <- fit.hals[[kk]][["hal.fit"]]$glmnet$beta@Dimnames[[1]]
            }
            return(out.beta)
        })

        if (length(derived.vars)>0) {
            message("FIX: still need to check that HAL handles derived.var correctly")
        }
        
        hal.vars <- unique(c(unlist(hal.vars.list)))
        
        if (length(A0.name)>0 & length(interaction.vars <- grep(":", hal.vars, value = TRUE))>0) {
            hal.vars.A0.dynamic <- unlist(sapply(N.vars, function(N.var) {
                tmp.vars <- #gsub("TRUE", "", grep(N.var, unlist(hal.vars.list), value = TRUE))
                    grep(N.var, unlist(hal.vars.list), value = TRUE)
                tmp.vars <- tmp.vars[tmp.vars %in% grep(A0.name, unlist(hal.vars.list), value = TRUE)]
                return(unique(tmp.vars[!tmp.vars %in% grep(paste0("T.", N.var), tmp.vars, value = TRUE)]))
            }))
        } else {
            hal.vars.A0.dynamic <- NULL
        }
        
        ## browser()
        
        hal.vars.dynamic <- setdiff(unlist(sapply(N.vars, function(N.var) {
            tmp.vars <- #gsub("TRUE", "", grep(N.var, unlist(hal.vars.list), value = TRUE))
                grep(N.var, unlist(hal.vars.list), value = TRUE)
            return(unique(tmp.vars[!tmp.vars %in% grep(paste0("T.", N.var), tmp.vars, value = TRUE)]))
        })), hal.vars.A0.dynamic)
        
        if (length(A0.name)>0) {
            hal.vars.A0 <- setdiff(
                grep(A0.name, hal.vars[!(hal.vars %in% grep("FALSE", hal.vars, value = TRUE))], value = TRUE),
                hal.vars.A0.dynamic)
        } else {
            hal.vars.A0 <- NULL
        }
        
        hal.vars.static <- setdiff(hal.vars, c(hal.vars.dynamic, hal.vars.A0, hal.vars.A0.dynamic,
                                               unique(unlist(sapply(N.vars, function(N.var) grep(N.var, unlist(hal.vars.list), value = TRUE))))))

        if (FALSE) { ## may want to change to TRUE
            flip.interactions <- function(X.mat, var.vec) {
                if (length(cols.diff <- setdiff(colnames(X.mat), var.vec)) > 0) {
                    if (length(cols.diff <- setdiff(gsub("TRUE|FALSE", "", cols.diff),
                                                    var.vec)) > 0) {
                        colnames(X.mat)[gsub("TRUE|FALSE", "", colnames(X.mat)) %in% cols.diff] <-
                            sapply(colnames(X.mat)[gsub("TRUE|FALSE", "", colnames(X.mat)) %in% cols.diff], function(x) {
                                if (!grepl(":", x)) return(x)
                                parts <- strsplit(x, ":", fixed = TRUE)[[1]]
                                paste(rev(parts), collapse = ":")
                            })
                    }
                }
                return(X.mat)
            }
        } else {
            flip.interactions <- function(X.mat, var.vec) {
                if (length(cols.diff <- setdiff(var.vec, colnames(X.mat))) > 0) {
                    cols.diff.rev <- sapply(cols.diff, function(x) {
                        if (!grepl(":", x)) return(x)
                        parts <- strsplit(x, ":", fixed = TRUE)[[1]]
                        paste(rev(parts), collapse = ":")
                    })
                    cols.diff.TF <- grep("TRUE|FALSE", cols.diff, value = TRUE)
                    if (length(cols.diff.TF)>0) {
                        idx <- match(cols.diff.rev, colnames(X.mat))
                        valid <- !is.na(idx)
                        colnames(X.mat)[idx[valid]] <- cols.diff[valid]
                    }
                    if (length(cols.diff.TF)<length(cols.diff)) {
                        idx <- match(cols.diff.rev, gsub("TRUE|FALSE", "", colnames(X.mat)))
                        valid <- !is.na(idx)
                        colnames(X.mat)[idx[valid]] <- cols.diff[valid]
                    }
                }
                return(X.mat)
            }
        }
        
        X.hal.static <- Matrix(
            model.matrix(formula(paste0("delta", "~-1+",
                                        paste(paste0("(", gsub("FALSE|TRUE", "", gsub(":", "):(", hal.vars.static)), ")"),
                                              collapse = "+"))), 
                         data=tmp.long), sparse=TRUE)
        
        X.hal.static <- flip.interactions(X.hal.static, hal.vars.static)
        
        if (length(A0.name)>0) {
            X.hal.A0 <- Matrix(
                model.matrix(formula(paste0("delta", "~-1+",
                                            paste(paste0("(", gsub("FALSE|TRUE", "", gsub(":", "):(", hal.vars.A0)), ")"),
                                                  collapse = "+"))), 
                             data=tmp.long), sparse=TRUE)
            X.hal.A0 <- flip.interactions(X.hal.A0, hal.vars.A0)
        } else {
            X.hal.A0 <- Matrix(0, nrow = nrow(tmp.long), ncol = 0, sparse = TRUE)
        }
        if (length(hal.vars.A0.dynamic)>0) {
            X.hal.A0.dynamic <- Matrix(
                model.matrix(formula(paste0("delta", "~-1+",
                                            paste(paste0("(", gsub("FALSE|TRUE", "", gsub(":", "):(", hal.vars.A0.dynamic)), ")"),
                                                  collapse = "+"))), 
                             data=tmp.long), sparse=TRUE)
            X.hal.A0.dynamic <- flip.interactions(X.hal.A0.dynamic, hal.vars.A0.dynamic)
        } else {
            X.hal.A0.dynamic <- Matrix(0, nrow = nrow(tmp.long), ncol = 0, sparse = TRUE)
        }
        
        X.hal.dynamic <- Matrix(
            model.matrix(formula(paste0("delta", "~-1+",
                                        paste(paste0("(", gsub("FALSE|TRUE", "",
                                                               gsub(":", "):(", hal.vars.dynamic)), ")"),
                                              collapse = "+"))), 
                         data=tmp.long), sparse=TRUE)
        X.hal.dynamic <- flip.interactions(X.hal.dynamic, hal.vars.dynamic)
        
        ###browser()
    }
    
    #--------------------------------
    #-- compute needed quantities for (initial) estimation and targeting

    for (fit.type.jj in 1:length(fit.types)) {
        if (fit.type.jj %in% at.risk.ids) {
            tmp.long[, (paste0("at.risk.", names(fit.types)[fit.type.jj])) :=
                           fit.types[[fit.type.jj]][["at.risk"]](tmp.long)]
        }
        if (fit.types[[fit.type.jj]]$fit == "hal") {
            hal.index <- (1:length(which.hal))[names(which.hal) == names(fit.types)[fit.type.jj]]
            ##delta.value <- fit.hals[[hal.index]][["delta.value"]]
            #### HERE (below), does not work with interaction?
            tmp.long[, (paste0("P.", names(fit.types)[fit.type.jj])) :=
                           exp(predict(fit.hals[[hal.index]][["hal.fit"]], cbind(X.hal.static, X.hal.dynamic, X.hal.A0.dynamic, X.hal.A0)[, hal.vars.list[[hal.index]]],
                                       newoffset=0, s=fit.hals[[hal.index]][["lambda.cv"]]))*risk.time]
        } else {        
            tmp.long[, (paste0("exp.", names(fit.types)[fit.type.jj])) := exp(predict(fit.cox.types[[fit.type.jj]]["fit.cox"][[1]], newdata=tmp.long, type="lp"))]
            tmp.long[, (paste0("P.", names(fit.types)[fit.type.jj])) :=
                           tmp.long[[paste0("dhazard.", names(fit.types)[fit.type.jj])]]*tmp.long[[paste0("exp.", names(fit.types)[fit.type.jj])]]]
        } 
        if (fit.type.jj == cens.process.id) {
            tmp.long[, surv.0 := exp(-cumsum(get(paste0("P.", names(fit.types)[fit.type.jj])))), by = "id"]
            tmp.long[, surv.0.1 := c(1, surv.0[-.N]), by = "id"]
        }
    }

    #--------------------------------    
    #-- to handle dependence on jumps in the past:

    if (browse) browser()

    if (length(state.names)>0) {
        grid.list <- lapply(state.names,
                            function(varname) 0:min(ifelse(varname %in% process.names, max.count, Inf),
                                                    max(unique(tmp.long[[varname]]))))
        depend.matrix <- as.data.table(do.call(expand.grid, grid.list))
        #names(depend.matrix) <- state.names
        setnames(depend.matrix, state.names)
        if (any.hal) {
            if (length(non.zero.state.vars)>0) {
                grid.list <- lapply(1:length(non.zero.state.vars), function(non.zero.jj) {
                    if ((non.zero.var <- names(non.zero.state.vars)[non.zero.jj]) %in% names(derived.vars)) {
                        return(0)#(unique(tmp.long[[non.zero.var]]))
                    } else if ((non.zero.var <- names(non.zero.state.vars)[non.zero.jj]) %in% setdiff(N.vars, T.vars)) {
                        return(0:max(as.numeric(gsub(paste0(non.zero.var, ">="), "", gsub(" ", "", non.zero.state.vars[[non.zero.jj]])))))
                    } else {
                        return(c(0, sapply(1:length(non.zero.state.vars[[non.zero.var]]), function(sss) {
                            return(as.numeric(gsub(paste0(paste0(non.zero.var, ">=")), "",
                                                   gsub(" ", "", non.zero.state.vars[[non.zero.var]][sss]))))
                        })))
                    }
                })
                depend.matrix <- as.data.table(do.call(expand.grid, grid.list))
                setnames(depend.matrix, names(non.zero.state.vars))
                #names(depend.matrix) <- names(non.zero.state.vars)
                if (length(T.vars <- grep("T.", names(non.zero.state.vars), value = TRUE))>0) {
                    for (T.var in T.vars) {
                        time.var <- strsplit(T.var, "\\.")[[1]][2]
                        time.value <- as.numeric(strsplit(T.var, "\\.")[[1]][3])
                        if (is.na(time.value)) time.value <- 1
                        depend.matrix <- depend.matrix[!(depend.matrix[[time.var]] <= time.value-1 & depend.matrix[[T.var]] > 0)]
                    }
                }

                if (FALSE & length(derived.vars)>0) {
                    for (derived.jj in 1:length(derived.vars)) {
                        derived.var <- names(derived.vars)[derived.jj]
                        derived.fun <- derived.vars[[derived.jj]]
                        if (derived.var %in% names(depend.matrix)) {
                            depend.matrix <- depend.matrix[, names(depend.matrix) != derived.var, with = FALSE]
                        }
                        print(colnames(depend.matrix))
                        depend.matrix[, derived.var := derived.fun(depend.matrix)]
                        if (any(sapply(derived.cols <- extract_dt_cols(derived.fun),
                                       function(varname) max.count < max(unique(tmp.long[[varname]]))))) {
                            depend.matrix.max <- depend.matrix[
                                apply(sapply(derived.cols, function(derived.col) depend.matrix[[derived.col]] == max(depend.matrix[[derived.col]])),
                                      1, all)]
                            unique.vals.max <- unique(tmp.long[apply(sapply(derived.cols, function(derived.col) tmp.long[[derived.col]] >= max(depend.matrix[[derived.col]])), 1, all)][["derived.var"]])
                            for (vals.max in setdiff(unique.vals.max, depend.matrix.max[["derived.var"]])) {
                                depend.matrix <-
                                    rbind(depend.matrix,
                                          depend.matrix[
                                              apply(sapply(derived.cols, function(derived.col) depend.matrix[[derived.col]] == max(depend.matrix[[derived.col]])),
                                                    1, all)][, derived.var := vals.max]) # <- (this may need be fixed for more general derived.vars)
                            }
                        }
                        vals <- unique(tmp.long[[derived.var]])
                        depend.matrix <- depend.matrix[derived.var %in% vals]
                        setnames(depend.matrix, "derived.var", derived.var)
                        print(colnames(depend.matrix))
                    }
                }
                for (non.zero.jj in (1:length(non.zero.state.vars))[names(non.zero.state.vars) %in% T.vars]) {
                    non.zero.var <- names(non.zero.state.vars)[non.zero.jj]
                    non.zero.var.states <- non.zero.state.vars[[non.zero.var]]
                    tmp.long[eval(parse(text = paste0("!(", non.zero.var.states[1], ")"))),
                    (paste0(non.zero.var, "_1")) := 0]
                    tmp.long[eval(parse(text = paste0("!(", gsub(non.zero.var, "time", non.zero.var.states[1]), ")"))),
                    (gsub("T.", "time.", non.zero.var)) := 0]
                    for (each.jj in 1:length(non.zero.var.states)) {
                        value.jj <- as.numeric(gsub(paste0(paste0(non.zero.var, ">=")), "",
                                                    gsub(" ", "", non.zero.state.vars[[non.zero.var]][each.jj])))
                        tmp.long[
                            eval(parse(text = non.zero.var.states[each.jj])),
                            (paste0(non.zero.var, "_1")) := value.jj]
                        tmp.long[
                            eval(parse(text = gsub(non.zero.var, "time", non.zero.var.states[each.jj]))),
                            (gsub("T.", "time.", non.zero.var)) := value.jj]
                    }
                    #tmp.long[[paste0(non.zero.var)]] <- tmp.long[[paste0(non.zero.var, "_1")]]
                    #tmp.long[[paste0(non.zero.var, "_1")]] <- NULL
                    tmp.long[, (non.zero.var) := get(paste0(non.zero.var, "_1"))]
                    tmp.long[, (paste0(non.zero.var, "_1")) := NULL]
                }
            }
        } else if (length(depend.time)>0) {
            for (time.jj in 1:length(depend.time)) {
                time.var <- names(depend.time)[time.jj]
                time.value <- as.numeric(strsplit(time.var, "\\.")[[1]][2])
                if (is.na(time.value)) time.value <- 1
                time.fun <- depend.time[[time.jj]]
                ##tmp.long[, (paste0("T.", time.var)) := min(time[c(get(time.var)[-1],get(time.var)[.N]) == 1]), by = "id"]
                time.var0 <- strsplit(time.var, "\\.")[[1]][1]
                depend.matrix <-
                    depend.matrix[!(get(time.var0) <= time.value-1 & get(paste0("T.", time.var)) >= 1)]
                ## check -> tmp.long[tstart >= get(paste0("T.", time.var)), (paste0("T.", time.var)) := time.fun(get(paste0("T.", time.var)))]
            }
        }
        if (length(derived.vars)>0) {
            for (derived.jj in 1:length(derived.vars)) {
                derived.var <- names(derived.vars)[derived.jj]
                derived.fun <- derived.vars[[derived.jj]]
                if (derived.var %in% names(depend.matrix)) {
                    depend.matrix <- depend.matrix[, names(depend.matrix) != derived.var, with = FALSE]
                }
                depend.matrix[, derived.var := derived.fun(depend.matrix)]
                if (any(sapply(derived.cols <- extract_dt_cols(derived.fun), function(varname) max.count < max(unique(tmp.long[[varname]]))))) {
                    depend.matrix.max <- depend.matrix[
                        apply(sapply(derived.cols, function(derived.col) depend.matrix[[derived.col]] == max(depend.matrix[[derived.col]])),
                              1, all)]
                    unique.vals.max <- unique(tmp.long[apply(sapply(derived.cols, function(derived.col) tmp.long[[derived.col]] >= max(depend.matrix[[derived.col]])), 1, all)][[derived.var]])
                    for (vals.max in setdiff(unique.vals.max, depend.matrix.max[["derived.var"]])) {
                        depend.matrix <-
                            rbind(depend.matrix,
                                  depend.matrix[
                                      apply(sapply(derived.cols, function(derived.col) depend.matrix[[derived.col]] == max(depend.matrix[[derived.col]])),
                                            1, all)][, derived.var := vals.max]) # <- (this may need be fixed for more general derived.vars)
                    }
                    ## setnames(depend.matrix, "derived.var", paste0(derived.var, " = ", paste0(derived.cols, collapse = ",")))
                    ## tmp.long[, At := At(tmp.long)]
                    ## setnames(tmp.long, derived.var, paste0("At = ", paste0(At.cols, collapse = ",")))
                }
                vals <- unique(tmp.long[[derived.var]])
                depend.matrix <- depend.matrix[derived.var %in% vals]
                setnames(depend.matrix, "derived.var", derived.var)
            }
        }
        depend.matrix[, state := 1:.N]
    } else {
        depend.matrix <- data.table(state = 1L)
    }

    if (prune.states) {
        ## this was added to handle the case that some jumps are not admissible/possible
        for (fit.type.jj in 1:length(fit.types)) {
            if (fit.type.jj %in% at.risk.ids) {
                depend.matrix[, (paste0("at.risk.", names(fit.types)[fit.type.jj])) :=
                                    fit.types[[fit.type.jj]][["at.risk"]](depend.matrix) |
                                                                         fit.types[[fit.type.jj]][["at.risk"]](depend.matrix-1)]
                if (all(depend.matrix[[paste0("at.risk.", names(fit.types)[fit.type.jj])]])) {
                    depend.matrix[[paste0("at.risk.", names(fit.types)[fit.type.jj])]] <- NULL
                }
            }
        }

        at.risk.names <- grep("at.risk.", names(depend.matrix), value = TRUE)

        if (length(at.risk.names)>0) {
            depend.matrix <- depend.matrix[rowSums(depend.matrix[, at.risk.names, with = FALSE])>0]
            depend.matrix <- depend.matrix[, !(names(depend.matrix) %in% at.risk.names), with = FALSE]
        }

        depend.matrix <- depend.matrix[order(state)]
        depend.matrix[, state := 1:.N]

    }

    if (length(a) == 0) { 
        intervene.A0 <- FALSE
        a <- 1
    } else {
        intervene.A0 <- TRUE
        setnames(tmp.long, A0.name, "A0.obs")
        ##tmp.long[["pi.A0.1"]] <- predict(fit.A0, newdata = tmp.long, type = "response")
    }

    #--------------------------------    
    #-- initial steps to compute clever weights:

    tmp.long[, C := cumsum(1*(time == time.obs & delta == 0 & time %in% dt[["time"]])), by = "id"]
    tmp.long[, C.1 := c(0, C[-.N]), by = "id"]
    tmp.long[, clever.weight := (C.1 == 0)/surv.0.1]
    
    #--------------------------------    
    #-- predict for different values of state and intervention on baseline treatment:
    
    tmp.long[, state := 0]
   
    for (aa in a) {

        if (intervene.A0) {
            tmp.long[, (A0.name) := aa]
            tmp.long[, (paste0("clever.weight.a", aa)) := clever.weight*
                           (A0.obs == aa)/((pi.A0.1^(A0.obs == 1)*(1-pi.A0.1)^(A0.obs == 0)))]

            if (any.hal) {
                X.hal.A0 <- Matrix(
                    model.matrix(formula(paste0("delta", "~-1+",
                                                paste(paste0("(", gsub("FALSE|TRUE", "", gsub(":", "):(", hal.vars.A0)), ")"),
                                                      collapse = "+"))), 
                                 data=tmp.long), sparse=TRUE)
                X.hal.A0 <- flip.interactions(X.hal.A0, hal.vars.A0)

                if (length(hal.vars.A0.dynamic)>0) {
                    X.hal.A0.dynamic <- Matrix(
                        model.matrix(formula(paste0("delta", "~-1+",
                                                    paste(paste0("(", gsub("FALSE|TRUE", "", gsub(":", "):(", hal.vars.A0.dynamic)), ")"),
                                                          collapse = "+"))), 
                                     data=tmp.long), sparse=TRUE)
                    X.hal.A0.dynamic <- flip.interactions(X.hal.A0.dynamic, hal.vars.A0.dynamic)
                }
                
                for (fit.type.jj in (1:length(fit.types))[-cens.process.id]) {
                    if (fit.types[[fit.type.jj]]$fit == "hal") {
                        hal.index <- (1:length(which.hal))[names(which.hal) == names(fit.types)[fit.type.jj]]
                        tmp.long[, (paste0("P.a", aa, ".", names(fit.types)[fit.type.jj])) :=
                                       exp(predict(fit.hals[[hal.index]][["hal.fit"]], cbind(X.hal.static, X.hal.dynamic, X.hal.A0.dynamic, X.hal.A0)[, hal.vars.list[[hal.index]]],
                                                   newoffset=0, s=fit.hals[[hal.index]][["lambda.cv"]]))*risk.time]
                    }
                }
            } else {        
                for (fit.type.jj in (1:length(fit.types))[-cens.process.id]) {
                    tmp.long[, (paste0("exp.a", aa, ".", names(fit.types)[fit.type.jj])) := exp(predict(fit.cox.types[[fit.type.jj]]["fit.cox"][[1]], newdata=tmp.long, type="lp"))]
                    tmp.long[, (paste0("P.a", aa, ".", names(fit.types)[fit.type.jj])) :=
                                   tmp.long[[paste0("dhazard.", names(fit.types)[fit.type.jj])]]*tmp.long[[paste0("exp.a", aa, ".", names(fit.types)[fit.type.jj])]]]
                }
            }
        }

        ###browser()

        if (length(depend.matrix[, unique(state)])>10) {
            if (length(depend.matrix[, unique(state)])<25) {
                print.state.jj <- floor(seq(1, length(depend.matrix[, unique(state)]),
                                            length = length(depend.matrix[, unique(state)])/3))
            } else{ 
                print.state.jj <- floor(seq(1, length(depend.matrix[, unique(state)]),
                                            length = length(depend.matrix[, unique(state)])/10))
            }
        } else {
            print.state.jj <- depend.matrix[, unique(state)]
        }

        tmp.long[, row_id := .I]
        #tmp.long.jj <-
        tmp.state <- tmp.long[, names(tmp.long) %in%
                                unique(
                                    c("state", "id", "row_id", "risk.time", "delta", "tstart", "time", "A0.obs",
                                      grep("dhazard.", names(tmp.long), value = TRUE),
                                      grep("at.risk", names(tmp.long), value = TRUE),
                                      "grid.time", varnames, N.vars)), with = FALSE]
        
        state.vec <- depend.matrix[, unique(state)]#sort(unique(depend.matrix[["state"]]))
        chunk.size <- ceiling(length(state.vec) / min(detectCores()-5, use.cores.prediction, use.cores))
        state.chunks <- split(state.vec, ceiling(seq_along(state.vec) / chunk.size))

        ## browser()
        tmp.list <- mclapply(
            state.chunks,
            function(state.chunk) {

                chunk.out <- #vector("list", length(state.chunk))
                    data.table(row_id = tmp.state[["row_id"]])

                chunk.out[["state"]] <- 0
                     
                for (state.jj in state.chunk) {

                    ## tmp.state <- copy(tmp.state)

                    if (verbose & (state.jj %in% print.state.jj))
                        print(paste0(state.jj, "/", max(depend.matrix[["state"]])))

                    ##print(state.jj)

                    ## tmp.state <- copy(tmp.long)

                    which.jj <- rep(TRUE, nrow(tmp.state))

                    for (varname in setdiff(names(depend.matrix), "state")) { #state.names

                        ##message("before assigning ", varname, " for state ", state.jj)
                        tmp.state[[varname]] <- depend.matrix[state == state.jj][[varname]]

                        if (!any.hal) #{
                            # which.jj <- which.jj*(tmp.long[[varname]] == tmp.state[[varname]])
                            #} else {
                            if (length(which.recurrent) > 0) {
                                if (varname %in% which.recurrent) {
                                    for (count.jj in 1:max.count) {
                                        tmp.state[, (paste0(varname, ".", count.jj)) := (tmp.state[[varname]] >= count.jj)]
                                    }
                                }
                            }
                        #   which.jj <- which.jj*(tmp.long[[varname]] >= tmp.state[[varname]])
                        which.jj <- which.jj & (tmp.long[[varname]] >= tmp.state[[varname]])

                    }

                    if (any.hal) {
                        X.hal.dynamic <- Matrix(
                            model.matrix(formula(paste0("delta", "~-1+",
                                                        paste(paste0("(", gsub("FALSE|TRUE", "",
                                                                               gsub(":", "):(", hal.vars.dynamic)), ")"),
                                                              collapse = "+"))), 
                                         data=tmp.state), sparse=TRUE)
                        X.hal.dynamic <- flip.interactions(X.hal.dynamic, hal.vars.dynamic)
                        if (length(hal.vars.A0.dynamic)>0) {
                            X.hal.A0.dynamic <- Matrix(
                                model.matrix(formula(paste0("delta", "~-1+",
                                                            paste(paste0("(", gsub("FALSE|TRUE", "",
                                                                                   gsub(":", "):(", hal.vars.A0.dynamic)), ")"),
                                                                  collapse = "+"))), 
                                             data=tmp.state), sparse=TRUE)
                            X.hal.A0.dynamic <- flip.interactions(X.hal.A0.dynamic, hal.vars.A0.dynamic)
                        }

                        X.full <- cbind(X.hal.static, X.hal.dynamic, X.hal.A0.dynamic, X.hal.A0)
                    }


                    for (fit.type.jj in (1:length(fit.types))[-cens.process.id]) {
                        if (intervene.A0) {
                            colname <- paste0("P.a", aa, ".", names(fit.types)[fit.type.jj], ".", state.jj)
                        } else {
                            colname <- paste0("P.", names(fit.types)[fit.type.jj], ".", state.jj)
                        }
                        if (fit.types[[fit.type.jj]]$fit == "hal") {
                            hal.index <- (1:length(which.hal))[names(which.hal) == names(fit.types)[fit.type.jj]]
                            ##message("fit.type.jj ", fit.type.jj, " for state ", state.jj)
                            chunk.out[[colname]] <- as.numeric(exp(predict(fit.hals[[hal.index]][["hal.fit"]], X.full[, hal.vars.list[[hal.index]]],
                                                                           newoffset=0, s=fit.hals[[hal.index]][["lambda.cv"]]))*tmp.state[["risk.time"]])
                        } else {
                            lp <- predict(fit.cox.types[[fit.type.jj]]["fit.cox"][[1]], newdata=tmp.state, type="lp")
                            chunk.out[[colname]] <- as.numeric(tmp.long[[paste0("dhazard.", names(fit.types)[fit.type.jj])]] * exp(lp))
                            ##tmp.state[, (paste0("exp.", names(fit.types)[fit.type.jj])) := exp()]
                            #tmp.long[, (paste0("P.", names(fit.types)[fit.type.jj], ".", state.jj)) :=
                            #               tmp.state[[paste0("dhazard.", names(fit.types)[fit.type.jj])]]*tmp.state[[paste0("exp.", names(fit.types)[fit.type.jj])]]]
                        }
                        if (fit.type.jj %in% at.risk.ids) { # <- FIX: may just want to remove again
                            chunk.out[[colname]] <- chunk.out[[colname]] *
                                fit.types[[fit.type.jj]][["at.risk"]](depend.matrix[state == state.jj])
                            #tmp.long[, (paste0("P.", names(fit.types)[fit.type.jj], ".", state.jj)) :=
                            #               tmp.long[[paste0("P.", names(fit.types)[fit.type.jj], ".", state.jj)]]*
                            #               fit.types[[fit.type.jj]][["at.risk"]](depend.matrix[state == state.jj])]
                        }
                        ## if (intervene.A0) {
                        ##     setnames(
                        ##         chunk.out,
                        ##         old = colname,
                        ##         new = paste0("P.a", aa, ".", names(fit.types)[fit.type.jj], ".", state.jj)
                        ##     )
                        ##     #names(chunk.out)[names(chunk.out) == colname] <-
                        ##     #    paste0("P.a", aa, ".", names(fit.types)[fit.type.jj], ".", state.jj)
                        ##     ## setnames(tmp.long, paste0("P.", names(fit.types)[fit.type.jj], ".", state.jj),
                        ##     ##        paste0("P.", "a", aa, ".", names(fit.types)[fit.type.jj], ".", state.jj))
                        ## }
                    }

                    ## outchunk.out[["state"]] <- 0
                    chunk.out[["state"]][which.jj == 1] <- state.jj

                    #tmp.long[
                    #    which.jj == 1,
                    #    state := state.jj]
                }

               
                
                return(chunk.out)
            },
            mc.cores = min(detectCores()-5, use.cores.prediction, use.cores))

        for (out in tmp.list) {

            # add P.* columns
            cols <- setdiff(names(out), c("row_id", "state"))
            tmp.long[, (cols) := out[, cols, with = FALSE]]

            # combine state
            tmp.long[, state := pmax(state, out$state)]
        }

    }

    names.P <- names(tmp.long)[substr(names(tmp.long), 1, 2) == "P."]

    for (name.P in names.P) {
        if (any(tmp.long[[name.P]]>1) | use.exponential) {
            if (!use.exponential) {
                if (verbose.exponential) print(paste0("transform ", name.P, " with 1-exp to avoid values >1"))
            }
            tmp.long[, (name.P) := 1-exp(-tmp.long[[name.P]])]
        }
    }

    if (length(derived.vars)>0) {
        for (derived.jj in 1:length(derived.vars)) {
            derived.var <- names(derived.vars)[derived.jj]
            derived.fun <- derived.vars[[derived.jj]]
            derived.cols <- extract_dt_cols(derived.fun)
            if (derived.var %in% names(depend.matrix)) {
                setnames(depend.matrix, derived.var, paste0(derived.var, " = ", paste0(derived.cols, collapse = ",")))
            }
            ## setnames(tmp.long, derived.var, paste0("At = ", paste0(At.cols, collapse = ",")))
        }
    }

    var.depend <- setdiff(names(depend.matrix), "state")
    process.depend <- var.depend[var.depend %in% process.names]

    ### browser()
     
    if (any.hal) { 
        if (length(T.vars <- grep("T.", names(non.zero.state.vars), value = TRUE))>0) {
            for (time.var in T.vars) {
                varname <- strsplit(time.var, "\\.")[[1]][2]
                count.value <- as.numeric(strsplit(time.var, "\\.")[[1]][3])
                if (is.na(count.value)) count.value <- 1
                time.values <- unique(c(0, tmp.long[[time.var2 <- gsub("T.", "time.", time.var)]]))
                depend.matrix[[time.var2]] <- Inf
                depend.matrix[depend.matrix[[varname]] == count.value-1][[time.var2]] <-
                    time.values[1]
                for (t2 in time.values[-1]) {
                    depend.matrix.t2 <- depend.matrix[depend.matrix[[varname]] == count.value-1 &
                                                      get(time.var2) == 0]
                    depend.matrix.t2[[time.var2]] <- t2
                    depend.matrix <- rbind(depend.matrix, depend.matrix.t2)[order(state, get(time.var2))]
                }
            }
        }
    } else if (length(depend.time)>0) {
        for (varname in process.depend) { # processess in states
            if (length(grep(varname, names(depend.time)))>0) {
                for (time.var in paste0("T.", grep(varname, names(depend.time), value = TRUE)## names(depend.time)
                                        )) {
                    count.value <- as.numeric(strsplit(time.var, "\\.")[[1]][3])
                    if (is.na(count.value)) count.value <- 1
                    time.values <- unique(tmp.long[[time.var2 <- gsub("T.", "time.", time.var)]])
                    depend.matrix[[time.var2]] <- Inf
                    depend.matrix[depend.matrix[[varname]] == count.value-1][[time.var2]] <-
                        time.values[1]
                    for (t2 in time.values[-1]) {
                        depend.matrix.t2 <- depend.matrix[depend.matrix[[varname]] == count.value-1 &
                                                          get(time.var2) == 0]
                        depend.matrix.t2[[time.var2]] <- t2
                        depend.matrix <- rbind(depend.matrix, depend.matrix.t2)[order(state, get(time.var2))]
                    }
                }
            }
        }
    }

    depend.matrix[, state.row.index := 1:.N]
    time.vars <- grep("time.", names(depend.matrix), value = TRUE)

    ###browser()

    for (state.jj in depend.matrix[, unique(state)]) {

        for (varname in process.depend) { # processess in states

            depend.state.jj <- depend.matrix[state == state.jj]

            var.value <- depend.state.jj[[varname]]
            var1.value <- pmin(var.value+1L, max(depend.matrix[[varname]]))

            depend.state.jj[[varname]] <- var1.value

            if (any.hal) { ## 
                if (length(T.varname.vars <-
                               grep(paste0("T.", varname), names(non.zero.state.vars), value = TRUE))>0) {
                    for (time.var in T.varname.vars) {
                        count.value <- as.numeric(strsplit(time.var, "\\.")[[1]][3])
                        if (is.na(count.value)) count.value <- 1
                        if (all(var.value == count.value-1)) {
                            depend.state.jj[[time.var]] <- depend.state.jj[[gsub("^T\\.", "time.", time.var)]]
                        }
                    }
                }
            } else if (length(depend.time)>0) {
                if (length(grep(varname, names(depend.time)))>0) {
                    for (time.var in grep(varname, paste0("T.", names(depend.time)), value = TRUE)) {
                        count.value <- as.numeric(strsplit(time.var, "\\.")[[1]][3])
                        if (is.na(count.value)) count.value <- 1
                        if (all(var.value == count.value-1)) {
                            depend.state.jj[[time.var]] <- depend.state.jj[[gsub("^T\\.", "time.", time.var)]]
                        }
                    }
                }
            }

            if (length(derived.vars)>0) { # <- this still needs checking 
                if (length(grep(varname, grep(names(derived.vars), names(depend.matrix), value = TRUE)))>0) {
                    for (derived.var in grep(varname, grep(names(derived.vars), names(depend.matrix), value = TRUE), value = TRUE)) {
                        depend.state.jj[[derived.var]] <- 1-depend.state.jj[[derived.var]]
                    }
                }
            }

            # # # # # #
            if (length(state.1 <- unique(
                           merge(depend.state.jj[, (names(depend.state.jj) %in% c("state.row.index", var.depend)), with = FALSE],
                                 unique(depend.matrix[, (names(depend.matrix) != "state.row.index" &
                                                         names(depend.state.jj) %in% c("state", var.depend)), with = FALSE]),
                                 by = var.depend))[order(state.row.index)][["state"]]) > 0) {
                depend.matrix[state == state.jj, (paste0("s.", varname)) := state.1]
            } else {
                depend.matrix[state == state.jj, (paste0("s.", varname)) := state.jj]
            }
        }

        if (length(time.vars)>0) {
            which.non.Inf <- time.vars[sapply(time.vars, function(time.var) (any(depend.matrix[state == state.jj][[time.var]]<Inf)))]
            if (length(which.non.Inf)>0) {
                ### NEED CORRECT ORDER HERE!
                #merge(tmp.long[state == state.jj],
                #      depend.matrix[state == state.jj][, c(which.non.Inf, "state.row.index"), with = FALSE],
                #      by = which.non.Inf)
                ## state.row.index
                ## THIS - below? 
                tmp.long[state == state.jj, state.row.index := 
                                                depend.matrix[state == state.jj][, c("state.row.index", which.non.Inf),
                                                                                 with = FALSE][tmp.long[state == state.jj], on = which.non.Inf][["state.row.index"]]]
            } else {
                tmp.long[state == state.jj, state.row.index :=
                                                depend.matrix[state == state.jj][["state.row.index"]]]
            }
        } else {
            tmp.long[state == state.jj, state.row.index :=
                                            depend.matrix[state == state.jj][["state.row.index"]]]
        }
    }

    depend.matrix <- depend.matrix[, !names(depend.matrix)[names(depend.matrix) == "state.row.index"], with = FALSE]
    
    return(list(
        tmp.long = tmp.long,
        depend.matrix = depend.matrix,
        process.names = process.names,
        process.deltas = process.deltas,
        process.types = process.types,
        cens.process.id = cens.process.id,
        at.risks = at.risks
    ))
}



######################################################################
### prepare.initial.R ends here


######################################################################
### faster version to get unique rows: 

library(digest)

hash_sparse_rows_dgC <- function(M) {
    stopifnot(inherits(M, "dgCMatrix"))

    p  <- M@p
    i  <- M@i
    x  <- M@x
    nr <- M@Dim[1]

    # storage: list of integer/values per row
    idx_list <- vector("list", nr)
    val_list <- vector("list", nr)

    # loop over columns, add entries to each row
    for (col in seq_len(M@Dim[2])) {
        start <- p[col] + 1
        end   <- p[col + 1]

        if (end >= start) {
            rows <- i[start:end] + 1        # row indices (1-based)
            vals <- x[start:end]            # values

            # append col index + value to each row's list
            for (k in seq_along(rows)) {
                r <- rows[k]
                idx_list[[r]] <- c(idx_list[[r]], col)
                val_list[[r]] <- c(val_list[[r]], vals[k])
            }
        }
    }

    # hash each row using only the sparse pattern
    out <- character(nr)
    for (r in seq_len(nr)) {
        out[r] <- digest::digest(
                              list(idx_list[[r]], val_list[[r]]),
                              algo = "xxhash64",
                              serialize = TRUE
                          )
    }

    out
}



######################################################################
### prepare.initial.R ends here
